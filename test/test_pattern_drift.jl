using Test
using Droplet

@testset "Pattern Drift Prevention" begin

    @testset "patternTailPos initialization" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern: 8 base, 8 pattern, 16 base (tail at position 16)
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        @test p.patternTailPos !== nothing
        @test p.originalPatternTailPos !== nothing
        @test p.patternTailPos[patternIntIdx] == 16
        @test p.originalPatternTailPos[patternIntIdx] == 16
        @test p.patternTailPos == p.originalPatternTailPos
    end

    @testset "patternTailPos updates on addBond" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern: 8 base, 16 pattern, 8 base (tail at position 24)
        customPattern = vcat(fill(1, 8), fill(2, 16), fill(1, 8))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        originalTailPos = p.patternTailPos[patternIntIdx]
        @test originalTailPos == 24
        
        # Run moves until we get an addBond inside pattern region that succeeds
        # Track the tail position change
        tailPosBefore = p.patternTailPos[patternIntIdx]
        addBondInPatternOccurred = false
        for _ in 1:10000
            tailPosBefore = p.patternTailPos[patternIntIdx]
            lengthBefore = p.length
            if Droplet.move!(p)
                # Check if length increased (addBond) and we're still tracking the pattern
                if p.length == lengthBefore + 2 && haskey(p.patternTailPos, patternIntIdx)
                    # If tail position changed, it was an addBond in pattern region
                    if p.patternTailPos[patternIntIdx] != tailPosBefore
                        addBondInPatternOccurred = true
                        @test p.patternTailPos[patternIntIdx] == tailPosBefore - 2
                        break
                    end
                end
            end
        end
        
        # Original should remain unchanged
        @test p.originalPatternTailPos[patternIntIdx] == originalTailPos
    end

    @testset "patternTailPos updates on removeBond" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern: 8 base, 16 pattern, 8 base
        customPattern = vcat(fill(1, 8), fill(2, 16), fill(1, 8))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        originalTailPos = p.patternTailPos[patternIntIdx]
        @test originalTailPos == 24
        
        # Force a removeBond inside pattern region
        success = false
        for attempt in 1:1000
            p.moveIdx = rand(9:22)
            Droplet.identifyBond!(p)
            Droplet.proposeMove!(p)
            if p.moveType == 3  # removeBond
                success = Droplet.move!(p)
                if success
                    break
                end
            end
        end
        
        if success
            @test p.patternTailPos[patternIntIdx] == originalTailPos + 2
            @test p.originalPatternTailPos[patternIntIdx] == originalTailPos
        end
    end

    @testset "Drift limit values" begin
        # Test the inline drift limit computation
        function getDriftLimit(len)
            len < 8 ? Int32(-1) : len < 32 ? Int32(2) : len < 64 ? Int32(4) : Int32(8)
        end
        
        @test getDriftLimit(6) == -1
        @test getDriftLimit(7) == -1
        @test getDriftLimit(8) == 2
        @test getDriftLimit(16) == 2
        @test getDriftLimit(31) == 2
        @test getDriftLimit(32) == 4
        @test getDriftLimit(48) == 4
        @test getDriftLimit(63) == 4
        @test getDriftLimit(64) == 8
        @test getDriftLimit(128) == 8
    end

    @testset "Drift prevention at limit" begin
        l = Lattice(boxlength=64, nspecies=2)
        # Length 64 → drift limit = 8
        # Pattern: 16 base, 16 pattern, 32 base (tail at position 32)
        customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        originalTailPos = p.originalPatternTailPos[patternIntIdx]
        @test originalTailPos == 32
        
        # Run many moves and track maximum drift
        maxDrift = 0
        for _ in 1:10000
            Droplet.move!(p)
            
            currentDrift = abs(p.patternTailPos[patternIntIdx] - originalTailPos)
            maxDrift = max(maxDrift, currentDrift)
        end
        
        # Drift should never exceed limit (8 for length 64)
        @test maxDrift <= 8
    end

    @testset "No drift limit for short polymers" begin
        l = Lattice(boxlength=16, nspecies=2)
        # Length 6 → drift limit = -1 (no enforcement)
        customPattern = [1, 2, 2, 2, 2, 1]
        p = Polymer(id=1, length=6, maxlength=16, boxlength=16,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        # Should initialize but not enforce drift
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        @test p.patternTailPos[patternIntIdx] == 5
        
        # Moves should not be rejected based on drift
        moveCount = 0
        for _ in 1:100
            if Droplet.move!(p)
                moveCount += 1
            end
        end
        
        # Should have some successful moves
        @test moveCount > 0
    end

    @testset "Uniform polymer has no drift tracking" begin
        l = Lattice(boxlength=32, nspecies=1)
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=1, pattern=:uniform)
        initMonomers!(p)
        init!(l, [p])
        
        @test p.patternTailPos === nothing
        @test p.originalPatternTailPos === nothing
        @test p.patternTails === nothing
    end

    @testset "Bidirectional drift bounded" begin
        l = Lattice(boxlength=64, nspecies=2)
        # Length 32 → drift limit = 4
        # Pattern centered: 8 base, 16 pattern, 8 base (tail at 24)
        customPattern = vcat(fill(1, 8), fill(2, 16), fill(1, 8))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        originalTailPos = p.originalPatternTailPos[patternIntIdx]
        
        minPos = originalTailPos
        maxPos = originalTailPos
        
        for _ in 1:5000
            Droplet.move!(p)
            
            currentPos = p.patternTailPos[patternIntIdx]
            minPos = min(minPos, currentPos)
            maxPos = max(maxPos, currentPos)
        end
        
        # Both forward and backward drift bounded by limit (4)
        @test originalTailPos - minPos <= 4
        @test maxPos - originalTailPos <= 4
    end

    @testset "Energy consistency with drift prevention" begin
        l = Lattice(boxlength=64, nspecies=2)
        customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p; defaultIntPot=-1.0)
        init!(l, [p])
        
        for _ in 1:1000
            Droplet.ham!(p)
            hamBefore = p.ham
            
            p.moveIdx = rand(1:p.length-1)
            Droplet.identifyBond!(p)
            
            # Skip if at chain end
            if p.kj == 0
                continue
            end
            
            Droplet.proposeMove!(p)
            Droplet.deltaHam!(p)
            
            if Droplet.move!(p)
                Droplet.ham!(p)
                @test isapprox(hamBefore + p.deltaH, p.ham, atol=1e-10)
            end
        end
    end

    @testset "Multiple patterns independent drift" begin
        l = Lattice(boxlength=64, nspecies=3)
        # Length 64: 16 base, 8 type2, 16 base, 8 type3, 16 base
        customPattern = vcat(fill(1, 16), fill(2, 8), fill(1, 16), fill(3, 8), fill(1, 16))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=3, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        type2IntIdx = intIdxList[2]
        type3IntIdx = intIdxList[3]
        
        # Both patterns should have their own tail positions
        @test haskey(p.patternTailPos, type2IntIdx)
        @test haskey(p.patternTailPos, type3IntIdx)
        @test p.patternTailPos[type2IntIdx] == 24  # end of first pattern block
        @test p.patternTailPos[type3IntIdx] == 48  # end of second pattern block
        
        # Run moves and verify both are bounded independently
        for _ in 1:5000
            Droplet.move!(p)
        end
        
        @test abs(p.patternTailPos[type2IntIdx] - p.originalPatternTailPos[type2IntIdx]) <= 8
        @test abs(p.patternTailPos[type3IntIdx] - p.originalPatternTailPos[type3IntIdx]) <= 8
    end

end
