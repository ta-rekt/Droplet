using Droplet
using Test

@testset "Droplet.jl" begin
    @testset "Init, identify" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        # polymer
        p = Poly(;
            length = 16,
            maxlength = 20,
            intPot = 1.0,
            chemPot = 1.0
        )

        # dilute monomers
        b = Bulk(;
            length = 256,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p, b])
        
        for i in 1:16
            @test p.latS[i,8,8] == 1
            @test l.lat[1,i,8,8] == 1
        end

        @test p.pos[1][4] == 0
        @test p.pos[5][4] == 4
        @test p.pos[5][5] == 6
        @test p.pos[16][5] == 0

        p.moveIdx = 5

        identifyBond!(p)

        @test p.ki == 5
        @test p.kj == 6
        @test p.kprev == 4
        @test p.knext == 7
        @test p.bi == [5,8,8]
        @test p.bj == [6,8,8]
        @test p.bprev == [4,8,8]
        @test p.bnext == [7,8,8]

        p.moveIdx = 1

        identifyBond!(p)

        @test p.ki == 1
        @test p.kj == 2
        @test p.kprev == 0

        p.moveIdx = 15

        identifyBond!(p)

        @test p.ki == 15
        @test p.kj == 16
        @test p.knext == 0
    end

    @testset "Movetypes" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        # polymer
        p = Poly(;
            length = 16,
            maxlength = 20,
            intPot = 1.0,
            chemPot = 1.0
        )

        # dilute monomers
        b = Bulk(;
            length = 256,
            intPot = 1.0,
            chemPot = 1.0
        )

        # initialize the system
        init!(l, [p, b])

        ### ### ### Adding a bond

        p.moveIdx = 5
        identifyBond!(p)
        p.moveType = 1  # bond addition
        p.moveVec = [0,0,1]
        p.newPosi = p.bi .+ p.moveVec
        p.newPosj = p.bj .+ p.moveVec
        
        # Original length
        origLength = p.length
        
        makeMove!(p)
        p.length = length(keys(p.pos))

        p.moveIdx = 3
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [1,1,1,1]

        p.moveIdx = 4
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [1,1,2,1]

        p.moveIdx = 5
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [2,2,1,1]

        p.moveIdx = 17
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [1,1,1,3]

        p.moveIdx = 18
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [2,2,1,1]

        
        p.moveIdx = 18
        identifyBond!(p)
        p.moveType = 1  # bond addition
        p.moveVec = [0,1,0]
        p.newPosi = p.bi .+ p.moveVec
        p.newPosj = p.bj .+ p.moveVec
        makeMove!(p)
        p.length = length(keys(p.pos))

        p.moveIdx = 18
        identifyBond!(p)
        moveTypes!(p)

        @test p.moveTypes == [1,2,1,2]
    end

    @testset "Hamiltonian (poly)" begin

        l = Lattice(;
            boxlength = 42,
            nspecies = 2
        )

        # polymer
        p = Poly(;
            length = 32,
            maxlength = 24,
            intPot = 1,
            chemPot = 1,
            initPos = [5,21,21]
        )

        slist = [p]

        # initialize lattice and species
        init!(l, slist)

        # initial energy
        hamInitial = p.ham

        #### iterating multiple moves
        nmoves = 0
        frame = 1
        nframes = 100
        step = 50
        deltaH = 0
        moveMade = false
        
        hams = zeros(nframes)
        accDeltaH = zeros(nframes)
        
        # first frame
        accDeltaH[frame] = hamInitial
        
        ham!(p)
        hams[frame] = hamInitial

        while frame <= nframes - 1
            
            if moveMade

                deltaH += p.deltaH

                if mod(nmoves, step) == 0
                
                    frame += 1

                    accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                    deltaH = 0
                    
                    ham!(p)
                    hams[frame] = p.ham

                    @test hams[frame] == accDeltaH[frame]
                end
            end
            
            moveMade = move!(p)
            nmoves += moveMade
        end
    end


    @testset "Hamiltonian (bulk)" begin

        l = Lattice(;
            boxlength = 32,
            nspecies = 1
        )

        # polymer
        p = Bulk(;
            length = 1000,
            intPot = -1,
            chemPot = 1
        )

        slist = [p]

        # initialize lattice and species
        init!(l, slist)

        # initial energy
        hamInitial = p.ham

        #### iterating multiple moves
        nmoves = 0
        frame = 1
        nframes = 100
        step = 50
        deltaH = 0
        moveMade = false
        
        hams = zeros(nframes)
        accDeltaH = zeros(nframes)
        
        # first frame
        accDeltaH[frame] = hamInitial
        
        ham!(p)
        hams[frame] = hamInitial

        while frame <= nframes - 1
            
            if moveMade

                deltaH += p.deltaH

                if mod(nmoves, step) == 0
                
                    frame += 1

                    accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                    deltaH = 0
                    
                    ham!(p)
                    hams[frame] = p.ham

                    @test hams[frame] == accDeltaH[frame]
                end
            end
            
            moveMade = move!(p)
            nmoves += moveMade
        end
    end

    @testset "Connected components (multiple)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 2,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 2
        )

        fb = FastBulk(;
            length = 2,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 2
        )

        init!(l, [b, fb])

        @test cc(b) == Set([Set([4,3]), Set([2,1])])
        @test cc(fb) == Set([Set([4,3]), Set([2,1])])
    end

    @testset "Connected components (single)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 20,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1
        )

        fb = FastBulk(;
            length = 20,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1
        )

        init!(l, [b, fb])

        fx = cc(fb)
        fy = Set{Set{Int32}}()
        fy0 = Set{Int32}(range(1,20))
        push!(fy, fy0)
    
        @test fx == fy

        x = cc(b)
        y = Set{Set{Int32}}()
        y0 = Set{Int32}(range(1,20))
        push!(y, y0)
    
        @test x == y
    end

    @testset "Connected components (line)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "line"
        )

        fb = FastBulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "line"
        )

        init!(l, [b, fb])

        x = cc(b)
        y = Set{Set{Int32}}()
        y0 = Set{Int32}(range(1,8))
        push!(y, y0)

        @test x == y

        fx = cc(fb)
        fy = Set{Set{Int32}}()
        fy0 = Set{Int32}(range(1,8))
        push!(fy, fy0)

        @test fx == fy
    end

    @testset "Box density (uniform)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 256,
            intPot = -1,
            chemPot = 1,
            initShape = "uniform",
            initParams = 2
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        @test boxDensity(b, [[Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength]]) == b.length / b.boxlength ^ 3
    end

    @testset "Box density (full box)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 27,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        y = b.length / 3 ^ 3
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y

        b = Bulk(;
            length = 26,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        y = 26/27
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y

        y = 7/27
        x = boxDensity(b, [[Int32(2), Int32(4)], 
                        [Int32(2), Int32(4)], 
                        [Int32(2), Int32(4)]])

        @test x == y

        y = 0.0
        x = boxDensity(b, [[Int32(3), Int32(5)], 
                        [Int32(3), Int32(5)], 
                        [Int32(3), Int32(5)]])

        @test x == y

        y = 0.0
        x = boxDensity(b, [[Int32(4), Int32(6)], 
                        [Int32(4), Int32(6)], 
                        [Int32(4), Int32(6)]])

        @test x == y


        b = Bulk(;
            length = 27,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        y = 1.0
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(0), Int32(2)]])
        @test x == y

        y = 2/3
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y
    end

end
