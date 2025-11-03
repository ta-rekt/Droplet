using Droplet
using Test

@testset "Droplet.jl" begin
    @testset "Init, identify" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        # polymer
        s1 = Species(;
            length = 16,
            maxlength = 20,
            boxlength = l.boxlength,
            intPot = 1.0,
            chemPot = 1.0
        )

        # dilute monomers
        s2 = Species(;
            length = 256,
            maxlength = 256,
            boxlength = l.boxlength,
            intPot = 1.0,
            chemPot = 1.0
        )

        initPoly!(s1)
        initBulk!(s2)
        initLattice!(l, [s1, s2])
        
        for i in 1:16
            @test s1.latS[i,8,8] == 1
            @test l.lat[1,i,8,8] == 1
        end

        @test s1.pos[1][4] == 0
        @test s1.pos[5][4] == 4
        @test s1.pos[5][5] == 6
        @test s1.pos[16][5] == 0

        s1.moveIdx = 5

        identifyBond!(s1)

        @test s1.ki == 5
        @test s1.kj == 6
        @test s1.kprev == 4
        @test s1.knext == 7
        @test s1.bi == [5,8,8]
        @test s1.bj == [6,8,8]
        @test s1.bprev == [4,8,8]
        @test s1.bnext == [7,8,8]

        s1.moveIdx = 1

        identifyBond!(s1)

        @test s1.ki == 1
        @test s1.kj == 2
        @test s1.kprev == 0

        s1.moveIdx = 15

        identifyBond!(s1)

        @test s1.ki == 15
        @test s1.kj == 16
        @test s1.knext == 0
    end

    @testset "Movetypes" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        # polymer
        s1 = Species(;
            length = 16,
            maxlength = 20,
            boxlength = l.boxlength,
            intPot = 1.0,
            chemPot = 1.0
        )

        # dilute monomers
        s2 = Species(;
            length = 256,
            maxlength = 256,
            boxlength = l.boxlength,
            intPot = 1.0,
            chemPot = 1.0
        )

        # initializing lattice and species
        initPoly!(s1)
        initBulk!(s2)
        initLattice!(l, [s1, s2])

        ### ### ### Adding a bond

        s1.moveIdx = 5
        identifyBond!(s1)
        s1.moveType = 1  # bond addition
        s1.moveVec = [0,0,1]
        s1.newPosi = s1.bi .+ s1.moveVec
        s1.newPosj = s1.bj .+ s1.moveVec
        
        # Original length
        origLength = s1.length
        
        makeMove!(s1)
        s1.length = length(keys(s1.pos))

        s1.moveIdx = 3
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [1,1,1,1]

        s1.moveIdx = 4
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [1,1,2,1]

        s1.moveIdx = 5
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [2,2,1,1]

        s1.moveIdx = 17
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [1,1,1,3]

        s1.moveIdx = 18
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [2,2,1,1]

        
        s1.moveIdx = 18
        identifyBond!(s1)
        s1.moveType = 1  # bond addition
        s1.moveVec = [0,1,0]
        s1.newPosi = s1.bi .+ s1.moveVec
        s1.newPosj = s1.bj .+ s1.moveVec
        makeMove!(s1)
        s1.length = length(keys(s1.pos))

        s1.moveIdx = 18
        identifyBond!(s1)
        moveTypes!(s1)

        @test s1.moveTypes == [1,2,1,2]
    end

    @testset "Hamiltonian" begin

        global l = Lattice(;
            boxlength = 42,
            nspecies = 2
        )

        # polymer
        global s1 = Species(;
            length = 32,
            maxlength = 24,
            boxlength = l.boxlength,
            intPot = 1,
            chemPot = 1,
            initPos = [5,21,21]
        )

        slist = [s1]

        # initialize lattice and species
        initPoly!(s1)
        initLattice!(l, slist)
        
        # initialize hamiltonians
        initHamPoly!(s1)
        hamInitial = s1.ham

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
        
        hamPoly!(s1)
        hams[frame] = hamInitial

        while frame <= nframes - 1
            
            if moveMade

                deltaH += s1.deltaH

                if mod(nmoves, step) == 0
                
                    frame += 1

                    accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                    deltaH = 0
                    
                    hamPoly!(s1)
                    hams[frame] = s1.ham

                    @test hams[frame] == accDeltaH[frame]
                end
            end
            
            moveMade = movePoly!(s1)
            nmoves += moveMade
        end
    end
end