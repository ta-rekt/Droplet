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
            @test l.lat[2,i,8,8] == 1
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

    @testset "Hamiltonian (multi-component bulk)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        a = Bulk(;
            length = 128,
            intPot = -1,
            intPotP = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        b = Bulk(;
            length = 128,
            intPot = -1,
            intPotP = -1,
            chemPot = 1,
            initShape = "uniform"
        )

        b.partner = a
        a.partner = b

        slist = [a, b]

        # initialize lattice and species
        init!(l, slist)

        #### iterating multiple moves
        nmoves = [0,0]
        frame = [1,1]
        frame_joint = 1
        nframes = 10000
        step = 1
        deltaH = [0,0]
        niter = 0
        moveMade = [false, false]

        hamsA = Vector{Float64}(undef, 0)
        hamsB = Vector{Float64}(undef, 0)
        hams = [hamsA, hamsB]
        accDeltaHA = Vector{Float64}(undef, 0)
        accDeltaHB = Vector{Float64}(undef, 0)
        accDeltaH = [accDeltaHA, accDeltaHB]
        
        # first frame, initial energy
        push!(hams[1], a.ham)
        push!(hams[2], b.ham)
        push!(accDeltaH[1], a.ham)
        push!(accDeltaH[2], b.ham)

        ham!(a)
        ham!(b)
        a.deltaH = 0
        b.deltaH = 0

        a.name = string(a.name, "_a")
        b.name = string(b.name, "_b")
        color = ["blue (A)", "orange (B)"]

        while frame_joint <= nframes - 2
            for (i, s) in enumerate(slist)
                if move!(s)
                    nmoves[i] += 1
                    deltaH[i] += s.deltaH
                    ip = i == 1 ? 2 : 1

                    # update deltaH for the partner
                    deltaH[ip] += s.partner.deltaHP
                    
                    push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                    deltaH[i] = 0
                    
                    ham!(s)
                    push!(hams[i], s.ham)

                    frame[i] += 1
                    frame_joint += 1
                end
            end
        end

        @test hams[1] == accDeltaH[1]
        @test hams[2] == accDeltaH[2]
    end


        @testset "Hamiltonian (multi-component fast bulk)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        a = FastBulk(;
            length = 128,
            intPot = -1,
            intPotP = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        b = FastBulk(;
            length = 128,
            intPot = -1,
            intPotP = -1,
            chemPot = 1,
            initShape = "uniform"
        )

        b.partner = a
        a.partner = b

        slist = [a, b]

        # initialize lattice and species
        init!(l, slist)

        #### iterating multiple moves
        nmoves = [0,0]
        frame = [1,1]
        frame_joint = 1
        nframes = 10000
        step = 1
        deltaH = [0,0]
        niter = 0
        moveMade = [false, false]

        hamsA = Vector{Float64}(undef, 0)
        hamsB = Vector{Float64}(undef, 0)
        hams = [hamsA, hamsB]
        accDeltaHA = Vector{Float64}(undef, 0)
        accDeltaHB = Vector{Float64}(undef, 0)
        accDeltaH = [accDeltaHA, accDeltaHB]
        
        # first frame, initial energy
        push!(hams[1], a.ham)
        push!(hams[2], b.ham)
        push!(accDeltaH[1], a.ham)
        push!(accDeltaH[2], b.ham)

        ham!(a)
        ham!(b)
        a.deltaH = 0
        b.deltaH = 0

        a.name = string(a.name, "_a")
        b.name = string(b.name, "_b")
        color = ["blue (A)", "orange (B)"]

        while frame_joint <= nframes - 2
            for (i, s) in enumerate(slist)
                if move!(s)
                    nmoves[i] += 1
                    deltaH[i] += s.deltaH
                    ip = i == 1 ? 2 : 1

                    # update deltaH for the partner
                    deltaH[ip] += s.partner.deltaHP
                    
                    push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                    deltaH[i] = 0
                    
                    ham!(s)
                    push!(hams[i], s.ham)

                    frame[i] += 1
                    frame_joint += 1
                end
            end
        end

        @test hams[1] == accDeltaH[1]
        @test hams[2] == accDeltaH[2]
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

        @test cc(b) == Vector([Set([4,3]), Set([2,1])])
        @test cc(fb) == Vector([Set([4,3]), Set([2,1])])
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
        fy = Vector{Set{Int32}}()
        fy0 = Set{Int32}(range(1,20))
        push!(fy, fy0)
    
        @test fx == fy

        x = cc(b)
        y = Vector{Set{Int32}}()
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
        y = Vector{Set{Int32}}()
        y0 = Set{Int32}(range(1,8))
        push!(y, y0)

        @test x == y

        fx = cc(fb)
        fy = Vector{Set{Int32}}()
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

    @testset "Split droplet" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5, 6]), Set([4, 3]), Set([7, 8]), Set([2, 1])]

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5]), Set([4]), Set([6]), Set([7]), 
                      Set([2]), Set([8]), Set([3]), Set([1])]

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        
        @test dlist == [Set([5, 6, 7, 8]), Set([4, 2, 3, 1])]

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])

        @test dlist == [Set([5, 4, 6, 7, 2, 8, 3, 1])]
    end


    @testset "Droplet bbox" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[6, 9], [6, 9], [6, 9]]

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,3],[1,3],[1,3]]

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,2,2]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,4],[1,3],[1,3]]
    end


    @testset "Fill/clear bbox" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        for i in range(1, 100)

            b = Bulk(;
                length = i,
                intPot = -1,
                chemPot = 1,
                initShape = "box",
                initParams = 1,
                initPos = [3,3,3]
            )

            slist = [b]

            # initialize lattice and species
            init!(l, slist)

            ccs = cc(b)

            bbox = dropletBbox(b, ccs[1])

            fillBbox!(b, ccs[1], bbox)

            @test length(collect(keys(b.ghost.pos))) == vol(bbox) - length(ccs[1])
            @test sum(map(x -> x > 0, b.ghost.latS)) == vol(bbox) - length(ccs[1])

            clearBbox!(b, bbox)

            @test iszero(b.ghost.latS)
            @test isempty(b.ghost.pos)
        end
    end

    @testset "Droplet volume" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 64,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
             [3,2,2],[3,2,3],[3,3,2],[3,3,3],[1,1,1],[2,1,2]]  # poke a hole in it
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 54

        b = Bulk(;
            length = 64,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
             [3,2,2],[3,2,3],[3,3,2],[3,3,3],[1,1,1]]  # no hole: outer volume preserved
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 63


        b = Bulk(;
            length = 64,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
                 [3,2,2],[3,2,3],[3,3,2],[3,3,3]]
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 64

        b = Bulk(;
            length = 24,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = Bulk(;
            length = 29,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = Bulk(;
            length = 2,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))

        b = Bulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))
    end

    @testset "Movable indices" begin
        l = Lattice(;
            boxlength = 5,
            nspecies = 1
        )

        b = FastBulk(;
            length = 27,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        for i in range(1, 27)
            if i == 14
                @test !(i in b.movableIdxs)
            else
                @test i in b.movableIdxs
            end
        end

        # exportSpecies!(slist, string(outpath, "_frame$(1)"))

        ########### corner ###########

        b.moveIdx = 1

        # find coordinates and keys of selected bond
        identifyBond!(b)
        @test b.bi == [2,2,2]

        nnIdxs!(b, b.bi)
        @test b.freeNNs == Set([2,4,6])

        for i in range(1, 3)
            b.tempVec[i] = b.pos[b.moveIdx][i]
        end
        
        # test each move direction
        b.moveType = 2
        axis = div(b.moveType + 1, 2)
        dir = - mod(b.moveType + 1, 2) * 2 + 1
        
        b.tempVec[axis] += dir
        @test b.tempVec == [1,2,2]

        # clip to within bounding box
        boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
        @test !overlap(b)

        makeMove!(b)
        @test b.bi == [1,2,2]


        ########### edge ###########

        b = FastBulk(;
            length = 27,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        b.moveIdx = 2

        # find coordinates and keys of selected bond
        identifyBond!(b)
        @test b.bi == [2,2,3]

        nnIdxs!(b, b.bi)
        @test b.freeNNs == Set([2,4])

        for i in range(1, 3)
            b.tempVec[i] = b.pos[b.moveIdx][i]
        end
        
        # test each move direction
        b.moveType = 4
        axis = div(b.moveType + 1, 2)
        dir = - mod(b.moveType + 1, 2) * 2 + 1
        
        b.tempVec[axis] += dir
        @test b.tempVec == [2,1,3]

        # clip to within bounding box
        boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
        @test !overlap(b)

        makeMove!(b)
        @test b.bi == [2,1,3]


        ########### center ###########

        b = FastBulk(;
            length = 27,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        b.moveIdx = 14

        # find coordinates and keys of selected bond
        identifyBond!(b)
        @test b.bi == [3,3,3]

        nnIdxs!(b, b.bi)
        @test isempty(b.freeNNs)


        ########### across boundary ###########

        b = FastBulk(;
            length = 8,
            intPot = -1,
            chemPot = 1,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # initialize lattice and species
        init!(l, slist)

        b.moveIdx = 8

        # find coordinates and keys of selected bond
        identifyBond!(b)
        @test b.bi == [1,1,1]

        nnIdxs!(b, b.bi)
        @test b.freeNNs == Set([1,3,5])

        for i in range(1, 3)
            b.tempVec[i] = b.pos[b.moveIdx][i]
        end
        
        # test each move direction
        b.moveType = 2
        axis = div(b.moveType + 1, 2)
        dir = - mod(b.moveType + 1, 2) * 2 + 1
        
        b.tempVec[axis] += dir

        # clip to within bounding box
        boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
        @test overlap(b)
    end
end
