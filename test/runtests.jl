using Droplet
using Test

# MARK: -Droplet
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
            length = 256
        )
        
        # Set interaction potential after init
        intPot!(b, 1.0)

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

# MARK: -Droplet
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
            length = 256
        )

        # Set interaction potential
        intPot!(b, 1.0)

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

# MARK: -Ham poly
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
                end
            end
            
            moveMade = move!(p)
            nmoves += moveMade
        end

        @test hams == accDeltaH
    end

# MARK: -Ham bulk
    @testset "Hamiltonian (bulk)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        # bulk species
        p = Bulk(;
            length = 500
        )

        slist = [p]

        # Set interaction potential after init
        intPot!(p, -1.0)

        # initialize lattice and species
        init!(l, slist)

        # initial energy
        hamInitial = p.ham

        #### iterating multiple moves
        nmoves = 0
        frame = 1
        nframes = 1000
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
                end
            end
            
            moveMade = move!(p)
            nmoves += moveMade
        end

        @test hams == accDeltaH
    end

# MARK: -Ham fastbulk
    @testset "Hamiltonian (fast bulk)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        # bulk species
        p = FastBulk(;
            length = 500
        )

        slist = [p]

        # Set interaction potential after init
        intPot!(p, -1.0)

        # initialize lattice and species
        init!(l, slist)

        # initial energy
        hamInitial = p.ham

        #### iterating multiple moves
        nmoves = 0
        frame = 1
        nframes = 1000
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

                end
            end
            
            moveMade = move!(p)
            nmoves += moveMade
        end
        
        @test hams == accDeltaH
    end

# MARK: -Ham multi
    @testset "Hamiltonian (multi-component bulk)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        a = Bulk(;
            length = 128,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        b = Bulk(;
            length = 128,
            initShape = "uniform"
        )

        c = Bulk(;
            length = 128,
            initShape = "uniform"
        )

        # self and interaction potentials
        intPot!(a, -1.0)
        intPot!(b, -2.0)
        intPot!(c, -3.0)
        intPot!(a, b, -1.5)
        intPot!(a, c, -0.5)

        slist = [a, b, c]

        # initialize lattice and species
        init!(l, slist)

        #### iterating multiple moves
        nmoves = [0,0,0]
        frame = [1,1,1]
        frame_joint = 1
        nframes = 10000
        step = 1
        deltaH = [0.0,0.0,0.0]
        niter = 0
        moveMade = [false, false, false]

        hamsA = Vector{Float64}(undef, 0)
        hamsB = Vector{Float64}(undef, 0)
        hamsC = Vector{Float64}(undef, 0)
        hams_ = [hamsA, hamsB, hamsC]
        hams = []
        accDeltaHA = Vector{Float64}(undef, 0)
        accDeltaHB = Vector{Float64}(undef, 0)
        accDeltaHC = Vector{Float64}(undef, 0)
        accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
        accDeltaH = []
        
        # first frame, initial energy
        push!(hams_[1], a.ham)
        push!(hams_[2], b.ham)
        push!(hams_[3], c.ham)
        push!(accDeltaH_[1], a.ham)
        push!(accDeltaH_[2], b.ham)
        push!(accDeltaH_[3], c.ham)

        for (sp, idx) in a.intSpecies
            push!(hams, hams_[idx])
            push!(accDeltaH, accDeltaH_[idx])
        end

        ham!(a)
        ham!(b)
        ham!(c)
        a.deltaH = 0.0
        b.deltaH = 0.0
        c.deltaH = 0.0

        while frame_joint <= nframes - 2
            for (s, i) in a.intSpecies
                if move!(s)
                    nmoves[i] += 1
                    deltaH[i] += s.deltaH

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end

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
        @test hams[3] == accDeltaH[3]
    end

# MARK: -Ham multi fast
    @testset "Hamiltonian (multi-component fast bulk)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        a = FastBulk(;
            length = 128,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        b = FastBulk(;
            length = 128,
            initShape = "uniform"
        )

        c = FastBulk(;
            length = 128,
            initShape = "uniform"
        )

        # self and interaction potentials
        intPot!(a, -1.0)
        intPot!(b, -2.0)
        intPot!(c, -3.0)
        intPot!(a, b, -1.5)
        intPot!(a, c, -0.5)

        slist = [a, b, c]

        # initialize lattice and species
        init!(l, slist)

        #### iterating multiple moves
        nmoves = [0,0,0]
        frame = [1,1,1]
        frame_joint = 1
        nframes = 10000
        step = 1
        deltaH = [0.0,0.0,0.0]
        niter = 0
        moveMade = [false, false, false]

        hamsA = Vector{Float64}(undef, 0)
        hamsB = Vector{Float64}(undef, 0)
        hamsC = Vector{Float64}(undef, 0)
        hams_ = [hamsA, hamsB, hamsC]
        hams = []
        accDeltaHA = Vector{Float64}(undef, 0)
        accDeltaHB = Vector{Float64}(undef, 0)
        accDeltaHC = Vector{Float64}(undef, 0)
        accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
        accDeltaH = []
        
        # first frame, initial energy
        push!(hams_[1], a.ham)
        push!(hams_[2], b.ham)
        push!(hams_[3], c.ham)
        push!(accDeltaH_[1], a.ham)
        push!(accDeltaH_[2], b.ham)
        push!(accDeltaH_[3], c.ham)

        for (sp, idx) in a.intSpecies
            push!(hams, hams_[idx])
            push!(accDeltaH, accDeltaH_[idx])
        end

        ham!(a)
        ham!(b)
        ham!(c)
        a.deltaH = 0.0
        b.deltaH = 0.0
        c.deltaH = 0.0

        while frame_joint <= nframes - 2
            for (s, i) in a.intSpecies
                if move!(s)
                    nmoves[i] += 1
                    deltaH[i] += s.deltaH

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end

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
        @test hams[3] == accDeltaH[3]
    end

# MARK: -CC multi
    @testset "Connected components (multiple)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 2,
            initShape = "box",
            initParams = 2
        )

        fb = FastBulk(;
            length = 2,
            initShape = "box",
            initParams = 2
        )

        # Set interaction potentials after init
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

        init!(l, [b, fb])

        @test cc(b) == Vector([Set([4,3]), Set([2,1])])
        @test cc(fb) == Vector([Set([4,3]), Set([2,1])])
    end

# MARK: -CC single
    @testset "Connected components (single)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 20,
            initShape = "box",
            initParams = 1
        )

        fb = FastBulk(;
            length = 20,
            initShape = "box",
            initParams = 1
        )

        # Set interaction potentials after init
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

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

# MARK: -CC line
    @testset "Connected components (line)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        b = Bulk(;
            length = 8,
            initShape = "line"
        )

        fb = FastBulk(;
            length = 8,
            initShape = "line"
        )

        # Set interaction potentials 
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

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

# MARK: -Density uniform
    @testset "Box density (uniform)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 256,
            initShape = "uniform",
            initParams = 2
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        @test boxDensity(b, [[Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength]]) == b.length / b.boxlength ^ 3
    end

# MARK: -Density full
    @testset "Box density (full box)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 27,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        y = b.length / 3 ^ 3
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y

        b = Bulk(;
            length = 26,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential
        intPot!(b, -1.0)

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
            initShape = "box",
            initParams = 1,
            initPos = [2,2,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

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

# MARK: -Split droplet
    @testset "Split droplet" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5, 6]), Set([4, 3]), Set([7, 8]), Set([2, 1])]

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)
    
        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5]), Set([4]), Set([6]), Set([7]), 
                      Set([2]), Set([8]), Set([3]), Set([1])]

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        
        @test dlist == [Set([5, 6, 7, 8]), Set([4, 2, 3, 1])]

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])

        @test dlist == [Set([5, 4, 6, 7, 2, 8, 3, 1])]
    end

# MARK: -Droplet bbox
    @testset "Droplet bbox" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[6, 9], [6, 9], [6, 9]]

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,3],[1,3],[1,3]]

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [3,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,4],[1,3],[1,3]]
    end

# MARK: -Fill/clear
    @testset "Fill/clear bbox" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        for i in range(1, 100)

            b = Bulk(;
                length = i,
                initShape = "box",
                initParams = 1,
                initPos = [3,3,3]
            )

            slist = [b]

            # Set interaction potential 
            intPot!(b, -1.0)

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
# MARK: -Droplet volume
    @testset "Droplet volume" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = Bulk(;
            length = 64,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

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
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

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
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

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
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = Bulk(;
            length = 29,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = Bulk(;
            length = 2,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))

        b = Bulk(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))
    end
# # MARK: -Movable indices
#     @testset "Movable indices" begin
#         l = Lattice(;
#             boxlength = 5,
#             nspecies = 1
#         )

#         b = FastBulk(;
#             length = 27,
#             initShape = "box",
#             initParams = 1,
#             initPos = [3,3,3]
#         )

#         slist = [b]

#         # Set interaction potential 
#         intPot!(b, -1.0)

#         # initialize lattice and species
#         init!(l, slist)

#         for i in range(1, 27)
#             if i == 14
#                 @test !(i in b.movableIdxs)
#             else
#                 @test i in b.movableIdxs
#             end
#         end

#         ########### corner ###########

#         b.moveIdx = 1

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [2,2,2]

#         nnIdxs!(b, b.bi)
#         @test b.freeNNs == Set([2,4,6])

#         for i in range(1, 3)
#             b.tempVec[i] = b.pos[b.moveIdx][i]
#         end
        
#         # test each move direction
#         b.moveType = 2
#         axis = div(b.moveType + 1, 2)
#         dir = - mod(b.moveType + 1, 2) * 2 + 1
        
#         b.tempVec[axis] += dir
#         @test b.tempVec == [1,2,2]

#         # clip to within bounding box
#         boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
#         @test !overlap(b)

#         makeMove!(b)
#         @test b.bi == [1,2,2]


#         ########### edge ###########

#         b = FastBulk(;
#             length = 27,
#             initShape = "box",
#             initParams = 1,
#             initPos = [3,3,3]
#         )

#         slist = [b]

#         # initialize lattice and species
#         init!(l, slist)

#         b.moveIdx = 2

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [2,2,3]

#         nnIdxs!(b, b.bi)
#         @test b.freeNNs == Set([2,4])

#         for i in range(1, 3)
#             b.tempVec[i] = b.pos[b.moveIdx][i]
#         end
        
#         # test each move direction
#         b.moveType = 4
#         axis = div(b.moveType + 1, 2)
#         dir = - mod(b.moveType + 1, 2) * 2 + 1
        
#         b.tempVec[axis] += dir
#         @test b.tempVec == [2,1,3]

#         # clip to within bounding box
#         boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
#         @test !overlap(b)

#         makeMove!(b)
#         @test b.bi == [2,1,3]


#         ########### center ###########

#         b = FastBulk(;
#             length = 27,
#             initShape = "box",
#             initParams = 1,
#             initPos = [3,3,3]
#         )

#         slist = [b]

#         # initialize lattice and species
#         init!(l, slist)

#         b.moveIdx = 14

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [3,3,3]

#         nnIdxs!(b, b.bi)
#         @test isempty(b.freeNNs)


#         ########### across boundary ###########

#         b = FastBulk(;
#             length = 8,
#             initShape = "box",
#             initParams = 1,
#             initPos = [1,1,1]
#         )

#         slist = [b]

#         # initialize lattice and species
#         init!(l, slist)

#         b.moveIdx = 8

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [1,1,1]

#         nnIdxs!(b, b.bi)
#         @test b.freeNNs == Set([1,3,5])

#         for i in range(1, 3)
#             b.tempVec[i] = b.pos[b.moveIdx][i]
#         end
        
#         # test each move direction
#         b.moveType = 2
#         axis = div(b.moveType + 1, 2)
#         dir = - mod(b.moveType + 1, 2) * 2 + 1
        
#         b.tempVec[axis] += dir

#         # clip to within bounding box
#         boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
#         @test overlap(b)
#     end
# # MARK: -Hastings factor
#     @testset "Hastings factor" begin
#         ############ full box ############

#         l = Lattice(;
#             boxlength = 4,
#             nspecies = 2
#         )

#         b = FastBulk(;
#             length = 64,
#             initShape = "box",
#             initParams = 1,
#             initPos = [4,4,4]
#         )

#         slist = [b]

#         # Set interaction potential 
#         intPot!(b, -1.0)

#         # initialize lattice and species
#         init!(l, slist)
        
#         @test isempty(b.movableIdxs)
        
#         # surgically remove one bead
#         remKey = 17
#         b.latS[b.pos[remKey][1], b.pos[remKey][2], b.pos[remKey][3]] = 0
#         delete!(b.pos, remKey)

#         # re-initialize movable set
#         initMovableSet!(b)
#         @test b.movableIdxs == Set([20, 18, 29, 21, 1, 33])

#         b.moveIdx = 20

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [3,2,1]

#         # sample a move uniformly in each direction
#         proposeMove!(b)
#         @test b.hastingsFactor == 1.0

#         @test !overlap(b)

#         deltaHam!(b)
#         @test b.deltaH == 0.0

#         if acceptMove(b)
#             makeMove!(b)
#         end

#         @test b.movableIdxs == Set([20, 19, 32, 24, 36, 4])


#         ############ droplet boundary ############

#         l = Lattice(;
#             boxlength = 4,
#             nspecies = 2
#         )

#         b = FastBulk(;
#             length = 7,
#             initShape = "box",
#             initParams = 1,
#             initPos = [3,3,3]
#         )

#         slist = [b]

#         # Set interaction potential 
#         intPot!(b, -1.0)

#         # initialize lattice and species
#         init!(l, slist)
        
#         @test b.movableIdxs == Set(range(1, 7))

#         b.moveIdx = 6

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [3,2,3]

#         # sample a move uniformly in each direction
#         nnIdxs!(b, b.bi)
#         @test b.nnIdxs == [0,2,0,0,0,5]
#         @test b.freeNNs == Set([1,3,4,5])

#         for i in range(1, 3)
#             b.tempVec[i] = b.pos[b.moveIdx][i]
#         end

#         # propose a move uniformly among the free sites
#         b.moveType = 3
#         axis = div(b.moveType + 1, 2)
#         dir = - mod(b.moveType + 1, 2) * 2 + 1
#         b.tempVec[axis] += dir

#         @test b.tempVec == [3,3,3]

#         # clip to within bounding box
#         boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)

#         hastingsFactor!(b)
#         @test b.hastingsFactor == 1.0

#         deltaHam!(b)
#         @test b.deltaH == 0.0

#         makeMove!(b)

#         @test b.movableIdxs == Set(range(1, 7))


#         ############ droplet boundary 2 ############

#         l = Lattice(;
#             boxlength = 4,
#             nspecies = 2
#         )

#         b = FastBulk(;
#             length = 7,
#             initShape = "box",
#             initParams = 1,
#             initPos = [3,3,3]
#         )

#         slist = [b]

#         # Set interaction potential 
#         intPot!(b, -1.0)

#         # initialize lattice and species
#         init!(l, slist)
        
#         @test b.movableIdxs == Set(range(1, 7))

#         b.moveIdx = 6

#         # find coordinates and keys of selected bond
#         identifyBond!(b)
#         @test b.bi == [3,2,3]

#         # sample a move uniformly in each direction
#         nnIdxs!(b, b.bi)
#         @test b.nnIdxs == [0,2,0,0,0,5]
#         @test b.freeNNs == Set([1,3,4,5])

#         for i in range(1, 3)
#             b.tempVec[i] = b.pos[b.moveIdx][i]
#         end

#         # propose a move uniformly among the free sites
#         b.moveType = 4
#         axis = div(b.moveType + 1, 2)
#         dir = - mod(b.moveType + 1, 2) * 2 + 1
#         b.tempVec[axis] += dir

#         @test b.tempVec == [3,1,3]

#         # clip to within bounding box
#         boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)

#         hastingsFactor!(b)
#         @test b.hastingsFactor == 2/3

#         deltaHam!(b)
#         @test b.deltaH == 2.0

#         makeMove!(b)

#         @test b.movableIdxs == Set(range(1, 7))
#     end

# MARK: -intPot self
    @testset "intPot! (self)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 1
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Test setting self interaction potential
        intPot!(a, -2.0)

        init!(l, [a])
        @test a.intPot == -2.0

        # Test with Bulk type as well
        b = Bulk(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        init!(l, [b])
        
        intPot!(b, -4.0)
        @test b.intPot == -4.0

        intIdxs = map(x -> x[2], b.intSpecies)

        @test b.intIdx in intIdxs
    end

# MARK: -intPot two species
    @testset "intPot! (two species)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set cross-interaction potential
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        # Both species should have interaction structures
        @test a.intSpecies !== nothing
        @test b.intSpecies !== nothing
        @test a.intMat !== nothing
        @test b.intMat !== nothing
        @test a.intIdx !== nothing
        @test b.intIdx !== nothing

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB

        # They should share the same matrix
        @test a.intMat === b.intMat

        # Matrix should be 2x2
        @test size(a.intMat) == (2, 2)

        # Matrix should be symmetric
        @test a.intMat[a.intIdx, b.intIdx] == -3.0
        @test a.intMat[b.intIdx, a.intIdx] == -3.0

        # Diagonal entries should equal self intPot
        @test a.intMat[a.intIdx, a.intIdx] == a.intPot
        @test a.intMat[b.intIdx, b.intIdx] == b.intPot
    end

# MARK: -intPot three species
    @testset "intPot! (three species)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 3
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulk(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set pairwise interactions
        intPot!(a, b, -4.0)
        intPot!(a, c, -5.0)
        intPot!(b, c, -6.0)

        init!(l, [a, b, c])

        # All three should share the same matrix
        @test a.intMat === b.intMat
        @test b.intMat === c.intMat

        # Matrix should be 3x3
        @test size(a.intMat) == (3, 3)

        # Check all pairwise interactions (symmetric)
        @test a.intMat[a.intIdx, b.intIdx] == -4.0
        @test a.intMat[b.intIdx, a.intIdx] == -4.0
        @test a.intMat[a.intIdx, c.intIdx] == -5.0
        @test a.intMat[c.intIdx, a.intIdx] == -5.0
        @test a.intMat[b.intIdx, c.intIdx] == -6.0
        @test a.intMat[c.intIdx, b.intIdx] == -6.0

        # Diagonal entries should equal self intPot
        @test a.intMat[a.intIdx, a.intIdx] == a.intPot
        @test a.intMat[b.intIdx, b.intIdx] == b.intPot
        @test a.intMat[c.intIdx, c.intIdx] == c.intPot

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB
    end

# MARK: -intPot update
    @testset "intPot! (update existing)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set initial interaction
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        @test a.intMat[a.intIdx, b.intIdx] == -3.0

        # Update interaction - should overwrite
        intPot!(a, b, -5.0)
        @test a.intMat[a.intIdx, b.intIdx] == -5.0
        @test a.intMat[b.intIdx, a.intIdx] == -5.0

        # Update in reverse order
        intPot!(b, a, -7.0)
        @test a.intMat[a.intIdx, b.intIdx] == -7.0
        @test a.intMat[b.intIdx, a.intIdx] == -7.0
    end

# MARK: -intPot self after pair
    @testset "intPot! (self after pair interaction)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set pairwise interaction first
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        # Now update self interaction - should update diagonal
        intPot!(a, -10.0)
        @test a.intPot == -10.0
        @test a.intMat[a.intIdx, a.intIdx] == -10.0

        # Cross interaction should be unchanged
        @test a.intMat[a.intIdx, b.intIdx] == -3.0
    end

# MARK: -getIntPot
    @testset "getIntPot" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 3
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulk(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Before setting interactions, getIntPot should return nothing
        @test getIntPot(a, b) === nothing
        
        # Set interactions
        intPot!(a, b, -4.0)
        intPot!(a, c, -5.0)
        intPot!(b, c, -6.0)
        
        init!(l, [a, b, c])

        # Now getIntPot should return correct values
        @test getIntPot(a, b) == -4.0
        @test getIntPot(b, a) == -4.0  # symmetric
        @test getIntPot(a, c) == -5.0
        @test getIntPot(c, a) == -5.0
        @test getIntPot(b, c) == -6.0
        @test getIntPot(c, b) == -6.0

        # Self interactions
        @test getIntPot(a, a) == a.intPot
        @test getIntPot(b, b) == b.intPot
        @test getIntPot(c, c) == c.intPot
    end

# MARK: -intPot chain merge
    @testset "intPot! (chain merge)" begin
        # Test merging two separate interaction networks
        l = Lattice(;
            boxlength = 8,
            nspecies = 4
        )

        a = FastBulk(;
            length = 5,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 5,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulk(;
            length = 5,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        d = FastBulk(;
            length = 5,
            intPot = -4.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Create two separate networks: (a,b) and (c,d)
        intPot!(a, b, -10.0)
        intPot!(c, d, -20.0)
        
        init!(l, [a, b, c, d])

        # Verify separate networks
        @test a.intMat === b.intMat
        @test c.intMat === d.intMat
        @test a.intMat !== c.intMat

        # Now connect the networks
        intPot!(b, c, -30.0)

        # All should now share the same matrix
        @test a.intMat === b.intMat
        @test b.intMat === c.intMat
        @test c.intMat === d.intMat

        # Matrix should be 4x4
        @test size(a.intMat) == (4, 4)

        # All previous interactions should be preserved
        @test getIntPot(a, b) == -10.0
        @test getIntPot(c, d) == -20.0
        @test getIntPot(b, c) == -30.0

        # Diagonal entries preserved
        @test a.intMat[a.intIdx, a.intIdx] == -1.0
        @test a.intMat[b.intIdx, b.intIdx] == -2.0
        @test a.intMat[c.intIdx, c.intIdx] == -3.0
        @test a.intMat[d.intIdx, d.intIdx] == -4.0

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)
        intIdxsC = map(x -> x[2], c.intSpecies)
        intIdxsD = map(x -> x[2], d.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB
        @test c.intIdx in intIdxsC
        @test d.intIdx in intIdxsD
    end

# MARK: -intPot with Bulk
    @testset "intPot! (mixed Bulk and FastBulk)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = Bulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set cross-interaction between Bulk and FastBulk
        intPot!(a, b, -5.0)
        
        init!(l, [a, b])

        # Both should have interaction structures
        @test a.intMat !== nothing
        @test b.intMat !== nothing
        @test a.intMat === b.intMat

        # Check values
        @test getIntPot(a, b) == -5.0
        @test getIntPot(b, a) == -5.0
    end

# MARK: -intPot partner
    @testset "intPot! (partner update)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulk(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set partners before init
        a.partner = b
        b.partner = a
        
        # Set interaction potential between a and b
        intPot!(a, b, -5.0)

        init!(l, [a, b])

        # Check interaction is set in intMat
        @test a.intMat[a.intIdx, b.intIdx] == -5.0
        @test b.intMat[b.intIdx, a.intIdx] == -5.0
    end

# MARK: -intPot same species
    @testset "intPot! (same species twice)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 1
        )

        a = FastBulk(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )
        
        # Set interaction with itself
        intPot!(a, a, -5.0)

        init!(l, [a])

        # Should update diagonal
        @test a.intMat !== nothing
        @test size(a.intMat) == (1, 1)
        @test a.intMat[1, 1] == -5.0
    end


# MARK: -Subspecies Tests
    @testset "Subspecies initialization" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        p = Poly(;
            length = 16,
            maxlength = 20,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        
        # Initialize with 2 subspecies, uniform pattern (all beads = subspecies 1)
        initSubspecies!(p, Int32(2), pattern=:uniform)
        
        @test p.nSubspecies == 2
        @test length(p.subspecies) == 2
        @test length(p.beadTypes) == p.maxlength
        @test all(p.beadTypes .== 1)
        
        # Each subspecies should inherit polymer's chemPot and intPot
        @test p.subspecies[1].chemPot == p.chemPot
        @test p.subspecies[1].intPot == p.intPot
        @test p.subspecies[2].chemPot == p.chemPot
        @test p.subspecies[2].intPot == p.intPot
    end

    @testset "Subspecies patterns" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        # Test alternating pattern
        p1 = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )
        init!(l, [p1])
        initSubspecies!(p1, Int32(2), pattern=:alternating)
        
        @test p1.beadTypes[1] == 1
        @test p1.beadTypes[2] == 2
        @test p1.beadTypes[3] == 1
        @test p1.beadTypes[4] == 2

        # Test block pattern
        p2 = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )
        init!(l, [p2])
        initSubspecies!(p2, Int32(2), pattern=:block)
        
        # First half should be subspecies 1, second half subspecies 2
        @test p2.beadTypes[1] == 1
        @test p2.beadTypes[5] == 1
        @test p2.beadTypes[6] == 2
        @test p2.beadTypes[10] == 2
    end

    @testset "Subspecies intPot!" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        p = Poly(;
            length = 16,
            maxlength = 20,
            intPot = 1.0,
            chemPot = 1.0
        )

        b = Bulk(;
            length = 100,
            intPot = 1.0,
            chemPot = 1.0
        )
        intPot!(b, 1.0)

        init!(l, [p, b])
        
        # Initialize subspecies
        initSubspecies!(p, Int32(2), pattern=:alternating)
        
        # Set different intPot for each subspecies
        intPot!(p, 1, 2.0)  # subspecies 1 self-interaction
        intPot!(p, 2, 3.0)  # subspecies 2 self-interaction
        
        @test p.subspecies[1].intPot == 2.0
        @test p.subspecies[2].intPot == 3.0
        
        # Set cross-subspecies interaction (uses shared intMat on subspecies)
        intPot!(p, 1, 2, 1.5)  # interaction between subspecies 1 and 2
        
        @test getIntPot(p, 1, 2) == 1.5
        @test getIntPot(p, 2, 1) == 1.5  # symmetric
        
        # Verify subspecies share the same intMat
        @test p.subspecies[1].intMat === p.subspecies[2].intMat
        
        # Set interaction between subspecies and bulk species
        intPot!(p, 1, b, -1.0)  # subspecies 1 with bulk
        intPot!(p, 2, b, -2.0)  # subspecies 2 with bulk
        
        @test getIntPot(p, 1, b) == -1.0
        @test getIntPot(p, 2, b) == -2.0
        
        # Test bulkIntPot structure
        @test p.bulkIntPot !== nothing
        @test p.bulkSpecies !== nothing
        @test b in p.bulkSpecies
    end

    @testset "Subspecies getSubspeciesIdx" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        p = Poly(;
            length = 8,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        initSubspecies!(p, Int32(2), pattern=:alternating)
        
        # Beads 1,3,5,7 should be subspecies 1
        # Beads 2,4,6,8 should be subspecies 2
        @test getSubspeciesIdx(p, Int32(1)) == 1  # head
        @test getSubspeciesIdx(p, Int32(2)) == 2  # second bead
        @test getSubspeciesIdx(p, Int32(3)) == 1  # third bead
    end

    @testset "Subspecies ham! backward compatibility" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        # Test without subspecies (backward compatible)
        p = Poly(;
            length = 16,
            maxlength = 20,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        ham!(p)
        legacy_ham = p.ham
        
        # Now initialize subspecies with uniform (all same) and same intPot
        initSubspecies!(p, Int32(1), pattern=:uniform)
        ham!(p)
        subspecies_ham = p.ham
        
        # Should give same result
        @test legacy_ham ≈ subspecies_ham
    end

    @testset "Subspecies ham! different potentials" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        p = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        
        # Initialize with 2 subspecies, block pattern
        initSubspecies!(p, Int32(2), pattern=:block)
        
        # Set different potentials for each subspecies
        intPot!(p, 1, 2.0)  # higher interaction for first half
        intPot!(p, 2, 0.5)  # lower interaction for second half
        
        # Calculate Hamiltonian
        ham!(p)
        
        # The result should be different from uniform potential case
        p2 = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )
        init!(l, [p2])
        ham!(p2)
        
        # Different potentials should give different energy
        @test p.ham != p2.ham
    end

    @testset "Subspecies polymer-bulk interaction" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        p = Poly(;
            length = 6,
            maxlength = 8,
            intPot = 1.0,
            chemPot = 1.0
        )

        b = Bulk(;
            length = 50,
            intPot = 1.0,
            chemPot = 1.0
        )
        intPot!(b, 1.0)

        init!(l, [p, b])
        
        # Initialize subspecies
        initSubspecies!(p, Int32(2), pattern=:alternating)
        
        # Set different bulk interactions for each subspecies
        intPot!(p, 1, b, -1.0)  # subspecies 1 attracted to bulk
        intPot!(p, 2, b, 1.0)   # subspecies 2 repelled by bulk
        
        # Calculate Hamiltonian - should account for bead-specific interactions
        ham!(p)
        
        @test p.ham !== nothing
        @test p.bulkIntPot !== nothing
        
        # Verify interactions are stored correctly
        sub1_intIdx = p.subspecies[1].intIdx
        sub2_intIdx = p.subspecies[2].intIdx
        @test p.bulkIntPot[(sub1_intIdx, b.id)] == -1.0
        @test p.bulkIntPot[(sub2_intIdx, b.id)] == 1.0
    end

    @testset "Subspecies intPot! all subspecies convenience" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        p = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )

        b = Bulk(;
            length = 50,
            intPot = 1.0,
            chemPot = 1.0
        )
        intPot!(b, 1.0)

        init!(l, [p, b])
        initSubspecies!(p, Int32(3), pattern=:block)
        
        # Set interaction for ALL subspecies at once
        intPot!(p, b, -2.0)
        
        # All subspecies should have same interaction
        @test getIntPot(p, 1, b) == -2.0
        @test getIntPot(p, 2, b) == -2.0
        @test getIntPot(p, 3, b) == -2.0
        
        # Partner should be set
        @test p.partner === b
    end

    @testset "Subspecies custom pattern" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        p = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        
        # Create a custom pattern: first 3 beads subspecies 1, next 4 subspecies 2, last 3 subspecies 1
        customPattern = Int32[1, 1, 1, 2, 2, 2, 2, 1, 1, 1]
        initSubspecies!(p, Int32(2), pattern=:custom, customPattern=customPattern)
        
        @test p.beadTypes == customPattern
        @test p.subspecies[1].intMat === p.subspecies[2].intMat  # shared matrix
    end

    @testset "Subspecies intMat synchronization" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        p = Poly(;
            length = 10,
            maxlength = 10,
            intPot = 1.0,
            chemPot = 1.0
        )

        init!(l, [p])
        initSubspecies!(p, Int32(3), pattern=:block)
        
        # All subspecies should share the same intMat and intSpecies
        @test p.subspecies[1].intMat === p.subspecies[2].intMat
        @test p.subspecies[2].intMat === p.subspecies[3].intMat
        @test p.subspecies[1].intSpecies === p.subspecies[2].intSpecies
        
        # Setting interaction should update the shared matrix
        intPot!(p, 1, 2, 0.5)
        intPot!(p, 2, 3, 0.7)
        
        # Access via any subspecies should give same result
        @test p.subspecies[1].intMat[p.subspecies[1].intIdx, p.subspecies[2].intIdx] == 0.5
        @test p.subspecies[2].intMat[p.subspecies[2].intIdx, p.subspecies[3].intIdx] == 0.7
    end
end
