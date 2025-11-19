# returns the center of mass of all (or a subset of) particles of a species
export centerMass
function centerMass(s::AbstractSpecies, 
                    subset::Union{Nothing, Set{Int32}} = nothing)
    cm = zeros(3)
    subset = subset === nothing ? eachindex(s.pos) : subset

    # compute average of all positions in pos
    for i in subset
        cm = cm .+ view(s.pos[i], 1:3)
    end

    return cm ./ length(subset)
end


# radius of gyration, mean squared distance between each bead and
# the center of mass
export radius
function radius(s::AbstractSpecies)
    cm = centerMass(s)
    rg = 0

    for i in eachindex(s.pos)
        rg += ((cm[1] .- s.pos[i][1])^2
             + (cm[2] .- s.pos[i][2])^2
             + (cm[3] .- s.pos[i][3])^2)
    end

    return sqrt(rg / s.length)
end


# gyration radius of a line
export radiusLine
function radiusLine(length)
    n = div(length, 2)
    return sqrt((n - 1) / 2)
end


# gyration radius of a sphere
export radiusSphere
function radiusSphere(length)
    r = (3*length/(4*pi))^(1/3)
    return sqrt(3/5)*r
end


# density of particles of species s inside a bounding box
# bbox[1] = [xl, xh], bbox[2] = [yl, yh], bbox[3] = [zl, zh]
export boxDensity
function boxDensity(s::AbstractSpecies, bbox::Vector{Vector{Int32}})
    count = 0
    pos = zeros(Int32, 3)
    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    vol = abs(xh - xl + 1) * abs(yh - yl + 1) * abs(zh - zl + 1)

    for x in range(xl, xh)
        for y in range(yl, yh)
            for z in range(zl, zh)
                boxmod!([x,y,z], pos, Int32(1), s.boxlength)
                count += (s.latS[pos[1], pos[2], pos[3]] > 0)
            end
        end
    end

    return count / vol
end


export dropletDensity
function dropletDensity(s::AbstractSpecies, ccs::Set{Set{Int32}})
    
    dens = Vector{Float64}(undef, 0)

    for (i, cc) in enumerate(ccs)
        # center of mass of each droplet
        cm = centerMass(s, cc)

        # compute bbox given droplet size
        b = max(1, 0.25 * (length(cc))^(1/3))
        bbox = [[floor(Int32, cm[1] - 0.5*b), floor(Int32, cm[1] + 0.5*b)],
                [floor(Int32, cm[2] - 0.5*b), floor(Int32, cm[2] + 0.5*b)],
                [floor(Int32, cm[3] - 0.5*b), floor(Int32, cm[3] + 0.5*b)]]

        # droplet density
        push!(dens, boxDensity(s, bbox))
    end

    return dens
end


# density of dilute medium 
export diluteDensity
function diluteDensity(s::AbstractBulk)

    # 8 bounding boxes, one at each corner of the lattice
    l = s.boxlength
    b = div(l, 8)
    bboxes = zeros(Int32, 8, 3, 2)
    bboxes[1] = [[1, b],[1, b],[1, b]]
    bboxes[2] = [[1, b],[1, b],[l - b + 1, l]]
    bboxes[3] = [[1, b],[l - b + 1, l],[1, b]]
    bboxes[4] = [[1, b],[l - b + 1, l],[l - b + 1, l]]
    bboxes[5] = [[l - b + 1, l],[1, b],[1, b]]
    bboxes[6] = [[l - b + 1, l],[1, b],[l - b + 1, l]]
    bboxes[7] = [[l - b + 1, l],[l - b + 1, l],[1, b]]
    bboxes[8] = [[l - b + 1, l],[l - b + 1, l],[l - b + 1, l]]

    # return the densities inside each of the 8 bounding boxes
    dens = Dict{Int32, Float64}()

    for i in range(1, 8)
        push!(dens, (i => boxDensity(s, bboxes[i])))
    end

    # select 4 boxes with the lowest density
    for i in range(1, 4)
        pop!(dens, argmin(dens))
    end

    # return average density of those 4 boxes
    return sum(values(dens)) / length(dens)
end


# depth-first search: takes a bead index bi and marks all beads
# traversed starting from bi (adds their indices to the marked list)
export dfs
function dfs(s::AbstractBulk, bi::Int32, marked::Set{Int32})
    if !(bi in marked)
        push!(marked, bi)
    end

    nnIdxs!(s, s.pos[bi][range(1,3)])
    nnIdxs = deepcopy(s.nnIdxs)

    for bj in nnIdxs
        if bj != 0
            if !(bj in marked)
                marked = dfs(s, bj, marked)
            end
        end
    end

    return marked
end


# connected components
export cc
function cc(s::AbstractBulk)

    ccs = Set{Set{Int32}}()
    marked = Set{Int32}(keys(s.pos))  # all keys
    marked_bi = Set{Int32}()

    while !isempty(marked)
        bi = pop!(marked)
        marked_bi = dfs(s, bi, marked_bi)

        for bj in marked_bi
            if bi != bj
                pop!(marked, bj)
            end
        end

        push!(ccs, marked_bi)
        marked_bi = Set{Int32}()
    end

    return ccs
end
