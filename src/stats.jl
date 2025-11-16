# returns the center of mass of all (or a subset of) particles of a species
export centerMass
function centerMass(s::AbstractSpecies, 
                    subset::Union{Nothing, Vector{Int32}} = nothing)
    cm = zeros(3)
    subset = subset === nothing ? eachindex(s.pos) : subset

    # compute average of all positions in pos
    for i in subset
        cm = cm .+ view(s.pos[i], 1:3)
    end

    return cm ./ s.length
end

# radius of gyration, mean squared distance between each bead and
# the center of mass
export radius
function radius(s::AbstractSpecies)
    cm = centerMass(s)
    rg2 = 0

    for i in eachindex(s.pos)
        rg2 += ((cm[1] .- s.pos[i][1])^2
             + (cm[2] .- s.pos[i][2])^2
             + (cm[3] .- s.pos[i][3])^2)
    end

    return sqrt(rg2 / s.length)
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
export density
function density(s::AbstractSpecies, bbox::Vector{Vector{Int32}})
    count = 0
    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    vol = abs(xh - xl) * abs(yh - yl) * abs(zh - zl)

    for x in range(xl, xh)
        for y in range(yl, yh)
            for z in range(zl, zh)
                count += s.latS[x,y,z]
            end
        end
    end

    return count / vol
end


# Depth-first search: takes a bead index bi and marks all beads
# traversed starting from bi (adds their indices to the marked list)
export dfs
function dfs(s::AbstractFastBulk, bi::Int32, marked::Set{Int32})
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

# Connected components in the lattice
export cc
function cc(s::AbstractFastBulk)

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
