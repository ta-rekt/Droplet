# returns the center of mass of a species
export centerMass
function centerMass(s::AbstractSpecies)
    cm = zeros(3)

    # compute average of all positions in pos
    for i in eachindex(s.pos)
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


# Depth-first search: takes a bead index bi and marks all beads
# traversed starting from bi (adds their indices to the marked list)
export dfs
function dfs(s::AbstractFastBulk, bi::Int32, marked::Vector{Int32})
    if bi not in marked
        push!(marked, bi)
    end

    nnIdxs!(s, s.pos[bi])

    for bj in s.nnIdxs
        if bj != 0
            if bj not in marked
                marked = dfs(s, bj, marked)
            end
        end
    end

    return marked
end

# Connected components in the lattice
export cc
function cc(s::AbstractFastBulk)

    ccs = []
    marked = Set(range(1, s.length))  # all keys
    marked_bi = []

    while !isempty(marked)
        bi = marked[1]
        marked_bi = dfs(s, bi, marked_bi)

        for bj in marked_bi
            pop!(marked, bj)
        end

        push!(ccs, marked_bi)
    end

    return ccs
end
