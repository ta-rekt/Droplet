############# Generic utilities #############

# ensures positions are periodic within the lattice's bounding box
export boxmod!
function boxmod!(v::Union{SubArray{Int32}, Vector{Int32}}, 
                 w::Union{SubArray{Int32}, Vector{Int32}}, 
                 lb::Int32, ub::Int32)

    for i in eachindex(v)
        if v[i] > ub
            w[i] = v[i] - (ub - lb + 1)
        elseif v[i] < lb
            w[i] = v[i] + (ub - lb + 1)
        else
            w[i] = v[i]
        end
    end
    
    return w
end

# Number of nearest neighbors at a given position
export nn
function nn(s::AbstractSpecies, pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0

    for i in range(1, 3)
        pos[i] += 1
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += (s.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]] > 0)

        pos[i] -= 2
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += (s.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]] > 0)

        pos[i] += 1
    end

    return n
end

export coords_equal
@inline function coords_equal(a, b)
    return a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
end


############# FastBulk utilities #############

# Set of nearest neighbors of a given position
export nnIdxs
function nnIdxs(b::AbstractFastBulk, pos::Union{SubArray{Int32}, Vector{Int32}})
    k = 1

    for i in range(1, 3)
        pos[i] += 1
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        k += 1

        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end

        pos[i] -= 2
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        k += 1

        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end

        pos[i] += 1
    end

    return b
end