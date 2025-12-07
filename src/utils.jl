############# Generic utilities #############
# MARK: -Generic

export coords_equal
@inline function coords_equal(a, b)
    return a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
end


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


# volume of a bounding box
export vol
function vol(bbox::Vector{Vector{Int32}})

    @assert length(bbox) == 3
    @assert length(bbox[1]) == 2 && length(bbox[2]) == 2 && length(bbox[3]) == 2

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    
    return abs(xh - xl + 1) * abs(yh - yl + 1) * abs(zh - zl + 1)
end


# returns true if at least one element of x is adjacent to the bbox
export onEdge
function onEdge(s::AbstractSpecies, x::Set{Int32}, bbox::Vector{Vector{Int32}})

    @assert length(bbox) == 3
    @assert length(bbox[1]) == 2 && length(bbox[2]) == 2 && length(bbox[3]) == 2

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    
    for p in x
        if (s.pos[p][1] == xl || s.pos[p][1] == xh
            || s.pos[p][2] == yl || s.pos[p][2] == yh
            || s.pos[p][3] == zl || s.pos[p][3] == zh)

            return true
        end
    end

    return false
end


# Number of same-species nearest neighbors at a given position
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

    # clear temporary containers
    clamp!(s.tempVec, 0, 0)

    return n
end


# Number of nearest neighbors between species s and p at a given position
export nnp
function nnp(s::AbstractSpecies, p::AbstractSpecies, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0

    for i in range(1, 3)
        pos[i] += 1
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += (p.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]] > 0)

        pos[i] -= 2
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += (p.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]] > 0)

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(s.tempVec, 0, 0)

    return n
end


# Set of nearest neighbors of a given position
# WARNING: will break if called on a Poly TODO: fix
export nnIdxs!
function nnIdxs!(b::AbstractSpecies, 
                 pos::Union{SubArray{Int32}, Vector{Int32}})
    k = 1

    # clear working containers
    clamp!(b.tempVec, 0, 0)
    clamp!(b.nnIdxs, 0, 0)
    empty!(b.freeNNs)

    for i in range(1, 3)
        pos[i] += 1
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        
        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end
        k += 1

        pos[i] -= 2
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        
        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end
        k += 1

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(b.tempVec, 0, 0)

    return b
end


############# Bulk utilities #############
# MARK: -Fast/Bulk

# splits droplets that sit across a lattice boundary
# WARNING: assumes droplet max diameter < boxlength / 2, 
export splitDroplet!
function splitDroplet!(s::AbstractBulk, droplet::Set{Int32})

    droplet_vecs = []
    droplet_vec = collect(droplet)
    push!(droplet_vecs, [droplet_vec[1]])
    
    for (i, p) in enumerate(droplet_vec)
        if i > 1
            belonged = false
            for (j, d) in enumerate(droplet_vecs)
                if (abs(s.pos[p][1] - s.pos[d[1]][1]) < s.boxlength / 2 &&
                    abs(s.pos[p][2] - s.pos[d[1]][2]) < s.boxlength / 2 &&
                    abs(s.pos[p][3] - s.pos[d[1]][3]) < s.boxlength / 2)

                    push!(droplet_vecs[j], p)
                    belonged = true
                end
            end

            if !belonged
                push!(droplet_vecs, [p])
            end
        end
    end

    return [Set(d) for d in droplet_vecs]
end


# finds the bounding box of a droplet, adds 1 along each axis for the empty
# space around the droplet to be connected
export dropletBbox
function dropletBbox(s::AbstractBulk, droplet::Set{Int32})

    # start with position of an arbitrary particle in the droplet
    key0 = collect(droplet)[1]
    xb = [s.pos[key0][1], s.pos[key0][1]]
    yb = [s.pos[key0][2], s.pos[key0][2]]
    zb = [s.pos[key0][3], s.pos[key0][3]]

    for (i, key) in enumerate(droplet)
        # compute lower bounds
        xb[1] = s.pos[key][1] < xb[1] ? s.pos[key][1] : xb[1]
        yb[1] = s.pos[key][2] < yb[1] ? s.pos[key][2] : yb[1]
        zb[1] = s.pos[key][3] < zb[1] ? s.pos[key][3] : zb[1]

        # compute upper bounds
        xb[2] = s.pos[key][1] > xb[2] ? s.pos[key][1] : xb[2]
        yb[2] = s.pos[key][2] > yb[2] ? s.pos[key][2] : yb[2]
        zb[2] = s.pos[key][3] > zb[2] ? s.pos[key][3] : zb[2]
    end

    # add 1 and clip to lattice bounds
    xb[1] = max(xb[1] - 1, 1)
    yb[1] = max(yb[1] - 1, 1)
    zb[1] = max(zb[1] - 1, 1)

    xb[2] = min(xb[2] + 1, s.boxlength)
    yb[2] = min(yb[2] + 1, s.boxlength)
    zb[2] = min(zb[2] + 1, s.boxlength)

    return [xb, yb, zb]
end


# fills up a box around a droplet in the co-lattice
export fillBbox!
function fillBbox!(s::AbstractBulk, droplet::Set{Int32}, bbox::Vector{Vector{Int32}})

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]

    # ghost species
    gs = s.ghost
    key = 0

    # loop across the bbox, if site not part of the droplet, add it to the ghost lattice
    for x in range(xl, xh)
        for y in range(yl, yh)
            for z in range(zl, zh)
                if !(s.latS[x, y, z] in droplet)
                    key += 1
                    gs.latS[x, y, z] = key
                    gs.pos[key] = zeros(5)
                    gs.pos[key] = [x, y, z, key - 1, key + 1]
                end
            end
        end
    end

    gs.pos[1][4] = 0
    gs.pos[key][5] = 0

    return s
end


# clears a box around a droplet in the ghost lattice
export clearBbox!
function clearBbox!(s::AbstractBulk, bbox::Vector{Vector{Int32}})

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]

    # ghost species
    gs = s.ghost

    # loop across x,y,z in the bbox, if empty, add it to the co-lattice
    for x in range(xl, xh)
        for y in range(yl, yh)
            for z in range(zl, zh)
                gs.latS[x, y, z]  = 0
            end
        end
    end

    gs.pos = OrderedDict{Int32, Vector{Int32}}()

    return s
end