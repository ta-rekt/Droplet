export extractCoords!, centerMass, radius, radiusLine, radiusSphere, boxDensity, 
       axialDensity , radialDensity, extractCoords!, dfs, cc, dropletVolume, 
       dropletDensities, diluteDensity, nspeciesBound, get_pattern_centers,
       patterns_in_contact, compute_contact_matrix!

############### Generic #############
# MARK: -Generic

# extracts the coordinates of a polymer's beads
function extractCoords!(b::AbstractSpecies, coords::Matrix)

    for (i, k) in enumerate(collect(keys(b.pos)))
        coords[i, 1:3] = b.pos[k][1:3]         # next position
        coords[i, 4] = b.latS[b.pos[k][1], b.pos[k][2], b.pos[k][3]]   # bead index
    end

    return coords
end

# returns the center of mass of all (or a subset of) particles of a species
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
function radiusLine(length)
    n = div(length, 2)
    return sqrt((n - 1) / 2)
end


# gyration radius of a sphere
function radiusSphere(length)
    r = (3*length/(4*pi))^(1/3)
    return sqrt(3/5)*r
end


# density of particles of species s inside a bounding box
# bbox[1] = [xl, xh], bbox[2] = [yl, yh], bbox[3] = [zl, zh]
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

# axial density along a parallelepiped with a given thickness
function axialDensity(s::AbstractSpecies, axis::String, 
                      pos::Vector{Int32}, thick::Int32)

end

# radial density along a disk with given thickness 
function radialDensity(s::AbstractSpecies, axis::String, 
                      pos::Vector{Int32}, thick::Int32)

end

########### Poly stats ###########
# MARK: -Poly

# extracts the coordinates of a polymer's beads
function extractCoords!(p::Polymer, coords::Matrix)
    key = p.head           # initial bead
    k = 1
    
    for (i, k) in enumerate(collect(keys(p.pos)))
        if key != 0
            coords[i, 1:3] = p.pos[key][1:3] # next position
            coords[i, 4] = p.pos[key][6]  # monomer type
            key = p.pos[key][5]    # next bead
        end
    end

    return coords
end

# helper function to calculate pattern center positions
# patterns: vector of tuples [(start1, end1), (start2, end2), ...]
function get_pattern_centers(poly::Polymer, patterns::Vector{Tuple{Int,Int}})
    if length(patterns) != 2
        error("get_pattern_centers requires exactly 2 patterns, got $(length(patterns))")
    end
    
    # check if polymer is long enough for patterns
    max_bead = maximum([p[2] for p in patterns])
    if poly.length < max_bead
        return nothing, nothing, Inf
    end
    
    centers = [zeros(3) for _ in patterns]
    pattern_lengths = [p[2] - p[1] + 1 for p in patterns]
    
    # iterate through polymer linked list
    key = poly.head
    bead_num = 1
    
    for i in 1:poly.length
        pos = poly.pos[key][1:3]
        
        # check which pattern(s) this bead belongs to
        for (pidx, (pstart, pend)) in enumerate(patterns)
            if bead_num >= pstart && bead_num <= pend
                centers[pidx] .+= pos
            end
        end
        
        key = poly.pos[key][5]  # next bead
        bead_num += 1
    end
    
    # normalize by pattern lengths
    for (i, plen) in enumerate(pattern_lengths)
        centers[i] ./= plen
    end
    
    # euclidean distance between the two centers
    dist = sqrt(sum((centers[1] .- centers[2]).^2))
    
    return centers[1], centers[2], dist
end

# helper function to check if patterns are in contact
# patterns: vector of tuples [(start1, end1), (start2, end2), ...]
function patterns_in_contact(poly::Polymer, patterns::Vector{Tuple{Int,Int}})
    if length(patterns) != 2
        error("patterns_in_contact requires exactly 2 patterns, got $(length(patterns))")
    end
    
    # check if polymer is long enough
    max_bead = maximum([p[2] for p in patterns])
    if poly.length < max_bead
        return false
    end
    
    # collect positions for each pattern
    pattern_positions = [[] for _ in patterns]
    
    key = poly.head
    bead_num = 1
    
    for i in 1:poly.length
        pos = poly.pos[key][1:3]
        
        for (pidx, (pstart, pend)) in enumerate(patterns)
            if bead_num >= pstart && bead_num <= pend
                push!(pattern_positions[pidx], copy(pos))
            end
        end
        
        key = poly.pos[key][5]
        bead_num += 1
    end
    
    # check all pairs of beads between pattern 1 and pattern 2
    for pos1 in pattern_positions[1]
        for pos2 in pattern_positions[2]
            # check if adjacent (within 1 lattice site in each dimension, periodic BC)
            dx = abs(pos1[1] - pos2[1])
            dy = abs(pos1[2] - pos2[2])
            dz = abs(pos1[3] - pos2[3])
            
            # account for periodic boundaries
            dx = min(dx, poly.boxlength - dx)
            dy = min(dy, poly.boxlength - dy)
            dz = min(dz, poly.boxlength - dz)
            
            if dx <= 1 && dy <= 1 && dz <= 1 && (dx + dy + dz) > 0
                return true
            end
        end
    end
    
    return false
end

########### BulkMonomer stats ###########
# MARK: -Fast/BulkMonomer

# depth-first search: takes a bead index bi and marks all beads
# traversed starting from bi (adds their indices to the marked list)
# WARNING: will break if called on Poly TODO: fix
function dfs(s::AbstractSpecies, bi::Int32)
    marked = Set{Int32}()
    stack_ = OrderedSet{Int32}()

    push!(marked, bi)
    push!(stack_, bi)

    while !isempty(stack_)
        bi = pop!(stack_)
        nnIdxs!(s, view(s.pos[bi], range(1,3)))

        for bj in s.nnIdxs
            if bj != 0 && !(bj in marked)
                push!(stack_, bj)
                push!(marked, bj)
            end
        end
    end

    return marked
end


# connected components
# WARNING: will break if called on Poly TODO: fix
function cc(s::AbstractSpecies)
    ccs = Vector{Set{Int32}}()
    unvisited = Set{Int32}(keys(s.pos))  # all keys

    while !isempty(unvisited)
        bi = pop!(unvisited)
        component = dfs(s, bi)

        for bj in component
            if bi != bj
                pop!(unvisited, bj)
            end
        end

        push!(ccs, component)
    end

    return ccs
end


# returns the exact volume of a droplet. Length of CCs in the empty space
# inside the bounding box, subtracts that from the volume of the box
function dropletVolume(s::AbstractBulkMonomer, droplet::Set{Int32})

    # bounding box around a droplet
    bbox = dropletBbox(s, droplet)

    # fill up the bounding box in the ghost lattice
    fillBbox!(s, droplet, bbox)

    # connected components adjacent to bbox walls around a droplet
    ccs = cc(s.ghost)
    length_ccs = sum(map(x -> onEdge(s.ghost, x, bbox) ? length(x) : 0, ccs))

    # cleans up sites in the ghost lattice
    clearBbox!(s, bbox)

    return vol(bbox) - length_ccs
end


# exact density of each droplet in a list 
function dropletDensities(s::AbstractBulkMonomer, droplets::Vector{Set{Int32}})
    densities = Vector{Float64}(undef, 0)

    for d in droplets
        push!(densities, length(d) / dropletVolume(s, d))
    end

    return densities
end


# (number of particles - total droplet size) / (volume - total droplet volume)
function diluteDensity(s::AbstractBulkMonomer, droplets::Vector{Set{Int32}})
    
    n = s.length
    ndrop = sum(map(x -> length(x), droplets))
    voldrop = 0

    for d in droplets
        voldrop += dropletVolume(s, d)
    end

    return (n - ndrop) / (s.boxlength^3 - voldrop)
end

# number of particles of a species bound to any of the binding sites
function nspeciesBound(sp::AbstractBulkMonomer, 
                       bsites::Vector{T}) where T <: AbstractSpecies
    n = 0
    for site in bsites
        for (_, val) in site.pos
            n += (sp.latS[val[1], val[2], val[3]] > 0)
        end
    end
    return n
end

# computes euclidean distance between two positions accounting for periodic boundaries
function distance_periodic(pos1::Vector, pos2::Vector, boxlength::Int32)
    dx = abs(pos1[1] - pos2[1])
    dy = abs(pos1[2] - pos2[2])
    dz = abs(pos1[3] - pos2[3])
    
    # account for periodic boundaries
    dx = min(dx, boxlength - dx)
    dy = min(dy, boxlength - dy)
    dz = min(dz, boxlength - dz)
    
    return sqrt(dx^2 + dy^2 + dz^2)
end

# accumulates contacts into a matrix by traversing polymer linked list
# contact_threshold: distance cutoff for considering beads in contact (default 2.0)
function compute_contact_matrix!(poly::Polymer, contact_matrix::Matrix{Int}, 
                                contact_threshold::Float64 = 2.0)
    # traverse polymer and compute all pairwise distances
    key_i = poly.head
    bead_i = 1
    
    for i in 1:poly.length
        pos_i = poly.pos[key_i][1:3]
        
        # check contacts with all subsequent beads
        key_j = poly.pos[key_i][5]  # next bead after i
        bead_j = bead_i + 1
        
        for j in bead_j:poly.length
            pos_j = poly.pos[key_j][1:3]
            
            dist = distance_periodic(pos_i, pos_j, Int32(poly.boxlength))
            
            if dist <= contact_threshold
                contact_matrix[bead_i, bead_j] += 1
                contact_matrix[bead_j, bead_i] += 1
            end
            
            key_j = poly.pos[key_j][5]
        end
        
        key_i = poly.pos[key_i][5]
        bead_i += 1
    end

    return contact_matrix
end