
export extractCoords!, centerMass, radius, radiusLine, radiusSphere, boxDensity, 
       axialDensity , radialDensity, extractCoords!, dfs, cc, dropletVolume, 
       dropletDensities, diluteDensity, nspeciesBound

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
function extractCoords!(p::AbstractPolymer, coords::Matrix)
    key = p.head           # initial bead
    k = 1
    
    for (i, k) in enumerate(collect(keys(p.pos)))
        if key != 0
            coords[i, 1:3] = p.pos[key][1:3] # next position
            coords[i, 4] = i
            key = p.pos[key][5]    # next bead
        end
    end

    return coords
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