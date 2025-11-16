# saves a lattice's state into a binary file
export save!
function save!(l::Lattice, filename::String)
    filename = string(filename, ".bin")
    Serialization.serialize(filename, l.lat)
end

# loads a lattice's state from a binary file
export load!
function load!(l::Lattice, filename::String)
    filename = string(filename, ".bin")
    l.lat = Serialization.deserialize(filename)

    return l
end

# extracts the coordinates of a polymer's beads
export extractCoords!
function extractCoords!(p::AbstractPoly, coords::Matrix)
    key = p.head           # initial bead
    k = 1
    
    for k in 1:p.length
        if key != 0
            coords[k, :] = p.pos[key][1:3] # next position
            key = p.pos[key][5]    # next bead
        end
    end

    return coords
end

# extracts the coordinates of a polymer's beads
export extractCoords!
function extractCoords!(b::AbstractBulk, coords::Matrix)

    for k in 1:b.length
        coords[k, :] = b.pos[k][1:3] # next position
    end

    return coords
end

# exports a list of species to JSON for python compatibility
# slist is a list of species s, export the bead coordinates s.pos in a 
# (s.length x 3)-shaped vector and associated metadata
export exportSpecies!
function exportSpecies!(slist::Vector{<:AbstractSpecies}, filename::String)
    filename = string(filename, ".json")
    data = Dict()
    
    for (i, si) in enumerate(slist)
        data[i] = Dict()

        # Extract coordinates (first 3 elements of each position vector) if possible
        coords = Array{Int32}(undef, 0, 3)
        if hasproperty(si, :length) && hasproperty(si, :pos)
            try
                coords = Matrix{Int32}(undef, Int(si.length), 3)
                extractCoords!(si, coords)
            catch
                # leave coords empty on failure
                coords = Array{Int32}(undef, 0, 3)
            end
        end
        data[i]["coords"] = coords

        # helper to return property value or nothing if missing
        issetprop = (obj, name) -> hasproperty(obj, name) ? getproperty(obj, name) : nothing

        data[i]["name"] = issetprop(si, :name)
        data[i]["maxlength"] = issetprop(si, :maxlength)
        data[i]["boxlength"] = issetprop(si, :boxlength)
        data[i]["chemPot"] = issetprop(si, :chemPot)
        data[i]["intPot"] = issetprop(si, :intPot)
        data[i]["initPos"] = issetprop(si, :initPos)
        data[i]["finalPos"] = issetprop(si, :finalPos)
    end
    
    open(filename, "w") do f
        JSON.print(f, data, 4)
    end
end

# exports the full state of a lattice
export export!
function export!(l::Lattice, filename::String)
    filename = string(filename, ".json")
    open(filename, "w") do f
        JSON.print(f, l.lat, 4)
    end
end