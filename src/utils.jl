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

# extracts the coordinates of a species' beads
export extractCoords!
function extractCoords!(s::AbstractSpecies, coords::Matrix)
    key = s.head           # initial bead
    k = 1
    
    for k in 1:s.length
        if key != 0
            coords[k, :] = s.pos[key][1:3] # next position
            key = s.pos[key][5]    # next bead
        end
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

        # Extract coordinates (first 3 elements of each position vector)
        coords = Matrix{Int32}(undef, si.length, 3)
        extractCoords!(si, coords)

        # @assert k == si.length begin
        #     "Tail is not N'th bead"
        # end

        data[i]["coords"] = coords
        data[i]["name"] = si.name
        data[i]["maxlength"] = si.maxlength
        data[i]["boxlength"] = si.boxlength
        data[i]["chemPot"] = si.chemPot
        data[i]["intPot"] = si.intPot
        data[i]["initPos"] = si.initPos
        data[i]["finalPos"] = si.finalPos
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