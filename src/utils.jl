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

# exports a lattice's state to JSON for python compatibility
export export!
function export!(l::Lattice, filename::String)
    filename = string(filename, ".json")
    open(filename, "w") do f
        JSON.print(f, l.lat, 4)
    end
end