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
        coords = Array{Int32}(undef, 0, 4)
        if hasproperty(si, :length) && hasproperty(si, :pos)
            try
                coords = Matrix{Int32}(undef, Int(si.length), 4)
                extractCoords!(si, coords)
            catch e
                println(e)
                # leave coords empty on failure
                coords = Array{Int32}(undef, 0, 4)
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

# saves species mutable state to a checkpoint file
export checkpoint!
function checkpoint!(filename::String, slist::Vector{<:AbstractSpecies})
    species = Dict{Symbol,Any}[]
    for s in slist
        d = Dict{Symbol,Any}(:id => s.id, :pos => s.pos, :length => s.length)
        if s isa Polymer
            d[:type] = :Polymer
            for f in (:head, :keys_vec, :pattern, :customPattern, :patternTails, :baseIntIdx)
                d[f] = getfield(s, f)
            end
        elseif s isa FastBulkMonomer
            d[:type] = :FastBulkMonomer
        elseif s isa BulkMonomer
            d[:type] = :BulkMonomer
        end
        push!(species, d)
    end
    Serialization.serialize("$filename.jls", Dict(:version => 1, :species => species))
end

# restores mutable state from a checkpoint dict onto a live species object
function _validateRestore(s::AbstractSpecies, d::Dict{Symbol,Any}, expected::Symbol)
    d[:type] == expected || error("Type mismatch: expected $expected, got $(d[:type])")
    d[:id] == s.id       || error("ID mismatch: expected $(s.id), got $(d[:id])")
end

function _restoreSpecies!(p::Polymer, d::Dict{Symbol,Any})
    _validateRestore(p, d, :Polymer)
    p.pos = d[:pos]; p.length = d[:length]; p.head = d[:head]
    p.keys_vec = d[:keys_vec]; p.pattern = d[:pattern]
    p.customPattern = d[:customPattern]; p.patternTails = d[:patternTails]
    p.baseIntIdx = d[:baseIntIdx]
    # scratch buffers (mirrors initSpecies!)
    p.tempVec = zeros(Int32, 3); p.bondVec = zeros(Int32, 3)
    p.moveVec = zeros(Int32, 3); p.newPosi = zeros(Int32, 3)
    p.newPosj = zeros(Int32, 3); p.moveTypes = zeros(Int32, 4)
    p.moveTypeProbs = [0, PKINK, 0]; p.nTypes = zeros(Int32, 3)
    p.hastingsFactor = 1.0; p.deltaHP = 0.0
end

function _restoreSpecies!(b::BulkMonomer, d::Dict{Symbol,Any})
    _validateRestore(b, d, :BulkMonomer)
    b.pos = d[:pos]; b.length = d[:length]
    b.freeNNs = Set{Int32}(); b.hastingsFactor = 1.0
    b.nnsOld = 0; b.nnsNew = 0
    b.newPos = zeros(Int32, 3); b.tempVec = zeros(Int32, 3)
    b.nnIdxs = zeros(Int32, 6); b.deltaH = 0.0; b.deltaHP = 0.0
    b.chemPot === nothing && (b.chemPot = 0.0)
    b.intPot === nothing && (b.intPot = 0.0)
end

function _restoreSpecies!(b::FastBulkMonomer, d::Dict{Symbol,Any})
    _validateRestore(b, d, :FastBulkMonomer)
    b.pos = d[:pos]; b.length = d[:length]
    b.movableIdxs = Set{Int32}(); b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0; b.nnsOld = 0; b.nnsNew = 0
    b.bi = zeros(Int32, 3); b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3); b.moveVec = zeros(Int32, 2)
    b.nnIdxs = zeros(Int32, 6); b.deltaH = 0.0; b.deltaHP = 0.0
    b.chemPot === nothing && (b.chemPot = 0.0)
    b.intPot === nothing && (b.intPot = 0.0)
end

# rebuilds lattice, latS views, ghost, movableIdxs, KEYCOUNT, and ham from restored species
function _rebuildLattice!(l::Lattice, slist::Vector{<:AbstractSpecies})
    # separate layered species from poly-monomers (same logic as init!)
    layered = AbstractSpecies[]
    parents = Set{AbstractPolymer}()
    for s in slist
        s isa AbstractPolyMonomer ? push!(parents, s.poly) : push!(layered, s)
    end
    all = vcat(layered, collect(parents))

    l.nspecies = length(all)
    l.lat = zeros(Int32, l.nspecies + 1, l.boxlength, l.boxlength, l.boxlength)
    l.speciesDir = Dict{Int32, Int32}()

    # ghost species
    gs = Species(;)
    gs.latS = view(l.lat, 1, :, :, :)
    gs.boxlength = l.boxlength
    initSpecies!(gs)

    # recompute KEYCOUNT from all pos dicts
    maxkey = Int32(0)
    for s in all
        for k in keys(s.pos)
            k > maxkey && (maxkey = k)
        end
    end
    KEYCOUNT[] = maxkey

    for (i, s) in enumerate(all)
        s.latS = view(l.lat, i + 1, :, :, :)
        s.ghost = gs
        s.boxlength = l.boxlength
        s.name = string(s.name, "_", ALPHABET[i])
        l.speciesDir[s.id] = i

        # populate lattice from pos
        if s isa AbstractBulkMonomer
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = j
            end
        elseif s isa AbstractPolymer && s.monomers !== nothing
            key = s.head
            while key != 0
                p = s.pos[key]
                s.latS[p[1], p[2], p[3]] = p[6]
                key = p[5]
            end
        else
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
            end
        end

        s isa FastBulkMonomer && initMovableSet!(s)
    end

    # recompute hamiltonians
    for s in all
        hasproperty(s, :ham) && ham!(s)
    end

    return l
end

# restores full simulation state from a checkpoint file
export restore!
function restore!(filename::String, l::Lattice, slist::Vector{<:AbstractSpecies})
    data = Serialization.deserialize("$filename.jls")
    saved = data[:species]
    length(saved) == length(slist) || error("Species count mismatch: checkpoint has $(length(saved)), got $(length(slist))")

    # match by id and restore each species
    saved_by_id = Dict(d[:id] => d for d in saved)
    for s in slist
        haskey(saved_by_id, s.id) || error("No checkpoint entry for species id=$(s.id)")
        _restoreSpecies!(s, saved_by_id[s.id])
    end

    _rebuildLattice!(l, slist)
    return l
end

# exports the full state of a lattice
export export!
function export!(l::Lattice, filename::String)
    filename = string(filename, ".json")
    open(filename, "w") do f
        JSON.print(f, l.lat, 4)
    end
end