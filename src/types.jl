const SPECIES_ID_COUNTER = Ref(1)

mutable struct Species
    id::Int
    name::String       # unique name of this species
    length::Int      # number of beads
    maxlength::Int    # max number of beads
    boxlength::Int     # box linear dimension

    pos::Array{Int32,2}      # xyz coordinates of all the beads
    initPos::Union{Nothing, Vector{Int}}    # anchor position on one side
    finalPos::Union{Nothing, Vector{Int}}     # anchor position on the other side
    auxPos1::Union{Nothing, Vector{Int}}    # placeholder position (e.g. to store next bond position)
    auxPos2::Union{Nothing, Vector{Int}}     # placeholder position (e.g. to store next bond position)
    bondVec::Union{Nothing, Vector{Int}}    # direction vector for a given bond
    
    # interaction energies
    # chemical potentials
    # number of monomer species, their density
end

export Species
function Species(; length, maxlength, boxlength)
    SPECIES_ID_COUNTER[] += 1
    return Species(SPECIES_ID_COUNTER[], 
                   "", 
                   length, 
                   maxlength, 
                   boxlength, 
                   zeros(Int32, length, 3), 
                   nothing,
                   nothing,
                   zeros(Int32, 3),
                   zeros(Int32, 3),
                   nothing)
end

export Lattice
Base.@kwdef mutable struct Lattice
    boxlength::Int     # lattice linear dimension 
    nspecies::Int      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing   # lattice array
    speciesIndex::Dict{Int,Int} = Dict{Int,Int}()
end