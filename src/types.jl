const SPECIES_ID_COUNTER = Ref(1)

mutable struct Species
    id::Int
    name::String       # unique name of this species
    length::Int      # number of beads
    maxlength::Int    # max number of beads
    boxlength::Int     # box linear dimension

    pos::OrderedDict{Int32,Vector{Int32}}      # xyz coordinates of all the beads
    initPos::Union{Nothing, Vector{Int32}}    # anchor position on one side
    finalPos::Union{Nothing, Vector{Int32}}     # anchor position on the other side
    moveVec::Union{Nothing, Vector{Int32}}     # move vector
    bondVec::Union{Nothing, Vector{Int32}}    # direction vector for a given bond
    moveIdx::Int            # bead index to move next
    moveType::Int           # type of move: bond addition (1), bond removal (2) or kink (3)
    newPosi::Vector{Int32}    # placeholder position (e.g. to store next bond position)
    newPosj::Vector{Int32}    # placeholder position (e.g. to store next bond position)
    
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
                   OrderedDict{Int32, Vector{Int32}}(), 
                   nothing,
                   nothing,
                   nothing,
                   nothing,
                   0,
                   0,
                   zeros(Int32, 3),
                   zeros(Int32, 3))
end

export Lattice
Base.@kwdef mutable struct Lattice
    boxlength::Int     # lattice linear dimension 
    nspecies::Int      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing   # lattice array
    speciesDir::Dict{Int,Int} = Dict{Int,Int}()   # id-indexed locations of each species in the lattice array
end