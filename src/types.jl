const SPECIES_ID_COUNTER = Ref(1)

mutable struct Species
    id::Int32
    name::String         # unique name of this species
    length::Int32        # number of beads
    maxlength::Int32     # max number of beads
    boxlength::Int32     # box linear dimension

    # position and initial, final bead positions
    pos::OrderedDict{Int32,Vector{Int32}} = nothing       # xyz coordinates of all the beads
    initPos::Union{Nothing, Vector{Int32}} = nothing      # anchor position on one side
    finalPos::Union{Nothing, Vector{Int32}} = nothing     # anchor position on the other side

    # addBond!, removeBond!, identifyBond! variables
    latS::Union{Nothing, AbstractArray{Int32}} = nothing
    ki::Union{Nothing, Int32} = nothing
    kj::Union{Nothing, Int32} = nothing
    kprev::Union{Nothing, Int32} = nothing
    knext::Union{Nothing, Int32} = nothing
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    bj::Union{Nothing, AbstractArray{Int32}} = nothing
    bprev::Union{Nothing, AbstractArray{Int32}} = nothing
    bnext::Union{Nothing, AbstractArray{Int32}} = nothing

    # moveTypes! variables
    bondVec::Union{Nothing, Vector{Int32}} = nothing         # direction vector for a given bond
    moveVec::Union{Nothing, Vector{Int32}} = nothing         # move vector
    moveIdx::Union{Nothing, Int32}  = nothing                # bead index to move next
    moveType::Union{Nothing, Int32} = nothing                # type of move: bond addition (1), bond removal (2) or kink (3)
    newPosi::Union{Nothing, Vector{Int32}} = nothing         # placeholder position (e.g. to store next bond position)
    newPosj::Union{Nothing, Vector{Int32}}  = nothing        # placeholder position (e.g. to store next bond position)
    moveTypes::Union{Nothing, Vector{Int32}} = nothing       # array of possible move types around a given bead

    # proposeMove! variables
    moveProbs::Union{Nothing, Vector{Float64}} = nothing     # move proposal probabilities
    
    # interaction energies
    # chemical potentials
    # number of monomer species, their density
end

export Species
function Species(; length, maxlength, boxlength)
    SPECIES_ID_COUNTER[] += 1
    return Species(; SPECIES_ID_COUNTER[], 
                   "", 
                   length, 
                   maxlength, 
                   boxlength)
end

export Lattice
Base.@kwdef mutable struct Lattice
    boxlength::Int32     # lattice linear dimension 
    nspecies::Int32      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing   # lattice array
    speciesDir::Dict{Int32,Int32} = Dict{Int32,Int32}()   # id-indexed locations of each species in the lattice array
end