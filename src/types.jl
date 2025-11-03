const KEYCOUNT = Ref(1)
const SENTINEL = Int32[-1, -1, -1]
const COORD_MOVES = [1,1,2,2]
const MOVE_ORIENTATIONS = [+1,-1,+1,-1]
const SPECIES_ID_COUNTER = Ref(1)
const PKINK = 0.8

abstract type AbstractSpecies end
abstract type AbstractPoly <: AbstractSpecies end
abstract type AbstractBulk <: AbstractSpecies end 

export Poly
Base.@kwdef mutable struct Poly <: AbstractPoly
    
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Int32        # number of beads
    maxlength::Int32     # max number of beads
    chemPot::Union{Nothing,Float64}   # monomer chemical potential
    intPot::Union{Nothing,Float64}    # monomer interaction potential
    
    # key vector, position and initial, final bead positions
    boxlength::Union{Nothing, Int32} = nothing    # box linear dimension
    keys_vec::Union{Nothing, Vector{Int32}} = nothing  # pre-allocated key vector
    pos::Union{Nothing, OrderedDict{Int32,Vector{Int32}}} = nothing      # xyz coordinates of all the beads
    initPos::Union{Nothing, Vector{Int32}} = nothing     # anchor position on one side
    finalPos::Union{Nothing, Vector{Int32}} = nothing   # anchor position on the other side

    # addBond!, removeBond!, identifyBond! variables
    latS::Union{Nothing, AbstractArray{Int32}} = nothing
    head::Union{Nothing, Int32} = nothing
    ki::Union{Nothing, Int32} = nothing
    kj::Union{Nothing, Int32} = nothing
    kprev::Union{Nothing, Int32} = nothing
    knext::Union{Nothing, Int32} = nothing
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    bj::Union{Nothing, AbstractArray{Int32}} = nothing
    bprev::Union{Nothing, AbstractArray{Int32}} = nothing
    bnext::Union{Nothing, AbstractArray{Int32}} = nothing
    bondVec::Union{Nothing, Vector{Int32}} = nothing       # direction vector for a given bond
    tempVec::Union{Nothing, Vector{Int32}} = nothing

    # moveTypes! variables
    moveVec::Union{Nothing, Vector{Int32}} = nothing       # move vector
    moveIdx::Union{Nothing, Int32}      = nothing          # bead index to move next
    moveType::Union{Nothing, Int32} = nothing              # type of move: bond addition (1), bond removal (2) or kink (3)
    newPosi::Union{Nothing, Vector{Int32}} = nothing       # placeholder position (e.g. to store next bond position)
    newPosj::Union{Nothing, Vector{Int32}} = nothing       # placeholder position (e.g. to store next bond position)
    moveTypes::Union{Nothing, Vector{Int32}} = nothing     # array of possible move types around a given bead

    # proposeMove! variables
    moveTypeProbs::Union{Nothing, Vector{Float64}} = nothing   # move proposal probabilities
    nTypes::Union{Nothing, Vector{Float64}} = nothing          # how many of each move type there are
    
    # energy
    ham::Union{Nothing, Float64} = nothing       # Hamiltonian
    deltaH::Union{Nothing, Float64} = nothing    # Hamiltonian increment of a given move
end


export Bulk
Base.@kwdef mutable struct Bulk <: AbstractBulk
    
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Int32        # number of beads
    chemPot::Union{Nothing,Float64}   # monomer chemical potential
    intPot::Union{Nothing,Float64}    # monomer interaction potential

    # key vector, position and initial, final bead positions
    boxlength::Union{Nothing, Int32} = nothing
    keys_vec::Union{Nothing, Vector{Int32}} = nothing  # pre-allocated key vector
    pos::Union{Nothing, OrderedDict{Int32,Vector{Int32}}} = nothing      # xyz coordinates of all the beads
    
    # lattice view
    latS::Union{Nothing, AbstractArray{Int32}} = nothing

    # moveTypes! variables
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    moveIdx::Union{Nothing, Int32}      = nothing 
    moveVec::Union{Nothing, Vector{Int32}} = nothing
    newPos::Union{Nothing, Vector{Int32}} = nothing
    moveType::Union{Nothing, Int32} = nothing     # array of possible move types around a given bead
    tempVec::Union{Nothing, Vector{Int32}} = nothing

    # proposeMove! variables
    moveTypeProbs::Union{Nothing, Vector{Float64}} = nothing   # move proposal probabilities
    nTypes::Union{Nothing, Vector{Float64}} = nothing          # how many of each move type there are
    
    # energy
    ham::Union{Nothing, Float64} = nothing       # Hamiltonian
    deltaH::Union{Nothing, Float64} = nothing    # Hamiltonian increment of a given move
  end


export Lattice
Base.@kwdef mutable struct Lattice
    boxlength::Int32     # lattice linear dimension 
    nspecies::Int32      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing  # lattice array
    speciesDir::Union{Nothing, Dict{Int32,Int32}} = nothing    # id-indexed locations of each species in the lattice array
end