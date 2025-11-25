const KEYCOUNT = Ref(1)
const SENTINEL = Int32[-1, -1, -1]
const COORD_MOVES = [1,1,2,2]
const MOVE_ORIENTATIONS = [+1,-1,+1,-1]
const SPECIES_ID_COUNTER = Ref(1)
const PKINK = 0.8

abstract type AbstractSpecies end
abstract type AbstractPoly <: AbstractSpecies end
abstract type AbstractBulk <: AbstractSpecies end 
abstract type AbstractFastBulk <: AbstractBulk end 


######## Generic species #######

export Species
Base.@kwdef mutable struct Species <: AbstractSpecies
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Union{Nothing, Int32} = nothing       # number of beads
    maxlength::Union{Nothing, Int32} = nothing   # max number of beads
    boxlength::Union{Nothing, Int32} = nothing    # box linear dimension
    pos::Union{Nothing, OrderedDict{Int32,Vector{Int32}}} = nothing   # xyz coordinates of all the beads
    head::Union{Nothing, Int32} = nothing
    latS::Union{Nothing, AbstractArray{Int32}} = nothing   # lattice
    tempVec::Union{Nothing, Vector{Int32}} = nothing
    nnIdxs::Union{Nothing, Vector{Int32}} = nothing
    freeNNs::Union{Nothing, Set{Int32}} = nothing
    
end

######## Polymer #######

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
    initShape::Union{Nothing, String} = nothing
    finalPos::Union{Nothing, Vector{Int32}} = nothing   # anchor position on the other side

    # lattice, ghost and partner species
    latS::Union{Nothing, AbstractArray{Int32}} = nothing   # lattice
    ghost::Union{Nothing, AbstractSpecies} = nothing
    partner::Union{Nothing, AbstractSpecies} = nothing
    intPotP::Union{Nothing, Float64} = nothing    # inter-species interaction potential
    
    # addBond!, removeBond!, identifyBond! variables
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
    hastingsFactor::Union{Nothing, Float64} = nothing
    
    # energy
    ham::Union{Nothing, Float64} = nothing       # Hamiltonian
    deltaH::Union{Nothing, Float64} = nothing    # Hamiltonian increment of a given move
end


######## Bulk #######

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
    initPos::Union{Nothing, Vector{Int32}} = nothing     # anchor position on one side
    initShape::Union{Nothing, String} = nothing
    initParams::Union{Nothing, Int32} = nothing

    # lattice, ghost and partner species
    latS::Union{Nothing, AbstractArray{Int32}} = nothing
    ghost::Union{Nothing, AbstractSpecies} = nothing
    partner::Union{Nothing, AbstractSpecies} = nothing
    intPotP::Union{Nothing, Float64} = nothing    # inter-species interaction potential

    # moveTypes! variables
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    moveIdx::Union{Nothing, Int32} = nothing 
    moveVec::Union{Nothing, Vector{Int32}} = nothing
    newPos::Union{Nothing, Vector{Int32}} = nothing
    moveType::Union{Nothing, Int32} = nothing     # array of possible move types around a given bead
    tempVec::Union{Nothing, Vector{Int32}} = nothing
    nnIdxs::Union{Nothing, Vector{Int32}} = nothing

    # proposeMove! variables
    movableIdxs::Union{Nothing, Vector{Int32}} = nothing  # coordinates of unconfined particles
    freeNNs::Union{Nothing, Set{Int32}} = nothing   # move proposal probabilities
    nnsOld::Union{Nothing, Int32} = nothing          # number of nearest neighbs of old site
    nnsNew::Union{Nothing, Int32} = nothing          # number of nearest neighbs of new site
    hastingsFactor::Union{Nothing, Float64} = nothing
    
    # energy
    ham::Union{Nothing, Float64} = nothing       # Hamiltonian
    deltaH::Union{Nothing, Float64} = nothing    # Hamiltonian increment of a given move
end


######## FastBulk #######

export FastBulk
Base.@kwdef mutable struct FastBulk <: AbstractFastBulk
    
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
    initPos::Union{Nothing, Vector{Int32}} = nothing     # anchor position on one side
    initShape::Union{Nothing, String} = nothing
    initParams::Union{Nothing, Int32} = nothing

    # lattice, ghost and partner species
    latS::Union{Nothing, AbstractArray{Int32}} = nothing
    ghost::Union{Nothing, AbstractSpecies} = nothing
    partner::Union{Nothing, AbstractSpecies} = nothing
    intPotP::Union{Nothing, Float64} = nothing    # inter-species interaction potential

    # moveTypes! variables
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    moveIdx::Union{Nothing, Int32} = nothing 
    moveVec::Union{Nothing, Vector{Int32}} = nothing
    newPos::Union{Nothing, Vector{Int32}} = nothing
    moveType::Union{Nothing, Int32} = nothing     # array of possible move types around a given bead
    tempVec::Union{Nothing, Vector{Int32}} = nothing
    nnIdxs::Union{Nothing, Vector{Int32}} = nothing

    # proposeMove! variables
    movableIdxs::Union{Nothing, Set{Int32}} = nothing  # coordinates of unconfined particles
    freeNNs::Union{Nothing, Set{Int32}} = nothing   # move proposal probabilities
    nnsOld::Union{Nothing, Int32} = nothing          # number of nearest neighbs of old site
    nnsNew::Union{Nothing, Int32} = nothing          # number of nearest neighbs of new site
    hastingsFactor::Union{Nothing, Float64} = nothing
    
    # energy
    ham::Union{Nothing, Float64} = nothing       # Hamiltonian
    deltaH::Union{Nothing, Float64} = nothing    # Hamiltonian increment of a given move
end


######## Lattice #######

export Lattice
Base.@kwdef mutable struct Lattice
    boxlength::Int32     # lattice linear dimension 
    nspecies::Int32      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing  # lattice array
    speciesDir::Union{Nothing, Dict{Int32,Int32}} = nothing    # id-indexed locations of each species in the lattice array
end
