const KEYCOUNT = Ref(1)
const SENTINEL = Int32[-1, -1, -1]
const COORD_MOVES = [1,1,2,2]
const MOVE_ORIENTATIONS = [+1,-1,+1,-1]
const SPECIES_ID_COUNTER = Ref(1)
const PKINK = 0.8
const ALPHABET = "abcdefghijklmnopqrstuvwxyz"

abstract type AbstractSpecies end
abstract type AbstractPolymer <: AbstractSpecies end
abstract type AbstractMonomer <: AbstractSpecies end
abstract type AbstractBulkMonomer <: AbstractMonomer end 
abstract type AbstractPolyMonomer <: AbstractMonomer end 
abstract type AbstractFastBulkMonomer <: AbstractBulkMonomer end 

export Species, Polymer, PolyMonomer, BulkMonomer, FastBulkMonomer, Lattice,
       KMCSystem, OneBurstDecoyModel, LLROptModel

######## Generic species #######
# MARK: -Generic

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
# MARK: -Polymer

# Forward declare to allow PolyMonomer to reference it

######## Monomer (for bead-dependent interactions) #######
# MARK: -PolyMonomer
Base.@kwdef mutable struct PolyMonomer <: AbstractPolyMonomer
    id::Int32 = begin
        SPECIES_ID_COUNTER[] += 1
        SPECIES_ID_COUNTER[]
    end
    
    name::String = ""    # unique name of this Monomer

    # host polymer
    poly::AbstractPolymer
    
    # interaction potentials
    chemPot::Union{Nothing, Float64} = nothing   # chemical potential for this Monomer
    intPot::Union{Nothing, Float64} = nothing    # self-interaction potential
    
    # multi-Monomer interaction fields (mirrors AbstractBulkMonomer pattern)
    intSpecies::Union{Nothing, Vector{Tuple{AbstractMonomer, Int32}}} = nothing  # list of (Monomer, canonical index) pairs
    intMat::Union{Nothing, Matrix{Float64}} = nothing  # interaction matrix (symmetric, diagonal = self intPot)
    intIdx::Union{Nothing, Int32} = nothing  # this Monomer' canonical index in intMat
  end


# Now define Polymer completely
Base.@kwdef mutable struct Polymer <: AbstractPolymer
    
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Int32        # number of beads
    maxlength::Int32     # max number of beads
    nMonomerTypes::Union{Nothing, Int32} = nothing  # number of distinct PolyMonomer
    
    # key vector, position and initial, final bead positions
    boxlength::Union{Nothing, Int32} = nothing    # box linear dimension
    keys_vec::Union{Nothing, Vector{Int32}} = nothing  # pre-allocated key vector
    pos::Union{Nothing, OrderedDict{Int32, Vector{Int32}}} = nothing      # xyz coordinates of all the beads
    initPos::Union{Nothing, Vector{Int32}} = nothing     # anchor position on one side
    initShape::Union{Nothing, String} = nothing
    finalPos::Union{Nothing, Vector{Int32}} = nothing   # anchor position on the other side
    
    # lattice, ghost and partner species
    latS::Union{Nothing, AbstractArray{Int32}} = nothing   # lattice
    ghost::Union{Nothing, AbstractSpecies} = nothing
    
    # Monomer for bead-dependent interactions
    # Dictionary keyed by intIdx for O(1) lookup; pos[key][6] stores the intIdx
    monomers::Union{Nothing, Dict{Int32, PolyMonomer}} = nothing
    pattern::Symbol=:uniform
    customPattern::Union{Nothing, Vector{Int}}=nothing
    patternTails::Union{Nothing, Dict{Int32, Int32}} = nothing  # intIdx -> tail key for non-base types
    baseIntIdx::Union{Nothing, Int32} = nothing                  # intIdx of base type (can grow/shrink)
    
    # addBond!, removeBond!, identifyBond! variables
    head::Union{Nothing, Int32} = nothing
    ki::Union{Nothing, Int32} = nothing
    kj::Union{Nothing, Int32} = nothing
    kprev::Union{Nothing, Int32} = nothing
    knext::Union{Nothing, Int32} = nothing
    bi::Union{Nothing, Vector{Int32}} = nothing      # pre-allocated bond coordinate vectors
    bj::Union{Nothing, Vector{Int32}} = nothing
    bprev::Union{Nothing, Vector{Int32}} = nothing
    bnext::Union{Nothing, Vector{Int32}} = nothing
    bondVec::Union{Nothing, Vector{Int32}} = nothing       # direction vector for a given bond
    tempVec::Union{Nothing, Vector{Int32}} = nothing

    # moveTypes! variables
    moveVec::Union{Nothing, Vector{Int32}} = nothing       # move vector
    moveIdx::Union{Nothing, Int32} = nothing          # bead index to move next
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
    deltaHP::Union{Nothing, Float64} = nothing   # Hamiltonian increment from partner move
end


######## BulkMonomer #######
# MARK: -BulkMonomer

Base.@kwdef mutable struct BulkMonomer <: AbstractBulkMonomer
    
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Int32        # number of beads

    chemPot::Union{Nothing,Float64} = nothing   # monomer chemical potential
    intPot::Union{Nothing,Float64} = nothing   # monomer interaction potential
    k0::Float64 = 1.0   # KMC base hop rate (sets simulation time unit)

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

    # multi-species interaction fields
    intSpecies::Union{Nothing, Vector{Tuple{AbstractMonomer, Int32}}} = nothing  # list of (species, canonical index) pairs
    intMat::Union{Nothing, Matrix{Float64}} = nothing  # interaction matrix (symmetric, diagonal = self intPot)
    intIdx::Union{Nothing, Int32} = nothing  # this species' canonical index in intMat

    # three-way coupling fields (BulkMonomer only)
    # couplingMat[i,j,k] = energy contribution when self (couplingIdx=i) overlaps context k and neighbors partner j
    # Symmetry: couplingMat[i,j,k] == couplingMat[i,k,j] (partner/context interchangeable)
    couplingSpecies::Union{Nothing, Vector{Tuple{AbstractBulkMonomer, Int32}}} = nothing  # list of (BulkMonomer, couplingIdx) pairs
    couplingMat::Union{Nothing, Array{Float64, 3}} = nothing  # m × n × n array (m = BulkMonomers, n = all species)
    couplingIdx::Union{Nothing, Int32} = nothing  # this BulkMonomer's index in couplingMat first dimension

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
    deltaHP::Union{Nothing, Float64} = nothing   # Hamiltonian increment from partner move
end


######## FastBulkMonomer #######
# MARK: -FastBulkMonomer

Base.@kwdef mutable struct FastBulkMonomer <: AbstractFastBulkMonomer
    
    id::Int32 = begin
      SPECIES_ID_COUNTER[] += 1
      SPECIES_ID_COUNTER[]
    end

    name::String = ""    # unique name of this species

    # required fields
    length::Int32        # number of beads

    chemPot::Union{Nothing,Float64} = nothing   # monomer chemical potential
    intPot::Union{Nothing,Float64} = nothing   # monomer interaction potential

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

    # multi-species interaction fields
    intSpecies::Union{Nothing, Vector{Tuple{AbstractMonomer, Int32}}} = nothing  # list of (species, canonical index) pairs
    intMat::Union{Nothing, Matrix{Float64}} = nothing  # interaction matrix (symmetric, diagonal = self intPot)
    intIdx::Union{Nothing, Int32} = nothing  # this species' canonical index in intMat

    # moveTypes! variables
    bi::Union{Nothing, AbstractArray{Int32}} = nothing
    moveIdx::Union{Nothing, Int32} = nothing 
    moveIdxNew::Union{Nothing, Int32} = nothing 
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
    deltaHP::Union{Nothing, Float64} = nothing   # Hamiltonian increment from partner move
end


######## Lattice #######

Base.@kwdef mutable struct Lattice
    boxlength::Int32     # lattice linear dimension
    nspecies::Int32      # number of species
    lat::Union{Nothing, Array{Int32,4}} = nothing  # lattice array
    speciesDir::Union{Nothing, Dict{Int32,Int32}} = nothing    # id-indexed locations of each species in the lattice array
end


######## KMCSystem #######
# MARK: -KMCSystem

# Kinetic Monte Carlo (Gillespie) state container.
# Each "slot" encodes a candidate hop (speciesIdx, particleId, dir ∈ 1..6).
# Slot id = species_offset[speciesIdx] + (particleId - 1) * 6 + dir.
Base.@kwdef mutable struct KMCSystem
    slist::AbstractVector             # species participating in KMC
    beta::Float64 = 1.0               # inverse temperature
    rates::Vector{Float64} = Float64[]
    fenwick::Vector{Float64} = Float64[]
    species_offset::Vector{Int} = Int[]
    nslots::Int = 0
    Q::Float64 = 0.0
    t_sim::Float64 = 0.0
    n_events::Int = 0
    src::Vector{Int32} = zeros(Int32, 3)
    dst::Vector{Int32} = zeros(Int32, 3)
    refresh_set::Set{Int} = Set{Int}()
    nbuf::Vector{Int32} = zeros(Int32, 3)
    ybuf::Vector{Int32} = zeros(Int32, 3)
end


######## Receptor models ########

# constant memory objects
Base.@kwdef mutable struct LLROptModel
    # input params
    weighted_edges
    input_edges
    c
    dc
    on_ind
    idist
    n
    k
    
    # build_mat
    Q_c = nothing; # to initialize
    Q_cdc = nothing; # to initialize

    # generate
    t = nothing; # to initialize
    s = nothing; # to initialize

    # get_residence_time
    s_in_on = nothing; # to initialize
    rt = nothing; # to initialize
    onrt = nothing; # to initialize
    offrt = nothing; # to initialize
    jc = nothing; # to initialize
    bc = nothing; # to initialize

    # indicators
    off_ind = nothing; # to initialize
    _idist0_on = nothing; # to initialize
    _idist0_off = nothing; # to initialize

    # ll_full
    ind_set = nothing; # to initialize
    t_set = nothing; # to initialize
    t_comp = nothing; # to initialize
    idist_set = nothing; # to initialize

    ll_c = nothing; # to initialize
    ll_cdc = nothing; # to initialize
    llr = nothing; # to initialize
    ir = nothing; # to initialize

end

Base.@kwdef mutable struct OneBurstDecoyModel
    circuitId
    nstates
    edges
    input_edges
    logrates
    cr
    cw
    gamma
    ind_set
    start_state
    hitting_set
    eps
    reltol
    maxevals

    # indicator array for set
    ind_comp = nothing;
    _set = nothing;
    _comp = nothing;

    # rate matrix
    weights = nothing;
    Q = nothing;

    # set
    D_set = nothing;
    W_set = nothing;
    W_set_inv = nothing;
    aux_set = nothing;
    P_set = nothing;
    Q_h_set = nothing;
    F_set = nothing;

    # comp
    D_comp = nothing;
    D_aux = nothing;
    W_comp = nothing;
    W_comp_inv = nothing;
    aux_comp = nothing;
    P_comp = nothing;
    Q_h_comp = nothing;
    F_comp = nothing;

    # eigenvalues
    eigvals_set = nothing;
    eigvals_comp = nothing;
    maxeigval_set = nothing;
    maxeigval_comp = nothing;

    # likelihood
    Poff = nothing;
    Pon = nothing;
    Pjoint = nothing;

    # log-likelihoods
    logPjoint = nothing;

end