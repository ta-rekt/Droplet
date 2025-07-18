export Poly
Base.@kwdef mutable struct Poly
    length      # number of beads
    maxlength    # max number of beads
    boxlength     # box linear dimension
    
    # interaction energies
    # chemical potentials
    # number of monomer species, their density
    
    pos = nothing          # xyz coordinates of all the beads
    init_pos = nothing;     # anchor position on one side
    final_pos = nothing;     # anchor position on the other side
end

export Lattice
Base.@kwdef mutable struct Lattice
    boxlength   # lattice linear dimension 
    nspecies    # number of species
    
    lat = nothing;   # lattice array
end