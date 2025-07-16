export Poly
Base.@kwdef mutable struct Poly

    pos        # xyz coordinates of all the beads
    length      # number of beads
    maxlength    # max number of beads

    init_pos = nothing;     # anchor position on one side
    final_pos = nothing;     # anchor position on the other side

end

export Lattice
Base.@kwdef mutable struct Lattice

    lat
    dim
    nspecies

end