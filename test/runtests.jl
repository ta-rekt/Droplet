using Droplet
using Test
using Random

@testset "Droplet.jl" begin
    include("test_init.jl")
    include("test_hamiltonians.jl")
    include("test_connected_components.jl")
    include("test_density_droplet.jl")
    include("test_interaction_potentials.jl")
    include("test_polymonomer.jl")
    include("test_heterogeneous_phase0.jl")
    include("test_heterogeneous_phase2.jl")
    include("test_heterogeneous_phase3.jl")
    include("test_heterogeneous_phase5.jl")
    include("test_heterogeneous_phase6.jl")
    include("test_updateMonomers.jl")
    include("test_pattern_stats.jl")
    include("test_centerMassPattern.jl")
    include("test_contact_matrix.jl")
    include("test_find_patterns.jl")
    include("test_grow_pattern.jl")
    include("test_checkpoint.jl")
    include("test_restore_species.jl")
    include("test_rebuild_lattice.jl")
    include("test_restore.jl")
    include("test_integration_checkpoint.jl")
    # Three-way interaction tests
    include("test_nnpq.jl")
    include("test_couplingPot.jl")
    include("test_ham_coupling.jl")
end
