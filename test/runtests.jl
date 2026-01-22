using Droplet
using Test

@testset "Droplet.jl" begin
    include("test_init.jl")
    include("test_hamiltonians.jl")
    include("test_connected_components.jl")
    include("test_density_droplet.jl")
    include("test_interaction_potentials.jl")
    include("test_polymonomer.jl")
end
