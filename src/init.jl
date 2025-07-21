# initializes bead locations of a long polymer
export initPoly!
function initPoly!(p::Poly)
    p.pos = zeros(Int32, p.length, 3)

    # line from center of one face to the other along x axis
    for i in 1:p.length
        v = view(p.pos, i, :)
        v[1] = i
        v[2] = p.init_pos === nothing ? div(p.length, 2) : p.init_pos[2]
        v[3] = p.init_pos === nothing ? div(p.length, 2) : p.init_pos[3]
    end

    return p
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initDilute!
function initDilute!(p::Poly)
    p.pos = rand(1:p.boxlength, p.length, 3)

    return p
end

# initializes lattice sites of each species from their bead locations
export initLattice!
function initLattice!(l::Lattice, polys::Vector{Poly})
    l.nspecies = max(l.nspecies, length(polys))   # update number of species
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))

    # fill out lattice sites for each species
    for (i, p) in enumerate(polys)
        lati = view(l.lat, i, :, :, :)

        for j in 1:p.length
            v = view(p.pos, j, :)
            lati[CartesianIndex((v[1], v[2], v[3]))] = 1
        end
    end

    return l
end