# returns the center of mass of a species
export centerMass
function centerMass(s::AbstractSpecies)
    cm = zeros(3)

    # compute average of all positions in pos
    for i in eachindex(s.pos)
        cm = cm .+ view(s.pos[i], 1:3)
    end

    return cm ./ s.length
end

# radius of gyration, mean squared distance btw each bead and
# the center of mass
export radius
function radius(s::AbstractSpecies)
    cm = centerMass(s)
    rg2 = 0

    for i in eachindex(s.pos)
        rg2 += ((cm[1] .- s.pos[i][1])^2
             + (cm[2] .- s.pos[i][2])^2
             + (cm[3] .- s.pos[i][3])^2)
    end

    return sqrt(rg2 / s.length)
end

# gyration radius of a line
export radiusLine
function radiusLine(length)
    n = div(length, 2)
    return sqrt((n - 1) / 2)
end

# gyration radius of a sphere
export radiusSphere
function radiusSphere(length)
    r = (3*length/(4*pi))^(1/3)
    return sqrt(3/5)*r
end

# connected components
# pseudocode:

# ccs = []
# DFS(1, [])

# function DFS(bi, marked)
#   if bi not in marked
#       push!(marked, bi)

#   allMarked = true
#
#   for bj in nearestNeighbs(bi)
#       if bj != 0
    #       if bj not in marked
    #           allMarked = false
    #           marked = DFS(bj, marked)
#   
#   if allMarked
#       return marked

export connComp
function connComp(s::AbstractSpecies)

    ccs = []
    

end