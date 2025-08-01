using Droplet
using Test

@testset "Droplet.jl" begin
    @testset "Move Operations" begin
        @testset "Basic Functionality" begin
            @testset "Bond Addition" begin
                l = Lattice(; boxlength = 16, nspecies = 1)
                s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
                # Create a simple straight polymer along x-axis
                s.pos[1] = [5, 5, 5, 0, 2]
                s.pos[2] = [6, 5, 5, 1, 3]
                s.pos[3] = [7, 5, 5, 2, 4]
                s.pos[4] = [8, 5, 5, 3, 5]
                s.pos[5] = [9, 5, 5, 4, 6]
                s.pos[6] = [10, 5, 5, 5, 0]
                
                initLattice!(l, [s])
                
                # Prepare for bond addition
                s.moveIdx = 3  # Choose middle bond (between positions 3 and 4)
                s.moveType = 2  # Bond addition
                s.moveVec = [0, 1, 0]  # Move in y direction
                s.newPosi = [7, 6, 5]
                s.newPosj = [8, 6, 5]
                
                # Original length
                origLength = s.length
                
                # Execute move
                makeMove!(l, s)
                
                # Tests
                @test s.length == origLength + 2  # Length increased by 2
                @test haskey(s.pos, 7)  # New positions exist
                @test haskey(s.pos, 8)
                @test s.pos[3][5] == 7  # Check connections were updated
                @test s.pos[7][4] == 3
                @test s.pos[7][5] == 8
                @test s.pos[8][4] == 7
                @test s.pos[8][5] == 4
                @test s.pos[4][4] == 8
                @test l.lat[1, 7, 6, 5] == 1  # New lattice sites are occupied
                @test l.lat[1, 8, 6, 5] == 1
            end
      
            @testset "Bond Removal" begin
                l = Lattice(; boxlength = 16, nspecies = 1)
                s = Species(; length = 8, maxlength = 10, boxlength = l.boxlength)
                
                # Create a polymer with a detour
                s.pos[1] = [5, 5, 5, 0, 2]
                s.pos[2] = [6, 5, 5, 1, 3]
                s.pos[3] = [6, 6, 5, 2, 4]  # These two form the detour we'll remove
                s.pos[4] = [7, 6, 5, 3, 5]  # These two form the detour we'll remove
                s.pos[5] = [7, 5, 5, 4, 6]
                s.pos[6] = [8, 5, 5, 5, 7]
                s.pos[7] = [9, 5, 5, 6, 8]
                s.pos[8] = [10, 5, 5, 7, 0]
                
                # Setup lattice accordingly
                initLattice!(l, [s])
                
                # Prepare for bond removal
                s.moveIdx = 3  # Choose the detour bond (between positions 3 and 4)
                s.moveType = 1  # Bond removal
                s.moveVec = [0, 1, 0]  # Direction doesn't matter for removal
                s.newPosi = s.pos[2][1:3]  # Connect 2 directly to 5
                s.newPosj = s.pos[5][1:3]
                
                # Original length
                origLength = s.length
                
                # Execute move
                makeMove!(l, s)
                
                # Tests
                @test s.length == origLength - 2  # Length decreased by 2
                @test !haskey(s.pos, 3)  # Positions were removed
                @test !haskey(s.pos, 4)
                @test s.pos[2][5] == 5  # Check connections were updated
                @test s.pos[5][4] == 2
                @test l.lat[1, 6, 6, 5] == 0  # Lattice sites are now empty
                @test l.lat[1, 7, 6, 5] == 0
            end
        end
    end
end      
            
#             @testset "Kink Bond" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple straight polymer along x-axis
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice accordingly
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Prepare for kink move on first position
#                 s.moveIdx = 3  # Choose middle bond (between positions 3 and 4)
#                 s.moveType = 3  # Kink bond
#                 s.moveVec = [0, 1, 0]  # Move in y direction
#                 s.newPosi = s.pos[2][1:3]  # Keep first position same
#                 s.newPosj = [8, 6, 5]  # Move second position in y direction
                
#                 # Original positions
#                 origPos3 = copy(s.pos[3][1:3])
#                 origPos4 = copy(s.pos[4][1:3])
                
#                 # Execute move
#                 makeMove!(l, s)
                
#                 # Tests
#                 @test s.length == 6  # Length remains the same
#                 @test s.pos[4][1:3] == [8, 6, 5]  # Position was updated
#                 @test l.lat[1, origPos4[1], origPos4[2], origPos4[3]] == 0  # Old site is empty
#                 @test l.lat[1, s.pos[4][1], s.pos[4][2], s.pos[4][3]] == 1  # New site is filled
#             end
#         end
        
#         @testset "Edge Cases" begin
#             @testset "Add at Polymer Ends" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 4, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 0]
                
#                 # Setup lattice
#                 for i in 1:4
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Try adding at beginning
#                 s.moveIdx = 1
#                 s.moveType = 2
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = [5, 6, 5]
#                 s.newPosj = [6, 6, 5]
#                 makeMove!(l, s)
                
#                 @test s.length == 6
#                 @test haskey(s.pos, 5)
#                 @test haskey(s.pos, 6)
#                 @test s.pos[1][5] == 5
#             end
            
#             @testset "Remove and Add Cancellation" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Store original configuration
#                 origConfig = Dict{Int32, Vector{Int32}}()
#                 for (k, v) in s.pos
#                     origConfig[k] = copy(v)
#                 end
#                 origLength = s.length
                
#                 # Add a bond
#                 s.moveIdx = 3
#                 s.moveType = 2
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = [7, 6, 5]
#                 s.newPosj = [8, 6, 5]
#                 makeMove!(l, s)
                
#                 # Now remove the same bond
#                 s.moveIdx = 7  # The index of the first new bead
#                 s.moveType = 1
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[3][1:3]
#                 s.newPosj = s.pos[4][1:3]
#                 makeMove!(l, s)
                
#                 # Verify we're back to the original state
#                 @test s.length == origLength
#                 for (k, v) in origConfig
#                     @test haskey(s.pos, k)
#                     @test s.pos[k] == v
#                 end
#             end
            
#             @testset "Kink Cancellation" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Store original configuration
#                 origConfig = Dict{Int32, Vector{Int32}}()
#                 for (k, v) in s.pos
#                     origConfig[k] = copy(v)
#                 end
                
#                 # Apply first kink
#                 s.moveIdx = 3
#                 s.moveType = 3
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[3][1:3]
#                 s.newPosj = [8, 6, 5]
#                 makeMove!(l, s)
                
#                 # Apply reverse kink
#                 s.moveIdx = 3
#                 s.moveType = 3
#                 s.moveVec = [0, -1, 0]
#                 s.newPosi = s.pos[3][1:3]
#                 s.newPosj = [8, 5, 5]
#                 makeMove!(l, s)
                
#                 # Verify we're back to the original state
#                 for (k, v) in origConfig
#                     @test haskey(s.pos, k)
#                     @test s.pos[k] == v
#                 end
#             end
            
#             @testset "Invalid Move Vector" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Move parallel to bond (should fail)
#                 s.moveIdx = 3
#                 s.moveType = 3
#                 s.moveVec = [1, 0, 0]  # Parallel to bond direction
#                 s.bondVec = [1, 0, 0]  # Will be calculated in makeMove!, but setting for test
#                 s.newPosi = [8, 5, 5]
#                 s.newPosj = [9, 5, 5]
                
#                 @test_throws AssertionError makeMove!(l, s)
#             end
            
#             @testset "Move Type Inconsistent with Actual Move" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Setup kink move but claim it's removal
#                 s.moveIdx = 3
#                 s.moveType = 1  # Claim removal
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[3][1:3]  # But setup for kink
#                 s.newPosj = [8, 6, 5]
                
#                 @test_throws AssertionError makeMove!(l, s)
#             end
            
#             @testset "Boundary Conditions" begin
#                 l = Lattice(; boxlength = 8, nspecies = 1)  # Smaller box for boundary testing
#                 s = Species(; length = 4, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create polymer at boundary
#                 s.pos[1] = [7, 5, 5, 0, 2]
#                 s.pos[2] = [8, 5, 5, 1, 3]  # At boundary
#                 s.pos[3] = [1, 5, 5, 2, 4]  # Wrapped around
#                 s.pos[4] = [2, 5, 5, 3, 0]
                
#                 # Setup lattice
#                 for i in 1:4
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Bond vector should wrap properly
#                 s.moveIdx = 2
#                 s.moveType = 3
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[2][1:3]
#                 s.newPosj = [1, 6, 5]
                
#                 # Calculate bond vector manually for comparison
#                 bondVec = [(s.pos[3][1] - s.pos[2][1]) % s.boxlength, 
#                            (s.pos[3][2] - s.pos[2][2]) % s.boxlength, 
#                            (s.pos[3][3] - s.pos[2][3]) % s.boxlength]
                
#                 makeMove!(l, s)
                
#                 @test s.bondVec == bondVec
#                 @test s.pos[3][1:3] == [1, 6, 5]
#             end

#             @testset "Polymer Too Short for Removal" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 4, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer with minimal length
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 0]
                
#                 # Setup lattice
#                 for i in 1:4
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Attempt to remove a bond (would make polymer too short)
#                 s.moveIdx = 2
#                 s.moveType = 1
#                 s.moveVec = [0, 0, 1]
#                 s.newPosi = s.pos[1][1:3]
#                 s.newPosj = s.pos[4][1:3]
                
#                 # This should fail or be handled gracefully if your implementation has a minimum length check
#                 # Uncomment if you have such checks:
#                 # @test_throws AssertionError makeMove!(l, s)
                
#                 # If no checks, at least verify the operation happened correctly
#                 makeMove!(l, s)
#                 @test s.length == 2
#                 @test !haskey(s.pos, 2)
#                 @test !haskey(s.pos, 3)
#             end

#             @testset "Maximum Length Constraint" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 8, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a polymer at almost max length
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 7]
#                 s.pos[7] = [11, 5, 5, 6, 8]
#                 s.pos[8] = [12, 5, 5, 7, 0]
                
#                 # Setup lattice
#                 for i in 1:8
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # First addition should work (brings to max length)
#                 s.moveIdx = 3
#                 s.moveType = 2
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = [7, 6, 5]
#                 s.newPosj = [8, 6, 5]
#                 makeMove!(l, s)
                
#                 @test s.length == 10
#                 @test haskey(s.pos, 9)
#                 @test haskey(s.pos, 10)
                
#                 # Second addition should fail due to maxlength constraint
#                 s.moveIdx = 5
#                 s.moveType = 2
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = [9, 6, 5]
#                 s.newPosj = [10, 6, 5]
                
#                 # If you have maxlength checks, uncomment:
#                 # @test_throws AssertionError makeMove!(l, s)
#             end

#             @testset "Occupied Lattice Site Rejection" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 6, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 5]
#                 s.pos[5] = [9, 5, 5, 4, 6]
#                 s.pos[6] = [10, 5, 5, 5, 0]
                
#                 # Setup lattice
#                 for i in 1:6
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Manually occupy a site where we want to add
#                 l.lat[1, 7, 6, 5] = 1
                
#                 # Attempt to add bond to occupied site
#                 s.moveIdx = 3
#                 s.moveType = 2
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = [7, 6, 5]  # Already occupied
#                 s.newPosj = [8, 6, 5]
                
#                 # If you have occupation checks, uncomment:
#                 # @test_throws AssertionError makeMove!(l, s)
#             end

#             @testset "Kink at Polymer Ends" begin
#                 l = Lattice(; boxlength = 16, nspecies = 1)
#                 s = Species(; length = 4, maxlength = 10, boxlength = l.boxlength)
                
#                 # Create a simple polymer
#                 s.pos[1] = [5, 5, 5, 0, 2]
#                 s.pos[2] = [6, 5, 5, 1, 3]
#                 s.pos[3] = [7, 5, 5, 2, 4]
#                 s.pos[4] = [8, 5, 5, 3, 0]
                
#                 # Setup lattice
#                 for i in 1:4
#                     l.lat[1, s.pos[i][1], s.pos[i][2], s.pos[i][3]] = 1
#                 end
                
#                 # Kink at beginning
#                 s.moveIdx = 1
#                 s.moveType = 3
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[1][1:3]
#                 s.newPosj = [6, 6, 5]
#                 makeMove!(l, s)
                
#                 @test s.pos[2][1:3] == [6, 6, 5]
#                 @test l.lat[1, 6, 5, 5] == 0
#                 @test l.lat[1, 6, 6, 5] == 1
                
#                 # Kink at end
#                 s.moveIdx = 3
#                 s.moveType = 3
#                 s.moveVec = [0, 1, 0]
#                 s.newPosi = s.pos[3][1:3]
#                 s.newPosj = [8, 6, 5]
#                 makeMove!(l, s)
                
#                 @test s.pos[4][1:3] == [8, 6, 5]
#                 @test l.lat[1, 8, 5, 5] == 0
#                 @test l.lat[1, 8, 6, 5] == 1
#             end
#         end
#     end
# end
