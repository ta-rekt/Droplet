# Overview of new heterogenous polymer feature

I would like to make polymers heterogenous, with multiple types of monomers with possibly different interaction potentials with other species in the lattice. 

I will describe to you how I'd like to go about this, and I want to iterate with you on this document to turn it into a detailed proposal of changes to the code and tests. The outcome of our iterations should be a detailed implementation plan that an agent can follow to implement and test all these changes. I want you to ask me questions if something's not clear to you, and I will help you narrow it down until everything is crystal clear before we finalize the plan.

## Sketch of implementation plan

Here is a sketch of how I'd like the heterogenous polymer to be implemented:

### Monomer types stored in `pos`
- monomer types should be stored in the Polymer.pos dictionary, no longer in Polymer.monomerIndex vector. for this, the entries of `pos` should become 6-element vectors instead of 5-element vectors, where the 6th element is the monomer type. all necessary changes in the code should be made to reflect this, and should be tested before moving on to the next change. monomerIndex structure should be removed from the polymer type.

### Base types and pattern types of monomers
- monomers in a polymer can have a base type, say 1, and additional 'pattern' types (2+). base types and pattern types are inserted into a pattern given at initialization, much like the initMonomers! function does already. here are example patterns that I will typically use to initialize polymers: 
    uniform: all 1s
    promoter: 1111112222111111
    promoter-enhancer: 1111112222111111333333111111
take note that the initialization mechanism should be changed to work with the p.pos structure. you cannot loop over it, you must follow the beads by accessing the 5th element of p.pos[key], which represents the key of the next bead, just like is done in ham! in the hamiltonian.jl file, line 38.


### Adding/removing monomers
- base types and pattern types behave differently with regard to adding and removing bonds. I would like the length of patterns (any type of pattern, that is, any sequence of monomers of type 2 and above) to be fixed, and the region of base type 1 to be variable in length. when adding/removing a bond into a region of the polymer of type 1, everything behaves as with the homogeneous polymer. However, when adding/removing a bond inside a pattern, the last two beads at the end of the pattern must change type (x -> 1 for addition, 1 -> x for removal) in order to keep the total number of beads of type x fixed. this can be done inside the add and remove functions in move.jl.


### Accessing pattern boundaries
- the keys of beads at the end of a pattern must be stored in a dictionary in the Polymer type, called patternTails. that dictionary maps patternType -> endKey. this makes it easy to access the keys of beads whose type should be switched at bond addition/removal inside a pattern. these keys should be updated after add/remove, and the update can be done inside the add/remove functions in move.jl. you should add a check somewhere at initialization that ensures there is only one pattern of any given type, otherwise this approach breaks.


### Reject move if adjacent beads of different type
- if a bead at a certain index is chosen for a move and the next bead has different type, the move should be rejected. this is to avoid dealing with addition/removal of bonds that straddle the boundaries of a pattern

### Account for energy change of type swaps
- if an add/remove happens inside a pattern and two beads change type, this must be reflected in the deltaHam! function in hamiltonian.jl. you can detect if an add/remove is done inside a pattern simply by checking the type of the bead currently being added/removed. you should account for the energy change of this bead swap for both self-interactions and partner interactions, analogous to how it's currently done in hamiltonian.jl.

# Objectives

## Understand the code, draft implementation plan, ask questions

You can start by taking a look at the code, understanding the key elements for our purposes. Then write a draft of an implementation plan in implementation_plan.md in the same directory as this file, and ask me questions in the chat. After each round of feedback I give you, you'll refine the implementation_plan.md.

## Building and refining context

If you learn important context about the code that might be useful for future agents and humans in carrying out feature implementations, tests, or general understanding, don't hesitate to add files in the context directory. For now, don't modify any code. Just add context and edit the implementation_plan.md file.

### Guidelines for context

- please be frugal with context, keep context bloat to a minimum. only include what is absolutely necessary to understand a particular implementation detail, some background information about code, styling guidelines, common errors, common patterns, domain-specific knowledge, tips for future optimization or generalization.

- keep things human-readable. don't overuse bullet-points, try to strike a balance between descriptive paragraphs and bullet lists.

- always include pointers and references to files and code so agents can navigate easily

- keep an index file to navigate context effectively. keep it up to date as necessary

- at the end of a feedback run, make a pass over the context files to make simplifications, rewrite things you think can be explained more clearly, and keep things up to date.