"""
    abstract type AbstractLattice{N} end

Abstract type of all lattice.
`N::Int` corresponds to the spatial dimension of the lattice, and N-dimensional space corresponds to N linearly independent basis vectors.

# Note
Any `Lattice` struct inherited from `AbstractLattice` follows the conventions:
    1. Has fields `L::Int`, `W::Int`, ..., indicating the length of a dimension.
       Here, `L::Int` corresponds to the main dirrection of the quasi-1D lattice, with open boundary condition(:OBC).
    2. Lattice site sorting along the Zigzag path.
    3. Some methods for obtaining all local sites or pairs have been provided.
"""
abstract type AbstractLattice{N} end