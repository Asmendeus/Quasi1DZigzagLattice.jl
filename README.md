# Quasi1DZigzagLattices

A package for MPS/MPO algorithms for generating Zigzag path on a lattice

## Installation

```
julia> ]
pkg> add https://github.com/Asmendeus/Quasi1DZigzagLattices.jl.git#main
```

## Guide

We first generate a two-dimensional lattice (as an example of a square lattice)

```julia
latt = SquaLatt(4, 4)
```

Get all its sites:

```julia
getAllSites(latt)
```

Get all its nearest neighbor pairs with periodic boundary condition:

```julia
getAllNNPairs(latt)
```

Get all its next nearest neighbor pairs with open boundary condition:

```julia
getAllNNNPairs(latt; boundary=:OBC)
```

Currently supported lattice types:

```
SquareLattice (SquaLatt)
TriangularLattice (TriaLatt)
HexagonalLattice (HexaLatt)
KagomeLattice (KagoLatt)
DiagonalSquareLattice (DiagSquaLatt)
DiagonalTriangularLattice (DiagTriaLatt)
DiagonalHexagonalLattice (DiagHexaLatt)
DiagonalKagomeLattice (DiagKagoLatt)
```
