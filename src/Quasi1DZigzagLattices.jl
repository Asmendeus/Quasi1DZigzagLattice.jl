module Quasi1DZigzagLattices

export AbstractLattice
include("AbstractLattice.jl")

export SquareLattice, SquaLatt, getSite, getAllSites, getAllNNPairs, getAllNNNPairs
include("2DLattice/SquareLattice.jl")
export DiagonalSquareLattice, DiagSquaLatt
include("2DLattice/DiagonalSquareLattice.jl")
export TriangularLattice, TriaLatt
include("2DLattice/TriangularLattice.jl")
export DiagonalTriangularLattice, DiagTriaLatt
include("2DLattice/DiagonalTriangularLattice.jl")
export HexagonalLattice, HoneycombLattice, HexaLatt
include("2DLattice/HexagonalLattice.jl")
export DiagonalHexagonalLattice, DiagHexaLatt
include("2DLattice/DiagonalHexagonalLattice.jl")
export KagomeLattice, KagoLatt
include("2DLattice/KagomeLattice.jl")
export DiagonalKagomeLattice, DiagKagoLatt
include("2DLattice/DiagonalKagomeLattice.jl")

end # module Quasi1DZigzagLattice