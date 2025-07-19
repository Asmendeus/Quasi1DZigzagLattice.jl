# Hexagonal lattice:
#         ╲   ╱   ╲   ╱   ╲   ╱   ╲   ╱
#           B       B       B       B
#           |       |       |       |
#           A       A       A       A
#         ╱   ╲   ╱   ╲   ╱   ╲   ╱
#       B       B       B       B
#       |       |       |       |
#       A       A       A       A
#     ╱   ╲   ╱   ╲   ╱   ╲   ╱   ╲
#   B       B       B       B
#   |       |       |       |
#   A       A       A       A

struct HexagonalLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function HexagonalLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional hexagonal lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const HoneycombLattice = HexagonalLattice
const HexaLatt = HexagonalLattice

Base.size(latt::HexagonalLattice) = (latt.L, latt.W)

function getSite(latt::HexagonalLattice, ci::Int, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return ci + 2 * ((lp-1) * latt.W + (wp-1))
end

function getAllSites(latt::HexagonalLattice)
    return collect(1:2*latt.L*latt.W)
end

function getAllNNPairs(latt::HexagonalLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::HexagonalLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l-1, w+1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l-1, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l-1, w+1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l-1, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end
