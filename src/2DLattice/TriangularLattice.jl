# Triangular lattice:
#          ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#         A — A — A — A — A — A —
#        ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#       A — A — A — A — A — A —
#      ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#     A — A — A — A — A — A —
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#   A — A — A — A — A — A —

struct TriangularLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function TriangularLattice(L::Int, W::Int)
        @assert L ≥ W
        return new(L, W)
    end
end
const TriaLatt = TriangularLattice

Base.size(latt::TriangularLattice) = (latt.L, latt.W)

function getSite(latt::TriangularLattice, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return (lp-1) * latt.W + wp
end

function getAllSites(latt::TriangularLattice)
    return collect(1:latt.L*latt.W)
end

function getAllNNPairs(latt::TriangularLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 2:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 2:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::TriangularLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+2)))
        end
        for l in 3:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-2, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
        end
        for l in 2:latt.L, w in 1:latt.W-2
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+2)))
        end
        for l in 3:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l-2, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end
