# Kagome lattice:
#              ╱     ╲ ╱     ╲ ╱     ╲ ╱
#             B       B       B       B
#            ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#           A — C — A — C — A — C — A — C —
#          ╱     ╲ ╱     ╲ ╱     ╲ ╱
#         B       B       B       B
#        ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#       A — C — A — C — A — C — A — C —
#      ╱     ╲ ╱     ╲ ╱     ╲ ╱
#     B       B       B       B
#    ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#   A — C — A — C — A — C — A — C —

struct KagomeLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function KagomeLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional Kagome lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const KagoLatt = KagomeLattice

Base.size(latt::KagomeLattice) = (latt.L, latt.W)

function getSite(latt::KagomeLattice, ci::Int, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    if ci == 1
        return 3 * (lp-1) * latt.W + 2 * wp - 1
    elseif ci == 2
        return 3 * (lp-1) * latt.W + 2 * wp
    elseif ci == 3
        return 3 * (lp-1) * latt.W + 2 * latt.W + wp
    else
        throw(ArgumentError("The site number `ci` in the cells should be in (1, 2, 3)"))
    end
end

function getAllSites(latt::KagomeLattice)
    return collect(1:3*latt.L*latt.W)
end

function getAllNNPairs(latt::KagomeLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w-1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w)))
        end
        for l in 1:latt.L-1, w in 2:latt.W
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w-1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::KagomeLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l+1, w-1)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w)))
        end
        for l in 2:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l-1, w+1)))
        end
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 2:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l+1, w-1)))
        end
        for l in 2:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l-1, w+1)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end
