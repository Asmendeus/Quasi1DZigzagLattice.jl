# Diagonal Triangular lattice:
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#     A — A — A — A — A — A —
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#   A — A — A — A — A — A —
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#     A — A — A — A — A — A —
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#   A — A — A — A — A — A —

struct DiagonalTriangularLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function DiagonalTriangularLattice(L::Int, W::Int)
        @assert L ≥ W
        return new(L, W)
    end
end
const DiagTriaLatt = DiagonalTriangularLattice

Base.size(latt::DiagonalTriangularLattice) = (latt.L, latt.W)

function getSite(latt::DiagonalTriangularLattice, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return (lp-1) * latt.W + wp
end

function getAllSites(latt::DiagonalTriangularLattice)
    return collect(1:latt.L*latt.W)
end

function getAllNNPairs(latt::DiagonalTriangularLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for w in 1:latt.W
            if isodd(w)
                for l in 1:latt.L
                    push!(paris, (getSite(latt, l, w), getSite(latt, l, w+1)))
                end
                for l in 2:latt.L
                    push!(paris, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for w in 1:latt.W-1
            if isodd(w)
                for l in 1:latt.L
                    push!(paris, (getSite(latt, l, w), getSite(latt, l, w+1)))
                end
                for l in 2:latt.L
                    push!(paris, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::DiagonalTriangularLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in latt.L, w in latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+2)))
        end
        for w in 1:latt.W
            if isodd(w)
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
                end
                for l in 3:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l-2, w+1)))
                end
            else
                for l in 1:latt.L-2
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+2, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for l in latt.L, w in latt.W-2
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+2)))
        end
        for w in 1:latt.W-1
            if isodd(w)
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
                end
                for l in 3:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l-2, w+1)))
                end
            else
                for l in 1:latt.L-2
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l+2, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, l, w), getSite(latt, l-1, w+1)))
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end
