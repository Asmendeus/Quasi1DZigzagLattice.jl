# Diagonal Square Lattice:
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#     A   A   A   A   A
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#   A   A   A   A   A
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#     A   A   A   A   A
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#   A   A   A   A   A
# Diagonal square lattice is the square lattice rotating 45°

struct DiagonalSquareLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function DiagonalSquareLattice(L::Int, W::Int)
        @assert L ≥ W
        return new(L, W)
    end
end
const DiagSquaLatt = DiagonalSquareLattice

function Base.size(latt::DiagonalSquareLattice) = (latt.L, latt.W)

function getSite(latt::DiagonalSquareLattice, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return (lp-1) * latt.W + wp
end

function getAllSites(latt::DiagonalSquareLattice)
    return collect(1:latt.L*latt.W)
end

function getAllNNPairs(latt::DiagonalSquareLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
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

function getAllNNNPairs(latt::DiagonalSquareLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+2)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-2
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+2)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end
