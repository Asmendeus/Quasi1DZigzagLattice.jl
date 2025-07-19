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
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional triangular lattice requires L ≥ W!"))
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

function Base.show(io::IO, latt::TriangularLattice)

    function println_mainline(w::Int)
        print(repeat(" ", Int((maxlen+1)/2)*(w-1)))
        for l in 1:latt.L
            site = string(getSite(latt, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = floor(Int, (maxlen-len)/2)
            print(repeat(" ", len_l) * site * repeat(" ", len_r) * "—")
        end
        println()
    end

    maxlen = length(string(latt.L*latt.W)) + 2
    iseven(maxlen) && (maxlen += 1)
    halflen = Int((maxlen+1)/2)

    println(io, "$(latt.L) × $(latt.W) TriangularLattice:")
    for w = latt.W:-1:1
        print(repeat(" ", halflen*w + (isodd(halflen) ? 1 : 0)) * "╱")
        println(repeat(repeat(" ", halflen-1) * "╲" * repeat(" ", halflen-1) * "╱", latt.L-1))
        println_mainline(w)
    end
    return nothing
end
