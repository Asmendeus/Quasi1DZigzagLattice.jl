# Square Lattice:
#   |    |    |    |
#   A —— A —— A —— A
#   |    |    |    |
#   A —— A —— A —— A
#   |    |    |    |
#   A —— A —— A —— A
#   |    |    |    |
#   A —— A —— A —— A

struct SquareLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function SquareLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional square lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const SquaLatt = SquareLattice

Base.size(latt::SquareLattice) = (latt.L, latt.W)

function getSite(latt::SquareLattice, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return (lp-1) * latt.W + wp
end

function getAllSites(latt::SquareLattice)
    return collect(1:latt.L*latt.W)
end

function getAllNNPairs(latt::SquareLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        latt.W > 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width square lattice with periodic boundary condition is singular in geometry!"
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w)))
        end
        for l in 1:latt.L, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l, w+1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::SquareLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        latt.W > 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width square lattice with periodic boundary condition is singular in geometry!"
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w-1)))
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W-1
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w+1)))
        end
        for l in 1:latt.L-1, w in 2:latt.W
            push!(pairs, (getSite(latt, l, w), getSite(latt, l+1, w-1)))
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function Base.show(io::IO, latt::SquareLattice)

    function println_mainline(w::Int)
        for l in 1:latt.L
            site = string(getSite(latt, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = floor(Int, (maxlen-len)/2)
            print(repeat(" ", len_l) * site * repeat(" ", len_r) * (l == latt.L ? "" : "—"))
        end
        println()
    end

    maxlen = length(string(latt.L*latt.W)) + 2
    iseven(maxlen) && (maxlen += 1)
    halflen = Int((maxlen-1)/2)

    println(io, "$(latt.L) × $(latt.W) SquareLattice:")
    for w = latt.W:-1:1
        println(io, repeat(repeat(" ", halflen) * "|" * repeat(" ", halflen + 1), latt.L))
        println_mainline(w)
    end
    return nothing
end
