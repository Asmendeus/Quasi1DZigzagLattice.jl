# Hexagonal lattice:
#                 ╱   ╲   ╱   ╲   ╱   ╲   ╱
#               B       B       B       B
#               |       |       |       |
#               A       A       A       A
#             ╱   ╲   ╱   ╲   ╱   ╲   ╱
#           B       B       B       B
#           |       |       |       |
#           A       A       A       A
#         ╱   ╲   ╱   ╲   ╱   ╲   ╱
#       B       B       B       B
#       |       |       |       |
#       A       A       A       A
#     ╱   ╲   ╱   ╲   ╱   ╲   ╱
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
        latt.W ≤ 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width hexagonal lattice with periodic boundary condition is singular in geometry!"
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
        latt.W ≤ 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width hexagonal lattice with periodic boundary condition is singular in geometry!"
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

function Base.show(io::IO, latt::HexagonalLattice)

    function println_mainline_up(w::Int)
        print(repeat(" ", Int((maxlen+1)/2)*(w-1)))
        for l in 1:latt.L
            site = string(getSite(latt, 2, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = floor(Int, (maxlen-len)/2) + 1
            print(repeat(" ", len_l) * site * repeat(" ", len_r))
        end
        println()
    end
    function println_mainline_dn(w::Int)
        print(repeat(" ", Int((maxlen+1)/2)*(w-1)))
        for l in 1:latt.L
            site = string(getSite(latt, 1, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = floor(Int, (maxlen-len)/2) + 1
            print(repeat(" ", len_l) * site * repeat(" ", len_r))
        end
        println()
    end

    maxlen = length(string(2*latt.L*latt.W)) + 2
    iseven(maxlen) && (maxlen += 1)
    halflen = Int((maxlen+1)/2)

    println(io, "$(latt.L) × $(latt.W) HexagonalLattice:")
    for w = latt.W:-1:1
        print(repeat(" ", halflen*w + (isodd(halflen) ? 1 : 0)) * "╱")
        println(repeat(repeat(" ", halflen-1) * "╲" * repeat(" ", halflen-1) * "╱", latt.L-1))
        println_mainline_up(w)
        print(repeat(" ", 1 + halflen*(w-1) + (isodd(halflen) ? 1 : 0)) * "|")
        println(repeat(repeat(" ", maxlen) * "|", latt.L-1))
        println_mainline_dn(w)
    end
    return nothing
end
