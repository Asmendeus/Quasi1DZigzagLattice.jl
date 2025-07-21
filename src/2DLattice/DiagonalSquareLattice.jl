# Diagonal Square Lattice:
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#     A   A   A   A   A
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#   A   A   A   A   A
#    ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲
#     A   A   A   A   A
#    ╱ ╲ ╱ ╲ ╱ ╲ ╱ ╲ ╱
#   A   A   A   A   A
# Diagonal square lattice is the square lattice rotating 45°

struct DiagonalSquareLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function DiagonalSquareLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional diagonal square lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const DiagSquaLatt = DiagonalSquareLattice

Base.size(latt::DiagonalSquareLattice) = (latt.L, latt.W)

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
        latt.W > 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width diagonal square lattice with periodic boundary condition is singular in geometry!"
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
        latt.W > 4 || @warn "W ≤ 4 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width diagonal square lattice with periodic boundary condition is singular in geometry!"
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

function Base.show(io::IO, latt::DiagonalSquareLattice)

    function println_mainline(w::Int)
        print(repeat(" ", Int((maxlen+1)/2)*(isodd(w) ? 0 : 1)))
        for l in 1:latt.L
            site = string(getSite(latt, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = floor(Int, (maxlen-len)/2) + 1
            print(repeat(" ", len_l) * site * repeat(" ", len_r))
        end
        println()
    end

    maxlen = length(string(latt.L*latt.W)) + 2
    iseven(maxlen) && (maxlen += 1)
    halflen = Int((maxlen+1)/2)

    println(io, "$(latt.L) × $(latt.W) DiagonalSquareLattice:")
    for w = latt.W:-1:1
        if isodd(w)
            print(repeat(" ", halflen + (isodd(halflen) ? 1 : 0)) * "╱")
            println(repeat(repeat(" ", halflen-1) * "╲" * repeat(" ", halflen-1) * "╱", latt.L-1))
        else
            print(repeat(" ", halflen + (isodd(halflen) ? 1 : 0)) * "╲")
            println(repeat(repeat(" ", halflen-1) * "╱" * repeat(" ", halflen-1) * "╲", latt.L-1))
        end
        println_mainline(w)
    end
    return nothing
end
