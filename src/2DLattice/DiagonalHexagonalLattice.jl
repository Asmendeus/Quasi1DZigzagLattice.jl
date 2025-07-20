# Diagonal Hexagonal lattice:
#     ╱   ╲   ╱   ╲   ╱   ╲   ╱
#   B       B       B       B
#   |       |       |       |
#   A       A       A       A
#     ╲   ╱   ╲   ╱   ╲   ╱   ╲
#       B       B       B       B
#       |       |       |       |
#       A       A       A       A
#     ╱   ╲   ╱   ╲   ╱   ╲   ╱
#   B       B       B       B
#   |       |       |       |
#   A       A       A       A

struct DiagonalHexagonalLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function DiagonalHexagonalLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional diagonal hexagonal lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const DiagHexaLatt = DiagonalHexagonalLattice

Base.size(latt::DiagonalHexagonalLattice) = (latt.L, latt.W)

function getSite(latt::DiagonalHexagonalLattice, ci::Int, l::Int, w::Int)
    lp = mod(l-1, latt.L) + 1
    wp = mod(w-1, latt.W) + 1
    return ci + 2 * ((lp-1) * latt.W + (wp-1))
end

function getAllSites(latt::DiagonalHexagonalLattice)
    return collect(1:2*latt.L*latt.W)
end

function getAllNNPairs(latt::DiagonalHexagonalLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L, w in 1:latt.W
            push!(paris, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
        end
        for w in 1:latt.W
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for l in 1:latt.L, w in 1:latt.W
            push!(paris, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
        end
        for w in 1:latt.W-1
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::DiagonalHexagonalLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w)))
        end
        for w in 1:latt.W
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l-1, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w)))
        end
        for w in 1:latt.W-1
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l-1, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l+1, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 2, l, w+1)))
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function Base.show(io::IO, latt::DiagonalHexagonalLattice)

    function println_mainline_up(w::Int)
        print(repeat(" ", Int((maxlen+1)/2)*(isodd(w) ? 0 : 1)))
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
        print(repeat(" ", Int((maxlen+1)/2)*(isodd(w) ? 0 : 1)))
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

    println(io, "$(latt.L) × $(latt.W) DiagonalHexagonalLattice:")
    for w = latt.W:-1:1
        if isodd(w)
            print(repeat(" ", halflen + (isodd(halflen) ? 1 : 0)) * "╱")
            println(repeat(repeat(" ", halflen-1) * "╲" * repeat(" ", halflen-1) * "╱", latt.L-1))
        else
            print(repeat(" ", halflen + (isodd(halflen) ? 1 : 0)) * "╲")
            println(repeat(repeat(" ", halflen-1) * "╱" * repeat(" ", halflen-1) * "╲", latt.L-1))
        end
        println_mainline_up(w)
        print(repeat(" ", 1 + halflen*(isodd(w) ? 0 : 1) + (isodd(halflen) ? 1 : 0)) * "|")
        println(repeat(repeat(" ", maxlen) * "|", latt.L-1))
        println_mainline_dn(w)
    end
    return nothing
end
