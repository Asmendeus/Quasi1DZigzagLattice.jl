# Diagonal Kagome lattice:
#        ╲ ╱     ╲ ╱     ╲ ╱     ╲
#         B       B       B       B
#        ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#       A — C — A — C — A — C — A — C
#      ╱     ╲ ╱     ╲ ╱     ╲ ╱
#     B       B       B       B
#    ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#   A — C — A — C — A — C — A — C
#        ╲ ╱     ╲ ╱     ╲ ╱     ╲
#         B       B       B       B
#        ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#       A — C — A — C — A — C — A — C
#      ╱     ╲ ╱     ╲ ╱     ╲ ╱
#     B       B       B       B
#    ╱ ╲     ╱ ╲     ╱ ╲     ╱ ╲
#   A — C — A — C — A — C — A — C

struct DiagonalKagomeLattice <: AbstractLattice{2}
    L::Int
    W::Int
    function DiagonalKagomeLattice(L::Int, W::Int)
        L ≥ W || throw(ArgumentError("Quasi-one-dimensional Kagome lattice requires L ≥ W!"))
        return new(L, W)
    end
end
const DiagKagoLatt = DiagonalKagomeLattice

Base.size(latt::DiagonalKagomeLattice) = (latt.L, latt.W)

function getSite(latt::DiagonalKagomeLattice, ci::Int, l::Int, w::Int)
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

function getAllSites(latt::DiagonalKagomeLattice)
    return collect(1:3*latt.L*latt.W)
end

function getAllNNPairs(latt::DiagonalKagomeLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        latt.W > 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width diagonal Kagome lattice with periodic boundary condition is singular in geometry!"
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w)))
        end
        for w in 1:latt.W
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for l in 1:latt.L, w in 1:latt.W
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l, w)))
            push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w)))
            push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w)))
        end
        for l in 1:latt.L-1, w in 1:latt.W
            push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w)))
        end
        for w in 1:latt.W-1
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l-1, w+1)))
                end
            else
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l+1, w+1)))
                end
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function getAllNNNPairs(latt::DiagonalKagomeLattice; boundary::Symbol=:PBC)
    pairs = Tuple{Int, Int}[]
    if boundary == :PBC
        latt.W > 2 || @warn "W ≤ 2 may lead to singular behavior with periodic boundary condition!"
        iseven(latt.W) || @warn "Odd-width diagonal Kagome lattice with periodic boundary condition is singular in geometry!"
        for w in 1:latt.W
            for l in 2:latt.L, w in 1:latt.W
                push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l-1, w)))
            end
            for l in 1:latt.L-1, w in 1:latt.W
                push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w)))
            end
            if isodd(w)
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
                end
                for l in 2:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l-1, w+1)))
                end
            else
                for l in 1:latt.L
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                    push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w+1)))
                end
                for l in 1:latt.L-1
                    push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w+1)))
                    push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l+1, w+1)))
                end
            end
        end
    elseif boundary == :OBC
        for w in 1:latt.W
            for l in 2:latt.L, w in 1:latt.W
                push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 2, l-1, w)))
            end
            for l in 1:latt.L-1, w in 1:latt.W
                push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 2, l+1, w)))
            end
            for w in 1:latt.W-1
                if isodd(w)
                    for l in 1:latt.L
                        push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l, w+1)))
                        push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l, w+1)))
                    end
                    for l in 2:latt.L
                        push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l-1, w+1)))
                        push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l-1, w+1)))
                    end
                else
                    for l in 1:latt.L
                        push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 1, l, w+1)))
                        push!(pairs, (getSite(latt, 1, l, w), getSite(latt, 3, l, w+1)))
                    end
                    for l in 1:latt.L-1
                        push!(pairs, (getSite(latt, 3, l, w), getSite(latt, 1, l+1, w+1)))
                        push!(pairs, (getSite(latt, 2, l, w), getSite(latt, 3, l+1, w+1)))
                    end
                end
            end
        end
    else
        throw(ArgumentError("Undefined behavior of boundary condition `:$(boundary)`"))
    end
    return pairs
end

function Base.show(io::IO, latt::DiagonalKagomeLattice)

    function println_mainline_up(w::Int)
        print(repeat(" ", Int((maxlen+1)*(mod(w-1, 2)+1/2))))
        for l in 1:latt.L
            site = string(getSite(latt, 2, l, w))
            len = length(site)
            len_l = ceil(Int, (maxlen-len)/2)
            len_r = 2*(maxlen+1) - len_l - length(site)
            print(repeat(" ", len_l) * site * repeat(" ", len_r))
        end
        println()
    end
    function println_mainline_dn(w::Int)
        print(repeat(" ", (maxlen+1)*(mod(w-1, 2))))
        for l in 1:latt.L
            site1 = string(getSite(latt, 1, l, w))
            len1 = length(site1)
            len1_l = floor(Int, (maxlen-len1)/2)
            len1_r = ceil(Int, (maxlen-len1)/2)
            site2 = string(getSite(latt, 3, l, w))
            len2 = length(site1)
            len2_l = floor(Int, (maxlen-len2)/2)
            len2_r = ceil(Int, (maxlen-len2)/2)
            print(repeat(" ", len1_l) * site1 * repeat(" ", len1_r) * "—")
            print(repeat(" ", len2_l) * site2 * repeat(" ", len2_r) * (l == latt.L ? "" : "—"))
        end
        println()
    end

    maxlen = length(string(3*latt.L*latt.W)) + 2
    iseven(maxlen) && (maxlen += 1)
    halflen = Int((maxlen+1)/2)

    println(io, "$(latt.L) × $(latt.W) DiagonalKagomeLattice:")
    for w = latt.W:-1:1
        if isodd(w)
            print(repeat(" ", (maxlen+1)*(mod(w-1, 2)) + 2*halflen + (isodd(halflen) ? 1 : 0)) * "╱")
            println(repeat(repeat(" ", 3*halflen-1) * "╲" * repeat(" ", halflen-1) * "╱", latt.L-1))
        else
            print(repeat(" ", (maxlen+1)*(mod(w-1, 2)) + halflen + (isodd(halflen) ? 1 : 0)) * "╲")
            println(repeat(repeat(" ", halflen-1) * "╱" * repeat(" ", maxlen+halflen) * "╲", latt.L-1))
        end
        println_mainline_up(w)
        print(repeat(" ", (maxlen+1)*(mod(w-1, 2)) + halflen + (isodd(halflen) ? 1 : 0)))
        println(repeat("╱" * repeat(" ", halflen-1) * "╲" * repeat(" ", maxlen+halflen), latt.L))
        println_mainline_dn(w)
    end
    return nothing
end
