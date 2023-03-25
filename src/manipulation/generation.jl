using Combinatorics
using Random
using RandomNumbers.MersenneTwisters
Random.seed!(1)

function PowerSet(L::Vector{Int})
    return collect(powerset(L))
end


function TotalPowerRelation(N::Int)
    L = collect(1:N)
    PS = PowerSet(L)
    PS = shuffle(PS)
    s = 0
    k = ceil(Int, 2^N * rand())
    Class = PS[s+1:s+k]
    PR = [Class]
    s = s + k
    while s < 2^N
        k = ceil(Int, (2^N - s) * rand())
        Class = PS[s+1:s+k]
        push!(PR, Class)
        s = s + k
    end
    return PR
end

