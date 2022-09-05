using Distributed
using BenchmarkTools

addprocs(4, exeflags=["--project"])

@everywhere begin
    
    using Distributed, DistributedArrays

    """Sum the sum of `H`'s values."""
    f(H) = sum(sum.(values(H)))

    """Sum `x` runs of `f(H)`."""
    function f(H, x)
        sum(f(H) for _ ∈ 1:x) 
    end

    """Sum `x` runs of `f(ArrH[1])`."""
    function f(ArrH::DArray, x)
        H = ArrH[1]
        sum(f(H) for _ ∈ 1:x) 
    end

end

n = 10_000

G = Dict([x => collect(1:x) for x ∈ 1:n])

expected = sum(workers() * f(G))
got1 = sum(
    pmap(workers()) do worker
        f(G, worker)
    end
)
ArrG = @DArray [G for _ ∈ 1:1];
got2 = sum(
    pmap(workers()) do worker
        f(ArrG, worker)
    end
)
@assert expected == got1 == got2


function t1()
    n = 10_000
    G = Dict([x => collect(1:x) for x ∈ 1:n])
    sum(
        pmap(workers()) do worker
            f(G, worker)
        end
    )
end


function t2()
    n = 10_000
    G = Dict([x => collect(1:x) for x ∈ 1:n])
    ArrG = @DArray [G for _ ∈ 1:1]
    sum(
        pmap(workers()) do worker
            f(ArrG, worker)
        end
    )
end



@benchmark t1()

@benchmark t2()







