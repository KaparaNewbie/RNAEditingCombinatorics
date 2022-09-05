# basic slack versiob

# pmaps (didn't actually used these llinks in here)
# https://discourse.julialang.org/t/pmap-in-function-module/25896
# https://stackoverflow.com/questions/44741667/julia-equivalent-of-python-multiprocessing-pool-map




# JULIA_PROJECT=. julia

# JULIA_PROJECT=. julia -t 20



# julia -t 16 
using Pkg; Pkg.activate(".")

using Distributed
using BenchmarkTools
n_workers = 2
threads_per_worker = Int(Threads.nthreads() / n_workers)

# addprocs(workers, exeflags=["--threads=$threads_per_worker", "--project=$(Base.active_project())"])
addprocs(n_workers, exeflags=["--threads=$threads_per_worker", "--project"])

# addprocs(workers; exeflags="--project")
# addprocs(workers; exeflags="--project", env=["JULIA_NUM_THREADS" => string(threads_per_worker)])

@everywhere using Transducers

# 1 + 2: dcollect + threads

function f1(n, x)
    A = rand(n, n)
    sum(dcollect(f2(A, rand()) for _ in 1:nworkers() * x))
end

@everywhere function f2(A, x)
    B = zeros(size(A, 2))
    Threads.@threads for i in 1:size(A, 2)
        B[i] = sum(@views A[:, i])
    end
    sum(B .+ x)
end

# 3 + 4: dcollect

function f3(n, x)
    A = rand(n, n)
    sum(dcollect(f4(A, rand()) for _ in 1:nworkers() * x))
end

@everywhere function f4(A, x)
    B = zeros(size(A, 2))
    for i in 1:size(A, 2)
        B[i] = sum(@views A[:, i])
    end
    sum(B .+ x)
end


# 5 + 6: pmap + threads

function f5(n, x)
    A = rand(n, n)
    sum(pmap(x -> f6(A, x), [rand() for _ in 1:nworkers() * x]))
end

@everywhere function f6(A, x)
    B = zeros(size(A, 2))
    Threads.@threads for i in 1:size(A, 2)
        B[i] = sum(@views A[:, i])
    end
    sum(B .+ x)
end


# 7 + 8: pmap

function f7(n, x)
    A = rand(n, n)
    sum(pmap(x -> f8(A, x), [rand() for _ in 1:nworkers() * x]))
end

@everywhere function f8(A, x)
    B = zeros(size(A, 2))
    for i in 1:size(A, 2)
        B[i] = sum(@views A[:, i])
    end
    sum(B .+ x)
end


n = 10_000
x = 10

f1(n, x)
f3(n, x)
f5(n, x)
f7(n, x)

@benchmark f1($n, $x)  # dcollect + threads
@benchmark f3($n, $x)  # dcollect
@benchmark f5($n, $x)  # pmap + threads
@benchmark f7($n, $x)  # pmap


@everywhere using Distributed

@everywhere B() = println("Distributed.myid = $(Distributed.myid()), Threads.nthreads = $(Threads.nthreads())")
pmap(x -> B(), [worker for worker in workers()])


function f1(n)  # use all nthreads()
    G = ... # build a graph G
    tcollect(f2(G) for _ in 1:nworkers())
end


function f2(G) # use nthreads() / nworkers() threads
    # do some random work on G
end