# Performance tips

The majority of the computation involved in running SMC is the computation of the mutation and log potential functions. The goal in the design of this package has been to allow users to specify these functions in ways that avoid unnecessary computations and memory allocations.

For example, user-defined scratch space structs are intended to allow users to avoid any memory allocations entirely. Running the algorithm in serial is possible with no allocations whatsoever. Example scratch space structs are defined in the multivariate linear Gaussian model and Lorenz96 models in [SMCExamples](https://github.com/awllee/SMCExamples.jl).

It is sometimes possible to store important per-particle information as part of the data associated with a particle. This is done, e.g., in the SMC sampler demo in [SMCExamples](https://github.com/awllee/SMCExamples.jl) to store two different log densities of the state, which can therefore be computed only once.

In general, use of [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) can provide dramatic improvements for multivariate particles, and is highly recommended for fixed-size vector components of particles.
