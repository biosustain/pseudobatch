# This file runs all julia scripts. By running them from the same file they are run in the same julia REPL which saves loading time.
using Pkg
print("Running Julia scripts \n")
print("Activating the library. This can take several minutes \n")
Pkg.activate("julia-env")
Pkg.instantiate()

print("Simulating standard fed-batch process \n")
include("standard_fed-batch_process.jl")

print("Simulating fed-batch with product inhibition \n")
include("fed-batch_with_product_inhibition.jl")

print("Simulating fed-batch process with multiple feeds\n")
include("multiple_step_feed_process.jl")

print("Simulating fed-batch process with volatile product")
include("fed-batch_volatile_product.jl")