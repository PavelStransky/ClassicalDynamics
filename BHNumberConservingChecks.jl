# BHNumberConservingChecks.jl
#
# The consistency checks of the Liouvillian code in one run, before any production calculation of
# number-conserving-BH-paper-TODO.md is trusted:
#
#   * Checks() of BHNumberConservingLiouvillian.jl - symmetries, Lindblad structure, the Z_L sector
#     decomposition, and (new) the reflection-parity split at η = 0 and SectorDimensions;
#   * SparseChecks() of BHNumberConservingLiouvillianSparse.jl - the orbit-basis sectors against the
#     momentum-basis ones, eigenvalue by eigenvalue;
#   * ModesChecks() of BHNumberConservingLiouvillianModes.jl - the conversion between sector vectors
#     and operators, and every right eigenvector as an eigen-operator of the full Liouvillian.
#
# The slicing validation against dense spectra takes minutes and is run separately:
#   julia BHNumberConservingLiouvillianSparse.jl validate 12 --w 16
#
# Usage:  julia BHNumberConservingChecks.jl

ENV["CD_NO_PLOTS"] = "true"
include(joinpath(@__DIR__, "BHNumberConservingLiouvillianSparse.jl"))
include(joinpath(@__DIR__, "BHNumberConservingLiouvillianModes.jl"))

println("BHNumberConservingLiouvillian.jl")
Checks()
println("\nBHNumberConservingLiouvillianSparse.jl")
SparseChecks()
println("\nBHNumberConservingLiouvillianModes.jl")
ModesChecks()

println("\nK13: actual sector dimensions (d_N^2 / L is not an integer in general)")
for L in (3, 4), N in (6, 8, 10, 12, 14, 16)
    r = SectorDimensions(L, N)
    println("  L = $L, N = $N: d_N = $(r.dN), d_N^2/L = $(round(r.dN^2 / L, digits = 2)), sectors m = 0 ... $(L - 1): $(r.sectors)")
end
