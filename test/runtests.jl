using CryspBrainSim
using Test

@testset "CryspBrainSim" begin
    # Smoke test: the package and its RecoCryspTools dependency load.
    # The endpoint-port self-tests (dev/PLAN.md build step 2) land here next.
    @test isdefined(CryspBrainSim, :RecoCryspTools)
end
