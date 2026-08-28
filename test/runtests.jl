# `Pkg.test()`'s entry point. Each suite runs in its own module so their
# top-level constants (every file has its own MASTER_SEED) never collide,
# and Aqua runs first: the package-quality checks a registry reviewer
# would run (method ambiguities, stale or uncapped dependencies, type
# piracy, undefined exports). TESTING.md says what each suite protects.
using Test, Aqua
using ActiveSetLookup

@testset "ActiveSetLookup" begin
    @testset "Aqua: package quality" begin
        Aqua.test_all(ActiveSetLookup)
    end
    @eval module TestPublicAPI
        include("test_public_api.jl")
    end
    @eval module TestCertificate
        include("test_certificate.jl")
    end
    @eval module TestGuide
        include("test_guide.jl")
    end
end
