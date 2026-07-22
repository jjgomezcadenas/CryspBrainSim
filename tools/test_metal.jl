#!/usr/bin/env julia
# Minimal standalone Metal.jl smoke test.
# Run directly from a Terminal, not through a sandboxed runner:
#   julia --project=. tools/test_metal.jl

using Metal

println("Metal.jl version: ", Base.pkgversion(Metal))
println("Metal.functional(): ", Metal.functional())

devices = Metal.devices()
println("Metal devices found: ", length(devices))
for (index, device) in enumerate(devices)
    println("  [$index] ", device.name)
end

isempty(devices) && error(
    "No Metal device is visible to this Julia process. " *
    "Run this script from a normal local Terminal session."
)
Metal.functional() || error("Metal.jl found a device but is not functional")

# Exercise device allocation, a GPU broadcast kernel, reduction, and copy-back.
host = collect(Float32, 1:1024)
gpu = MtlArray(host)
gpu .= 2f0 .* gpu .+ 1f0
synchronize()
result = Array(gpu)

expected = 2f0 .* host .+ 1f0
result == expected || error("Metal result differs from the CPU reference")

println("GPU array smoke test: PASS")
println("First / last result: $(result[1]) / $(result[end])")
