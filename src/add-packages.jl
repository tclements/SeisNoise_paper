using Pkg

# this script installs all the packages necessary to reproduce the paper
# except for GMT.jl, which requires a local copy of GMT
# add following packages
pkgs = ["AWSCore", "AWSS3", "DataFrames", "CSV", "Glob", "JLD2",
        "Plots", "CuArrays", "CUDAnative", "CUDAdrv", "BenchmarkTools",
        "PyCall", "Conda","StatsPlots"]
for p in pkgs
      Pkg.add(p)
end
Pkg.add(PackageSpec(name="SeisIO",rev="master"))
Pkg.add(PackageSpec(name="SeisNoise",rev="master"))

using Conda
Conda.add("numpy")
Conda.add("scipy")
Conda.add("obspy")
