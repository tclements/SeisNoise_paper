using Pkg

# add following packages
pkgs = ["AWSCore", "AWSS3", "DataFrames", "CSV", "Glob", "JLD2", "SeisIO",
        "Plots, "CuArrays", "BenchmarkTools"]
for p in pkgs
      Pkg.add(p)
end
Pkg.add(PackageSpec(name="SeisNoise",rev="GPU"))
