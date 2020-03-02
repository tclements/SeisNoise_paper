using Pkg

# add following packages
pkgs = ["AWSCore", "AWSS3", "DataFrames", "CSV", "DSP", "FFTW", "Glob", "JLD2", "SeisIO",
       "Interpolations", "GLM", "Plots"]
for p in pkgs
      Pkg.add(p)
end
Pkg.add(PackageSpec(name="SeisNoise",rev="GPU"))
