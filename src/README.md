## SeisNoise.jl Benchmarking / Figure creation scripts

Files in this repo: 

- `install_julia.sh`: install Julia version 1.4 
- `add-packages.jl`: install Julia and Python packages necessary to run the benchmarks 
- `benchmarking.jl`: is the main benchmarking script that creates Figure 2
- `get_zenodo.jl`: script to download data for further benchmarks from [Zenodo](https://zenodo.org/record/3823283)
- `single-LHZ.jl`: performs cross-correlation for the dataset with a single core 
- `CPU-cc.jl`: performs cross-correlation with many CPU cores, requires running `splitseis.jl` first 
- `CPU-cc.jl`: performs cross-correlation with on the GPU, requires running `splitseis.jl` first
- `splitseis.jl`: Splits year-long sac files from Zenodo into many smaller files
- `moveout.jl`: plots the cross-correlation moveout from Figure 4 
- `stationmap.jl`: plots the map from Figure 3 with GMT.jl
