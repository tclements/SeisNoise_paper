# SeisNoise.jl: Ambient Seismic Noise Cross-Correlation on the CPU and GPU in Julia

This github repo hosts the benchmarks for Clements and Denolle, 2020: SeisNoise.jl: Ambient Seismic Noise Cross-Correlation on the CPU and GPU in Julia, submitted to Seismological Research Letters. 

Benchmarking scripts for this paper are stored in the `/src` directory. 

StationXML files for the 188 stations used in this paper are stored in `/DATA/XML`.

Waveforms for the 188 stations used in this paper are available on [Zenodo](https://zenodo.org/record/3823283). 

# How to Run Benchmarks
1. Clone this github repository `git clone https://github.com/tclements/SeisNoise_paper.git`.
2. Install Julia 1.4 and required packages using `. SeisNoise_paper/src/install_julia.sh`.
3. Download required files from Zenodo `julia SeisNoise_paper/src/get_zenodo.jl`.
4. Run `julia SeisNoise_paper/src/benchmarking.jl` (Note: This requires an Nvidia GPU for true benchmarking). 



