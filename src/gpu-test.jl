# this file implements cross-correlation on the GPU using one day of files from
# the AWS SCEDC dataset

using SeisIO, SeisNoise, Dates, Glob, CuArrays, CUDAdrv, BenchmarkTools, Plots

fs = 100.
cc_len = 7200
cc_step = 1800
freqmin = 0.1
freqmax = 0.3
maxlag = 100.
DATA = expanduser("~/data/continuous_waveforms/2020/2020_001")
CORR = expanduser("~/CORR")
files = glob("*",DATA)
Nfiles = length(files)
channel = basename(files[1])[8:10]
GPUMEM,GPUAVAIL = Mem.info()

println("Cross-correlating $Nfiles $channel stations from the CI Network")
println("Working on a $(name(device())) with $(round(GPUMEM / 1e9,digits=2)) GB VRAM")

if !isdir(CORR)
    mkpath(CORR)
end

"""
  gpu_cc(files, cc_len, cc_step, freqmin, freqmax, maxlag)


Cross-correlation on the GPU!

#  Arguments
- `files`: List of mseed files for input to cross-correlation.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `freqmin::Float64`: Pass band low corner frequency.
- `freqmax::Float64`: Pass band high corner frequency.
- `maxlag::Int`: Number of data points in cross-correlation to save,
                 e.g. `maxlag = 2000` will save lag times = -2000/fs:2000/fs s.
"""
function gpu_cc(files,cc_len::Int,cc_step::Int,freqmin::Float64,freqmax::Float64,maxlag::AbstractFloat)
    S = SeisData()
    read_data!(S,"mseed",files)
    ungap!(S)
    phase_shift!(S)
    Rs = Array{RawData}(undef,S.n)
    Fs = Array{FFTData}(undef,S.n)
    N = Int(cc_len * S[1].fs)
    for ii = 1:S.n
        Rs[ii] = RawData(S[ii],cc_len,cc_step) |> gpu
        SeisNoise.detrend!(Rs[ii])
        SeisNoise.demean!(Rs[ii])
        SeisNoise.taper!(Rs[ii])
        SeisNoise.bandpass!(Rs[ii],freqmin,freqmax)
        Fs[ii] = compute_fft(Rs[ii]) |> cpu
        Rs[ii] = RawData()
    end
    for ii = 1:S.n - 1
        for jj = 2:S.n
            compute_cc(Fs[ii],Fs[jj],maxlag)
        end
    end
    return nothing
end

function cpu_cc(files,cc_len::Int,cc_step::Int,freqmin::Float64,freqmax::Float64,maxlag::AbstractFloat)
    S = SeisData()
    read_data!(S,"mseed",files)
    ungap!(S)
    phase_shift!(S)
    Rs = Array{RawData}(undef,S.n)
    Fs = Array{FFTData}(undef,S.n)
    for ii = 1:S.n
        Rs[ii] = RawData(S[ii],cc_len,cc_step)
        SeisNoise.detrend!(Rs[ii])
        SeisNoise.demean!(Rs[ii])
        SeisNoise.taper!(Rs[ii])
        SeisNoise.bandpass!(Rs[ii],freqmin,freqmax)
        Fs[ii] = compute_fft(Rs[ii])
    end

    for ii = 1:S.n - 1
        for jj = 2:S.n
            compute_cc(Fs[ii],Fs[jj],maxlag)
        end
    end
    return nothing
end

numpairs = (2:Nfiles) .* (1:(Nfiles-1)) .÷ 2
idealpairs = 2 .^ (0:15)
actualpairs = Array{Int}(undef,length(idealpairs))
pairsind = Array{Int}(undef,length(idealpairs))
numfiles = 2:Nfiles
GPUtimes = Array{Float32}(undef,length(idealpairs))
CPUtimes = Array{Float32}(undef,length(idealpairs))

for ii = 1:length(idealpairs)
    # find closest number of pairs
    pairfinder = numpairs .- idealpairs[ii]
    pairfinder[pairfinder .< 0] .= typemax(Int)
    ind = argmin(pairfinder)
    if ii == length(idealpairs)
        ind = Nfiles - 1
    end
    actualpairs[ii] = numpairs[ind]
    pairsind[ii] = ind
end

for ii = 1:length(idealpairs)
    # time number of pairs
    GPUtimes[ii] = CUDAdrv.@elapsed gpu_cc(files[1:numfiles[pairsind[ii]]],cc_len,cc_step,freqmin,freqmax,maxlag)
    println("$(actualpairs[ii]) pairs took: $(GPUtimes[ii]) seconds")

end

for ii = 1:length(idealpairs)
    # time number of pairs
    CPUtimes[ii] = Base.@elapsed cpu_cc(files[1:numfiles[pairsind[ii]]],cc_len,cc_step,freqmin,freqmax,maxlag)
    println("$(actualpairs[ii]) pairs took: $(CPUtimes[ii]) seconds")
end

scatter(actualpairs,GPUtimes,label="GPU",xlabel="Number of Pairs",xscale=:log10,
        ylabel="Time [s]",yscale=:log10,legend=:topleft)
scatter!(actualpairs,CPUtimes,label="CPU",xlabel="Number of Pairs",xscale=:log10,
        ylabel="Time [s]",yscale=:log10,legend=:topleft)
