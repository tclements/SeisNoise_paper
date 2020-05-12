# GPU cross-correlation using files split into Npersplit windows per file
# run splitseis.jl with Npersplit set to 236 before running this
# Note: this will only work if an Nvidia GPU is available
using SeisIO, SeisNoise, Dates, Glob, BenchmarkTools, Printf

Npersplit = 236
cc_len, cc_step = 2^15, 26768
freqmin = 0.003
freqmax = 0.005
maxlag = 12000.
splits = glob("*",joinpath(expanduser("~/Clements-Denolle-2020/DATA/SPLITSEIS"),string(Npersplit)))
SPLITCORR = joinpath(expanduser("~/Clements-Denolle-2020/CORR/SPLITSEIS"),string(Npersplit))
XML = expanduser("~/Clements-Denolle-2020/DATA/XML")

if !isdir(SPLITCORR)
    mkpath(SPLITCORR)
end

# this function reads a file, applies pre-processing, converts to RawData
#  and then returns an FFTData
function seis2fft(file::String,XML::String,cc_len::Int,cc_step::Int,
                  freqmin::Float64,freqmax::Float64)
    S = SeisData()
    read_data!(S,"sac",file)
    detrend!(S)
    taper!(S)
    filtfilt!(S,fl=0.001,fh=0.4,np=3)

    # get xmlfile
    net,sta,loc,chan = split(S[1].id,'.')
    xmlfile = joinpath(XML,net*'.'*sta*".xml")
    Resp = read_meta("sxml",xmlfile)
    translate_resp!(S,Resp)
    remove_resp!(S)
    R = RawData(S[1],cc_len,cc_step) |> gpu
    demean!(R)
    taper!(R)
    bandpass!(R,freqmin,freqmax,corners=3)
    return rfft(R)
end

# this function cross-correlates data for a "split"
function GPUALL(splits,cc_len,cc_step,freqmin,freqmax,maxlag,OUTDIR,XML)
    Nsplits = length(splits)
    for ii = 1:Nsplits
        println("Cross-correlating chunk $ii of $Nsplits ",now())
        files = glob("*",splits[ii])
        Nfiles = length(files)

        # I/O + fft
        Fs = Array{FFTData}(undef,Nfiles)
        for jj = 1:Nfiles
            Fs[jj] = seis2fft(files[jj],XML,cc_len,cc_step,freqmin,freqmax)
        end

        # cross-correlate
        for jj = 1:Nfiles -1
            for kk = jj+1:Nfiles
                C = compute_cc(Fs[jj],Fs[kk],maxlag)
                stack!(C,allstack=true)
                C = cpu(C)
                finalize(C.corr)
                save_corr(C,OUTDIR)
            end
            finalize(Fs[jj].fft)
        end
    end
    return nothing
end

# warm up the GPU
GPUALL(splits[1:1],cc_len,cc_step,freqmin,freqmax,maxlag,SPLITCORR,XML)
toremove = glob("*",SPLITCORR)
rm.(toremove)

tGPU = @elapsed begin
    GPUALL(splits,cc_len,cc_step,freqmin,freqmax,maxlag,SPLITCORR,XML)
end
println("GPU XCORR took $(@sprintf("%.2f",tGPU)) seconds")
