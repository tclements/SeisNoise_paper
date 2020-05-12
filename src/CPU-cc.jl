# CPU cross-correlation using files split into Npersplit windows per file
# run splitseis.jl with Npersplit set to 25 before running this
# this assumes a 48-core machine

using Distributed
addprocs(48,topology=:master_worker)
@everywhere begin
using SeisIO, SeisNoise, Dates, Glob
fs = 1.
cc_len, cc_step = 2^15, 26768
freqmin = 0.003
freqmax = 0.005
maxlag = 12000.
Npersplit = 25
splits = glob("*",joinpath(expanduser("~/SeisNoise_paper/DATA/SPLITSEIS"),string(Npersplit)))
Nsplits = length(splits)
SPLITCORR = joinpath(expanduser("~/SeisNoise_paper/CORR/SPLITSEIS"),string(Npersplit))
XML = expanduser("~/SeisNoise_paper/DATA/XML")
N = length(splits)

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
    filtfilt!(S,fl=freqmin,fh=freqmax,np=3)
    R = RawData(S[1],cc_len,cc_step)
    demean!(R)
    taper!(R)
    return rfft(R)
end

# this function cross-correlates data for a "split"
function splitcorr(SPLITDIR::String,XML::String,cc_len::Int,cc_step::Int,
                  freqmin::Float64,freqmax::Float64,maxlag::Float64,OUTDIR::String)
    files = glob("*",SPLITDIR)
    Nfiles = length(files)
    println("Cross-correlating chunk $(basename(SPLITDIR)) ",now())
    SPLITOUT = joinpath(OUTDIR,basename(SPLITDIR))
    if !isdir(SPLITOUT)
        mkpath(SPLITOUT)
    end

    # I/O + fft
    Fs = Array{FFTData}(undef,Nfiles)
    for ii = 1:Nfiles
        Fs[ii] = seis2fft(files[ii],XML,cc_len,cc_step,freqmin,freqmax)
    end

    # cross-correlate
    for ii = 1:Nfiles -1
        for jj = ii+1:Nfiles
            C = compute_cc(Fs[ii],Fs[jj],maxlag)
            stack!(C,allstack=true)
            save_corr(C,SPLITOUT)
        end
    end
    return nothing

end
end

# precompile
pmap(splitcorr,splits[1:1],[XML],[cc_len],[cc_step],
                  [freqmin],[freqmax],[maxlag],[SPLITCORR])
toremove = glob("*/*",SPLITCORR)
rm.(toremove)

println("Beginning XCORR with $(nprocs()-1) processors...")
tALL = @elapsed begin
pmap(splitcorr,splits,fill(XML,Nsplits),fill(cc_len,Nsplits),fill(cc_step,Nsplits),
     fill(freqmin,Nsplits),fill(freqmax,Nsplits),fill(maxlag,Nsplits),fill(SPLITCORR,Nsplits))
end

println("XCORR took $tALL seconds")
