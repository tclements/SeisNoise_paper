# CPU cross-correlation using a single core
# this was tested on a machine with 32 GB of memory
# This will run out of memory with 16 GB of memory

using SeisIO, SeisNoise, Dates, Glob
fs = 1.
cc_len, cc_step = 2^15, 26768
freqmin = 0.003
freqmax = 0.005
maxlag = 12000.
YEARSEIS = expanduser("~/Clements-Denolle-2020/DATA/YEARSEIS")
SINGLECORR = expanduser("~/Clements-Denolle-2020/SINGLECORR")
XML = expanduser("~/Clements-Denolle-2020/DATA/XML")
files = glob("*seisio",YEARSEIS)
N = length(files)

if !isdir(SINGLECORR)
    mkpath(SINGLECORR)
end

function seis2fft(file::String,XML::String,cc_len::Int,cc_step::Int)
    S = SeisData()
    println("Reading file ",basename(file[1:end-7]))
    read_data!(S,"sac",file)
    # detrend!(S)
    # taper!(S,t_max=400.)
    filtfilt!(S,fl=0.01,fh=0.4,np=3)

    # get xmlfile
    net,sta,loc,chan = split(S[1].id,'.')
    xmlfile = joinpath(XML,net*'.'*sta*".xml")
    Resp = read_meta("sxml",xmlfile)
    translate_resp!(S,Resp)
    remove_resp!(S)
    filtfilt!(S,fh=0.05,rt="Lowpass")
    R = RawData(S[1],cc_len,cc_step)
    demean!(R)
    taper!(R,max_length=1000.)
    return rfft(R)
end

function corrmapper(A::Array{FFTData,1},maxlag::Float64,OUTDIR::String)
        N = size(A,1)
        # copy the current FFT and correlate against all remaining
        for ii = 1:N-1
            for jj = ii+1:N
                mapcc(A[ii],A[jj],maxlag,OUTDIR)
            end
            A[ii] = FFTData()
        end
end

function mapcc(FFT1::FFTData,FFT2::FFTData,maxlag::Float64,OUTDIR::String)
    println("Correlation $(FFT1.name), $(FFT2.name)")
    C = compute_cc(FFT1,FFT2,maxlag)
    if !isnothing(C)
        abs_max!(C)
        stack!(C,allstack=true)
        save_corr(C,OUTDIR)
    else
        println("FAILED $(FFT1.name), $(FFT2.name)")
    end
    return nothing
end

# precompile
Fs = map(seis2fft,files[1:24],fill(XML,24),fill(cc_len,24),fill(cc_step,24))
corrmapper(Fs,maxlag,SINGLECORR)
toremove = glob("*",SINGLECORR)
rm.(toremove)

println("Beginning XCORR $(today())")
println("$(length(files)) total stations to process...")
println("$(length(files) * (length(files)-1) / 2) total correlation to process...")
tALL = @elapsed begin
Fs = map(seis2fft,files,fill(XML,N),fill(cc_len,N),fill(cc_step,N))
corrmapper(Fs,maxlag,SINGLECORR)
end

println("XCORR took $tALL seconds")
