using Distributed, Dates
t1 = now()
addprocs()

@everywhere begin

    using Dates, Glob, SeisIO, SeisNoise, Serialization

    # directories
    BASEDIR = "/media/four/"
    DATADIR = joinpath(BASEDIR,"LASSO/SAC/")
    FFTDIR = joinpath(BASEDIR,"LASSO/FFT/")
    XMLDIR = joinpath(BASEDIR,"LASSO/XML")
    files = glob("*.SAC",DATADIR)

    if !isdir(FFTDIR)
        mkpath(FFTDIR)
    end

    # parameters
    fs = 250.
    cc_len = 900.
    cc_step = 450.
    freqmin = 0.1
    freqmax = 20.

    # fast saving of FFT
    function serialize_fft(F::FFTData,FFTDIR::String)
        filename = joinpath(FFTDIR,F.name)
        if isa(F.fft,SeisNoise.AbstractGPUArray)
            F = F |> cpu
        end
        serialize(filename,F)
        return nothing
    end

    function seis2fft(
        file::String,
        XMLDIR::String,
        FFTDIR::String,
        fs::Real,
        cc_len::Real,
        cc_step::Real,
        freqmin::Float64,
        freqmax::Float64,
    )
        println("Reading file $(basename(file))")
        S = read_data("sac",file)
        S.fs .= round.(S.fs)
        detrend!(S)
        taper!(S)
        ungap!(S)
        filtfilt!(S,rt="Lowpass",fh=fs/2,np=3)
        resample!(S,fs=fs)

        # get xmlfile
        net,sta,loc,chan = split(S[1].id,'.')
        xmlfile = joinpath(XMLDIR,net*'.'*sta*".xml")
        Resp = read_meta("sxml",xmlfile)
        S.loc[1] = Resp.loc[1]
        translate_resp!(S,Resp[1].resp)
        remove_resp!(S)
        filtfilt!(S,fl=freqmin,fh=freqmax,np=2)
        R = RawData(S[1],cc_len,cc_step)
        demean!(R)
        taper!(R)
        F = rfft(R)
        serialize_fft(F,FFTDIR)
        return nothing
    end
end

t2 = now()
pmap((file -> seis2fft(file,XMLDIR,FFTDIR,fs,cc_len,cc_step,freqmin,freqmax)),files)
t3 = now()

# time for stuff
compute = Dates.canonicalize(Dates.CompoundPeriod(t3-t2))
total = Dates.canonicalize(Dates.CompoundPeriod(t3-t1))
println("Total computation took $total")
println("Computation took $compute")
