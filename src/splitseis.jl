# this script splits the yearlong seisio files into smaller files to
# optimize memory use per core
# For the test data, using cc_len = 2^15 seconds and cc_step = 26788 seconds,
# there are 1,177 time windows per year
# using Npersplit = 25 gives 1,177 \ 25 = 48 splits
using SeisIO, Glob, Dates

files = glob("*SAC",expanduser("~/SeisNoise_paper/DATA/YEARSAC/"))
OUTDIR = expanduser("~/SeisNoise_paper/DATA/SPLITSEIS")
cc_len, cc_step = 2^15, 26768
Npersplit = 25  # number of time-windows per split
starttime = DateTime(2019,1,1)

function splitsta(file::String,cc_len::Int,cc_step::Int,starttime::DateTime,Npersplit::Int,OUTDIR::String)

    println("Reading file ",basename(file))
    S = read_data("sac",file)
    demean!(S)
    OUTDIR = joinpath(OUTDIR,string(Npersplit))
    _,endtime  = u2d.(SeisIO.t_win(S[1].t,S[1].fs) .* 1e-6)

    if !isdir(OUTDIR)
        mkpath(OUTDIR)
    end

    # find number of splits1:
    Nsplits = convert(Int,ceil(((endtime - starttime).value / 1000 - cc_len) / cc_step / Npersplit))
    t1 = starttime
    splitNum = 1

    while t1 < endtime
        t2 = t1 + Second((Npersplit - 1) * cc_step + cc_len - 1)
        t2 = min(t2,endtime)
        Ssplit = sync(S,s=t1,t=t2,pad=true)

        splitstr = lpad(string(splitNum), Int(ceil(log10(Nsplits))), '0')
        SPLITDIR = joinpath(OUTDIR,"SPLIT" * splitstr)
        if !isdir(SPLITDIR)
            mkpath(SPLITDIR)
        end
        fileout = joinpath(SPLITDIR,Ssplit[1].id * ".sac")
        writesac(Ssplit,fname=fileout)

        # update
        t1 += Second(Npersplit * cc_step)
        splitNum += 1
    end
    return nothing
end


for file in files
    splitsta(file,cc_len,cc_step,starttime,Npersplit,OUTDIR)
end
