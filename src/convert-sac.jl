# convert MSEED to SAC
using SeisIO, Glob, Dates

BASEDIR = "/media/four/"
DATADIR = joinpath(BASEDIR,"LASSO/DL/")
SACDIR = joinpath(BASEDIR,"LASSO/SAC/")

if !isdir(SACDIR)
    mkpath(SACDIR)
end

files = glob("*",DATADIR)
cd(SACDIR)

for ii = 1:length(files)
    println("Writing file $(basename(files[ii])), $ii of $(length(files)) $(now())")
    S = read_data("mseed",files[ii])
    writesac(S)
end
