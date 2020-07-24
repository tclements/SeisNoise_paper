using CUDA, Dates, Glob, SeisNoise, Serialization
include("MPI-CC-impl.jl")

# files and such
BASEDIR = "/media/four/"
FFTDIR = joinpath(BASEDIR,"LASSO/FFT/")
CORRDIR = joinpath(BASEDIR,"LASSO/CORR/")
if !isdir(CORRDIR)
    mkpath(CORRDIR)
end
files = glob("*",FFTDIR)
stas = [split(basename(s),'.')[2] |> (y -> parse(Int,y)) for s in files]
ind = sortperm(stas)
files = files[ind]

# create directory structure for saving CORRs
for sta in stas
    mkpath(joinpath(CORRDIR,"2A.$sta..DPZ"))
end

# parameters
maxlag = 200.
Nper = 18

# warm up functions
splits = splitFFT(files,Nper)
splitCC(splits[1:3],maxlag,CORRDIR)
testfiles = glob("*/*",CORRDIR)
rm.(testfiles)

t1 = now()
## xcorr on the GPU minimizing I/O
splitCC(splits,maxlag,CORRDIR)
t2 = now()
total = Dates.canonicalize(Dates.CompoundPeriod(t2-t1))
println("Total computation took $total")
