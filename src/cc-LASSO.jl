using CUDA, Dates, Glob, SeisNoise, Serialization

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

function corr2serial(F::FFTData,arr::Array{FFTData,1},maxlag::Real,CORRDIR::String)
    N = length(arr)
    for ii = 1:N
        C = correlate(F,arr[ii],maxlag)
        stack!(C,allstack=true)
        outpath = joinpath(CORRDIR,F.name,C.name)
        serialize(outpath,cpu(C))
    end
    return nothing
end

function all2all(arr1::Array{FFTData,1},arr2::Array{FFTData,1},maxlag::Real,CORRDIR::String)
    N = length(arr1)
    for ii = 1:N
        corr2serial(arr1[ii],arr2,maxlag,CORRDIR)
    end
    return nothing
end

function one2all(arr::Array{FFTData,1},maxlag::Real,CORRDIR::String)
    N = length(arr)
    for ii = 1:N-1
        corr2serial(arr[ii],arr[ii+1:end],maxlag,CORRDIR)
    end
    return nothing
end

function freemem!(arr::Array{FFTData,1})
	N = length(arr)
	for ii = 1:N
		CUDA.unsafe_free!(arr[ii].fft)
	end
	return nothing
end

function loadFFT(files::Array{String,1})
    N = length(files)
    arr = Array{FFTData}(undef,N)
    for ii = 1:N
        arr[ii] = deserialize(files[ii])
    end
    return arr
end

function splitFFT(files::AbstractArray,Nper::Int)
    Nfiles = length(files)
    splits = ceil(Int,Nfiles/Nper)
    out = []
    for ii = 1:splits
        startind = (ii-1)*Nper + 1
        endind = min(startind + Nper - 1,Nfiles)
        push!(out,files[startind:endind])
    end
    return out
end

function splitCC(splits::AbstractArray,maxlag::Real,CORRDIR::String)

	# load first split on the GPU
	arr1 = loadFFT(splits[1]) .|> gpu
	N = length(splits)
    for ii = 1:N
        # xcorr all in first array
		println("Correlating split $ii $(now())")
        one2all(arr1,maxlag,CORRDIR)

        # load other arrays
		arr2 = loadFFT(splits[N]) .|> gpu
        for jj = N:-1:ii+1
			println("Correlating split $jj $(now())")

			data = Channel{AbstractArray}(1)
			# do compute & I/O at same time
			@sync begin
				@async all2all(arr1,arr2,maxlag,CORRDIR)
				@async begin
					loadind = jj > ii + 1 ? jj - 1 : jj
					put!(data,loadFFT(splits[loadind]))
				end
			end

			# transfer CPU arr to gpu
			freemem!(arr2)
			arr2 = take!(data) .|> gpu
        end

		# transfer last arr2 to arr1 on GPU
		freemem!(arr1)
		arr1 = deepcopy(arr2)
		freemem!(arr2)
    end
    return nothing
end

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
