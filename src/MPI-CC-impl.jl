using MPI, SeisNoise, CUDA, Serialization, Dates

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
		data = Channel{AbstractArray}(1)

		@sync begin
			@async one2all(arr1,maxlag,CORRDIR)
			@async put!(data,loadFFT(splits[N]))
		end

        # load other arrays
		arr2 = take!(data) .|> gpu
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
