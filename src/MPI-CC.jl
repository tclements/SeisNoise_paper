using CUDA, Glob, MPI, SeisNoise, Dates, Serialization

# initialize MPI with CUDA
# suggested to use same number of CPU as GPU
MPI.Init()
comm = MPI.COMM_WORLD
lcomm = MPI.Comm_split_type(comm, MPI.MPI_COMM_TYPE_SHARED, MPI.Comm_rank(comm))
CUDA.device!(MPI.Comm_rank(lcomm) % length(devices()))
rank = MPI.Comm_rank(comm)
sz = MPI.Comm_size(comm)

# directory structure
FFTDIR = ""
CORRDIR = ""

# parameters
maxlag = 200.
Nper = 73

# files and such
files = glob("*",FFTDIR)
stas = [split(basename(s),'.')[2] |> (y -> parse(Int,y)) for s in files]
ind = sortperm(stas)
files = files[ind]
Nfiles = length(files)
Ntot = sz * Nper
splits = ceil(Int,Nfiles / Ntot)

# make output directories
if rank == 0
    for sta in stas
         mkpath(joinpath(CORRDIR,"2A.$sta..DPZ"))
     end
 end

MPI.Barrier(comm)

include("MPI-CC-impl.jl")
function main(files,splits,Nper,CORRDIR)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    sz = MPI.Comm_size(comm)
    println("rank $rank")
    Nfiles = length(files)
	Ntot = sz * Nper

    # loop through each split
    for ii = 1:splits
        if ii == splits
            extra = Nfiles % Ntot
            Nlocal = ceil(Int,extra / sz)
        else
            Nlocal = Nper
        end

        # get files for this rank
        startind = (ii - 1) * Ntot + rank * Nlocal + 1
        endind = (ii - 1) * Ntot + (rank + 1) * Nlocal
        endind = min(endind,Nfiles)
        rankfiles = files[startind:endind]

        # load files on rank
        arr1 = loadFFT(rankfiles) .|> gpu

        # correlate all files on this rank
        data = Channel{AbstractArray}(1)
        @sync begin
            @async one2all(arr1,maxlag,CORRDIR)
            @async begin
                src = mod(rank-1,sz)
                startind = (ii - 1) * Ntot + src * Nlocal + 1
                endind = (ii - 1) * Ntot + (src + 1) * Nlocal
                endind = min(endind,Nfiles)
                srcfiles = files[startind:endind]
                put!(data,loadFFT(srcfiles))
            end
        end
        arr2 = take!(data) .|> gpu
        MPI.Barrier(comm)

        # loop through splits
        for jj = 2:sz-1
			data = Channel{AbstractArray}(1)
            @sync begin
                @async all2all(arr1,arr2,maxlag,CORRDIR)
                @async begin
                    src = mod(rank-jj,sz)
                    startind = (ii - 1) * Ntot + src * Nlocal + 1
                    endind = (ii - 1) * Ntot + (src + 1) * Nlocal
                    endind = min(endind,Nfiles)
                    srcfiles = files[startind:endind]
                    put!(data,loadFFT(srcfiles))
                end
            end
            # transfer CPU arr to gpu
			freemem!(arr2)
			arr2 = take!(data) .|> gpu
        end
		freemem!(arr2)

        # loop through each remaining split
        for jj = ii + 1:splits

            if jj == splits
                extra = Nfiles % Ntot
                Nlocal = ceil(Int,extra / sz)
            else
                Nlocal = Nper
            end

            # get files for this rank
            startind = (jj - 1) * Ntot + rank * Nlocal + 1
            endind = (jj - 1) * Ntot + (rank + 1) * Nlocal
            endind = min(endind,Nfiles)
            rankfiles = files[startind:endind]

			# load onto GPU
            arr2 = loadFFT(rankfiles) .|> gpu

            for kk = 2:sz-1
				data = Channel{AbstractArray}(1)
				@sync begin
	            	@async all2all(arr1,arr2,maxlag,CORRDIR)
					@async begin
	                    src = mod(rank-kk,sz)
	                    startind = (jj - 1) * Ntot + src * Nlocal + 1
	                    endind = (jj - 1) * Ntot + (src + 1) * Nlocal
	                    endind = min(endind,Nfiles)
	                    srcfiles = files[startind:endind]
	                    put!(data,loadFFT(srcfiles))
	                end
				end
				freemem!(arr2)
				arr2 = take!(data) .|> gpu

            end
			freemem!(arr2)
            MPI.Barrier(comm)
        end
        MPI.Barrier(comm)
    end
    return nothing
end

# run test first
main(files[1:16],1,2,CORRDR)
if rank == 0
	testfiles = glob("*/*",CORRDIR)
	rm.(testfiles)
end
MPI.Barrier(comm)
main(files,splits,Nper,CORRDIR)

MPI.Finalize()
