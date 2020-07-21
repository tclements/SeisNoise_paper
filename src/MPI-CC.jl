# using CUDA, Glob, MPI, SeisNoise
using MPI, CUDA

MPI.Init()
comm = MPI.COMM_WORLD
lcomm = MPI.Comm_split_type(comm, MPI.MPI_COMM_TYPE_SHARED, MPI.Comm_rank(comm))
# CUDA.device!(MPI.Comm_rank(lcomm) % length(devices()))

rank = MPI.Comm_rank(comm)
sz = MPI.Comm_size(comm)

FFTDIR = ""
CORRDIR = ""

# parameters
maxlag = 200.
Nper = 72

# files and such
# files = glob("*",FFTDIR)
# stas = [split(basename(s),'.')[2] |> (y -> parse(Int,y)) for s in files]
# ind = sortperm(stas)
# files = files[ind]
# files = ["2A.$ii..DPZ" for ii = 1:1825]
files = [rand(5) |> cu for ii = 1:1825]
Nfiles = length(files)
Ntot = sz * Nper
splits = ceil(Int,Nfiles / Ntot)

# if rank == 0
#     for sta in stas
#          mkpath(joinpath(CORRDIR,"2A.$sta..DPZ"))
#      end
#  end

MPI.Barrier(comm)

include("MPI-CC-impl.jl")
function main(files,splits,Nper)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    sz = MPI.Comm_size(comm)
    println("rank $rank")
    Ntot = length(files)

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
        rank_arr = loadFFT(rankfiles) .|> gpu

        # correlate all files on this rank
        # one2all(rank_arr,maxlag,CORRDIR)
        MPI.Barrier(comm)

        # source / receiver for this rank
        dst = mod(rank+1,sz)
        src = mod(rank-1,sz)

        # get_arr = dosend_recv(rank_arr,Nlocal,src,dst)
        get_arr = do_sendrecv(rank_arr,Nlocal,src,dst)
        MPI.Barrier(comm)
        all2all(rank_arr,get_arr,maxlag,CORRDIR)

        # loop through splits
        for jj = 2:sz-1
            get_arr = do_sendrecv(get_arr,Nlocal,src,dst)
            all2all(rank_arr,get_arr,maxlag,CORRDIR)
        end

        freemem!(get_arr)

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
            getfiles = files[startind:endind]

            get_arr = getfiles
            get_arr = loadFFT(getfiles) .|> gpu
            all2all(rank_arr,get_arr,maxlag,CORRDIR)
            for kk = 2:sz-1
                get_arr = do_sendrecv(get_arr,Nlocal,src,dst)
                all2all(rank_arr,get_arr,maxlag,CORRDIR)
            end
            MPI.Barrier(comm)
        end
        MPI.Barrier(comm)
    end
    return nothing
end
main(files,splits,Nper)

MPI.Finalize()
