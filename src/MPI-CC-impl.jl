using MPI

function do_sendrecv(send_arr::AbstractArray,Nper::Int, src::Int, dst::Int)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println("rank $rank Sending $dst")
    MPI.send(send_arr, dst, rank+32, comm)
    MPI.Barrier(comm)
    # receive Array from src
    recv_arr,stat = MPI.recv(src,  src+32, comm)

    println("rank $rank Receiving $src")
    MPI.Barrier(comm)
    return recv_arr
end
