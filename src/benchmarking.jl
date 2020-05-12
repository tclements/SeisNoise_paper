using SeisIO, SeisNoise, Dates, Glob, BenchmarkTools, StatsPlots, Statistics
using CuArrays, PyCall, Plots.PlotMeasures, Printf
## This is a test of single core Julia CPU performace vs GPU performance vs python
# for one year of cross-correlation between two LHZ stations
# numpy, scipy and obspy must be installed via Conda for this benchmark

fs = 1.
cc_len, cc_step = 2^15, 26768
freqmin = 0.003
freqmax = 0.005
maxlag = 12000.
file1 = expanduser("~/SeisNoise_paper/YEARSAC/2019.001.00.00.00.000.BK.CMB.00.LHZ.R.SAC")
file2 = expanduser("~/SeisNoise_paper/YEARSAC/2019.001.00.00.00.000.US.WVOR.00.LHZ.R.SAC")

# function for cross-correlation on the CPU
function CPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag)
    S1 = read_data("sac",file1)
    S2 = read_data("sac",file2)
    R1 = RawData(S1[1],cc_len,cc_step)
    R2 = RawData(S2[1],cc_len,cc_step)
    detrend!(R1)
    detrend!(R2)
    taper!(R1)
    taper!(R2)
    highpass!(R1,0.001)
    highpass!(R2,0.001)
    onebit!(R1)
    onebit!(R2)
    F1 = rfft(R1)
    F2 = rfft(R2)
    C = compute_cc(F1,F2,maxlag)
    stack!(C,allstack=true)
    clean_up!(C,freqmin,freqmax)
    return C
end

# function for cross-correlation on the CPU
function GPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag)
    S1 = read_data("sac",file1)
    S2 = read_data("sac",file2)
    R1 = RawData(S1[1],cc_len,cc_step) |> gpu
    R2 = RawData(S2[1],cc_len,cc_step) |> gpu
    detrend!(R1)
    detrend!(R2)
    taper!(R1)
    taper!(R2)
    highpass!(R1,0.001)
    highpass!(R2,0.001)
    onebit!(R1)
    onebit!(R2)
    F1 = rfft(R1)
    F2 = rfft(R2)
    C = compute_cc(F1,F2,maxlag)
    stack!(C,allstack=true)
    clean_up!(C,freqmin,freqmax)
    return C |> cpu
end


#### TEST SINGLE CORE SINGLE CC ####
CPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag)
GPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag)

# benchmark total compute time
println("Benchmarking Total CPU Time")
println("---------------------------")
a = @benchmarkable CPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag) seconds=60 time_tolerance=0.01
TOTALCPU = run(a)
display("text/plain", TOTALCPU)

println("\n\nBenchmarking Total GPU Time")
println("---------------------------")
b = @benchmarkable begin CuArrays.@sync GPU_CC(file1,file2,cc_len,cc_step,fs,freqmin,freqmax,maxlag) end seconds=60 time_tolerance=0.01
TOTALGPU = run(b)
display("text/plain", TOTALGPU)

# benchmark file read speed
println("\n\nBenchmarking SeisIO Read Speed")
println("------------------------------")
SEISIO = @benchmark read_data("sac",file1)
display("text/plain", SEISIO)
S1 = read_data("sac",file1);
S2 = read_data("sac",file2);

# benchmark file read speed
println("\n\nBenchmarking Obspy Read Speed")
println("-----------------------------")
obspy = pyimport("obspy")
OBSPYIO = @benchmark obspy.read(file1)
display("text/plain", OBSPYIO)

# benchmark RawData creation
println("\n\nBenchmarking RawData CPU")
println("------------------------")
RAWCPU = @benchmark RawData(S1[1],cc_len,cc_step)
display("text/plain", RAWCPU)
R1 = RawData(S1[1],cc_len,cc_step);
R2 = RawData(S2[1],cc_len,cc_step);

# benchmark detrending
println("\n\nBenchmarking detrending CPU")
println("---------------------------")
DETRENDCPU = @benchmark detrend!(R1)
display("text/plain", DETRENDCPU)

# benchmark tapering
println("\n\nBenchmarking tapering CPU")
println("-------------------------")
TAPERCPU = @benchmark taper!(R1,max_length=400.)
display("text/plain", TAPERCPU)

# benchmark filtering
println("\n\nBenchmarking filtering CPU")
println("-------------------------")
FILTERCPU = @benchmark highpass!(R1,0.001)
display("text/plain", FILTERCPU)

# benchmark onebit
println("\n\nBenchmarking time-normalization CPU")
println("-----------------------------------")
ONEBITCPU = @benchmark onebit!(R1)
display("text/plain", ONEBITCPU)

# benchmark Fourier transform CPU
println("\n\nBenchmarking rfft CPU")
println("---------------------")
RFFTCPU = @benchmark rfft(R1)
display("text/plain", RFFTCPU)
F1 = rfft(R1);
F2 = rfft(R2);

# benchmark cross-correlation CPU
println("\n\nBenchmarking cross-correlation CPU")
println("----------------------------------")
CCCPU = @benchmark compute_cc(F1,F2,maxlag)
display("text/plain", CCCPU)
C = compute_cc(F1,F2,maxlag)

# benchmark stacking CPU
println("\n\nBenchmarking stacking CPU")
println("-------------------------")
# this times just the stack and not the reallocation
STACKCPU = @benchmark mean(C.corr,dims=2)
display("text/plain", STACKCPU)
stack!(C,allstack=true)

# benchmark final clean-up
println("\n\nBenchmarking post-processing CPU")
println("--------------------------------")
# this times just the stack and not the reallocation
POSTCPU = @benchmark clean_up!(C,freqmin,freqmax)
display("text/plain", POSTCPU)

# benchmark RawData creation
println("\n\nBenchmarking RawData GPU")
println("------------------------")
c = @benchmarkable begin CuArrays.@sync RawData(S1[1],cc_len,cc_step) |> gpu end
RAWGPU = run(c)
display("text/plain", RAWGPU)
R1 = RawData(S1[1],cc_len,cc_step) |> gpu;
R2 = RawData(S2[1],cc_len,cc_step) |> gpu;

println("\n\nBenchmarking detrending GPU")
println("---------------------------")
d = @benchmarkable begin CuArrays.@sync detrend!(R1) end
DETRENDGPU = run(d)
display("text/plain", DETRENDGPU)

# benchmark tapering
println("\n\nBenchmarking tapering GPU")
println("-------------------------")
e = @benchmarkable begin CuArrays.@sync taper!(R1,max_length=400.) end
TAPERGPU = run(e)
display("text/plain", TAPERGPU)

# benchmark filtering
println("\n\nBenchmarking filtering GPU")
println("-------------------------")
f = @benchmarkable begin CuArrays.@sync highpass!(R1,0.001) end
FILTERGPU = run(f)
display("text/plain", FILTERGPU)

# benchmark onebit
println("\n\nBenchmarking time-normalization GPU")
println("-----------------------------------")
g = @benchmarkable begin CuArrays.@sync onebit!(R1) end
ONEBITGPU = run(g)
display("text/plain", ONEBITGPU)

# benchmark Fourier transform CPU
println("\n\nBenchmarking rfft GPU")
println("---------------------")
h = @benchmarkable begin CuArrays.@sync rfft(R1) end
RFFTGPU = run(h)
display("text/plain", RFFTGPU)
F1 = rfft(R1);
F2 = rfft(R2);

# benchmark cross-correlation CPU
println("\n\nBenchmarking cross-correlation GPU")
println("----------------------------------")
i = @benchmarkable begin CuArrays.@sync compute_cc(F1,F2,maxlag) |> cpu end
CCGPU = run(i)
display("text/plain", CCGPU)
C = compute_cc(F1,F2,maxlag)

# benchmark stacking CPU
println("\n\nBenchmarking stacking GPU")
println("-------------------------")
# this times just the stack and not the reallocation
j = @benchmarkable begin CuArrays.@sync mean(C.corr,dims=2) end
STACKGPU = run(j)
display("text/plain", STACKGPU)
stack!(C,allstack=true)

# benchmark final clean-up
println("\n\nBenchmarking post-processing GPU")
println("--------------------------------")
k = @benchmarkable begin CuArrays.@sync clean_up!(cpu(C),freqmin,freqmax) end
POSTGPU = run(k)
display("text/plain", POSTGPU)

## python benchmarks ##
# reload RawData
R1 = RawData(S1[1],cc_len,cc_step);
R2 = RawData(S2[1],cc_len,cc_step);
# send to numpy
Rnp1 = PyReverseDims(R1.x)
Rnp2 = PyReverseDims(R2.x)
# import scipy modules
np = pyimport("numpy")
signal = pyimport("scipy.signal")
spf = pyimport("scipy.fft")

println("\n\nBenchmarking detrending SciPy")
println("-----------------------------")
DETRENDSCI = @benchmark @pycall signal.detrend(Rnp1)::PyArray
display("text/plain", DETRENDSCI)

# benchmark tapering
println("\n\nBenchmarking tapering SciPy")
println("---------------------------")
py"""
from scipy import signal
import numpy as np
def taper(A, fs,max_length=400):
    window = np.int(fs * max_length)
    tape = signal.hann(max_length * 2)
    if A.ndim == 2:
        A[:,:window] *= tape[:window]
        A[:,-window:] *= tape[-window:]
    elif A.ndim == 1:
        A[:window] *= tape[:window]
        A[-window:] *= tape[-window:]

    return A
"""
TAPERSCI = @benchmark @pycall py"taper"(Rnp1,fs)::PyArray
display("text/plain", TAPERCPU)

# benchmark filtering
println("\n\nBenchmarking filtering SciPy")
println("----------------------------")
sos = @pycall signal.butter(4,0.001,btype="highpass",fs=fs,output="sos")::PyArray
FILTERSCI = @benchmark @pycall signal.sosfiltfilt(sos,Rnp1)::PyArray
display("text/plain", FILTERSCI)

# benchmark onebit
println("\n\nBenchmarking time-normalization SciPy")
println("-------------------------------------")
ONEBITSCI = @benchmark @pycall np.sign(Rnp1)::PyArray
display("text/plain", ONEBITSCI)

# benchmark Fourier transform CPU
println("\n\nBenchmarking rfft SciPy")
println("-----------------------")
RFFTSCI = @benchmark @pycall spf.rfft(Rnp1)::PyArray
display("text/plain", RFFTSCI)
Fsc1 = @pycall spf.rfft(Rnp1)::PyArray
Fsc2 = @pycall spf.rfft(Rnp2)::PyArray

# benchmark cross-correlation
println("\n\nBenchmarking cross-correlation SciPy")
println("------------------------------------")
py"""
from scipy.fft import rfft, irfft, fftshift
import numpy as np
# cross-correlation function
def xcorr(A,B,N,maxlag):
    corr = irfft(np.conj(A) * B)
    t =  np.hstack((np.arange(0,N/2),np.arange(-N/2,0)))
    ind = np.where(np.abs(t) <= maxlag)[0]
    newind = fftshift(ind)
    return corr[:,newind]
"""
N = np.shape(Rnp1)[2]
CCSCI = @benchmark @pycall py"xcorr"(Fsc1,Fsc2,N,maxlag)::PyArray
display("text/plain", CCSCI)
Xsci = @pycall py"xcorr"(Fsc1,Fsc2,N,maxlag)::PyArray

# benchmark stacking CPU
println("\n\nBenchmarking stacking SciPy")
println("---------------------------")
# this times just the stack and not the reallocation
STACKSCI = @benchmark @pycall np.mean(Xsci,axis=0)::PyArray
display("text/plain", STACKSCI)
Xstack = @pycall np.mean(Xsci,axis=0)::PyArray

# benchmark final clean-up
println("\n\nBenchmarking post-processing SciPy")
println("----------------------------------")
py"""
def clean_up(A,freqmin,freqmax,fs):
    A = signal.detrend(A, type='linear')
    A = taper(A,fs)
    sos = signal.butter(4,[freqmin,freqmax],btype="bandpass",fs=fs,output="sos")
    A = signal.sosfiltfilt(sos,A)
    return A
"""
POSTSCI = @benchmark @pycall py"clean_up"(Xstack,freqmin,freqmax,fs)::PyArray
display("text/plain", POSTSCI)

processes = [SEISIO,RAWGPU,DETRENDGPU,TAPERGPU,FILTERGPU,ONEBITGPU,RFFTGPU,CCGPU,
             STACKGPU,POSTGPU,SEISIO, RAWCPU,DETRENDCPU,TAPERCPU,FILTERCPU,ONEBITCPU,
             RFFTCPU,CCCPU,STACKCPU,POSTCPU]
SCIPY = [DETRENDSCI,TAPERSCI,FILTERSCI,ONEBITSCI,RFFTSCI,CCSCI,STACKSCI,POSTSCI]
ctg = repeat(["CUDA", "Julia","Python"], inner = 10)
colors = repeat(["white", "black","grey"], inner = 10)
nam = repeat(["I/O","RawData","Detrend","Taper","Filter","Time\n Normalization","RFFT","Xcorr","Stack","Post-\nProcess"], outer = 3)

times = [median(p).time * 1e-6 for p in processes]
scitimes = vcat([median(OBSPYIO).time,0], [median(p).time for p in SCIPY]) .* 1e-6
scitotal = sum(scitimes .* [1,1,2,2,2,2,2,1,1,1])
scitotal += (median(OBSPYIO).time + median(RAWCPU).time) * 1e-6 * 2
cputotal = median(TOTALCPU).time * 1e-6
gputotal = median(TOTALGPU).time * 1e-6
times = vcat(times,scitimes)

# final timing results
println("\n\nFinal Timing Results")
println("--------------------")
println("GPU    Total Time: $(@sprintf("%.2f",gputotal)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",cputotal)) ms")
println("Python Total Time: $(@sprintf("%.2f",scitotal)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",cputotal/gputotal))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",scitotal/gputotal))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",scitotal/cputotal))x")

println("\n\nI/O Timing Results")
println("------------------")
seisiogpu = median(SEISIO).time * 1e-6
seisiocpu = median(SEISIO).time * 1e-6
obspyio = median(OBSPYIO).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",seisiogpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",seisiocpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",obspyio)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",seisiocpu/seisiogpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",obspyio/seisiogpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",obspyio/seisiocpu))x")

println("\n\nDetrend Results")
println("---------------")
detrendgpu = median(DETRENDGPU).time * 1e-6
detrendcpu = median(DETRENDCPU).time * 1e-6
detrendsci = median(DETRENDSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",detrendgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",detrendcpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",detrendsci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",detrendcpu/detrendgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",detrendsci/detrendgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",detrendsci/detrendcpu))x")

println("\n\nTaper Results")
println("-------------")
tapergpu = median(TAPERGPU).time * 1e-6
tapercpu = median(TAPERCPU).time * 1e-6
tapersci = median(TAPERSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",tapergpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",tapercpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",tapersci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",tapercpu/tapergpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",tapersci/tapergpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",tapersci/tapercpu))x")

println("\n\nFilter Results")
println("--------------")
filtergpu = median(FILTERGPU).time * 1e-6
filtercpu = median(FILTERCPU).time * 1e-6
filtersci = median(FILTERSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",filtergpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",filtercpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",filtersci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",filtercpu/filtergpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",filtersci/filtergpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",filtersci/filtercpu))x")

println("\n\nAmplitude Normalization Results")
println("-------------------------------")
onebitgpu = median(ONEBITGPU).time * 1e-6
onebitcpu = median(ONEBITCPU).time * 1e-6
onebitsci = median(ONEBITSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",onebitgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",onebitcpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",onebitsci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",onebitcpu/onebitgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",onebitsci/onebitgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",onebitsci/onebitcpu))x")

println("\n\nrFFT Results")
println("------------")
rfftgpu = median(RFFTGPU).time * 1e-6
rfftcpu = median(RFFTCPU).time * 1e-6
rfftsci = median(RFFTSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",rfftgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",rfftcpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",rfftsci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",rfftcpu/rfftgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",rfftsci/rfftgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",rfftsci/rfftcpu))x")

println("\n\nCross-Correlation Results")
println("-------------------------")
ccgpu = median(CCGPU).time * 1e-6
cccpu = median(CCCPU).time * 1e-6
ccsci = median(CCSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",ccgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",cccpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",ccsci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",cccpu/ccgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",ccsci/ccgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",ccsci/cccpu))x")

println("\n\nStacking Results")
println("----------------")
stackgpu = median(STACKGPU).time * 1e-6
stackcpu = median(STACKCPU).time * 1e-6
stacksci = median(STACKSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",stackgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",stackcpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",stacksci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",stackcpu/stackgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",stacksci/stackgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",stacksci/stackcpu))x")

println("\n\nPost-Processing Results")
println("-----------------------")
postgpu = median(POSTGPU).time * 1e-6
postcpu = median(POSTCPU).time * 1e-6
postsci = median(POSTSCI).time * 1e-6
println("GPU    Total Time: $(@sprintf("%.2f",postgpu)) ms")
println("CPU    Total Time: $(@sprintf("%.2f",postcpu)) ms")
println("Python Total Time: $(@sprintf("%.2f",postsci)) ms")
println("GPU vs CPU    Speedup: $(@sprintf("%.2f",postcpu/postgpu))x")
println("GPU vs Python Speedup: $(@sprintf("%.2f",postsci/postgpu))x")
println("CPU vs Python Speedup: $(@sprintf("%.2f",postsci/postcpu))x")

groupedbar(nam, times , group = ctg,xlabel = "Processing", ylabel = "Time [ms]",
        bar_width = 0.67, lw = 0, framestyle = :box,color=colors,legend=:best,
        size=(800,500),left_margin=10mm,bottom_margin=10mm,dpi=1000)
png(expanduser("~/SeisNoise_paper/CPU-GPU-BENCHMARK"))
