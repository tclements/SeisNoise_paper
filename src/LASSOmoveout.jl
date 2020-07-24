using Glob, PyCall, SeisNoise, Serialization
# read files
CORRDIR = "/media/four/LASSO/CORR/"
# do them all
files = glob("*/*",CORRDIR)
bins = collect(0.:0.05:55)
bins1 = collect(0.05:0.05:55.05)
counts = zeros(length(bins))
freqmin = 1.
freqmax = 2.
C = deserialize(files[1])
Cmat = zeros(eltype(Cs[1].corr),size(C.corr,1),length(bins))
for ii = 1:length(files)
    C = deserialize(files[ii])
    ind = findall((C.dist .>= bins) .& (C.dist .<= bins1))[1]
    clean_up!(C,freqmin,freqmax)
    abs_max!(C)
    Cmat[:,ind] .+= C.corr[:]
    counts[ind] += 1
    println("Reading file $ii")
end

maxdist = 25.
ind = findall((bins .> 0.15) .& (bins .<= maxdist))
Cmat = Cmat[:,ind]
counts = counts[ind]
bins = bins[ind]
Cmat ./= counts'
lags = -C.maxlag:1/C.fs:C.maxlag
maxlag = 15.
lagind = findall(abs.(lags) .<= maxlag)
Cmat = Cmat[lagind,:]
abs_max!(Cmat)
# plot envolopes at 5 km/s and 1.3 km/s
heatmap(
    lags[lagind],
    bins,
    Cmat',
    c=:balance,
    xlabel="Lag [s]",
    ylabel="Inter-station Distance [km]",
    legend = :none,
    dpi=500,
)
plot!([0.,5.],[0.1,25],line=(1.5,:dash),color=:black,label="")
plot!([0.,-5.],[0.1,25],line=(1.5,:dash),color=:black,label="")
plot!([0.,15.],[0.03,22.5],line=(1.5,:dash),color=:black,label="")
plot!([0.,-15.],[0.03,22.5],line=(1.5,:dash),color=:black,label="")
annotate!(2.5,17, Plots.text("5 km/s", 14, :dark, rotation = 78 ))
annotate!(8,10, Plots.text("1.5 km/s", 14, :dark, rotation = 48 ))
ylims!((minimum(bins),maximum(bins)))
xlims!((-maxlag,maxlag))
png(expanduser("~/SeisNoise_paper/LASSO-moveout"))
