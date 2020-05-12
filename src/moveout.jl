using SeisNoise, Glob, Plots, Statistics, Plots.PlotMeasures
# plot all correlations
CORRDIR = expanduser("~/SEISNOISE-TEST/SINGLECORR")
files = glob("*",CORRDIR)
freqmin = 1 / 50
freqmax = 1 / 30
max_length = 1000.

# load first file and allocate
C = load_corr(files[10000],"ZZ")
t = -C.maxlag:C.maxlag
corr = Array{eltype(C.corr)}(undef,size(C.corr,1),length(files))
dist = zeros(Float64,length(files))
azi = zeros(Float64,length(files))
baz = zeros(Float64,length(files))
lat = zeros(Float64,length(files))
lon = zeros(Float64,length(files))

# loop through and plot
for ii = 1:length(files)
    C = load_corr(files[ii],"ZZ")
    clean_up!(C,freqmin,freqmax,max_length=max_length)
    println(ii," ",C.name)
    corr[:,ii] .= C.corr[:,1]
    dist[ii] = C.dist
    azi[ii] = C.azi
    baz[ii] = C.baz
    lon[ii] = C.loc.lon
    lat[ii] = C.loc.lat
end

# sort by cross-correlation by distance
distind = sortperm(dist)
corr .= corr[:,distind]
dist .= dist[distind]
azi .= azi[distind]
baz .= baz[distind]
lon .= lon[distind]
lat .= lat[distind]

### plotting ###
# just find the highest snr traces
maxamp = maximum(abs.(corr),dims=1)
maxmad = mad(corr)
snr = (maxamp ./ maxmad)[:]
nonan = findall(.!isnan.(snr[:]))
snr = snr[nonan]
corr = corr[:,nonan]
dist = dist[nonan]
azi = azi[nonan]
baz = baz[nonan]
lon = lon[nonan]
lat = lat[nonan]

scaling = 10000 # scale waveforms by 1000
lw = 0.5  # linewidth
al = 0.1  # alpha
mindist = 100. # minimum distance in km for plot
maxdist = 10000.
distave = 10.
dists = collect(mindist:distave:maxdist)
maxt = 3000
tind = findall(abs.(t) .< maxt)
distplot = zeros(eltype(C.corr),length(tind),length(dists))
for ii = 1:length(dists)
    cind = findall((dist .> dists[ii]) .& (dist .< dists[ii] + distave))
    c = mean(corr[tind,cind] ./ (maximum(abs.(corr[tind,cind]),dims=1) .+ 1),dims=2)
    distplot[:,ii] = c ./ maximum(abs.(c))
end
heatmap(dists,t[tind],distplot,color=:RdBu, xlabel="Interstation Distance [km]",ylabel="Lag [s]")

p = plot()
for ii = 1:length(dists)
    cind = findall((dist .> dists[ii]) .& (dist .< dists[ii] + distave))
    c = mean(corr[tind,cind] ./ (maximum(abs.(corr[tind,cind]),dims=1) .+ 1),dims=2)
    # c = corr[tind,snrind[ii]] ./ maximum(abs.(corr[tind,snrind[ii]]),dims=1)
    plot!(p,t[tind],dists[ii] .+ c .* scaling,color=:black,linewidth=lw,alpha=al,
         label="",xlabel = "Lag [s]",ylabel="Interstation Distance [km]",
         xlims=(-maxt,maxt),ylims=(mindist,maxdist),dpi=500,left_margin=-60px,
                bottom_margin=-20px,right_margin=25px)
end
p
png(expanduser("~/SEISNOISE-TEST/LHZ-moveout"))
