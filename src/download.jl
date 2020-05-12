using Dates, CSV, DataFrames, SeisIO, SeisNoise, Glob

# This file downloads data in a cleanway from FDSN
startdate = DateTime(2019,1,1)
enddate = DateTime(2019,12,31)
dfmt = Dates.DateFormat("yyyy.ddd")
date_range = startdate:Day(1):enddate
N2download = 187
STATIONOUT = expanduser("~/SEISNOISE-TEST/DATA/BESTSTA")
YEARSEIS = expanduser("~/SEISNOISE-TEST/DATA/YEARSEIS")
stationcsv = expanduser("~/SEISNOISE-TEST/DATA/gmap-stations.csv")
df = CSV.read(stationcsv,delim="|")

if !isdir(STATIONOUT)
    mkpath(STATIONOUT)
end

if !isdir(YEARSEIS)
    mkpath(YEARSEIS)
end

# subset by time
df = df[df[:StartTime] .< startdate,:]
df = df[df[:EndTime] .> enddate,:]

# get information about channels
channels = Array{String}(undef,size(df,1))
netsta = Array{String}(undef,size(df,1))
for ii = 1:size(df,1)
    channels[ii] = join([df[ii,:Network],df[ii,:Station],df[ii,"Location"],"LHZ"],'.')
    netsta[ii] = join([df[ii,:Network],df[ii,:Station]],'.')
end
df[:NetSta] = netsta
STA = FDSNsta(channels,s=startdate,t=enddate)

# remove stations from df that are not in STA
goodnetsta = Array{String}(undef,STA.n)
for ii = 1:STA.n
    nslc = STA[ii].id
    goodnetsta[ii] = join(split(nslc,'.')[1:2],'.')
end

netstaind = findall(in(goodnetsta),df[!,:NetSta])
df = df[netstaind,:]
channels = channels[netstaind]

# calculate distance between all stations
Ngood = length(channels)
dists = zeros(Ngood,Ngood)
for ii = 1:Ngood-1
    for jj = ii+1:Ngood
        dists[ii,jj] = SeisNoise.surface_distance(df[ii,:Longitude],df[ii,:Latitude],
                                                  df[jj,:Longitude],df[jj,:Latitude],6378.137)
    end
end
for ii = 1:Ngood
    dists[ii,ii] = Inf
end

# select best stations for download
# start with IU & II networks
beststa = Array{String}(undef,length(channels))
bestind = zeros(Int,length(channels))
IU = findall(df[:Network] .== "IU")
beststa[1:length(IU)] = channels[IU]
bestind[1:length(IU)] = IU
II = findall(df[:Network] .== "II")
beststa[length(IU)+1:length(IU)+length(II)] = channels[II]
bestind[length(IU)+1:length(IU)+length(II)] = II

# now find stations which are the furthest from current stations
for ii = length(IU)+length(II) + 1 : length(beststa)
    outchannels = setdiff(1:length(channels),bestind[1:ii-1])
    outdist = zeros(length(outchannels))
    for jj = 1:length(outchannels)
        outdist[jj] = minimum(dists[bestind[1:ii-1],outchannels[jj]])
    end
    newbest = argmax(outdist)
    println(maximum(outdist))
    beststa[ii] = channels[outchannels[newbest]]
    bestind[ii] = outchannels[newbest]
end

# download in order
stationlist = channels[bestind]
outfile = expanduser("~/SEISNOISE-TEST/DATA/channelorder.txt")
open(outfile, "w") do f
  for i in stationlist
    println(f, i)
  end
end

# download with obspy due to error with SeisIO mass download
# obspy file is download.py

# find number of files in each dir -> want 365
dirs = glob("*",STATIONOUT)
numfiles = zeros(Int,length(dirs))
hasfirst = Array{Bool}(undef,length(dirs))
for ii = 1:length(numfiles)
    files = readdir(dirs[ii])
    numfiles[ii] = length(files)
    if length(files) > 0
        hasfirst[ii] = occursin("20190101",files[1])
    else
        hasfirst[ii] = false
    end
end

# find directories with > 360 files
ind = findall((numfiles .>= 359) .& hasfirst)

# function to read all files in
function readall(INDIR::String)
    S = SeisData()
    read_data!(S,"mseed",INDIR)
    demean!(S)
    ungap!(S)
    SeisIO.detrend!(S)
    taper!(S,t_max=1000.)
    phase_shift!(S[1], ϕshift=true)

    # check number of seconds
    demean!(S)
    sync!(S,s=startdate)
    return S
end

# convert to yearseis
for ii = 1:188
    netsta = join(split(dirs[ind[ii]],'.')[1:2],'.')
    mseeddir = joinpath(STATIONOUT,netsta)
    S = readall(mseeddir)
    fileout = joinpath(YEARSEIS,S.id[1] * ".seisio")
    println("Writing file $fileout")
    wseis(fileout,S)
end
