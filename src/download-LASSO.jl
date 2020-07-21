# download data from the LASSO array
using Dates, SeisIO, Printf, Glob

# download parameters
BASEDIR = "/media/four/"
DATADIR = joinpath(BASEDIR,"LASSO/DL/")
XMLDIR = joinpath(BASEDIR,"LASSO/XML")
s = DateTime(2016,5,1)
t = s + Day(1)
chans = ["2A.$ii..DPZ" for ii = 1:1850] # highest channel number

if !isdir(DATADIR)
    mkpath(DATADIR)
end

if !isdir(XMLDIR)
    mkpath(XMLDIR)
end

# download stations stations
cd(DATADIR)
for chan in chans
    try
        println("Downloading station $(@sprintf("%12s", chan)) $(now())")
        S = get_data("PH5",chan,s=s,t=t,w=true,autoname=true)
    catch err
        println(err)
    end
end

# check stations
rm("FDSNsta.xml")
files = glob("*",DATADIR)
stas = [split(s,'.')[8] |> (y -> parse(Int,y)) for s in files]
sort!(stas)

# grab station XML
for sta in stas
    println("Downloading station 2A.$sta $(now())")
    FDSNsta("2A.$sta..DPZ",s=s,t=t,xf=joinpath(XMLDIR,"2A.$sta.xml"),src="IRISPH5")
end
