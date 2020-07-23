# this script produces Figure 4

using CSV, DataFrames, DelimitedFiles, Glob, GMT, Statistics
stats = expanduser("~/SeisNoise_paper/DATA/gmap-stations.txt")
df = DataFrame!(CSV.File(stats,delim="|"))
lat,lon = df[:,:Latitude], df[:,:Longitude]
latmin, latmax = minimum(lat), maximum(lat)
lonmin, lonmax = minimum(lon), maximum(lon)
latpad = 0.025
lonpad = 0.05

# create figure
coast(
    region=(lonmin-lonpad,lonmax+lonpad,latmin-latpad,latmax+latpad),
    proj=:Mercator,
    res=:full,
    land=:white,
    water=:green,
    # scale=0.03,
    par=(:MAP_FRAME_TYPE,"fancy+"),
    map_scale="jBR+c$(mean(lat))+w10k+f+o2.5/2.5+at+u",
    figsize=15,
)
scatter!(lon,lat,fill=:gold,marker=:triangle,markeredgecolor=:black,markersize=0.2)
basemap!(inset=(anchor=:TL, width=0.25, offset=(1, 5), save="xx000"))
t = readdlm("xx000")
coast!(region=[-107 -91 26 41], proj=:merc,
           land=:lightgray, area=5000, shore=:faint,
           x_off=t[1], y_off=t[2],
           N="a",
           figsize=4,
           frame=:bare,
)
scatter!(
     [mean(lon)],
     [mean(lat)],
     marker=:star,
     markeredgecolor=:black,
     fill=:gold,
     markersize=0.5,
     savefig=expanduser("~/SeisNoise_paper/LASSOmap.png"),
)
