# this script produces Figure 3
# GMT was installed using miniconda3
# If GMT is working,

using SeisNoise, Libdl, Glob, GMT
push!(Libdl.DL_LOAD_PATH, expanduser("~/miniconda3/lib/"))

CORRDIR = expanduser("~/SeisNoise_paper/SINGLECORR")
files = glob("*",CORRDIR)

# load first file and allocate
lat = zeros(Float64,188)
lon = zeros(Float64,188)

C = load_corr(files[1],"ZZ")
lon[1] = C.loc.lon
lat[1] = C.loc.lat
# loop through and plot
for ii = 2:188
    C = load_corr(files[ii],"ZZ")
    loc = get_loc(C.loc,C.azi,C.dist)
    lon[ii] = loc.lon
    lat[ii] = loc.lat
end

# create figure
coast(region=:d, proj=:Winkel, frame=:g, res=:low, area=10000, land=:seashell4,
                 water=:skyblue, figsize=18, shore=:black)
scatter!(lon,lat,fill=:gold,marker=:triangle,markeredgecolor=:black,fmt=:png,
         savefig=expanduser("~/SeisNoise_paper/stationmap.png"))
