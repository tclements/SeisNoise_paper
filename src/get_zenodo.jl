# this files downloads files for reproducing from Zenodo
using HTTP, DelimitedFiles
# this script downloads file

zenodo = expanduser("~/SeisNoise_paper/DATA/zenodo.txt")
files = readdlm(zenodo,String)[:]
zenodourl = "https://zenodo.org/record/3823283/files/"

# make output directory
if !isdir(expanduser("~/SeisNoise_paper/DATA/YEARSAC"))
    mkpath(expanduser("~/SeisNoise_paper/DATA/YEARSAC"))
end

# download files from zenodo
for (ii,f) in enumerate(files)
    dlpath = joinpath(zenodourl,f * "?download=1")
    outpath = joinpath(expanduser("~/SeisNoise_paper/DATA/YEARSAC"), f)
    println("Downloading file $f, $ii of $(length(files))")
    HTTP.download(dlpath,outpath)
end
