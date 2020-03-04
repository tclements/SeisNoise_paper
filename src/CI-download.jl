# This will download data all data in the scedc-pds bucket on 2020-01-01
# for all HHZ stations in the CI network.
# Files will be downloaded to the ~/data directory.

using Distributed
addprocs()
@everywhere begin
    using SeisNoise, Dates
    startdate = Date(2020,1,1)
    network = "CI"
    channel = "HHZ"
    OUTDIR = "~/data"
end

scedctransfer(OUTDIR, startdate, network=network,channel=channel)
