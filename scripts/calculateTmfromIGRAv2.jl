using DrWatson
@quickactivate "CompareTmTs"

using ParseIGRAv2

include(srcdir("igrav2.jl"))

for ID in stationlist()
    IGRA2Tm(read(station(ID),path=datadir(),derived=true))
end