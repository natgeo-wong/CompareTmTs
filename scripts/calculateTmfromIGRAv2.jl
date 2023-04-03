using DrWatson
@quickactivate "CompareTmTs"

using ParseIGRAv2

include(srcdir("igrav2.jl"))

for ID in stationlist(derived=true)
    IGRA2Tm(read(station(ID,derived=true),path=datadir()))
end