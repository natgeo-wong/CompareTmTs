using DrWatson
@quickactivate "CompareTmTs"

using ParseIGRAv2

for ID in stationlist()
    download(station(ID),path=datadir())
end

for ID in stationlist(derived=true)
    download(station(ID,derived=true),path=datadir())
end