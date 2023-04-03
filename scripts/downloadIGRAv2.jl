using DrWatson
@quickactivate "CompareTmTs"

using ParseIGRAv2

for ID in stationlist(derived=true)
    stn = station(ID,derived=true)
    try
        download(stn,path=datadir())
        extract(stn,path=datadir())
    catch
    end
end