using DrWatson
@quickactivate "CompareTmTs"

using ParseIGRAv2

for ID in stationlist()
    stn = station(ID)
    try
        download(stn,path=datadir(),derived=true)
        extract(stn,path=datadir(),derived=true)
    catch
    end
end