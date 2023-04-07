using DrWatson
@quickactivate "CompareTmTs"

using Dates
using NCDatasets
using TmPi
using ParseIGRAv2

include(srcdir("igrav2.jl"))

dtvec = collect(DateTime(1979,1,1,0,0,0) : Hour(1) : DateTime(2022,1,1,0,0,0))
pop!(dtvec)
stations = stationinfodata(derived=true)
lon = stations[:,4]; ilon = zeros(length(lon))
lat = stations[:,3]; ilat = zeros(length(lat))

lsd = getLandSea(ERA5Region("GLB",gres=0.25))
for istn = 1 : length(lon)
    ilon[istn] = argmin(abs.(lon[istn].-lsd.lon))
    ilat[istn] = argmin(abs.(lat[istn].-lsd.lat))
end

TmMat = zeros(length(dtvec),length(stations)) * NaN
Tmtmp = zeros(length(lsd.lon),length(lsd.lat),24*31)

for dt in Date(1979,1) : Month(1) : Date(2021,12)

    tmds = readTm(dt,path=datadir())
    iiTm = @view Tmtmp[:,:,24*daysinmonth(dt)]
    NCDatasets.load!(tmds["Tm"].var,iiTm,:,:,:)
    close(tmds)

    for istn = 1 : length(lon)

        ibeg = Dates.value(dt-Date(1979)) * 24 + 1
        iend = Dates.value(dt+Month(1)-Date(1979)) * 24

        stnTm = @view TmMat[ibeg:iend,istn]
        @view @. TmMat[ibeg:iend,istn] = iiTm[ilon[istn],ilat[istn],:]

    end

end

fnc = datadir("IGRAv2","era5Tm.nc")
if isfile(fnc) rm(fnc,force=true) end
ds = NCDataset(fnc,"c")

ds.dim["station"] = length(lon)
ds.dim["time"] = length(dtvec)
nctime = defVar(ds,"time",Int32,("time",),attrib = Dict(
    "units"     => "hours since $(dtvec[1]) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))
ncTm = defVar(ds,"Tm",Float64,("station","time",),attrib = Dict(
    "units"     => "K",
))
nctime.var[:] = collect(1:length(dtvec)) .- 1
ncTm.var[:] = TmMat
close(ds)