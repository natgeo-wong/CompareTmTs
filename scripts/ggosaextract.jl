using Dates
using DrWatson
using Glob
using Logging
using NCDatasets

function ggosaextract(yr :: Int)

    flist = glob("*.h*",datadir("ggosa","rawtext","t$(yr)"))
    data = zeros(145*91,daysinyear(yr)*4)
    ii = 0

    for fID in flist

        @info "$(now()) - CompareTmTs - $fID"
        ii   += 1
        tdata = read(fID,String)
        tdata = replace(tdata,"\n"=>" ")
        tdata = split(tdata," ")
        tdata = tdata[.!(tdata.=="")][7:end]
        tdata = parse.(Float64,tdata)
        data[:,ii] .= tdata

    end

    data = reshape(data,145,91,:)[1:144,:,:]

    fnc = datadir("ggosa","ggosa-$(yr).nc")
    if isfile(fnc); rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = 144
    ds.dim["latitude"]  = 91
    ds.dim["time"]      = size(data,3)

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctime = defVar(ds,"time",Int32,("time",),attrib = Dict(
        "units"     => "hours since $(Date(yr)) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    ncvar = defVar(ds,"Tm",Float64,("longitude","latitude","time"),attrib = Dict(
        "long_name"     => "water_vapour_weighted_mean_temperature",
        "full_name"     => "Water Vapour Weighted Mean Temperature",
        "units"         => "K",
    ))

    nclon[:]  = collect(0 : 2.5 : 360)[1:(end-1)]
    nclat[:]  = collect(-90 : 2.0 : 90)
    nctime[:] = (collect(1:size(data,3)) .- 1) * 6
    ncvar[:]  = data

    close(ds)

end