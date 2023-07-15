using Dates
using Logging
using NCDatasets
using Statistics
using Trapz

function merraTm(;
    start :: Date,
    stop  :: Date
)

    tatmp = zeros(Float32,576,361,72)
    qvtmp = zeros(Float32,576,361,72)
    patmp = zeros(Float32,576,361,72)
    pstmp = zeros(Float32,576,361,3)
    tstmp = zeros(Float32,576,361,3)
    qstmp = zeros(Float32,576,361,3)

    ind = zeros(Bool,74)
    bot = zeros(Float64,74)
    ita = zeros(Float64,74)
    iqv = zeros(Float64,74)
    ipa = zeros(Float64,74)
    Tm = zeros(576,361,8)

    for idt in start : Day(1) : stop

        @info "$(now()) - Calculating MERRA-2 Tm for $idt ..."

        dtstr =  Dates.format(idt,dateformat"yyyymmdd")

        for ihr = 1 : 8

            ds  = NCDataset(datadir("merra2pre","merra2pre-$(dtstr).nc"))

            NCDatasets.load!(ds["T"].var,tatmp,:,:,:,ihr)
            NCDatasets.load!(ds["QV"].var,qvtmp,:,:,:,ihr)
            NCDatasets.load!(ds["PL"].var,patmp,:,:,:,ihr)

            close(ds)

            ds  = NCDataset(datadir("merra2sfc","merra2sfc-$(dtstr).nc"))

            NCDatasets.load!(ds["T2M"].var,tstmp,:,:,(ihr-1)*3 .+ 1:3)
            NCDatasets.load!(ds["QV2M"].var,qstmp,:,:,(ihr-1)*3 .+ 1:3)
            NCDatasets.load!(ds["PS"].var,pstmp,:,:,(ihr-1)*3 .+ 1:3)

            close(ds)

            for ilat = 1 : 361, ilon = 1 : 576
                
                ips      = mean(view(pstmp,ilon,ilat,:))
                ipa[end] = ips
                ita[end] = mean(view(tstmp,ilon,ilat,:))
                iqv[end] = mean(view(qstmp,ilon,ilat,:))


                for ip = 2 : 73
                    ipa[ip] = patmp[ilon,ilat,ip-1]
                    ita[ip] = tatmp[ilon,ilat,ip-1]
                    iqv[ip] = qvtmp[ilon,ilat,ip-1]
                end

                for ip = 2 : 74
                    bot[ip] = iqv[ip] / ita[ip]
                end

                for ip = 1 : 74
                    ind[ip] = ipa[ip] < ips
                end
                ind[end] = true

                top = @view iqv[ind]
                btm = @view bot[ind]
                ipp = @view ipa[ind]

                Tm[ilon,ilat,ihr] = trapz(ipp,top) / trapz(ipp,btm)

            end

        end

        mkpath(datadir("merra2Tm"))
        fnc = datadir("merra2Tm","merra2Tm-$(dtstr).nc")
        if isfile(fnc); rm(fnc,force=true) end
        ds = NCDataset(fnc,"c",attrib = Dict(
            "Conventions" => "CF-1.6",
        ))

        ds.dim["longitude"] = 576
        ds.dim["latitude"]  = 361
        ds.dim["time"] = 8

        nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
            "units"     => "degrees_east",
            "long_name" => "longitude",
        ))

        nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
            "units"     => "degrees_north",
            "long_name" => "latitude",
        ))

        nctime = defVar(ds,"time",Int32,("time",),attrib = Dict(
            "units"     => "hours since $(dt) 01:30:00.0",
            "long_name" => "time",
            "calendar"  => "gregorian",
        ))

        ncvar = defVar(ds,"Tm",Float64,("longitude","latitude","time"),attrib = Dict(
            "long_name"     => "water_vapour_weighted_mean_temperature",
            "full_name"     => "Water Vapour Weighted Mean Temperature",
            "units"         => "K",
        ))

        nclon[:] = collect(0:0.625:360)[1:(end-1)]
        nclat[:] = collect(-90:0.5:90)
        nctime[:] = collect(0:7) * 3
        ncvar[:] = Tm

        close(ds)

    end

end