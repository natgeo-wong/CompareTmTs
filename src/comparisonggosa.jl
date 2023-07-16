using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis
using NCDatasets
using StatsBase

function comparisonggosa(
    e5ds :: ERA5Hourly
)

    egeo  = ERA5Region("GLB",resolution=0.25)
    lsd   = getLandSea(e5ds,egeo)
    dtbeg = e5ds.start
    dtend = e5ds.stop

    evar_Tm = SingleVariable("Tm")

    lon = lsd.lon; glon = length(lon); lon = lon[1:10:end]; nlon = length(lon)
    lat = lsd.lat; glat = length(lat); lat = lat[1:8:end];  nlat = length(lat)

    e5tmp_int = zeros(Int16,nlon,nlat,31*24)
    e5tmp_24h = zeros(nlon,nlat,31*24)
    e5tmp_flt = zeros(nlon,nlat,366*6)
    gatmp_flt = zeros(nlon,nlat,366*6)
    e5gad_flt = zeros(nlon,nlat)
    e5gae_flt = zeros(nlon,nlat)
    e5gan_flt = zeros(nlon,nlat)

    for idt = dtbeg : Year(1) : dtend

        iyr = year(idt)
        dtb = Date(iyr,1)
        
        for imo = 1 : 12

            dt = Date(iyr,imo)

            nhr = daysinmonth(dt) * 24
            n6h = daysinmonth(dt) * 4
            ibeg = Dates.value(dt-dtb) * 4 + 1
            iend = Dates.value(dt-dtb) * 4 + nhr

            e5tmp = @view e5tmp_int[:,:,1:nhr]
            e524h = @view e5tmp_24h[:,:,1:nhr]
            e5dt  = @view e5tmp_flt[:,:,ibeg:iend]

            @info "$(Dates.now()) - CompareTmTs - Loading ERA5 Reanalysis Data for $dt ..."

            eds = read(e5ds,evar_Tm,egeo,dt,quiet=true)

            sc = eds["Tm"].attrib["scale_factor"]
            of = eds["Tm"].attrib["add_offset"]
            mv = eds["Tm"].attrib["missing_value"]
            fv = eds["Tm"].attrib["_FillValue"]
            NCDatasets.load!(eds["Tm"].var,e5tmp,1:10:glon,1:8:glat,:)
            ERA5Reanalysis.int2real!(e524h,e5tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            for i6h = 1 : n6h, ilat = 1 : nlat, ilon = 1 : nlon
                e5dt[ilon,ilat,i6h] = mean(view(e524h,ilon,ilat,(1:6).+(i6h-1)*6))
            end

            close(eds)

        end

        nhr = daysinyear(idt) * 4
        gadt  = @view gatmp_flt[:,:,1:nhr]

        @info "$(Dates.now()) - CompareTmTs - Loading GGOS Atmosphere Data $idt ..."

        gds = NCDataset(datadir("ggosa","ggosa-$(year(idt)).nc"))
        NCDatasets.load!(gds["Tm"].var,gadt,:,:,:)
        close(gds)

        for ihr = 1 : nhr, ilat = 1 : nlat, ilon = 1 : nlon
            e5gad_flt[ilon,ilat] +=  e5tmp_flt[ilon,ilat,ihr] - gatmp_flt[ilon,ilat,ihr]
            e5gae_flt[ilon,ilat] += (e5tmp_flt[ilon,ilat,ihr] - gatmp_flt[ilon,ilat,ihr])^2
            e5gan_flt[ilon,ilat] += 1
        end

    end

    for ilat = 1 : nlat, ilon = 1 : nlon
        e5gad_flt[ilon,ilat] /= e5gan_flt[ilon,ilat]
        e5gae_flt[ilon,ilat] /= e5gan_flt[ilon,ilat]
        e5gae_flt[ilon,ilat]  = sqrt(e5gae_flt[ilon,ilat])
    end

    fnc = datadir("e5gadiff.nc")
    if isfile(fnc) rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ncΔ = defVar(ds,"Δ",Float64,("longitude","latitude"),attrib = Dict(
        "units"     => "K",
        "full_name" => "Difference between ERA5 and GGOSA"
    ))
    ncσ = defVar(ds,"σ",Float64,("longitude","latitude"),attrib = Dict(
        "units"     => "K",
        "full_name" => "Standard Deviation between ERA5 and GGOSA"
    ))
    ncΔ.var[:] = e5gad_flt
    ncσ.var[:] = e5gae_flt
    close(ds)

end