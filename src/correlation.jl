using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis
using NCDatasets
using StatsBase

function correlationTmTs(
    e5ds :: ERA5Dataset,
    egeo :: ERA5Region
)

    lsd = getLandSea(e5ds,egeo)
    dtbeg = e5ds.start
    dtend = e5ds.stop
    ndt   = Dates.value((dtend + Day(1)) - dtbeg) * 24

    evar_Ts = SingleVariable("t2m")
    evar_Tm = SingleVariable("Tm")

    nlon = Int(length(lsd.lon)/4);       nnlon = length(lsd.lon)
    nlat = Int(ceil(length(lsd.lat)/4)); nnlat = length(lsd.lat)
    tstmp_int = zeros(Int16,nlon,nlat,31*24)
    tmtmp_int = zeros(Int16,nlon,nlat,31*24)
    tstmp_flt = zeros(nlon,nlat,ndt)
    tmtmp_flt = zeros(nlon,nlat,ndt)
    tstm_corr = zeros(nlon,nlat)

    for idt = dtbeg : Month(1) : dtend

        nhr = daysinmonth(idt) * 24
        ibeg = Dates.value(idt-dtbeg) * 24 + 1
        iend = Dates.value(idt-dtbeg) * 24 + nhr

        tstmp = @view tstmp_int[:,:,1:nhr]
        tmtmp = @view tmtmp_int[:,:,1:nhr]
        tsdt  = @view tstmp_flt[:,:,ibeg:iend]
        tmdt  = @view tmtmp_flt[:,:,ibeg:iend]

        @info "$(Dates.now()) - CompareTmTs - Loading data for $idt ..."

        tsds = read(e5ds,evar_Ts,egeo,idt,quiet=true)
        tmds = read(e5ds,evar_Tm,egeo,idt,quiet=true)

        sc = tsds["t2m"].attrib["scale_factor"]
        of = tsds["t2m"].attrib["add_offset"]
        mv = tsds["t2m"].attrib["missing_value"]
        fv = tsds["t2m"].attrib["_FillValue"]
        NCDatasets.load!(tsds["t2m"].var,tstmp,1:4:nnlon,1:4:nnlat,:)
        ERA5Reanalysis.int2real!(tsdt,tstmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = tmds["Tm"].attrib["scale_factor"]
        of = tmds["Tm"].attrib["add_offset"]
        mv = tmds["Tm"].attrib["missing_value"]
        fv = tmds["Tm"].attrib["_FillValue"]
        NCDatasets.load!(tmds["Tm"].var,tmtmp,1:4:nnlon,1:4:nnlat,:)
        ERA5Reanalysis.int2real!(tmdt,tmtmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        close(tsds)
        close(tmds)

    end

    for ilat = 1 : nlat, ilon = 1 : nlon
        tstmp_ilon = view(tstmp_flt,ilon,ilat,:)
        tmtmp_ilon = view(tmtmp_flt,ilon,ilat,:)
        tstm_corr[ilon,ilat] = cor(tstmp_ilon,tmtmp_ilon)
    end

    fnc = datadir("tmtscorr_$(egeo.ID).nc")
    if isfile(fnc) rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = Int(length(lsd.lon)/4)
    ds.dim["latitude"]  = Int(ceil(length(lsd.lat)/4))
    ncrho = defVar(ds,"rho",Float64,("longitude","latitude"),attrib = Dict(
        "units"     => "0-1",
        "full_name" => "Correlation between Ts and Tm"
    ))
    ncrho.var[:] = tstm_corr
    close(ds)

end

function compilecorrelationTmTs()

    corr_TsTm = zeros(360,181)

    for ilat = 1 : 18
        ibeg = (ilat-1) * 10 + 1
        iend =  ilat    * 10 + 1
        icorr = @view corr_TsTm[:,ibeg:iend,:]
        ds = NCDataset(datadir("tmtscorr_LAT$ilat.nc"))
        icorr .= ds["rho"][:]
        close(ds)
    end

    fnc = datadir("tmtscorr.nc")
    if isfile(fnc) rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = size(corr_TsTm,1)
    ds.dim["latitude"]  = size(corr_TsTm,2)
    ncrho = defVar(ds,"rho",Float64,("longitude","latitude"),attrib = Dict(
        "units"     => "0-1",
        "full_name" => "Correlation between Ts and Tm"
    ))
    ncrho.var[:] = corr_TsTm
    close(ds)

end