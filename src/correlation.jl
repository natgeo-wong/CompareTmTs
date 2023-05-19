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

    tstmp_int = zeros(Int16,length(lsd.lon),31*24)
    tmtmp_int = zeros(Int16,length(lsd.lon),31*24)
    tstmp_flt = zeros(length(lsd.lon),ndt)
    tmtmp_flt = zeros(length(lsd.lon),ndt)
    tstm_corr = zeros(length(lsd.lon),length(lsd.lat))

    for ilat = 1 : length(lsd.lat)

        @info "$(Dates.now()) - CompareTmTs - Performing correlation for latitude $ilat of $(length(lsd.lat))"

        for idt = dtbeg : Month(1) : dtend

            nhr = daysinmonth(idt) * 24
            ibeg = Dates.value(idt-dtbeg) * 24 + 1
            iend = Dates.value(idt-dtbeg) * 24 + nhr

            tstmp = @view tstmp_int[:,1:nhr]
            tmtmp = @view tmtmp_int[:,1:nhr]
            tsdt  = @view tstmp_flt[:,ibeg:iend]
            tmdt  = @view tmtmp_flt[:,ibeg:iend]

            tsds = read(e5ds,evar_Ts,egeo,idt,quiet=true)
            tmds = read(e5ds,evar_Tm,egeo,idt,quiet=true)

            sc = tsds["t2m"].attrib["scale_factor"]
            of = tsds["t2m"].attrib["add_offset"]
            mv = tsds["t2m"].attrib["missing_value"]
            fv = tsds["t2m"].attrib["_FillValue"]
            NCDatasets.load!(tsds["t2m"].var,tstmp,:,ilat,:)
            ERA5Reanalysis.int2real!(tsdt,tstmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            sc = tmds["Tm"].attrib["scale_factor"]
            of = tmds["Tm"].attrib["add_offset"]
            mv = tmds["Tm"].attrib["missing_value"]
            fv = tmds["Tm"].attrib["_FillValue"]
            NCDatasets.load!(tmds["Tm"].var,tmtmp,:,ilat,:)
            ERA5Reanalysis.int2real!(tmdt,tmtmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            close(tsds)
            close(tmds)

        end

        for ilon = 1 : length(lsd.lon)

            tstmp_ilon = view(tstmp_flt,ilon,:)
            tmtmp_ilon = view(tmtmp_flt,ilon,:)
            tstm_corr[ilon,ilat] = cor(tstmp_ilon,tmtmp_ilon)

        end

    end

    fnc = datadir("tmtscorr_$(egeo.geoID).nc")
    if isfile(fnc) rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)
    ncrho = defVar(ds,"rho",Float64,("longitude","latitude"),attrib = Dict(
        "units"     => "0-1",
        "full_name" => "Correlation between Ts and Tm"
    ))
    ncrho.var[:] = tstm_corr
    close(ds)

end

function compilecorrelationTmTs()

    corr_TsTm = zeros(1440,721)

    for ilat = 1 : 18
        ibeg = (ilat-1) * 40 + 1
        iend =  ilat    * 40 + 1
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