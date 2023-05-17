using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis
using NCDatasets
using StatsBase

function correlationTmTs(
    e5ds :: ERA5Dataset,
)

    lsd = getLandSea(e5ds,egeo)
    dtbeg = e5ds.start
    dtend = e5ds.stop
    ndt   = ((dtend + Month(1)) - dtbeg) * 24

    evar_Ts = SingleVariable("t2m")
    evar_Tm = SingleVariable("Tm")

    tstmp_int = zeros(Int16,length(lsd.lon),31*24)
    tmtmp_int = zeros(Int16,length(lsd.lon),31*24)
    tstmp_flt = zeros(length(lsd.lon),ndt)
    tmtmp_flt = zeros(length(lsd.lon),ndt)
    tstm_corr = zeros(length(lsd.lon),length(lsd.lat))

    for ilat = 1 : length(lsd.lat)

        for idt = dtbeg : Month(1) : dtend

            nhr = daysinmonth(idt) * 24
            ibeg = Dates.value(idt-dtbeg) * 24 + 1
            iend = Dates.value(idt-dtbeg) * 24 + nhr

            tstmp = view(tstmp_int,:,1:nhr)
            tmtmp = view(tmtmp_int,:,1:nhr)
            tsdt  = view(tstmp_flt,:,ibeg:iend)
            tmdt  = view(tmtmp_flt,:,ibeg:iend)

            tsds = read(e5ds,evar_Ts,egeo,idt)
            tmds = read(e5ds,evar_Tm,egeo,idt)

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
            tstmp_ilat = view(tstmp_flt,ilat,:)
            tstm_corr[ilon,ilat] = corr(tstmp_ilon,tstmp_ilat)

        end

    end

    fnc = datadir("tmtscorr.nc")
    if isfile(fnc) rm(fnc,force=true) end
    ds = NCDataset(fnc,"c")

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)
    ncrho = defVar(ds,"rho",Float64,("longitude","latitude","phase"),attrib = Dict(
        "units"     => "0-1",
        "full_name" => "Correlation between Ts and Tm"
    ))
    ncrho.var[:] = tstm_corr
    close(ds)

end

