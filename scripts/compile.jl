using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis
using NCDatasets
using StatsBase

function compile(
    e5ds :: ERA5Dataset,
    evar :: SingleLevel,
    egeo :: ERA5Region
)

    lsd = getLandSea(egeo,path=datadir("emask"))
    yrbeg = year(e5ds.start)
    yrend = year(e5ds.stop)

    ayr  = zeros(length(lsd.lon),length(lsd.lat))
    ast  = zeros(length(lsd.lon),length(lsd.lat),yrend+1-yrbeg)
    amo  = zeros(length(lsd.lon),length(lsd.lat),12)
    aitr = zeros(length(lsd.lon),length(lsd.lat),12)
    dhr  = zeros(length(lsd.lon),length(lsd.lat),24)

    iyr = 0

    for yr = yrbeg : yrend

        iyr += 1

        yds = read(e5ds,evar,egeo,Date(yr),analysis=true)
        ayr  += yds["domain_yearly_mean_climatology"][:]
        ast[:,:,iyr] = yds["domain_yearly_mean_climatology"][:]
        amo  += yds["domain_monthly_mean_climatology"][:]
        dhr  += yds["domain_yearly_mean_hourly"][:]
        aitr += yds["domain_monthly_maximum_climatology"][:] -
                        yds["domain_monthly_minimum_climatology"][:]
        close(yds)
        
    end

    ayr /= (yrend+1-yrbeg)
    amo /= (yrend+1-yrbeg)
    aitr/= (yrend+1-yrbeg)
    dhr /= (yrend+1-yrbeg)

    aitr = dropdims(mean(aitr,dims=3),dims=3)
    amo  = dropdims(maximum(amo,dims=3)  .- minimum(amo,dims=3),dims=3)
    dhr  = dropdims(maximum(dhr,dims=3)  .- minimum(dhr,dims=3),dims=3)

    lrg = ones(yrend+1-yrbeg,2)
    lrg[:,2] .= collect(1979:2021)
    avar = zeros(length(lsd.lon),length(lsd.lat))
    for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
        iitm = @view ast[ilon,ilat,:]
        iitm = iitm .- lrg * (lrg \ iitm)
        avar[ilon,ilat] = std(iitm) * 2
    end

    fnc = datadir("compile_$(evar.varID).nc"); if isfile(fnc); rm(fnc,force=true) end
    
    ds = NCDataset(fTm,"c",attrib = Dict(
        "history"     => "Created on $(Dates.now())",
        "comments"    => "Created by the CompareTmTs project, compiled Tm statistics from 1979 to 2021"
    ))

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncμ = defVar(ds,"mean",Float64,("longitude","latitude"))
    ncδd = defVar(ds,"variability_diurnal",Float64,("longitude","latitude"))
    ncδi = defVar(ds,"variability_intraseasonal",Float64,("longitude","latitude"))
    ncδs = defVar(ds,"variability_seasonal",Float64,("longitude","latitude"))
    ncδa = defVar(ds,"variability_interannual",Float64,("longitude","latitude"))
    
    nclon[:] = lsd.lon
    nclat[:] = lsd.lat
    ncμ[:]   = ayr
    ncδd[:]  = dhr
    ncδi[:]  = aitr
    ncδs[:]  = amo
    ncδa[:]  = avar

    close(ds)

end

compile(
    ERA5Hourly(start=Date(1979),stop=Date(2021),path=datadir()),
    SingleVariable("Tm"),
    ERA5Region("GLB",gres=0.25)
)