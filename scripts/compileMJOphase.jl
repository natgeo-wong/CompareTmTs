using DrWatson
@quickactivate "CompareTmTs"

using DelimitedFiles
using ERA5Reanalysis

data = readdlm(datadir("PCs_pp_rot.txt"),skipstart=1,',')

date = Date.(data[:,1])
yr   = year.(date)
mo   = month.(date)
pc1  = data[:,2]
pc2  = data[:,3]

function pc2mjophase(pc1::Real,pc2::Real)

    if sqrt(pc1^2+pc2^2) >= 1

        if abs(pc1) >= abs(pc2)
            if pc1 >= 0
                if pc2 >= 0
                    return 3
                else
                    return 2
                end
            else
                if pc2 >= 0
                    return 6
                else
                    return 7
                end
            end
        else
            if pc1 >= 0
                if pc2 >= 0
                    return 4
                else
                    return 1
                end
            else
                if pc2 >= 0
                    return 5
                else
                    return 8
                end
            end
        end

    else

        return NaN
        
    end
        

end

phase = pc2mjophase.(pc1,pc2)

e5ds = ERA5Daily(start=Date(1979),stop=Date(2018,7),path=datadir())
egeo = ERA5Region(GeoRegion("GLB"),gres=0.25)
evar = SingleVariable("Tm")
lsd  = getLandSea(e5ds,egeo)
nlon = length(lsd.lon)
nlat = length(lsd.lat)

dtvec = e5ds.start : Month(1) : e5ds.stop
nmjo  = zeros(nlon,nlat,8)
pnum  = zeros(8)

for idt in dtvec

    iyr = yr .== year(idt)
    imo = mo .== month(idt)
    ii  = iyr .& imo
    iip = phase[ii]

    ds = read(e5ds,evar,egeo,idt)
    Tm = nomissing(ds["Tm"][:],NaN)
    close(ds)
    
    for iphase = 1 : 8

        iphase_ii = iip .== iphase
        if !iszero(sum(iphase_ii))

            nmjo[:,:,iphase] += sum(Tm[:,:,iphase_ii],dims=3)
            pnum[iphase] += sum(iphase_ii)

        end

    end

    Tm = nothing

end

pnum = reshape(pnum,1,1,8)
nmjo_phase = nmjo ./ pnum
nmjo_avg   = dropdims(sum(nmjo,dims=3),dims=3) / sum(pnum)

fnc = datadir("mjophase.nc")
if isfile(fnc) rm(fnc,force=true) end
ds = NCDataset(fnc,"c")

ds.dim["longitude"] = nlon
ds.dim["latitude"]  = nlat
ds.dim["phase"]     = 8
ncTm = defVar(ds,"Tm",Float64,("longitude","latitude","phase"),attrib = Dict(
    "units"     => "K",
    "full_name" => "Tm compiled by MJO Phase"
))
ncavg = defVar(ds,"Tm_avg",Float64,("longitude","latitude","phase"),attrib = Dict(
    "units"     => "K",
    "full_name" => "Average Tm"
))
ncTm.var[:]  = nmjo_phase .- nmjo_avg
ncavg.var[:] = nmjo_avg
close(ds)