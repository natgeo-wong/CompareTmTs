using DrWatson
@quickactivate "CompareTmTs"

using Dates
using ParseIGRAv2

include(srcdir("igrav2.jl"))

dtvec = collect(DateTime(1979,1,1,0,0,0) : Hour(12) : DateTime(2022,1,1,0,0,0))
pop!(dtvec)
stations = stationlist(derived=true)

TmMat = zeros(length(dtvec),length(stations)) * NaN

for istn in 1 : length(stations)
    date,Tm = readTmfromIGRAv2(station(stations[istn],derived=true))
    ndt = length(date)
    for idt in 1 : ndt

        ii = findfirst(dtvec.==date[idt])
        if !isnan(Tm[idt]) && !isnothing(ii)
            TmMat[ii,istn] = Tm[idt]
        end

    end
end