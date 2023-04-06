using DelimitedFiles
using NumericalIntegration
using Statistics

function calcTd2e(T::Real)

    if T >= 273.16
    	return 611.21 * exp(17.502 * (Td-273.16) / (Td-32.19))
    elseif T <= 250.16
    	return 611.21 * exp(22.587 * (Td-273.16) / (Td+0.7))
    else
        α  = ((T - 250.16) / (273.16 - 250.16))^2
        ei = 611.21 * exp(22.587 * (T-273.16) / (T+0.7))
        ew = 611.21 * exp(17.502 * (T-273.16) / (T-32.19))
        return α * ew + (1-α) * ei
    end

end

calce2q(e::Real,p::Real) = e * 0.621981 / (p - 0.378019 * e)

function IGRA2Tm(station :: IGRAv2DataRaw)

    pres = profile_pressure(station)
    tair = profile_temperature(station) .+ 273.15
    rhum = profile_relativehumidity(station) ./ 100

    ntime = length(station.nlevels)
    tm = zeros(ntime)
    qair = zeros(maximum(station.nlevels)+1)
    btm  = zeros(maximum(station.nlevels)+1)
    ipres = zeros(maximum(station.nlevels)+1)
    itair = zeros(maximum(station.nlevels)+1)
    irhum = zeros(maximum(station.nlevels)+1)

    for it = 1 : ntime

        @views @. ipres[1:station.nlevels[it]] = pres[1:station.nlevels[it],it]
        @views @. itair[1:station.nlevels[it]] = tair[1:station.nlevels[it],it]
        @views @. irhum[1:station.nlevels[it]] = rhum[1:station.nlevels[it],it]

        ipres[(station.nlevels[it])+1] = 0
        itair[(station.nlevels[it])+1] = 0
        irhum[(station.nlevels[it])+1] = 0

        ix = @view ipres[1:(station.nlevels[it]+1)]

        for ilvl = 1 : station.nlevels[it]
            qair[ilvl] = calce2q(calcTd2e(itair[ilvl]),ipres[ilvl]) * irhum
            btm[ilvl] = qair[ilvl] ./ itair[ilvl]
        end
        iqair = @view qair[1:(station.nlevels[it]+1)]; iqair[end] = 0
        ibtm  = @view btm[1:(station.nlevels[it]+1)];  ibtm[end] = 0

        ind = .!isnan.(iqair) .& .!isnan.(ibtm)
        
        if sum(ind) > 3 # Sounding must have more than 3 valid data points of T and q
            tm[it] = integrate(view(ix,ind),view(iqair,ind)) /
                     integrate(view(ix,ind),view(ibtm ,ind))
        else
            tm[it] = NaN
        end

    end

    tmfol = joinpath(station.path,"IGRAv2","Tm","raw")
    if !isdir(tmfol); mkpath(tmfol) end
    tmID = joinpath(tmfol,"$(station.ID).txt")
    open(tmID, "w") do io
        writedlm(io, [station.dates tm], ',')
    end

end

function IGRA2Tm(station :: IGRAv2DataDerived)

    pres = profile_pressure(station)
    tair = profile_temperature(station)
    vapp = profile_vapourpressure(station)

    ntime = length(station.nlevels)
    tm = zeros(ntime)
    qair = zeros(maximum(station.nlevels)+1)
    btm  = zeros(maximum(station.nlevels)+1)
    ipres = zeros(maximum(station.nlevels)+1)
    itair = zeros(maximum(station.nlevels)+1)
    ivapp = zeros(maximum(station.nlevels)+1)

    for it = 1 : ntime

        @views @. ipres[1:station.nlevels[it]] = pres[1:station.nlevels[it],it]
        @views @. itair[1:station.nlevels[it]] = tair[1:station.nlevels[it],it]
        @views @. ivapp[1:station.nlevels[it]] = vapp[1:station.nlevels[it],it]

        ipres[(station.nlevels[it])+1] = 0
        itair[(station.nlevels[it])+1] = 0
        ivapp[(station.nlevels[it])+1] = 0

        ix = @view ipres[1:(station.nlevels[it]+1)]

        for ilvl = 1 : station.nlevels[it]
            qair[ilvl] = ivapp[ilvl] * 0.621981 / ipres[ilvl]
            btm[ilvl] = qair[ilvl] ./ itair[ilvl]
        end
        iqair = @view qair[1:(station.nlevels[it]+1)]; iqair[end] = 0
        ibtm  = @view btm[1:(station.nlevels[it]+1)];  ibtm[end] = 0

        ind = .!isnan.(iqair) .& .!isnan.(ibtm)
        
        if sum(ind) > 3 # Sounding must have at least 3 valid data points of T and q
            tm[it] = integrate(view(ix,ind),view(iqair,ind)) /
                     integrate(view(ix,ind),view(ibtm ,ind))
        else
            tm[it] = NaN
        end

    end

    tmfol = joinpath(station.path,"IGRAv2","Tm","derived")
    if !isdir(tmfol); mkpath(tmfol) end
    tmID = joinpath(tmfol,"$(station.ID).txt")
    open(tmID, "w") do io
        writedlm(io, [station.dates tm], ',')
    end

end

function readTmfromIGRAv2(station :: IGRAv2StationDerived)

    fID = datadir("IGRAv2","Tm","$(station.ID).txt")
    if isfile(fID)

        data = readdlm(fID,',')
        return DateTime.(data[:,1]), data[:,2]

    else

         error("$(now()) - CompareTmTs - $(station.ID) does not have a derived file to calculate Tm data for")

    end

end