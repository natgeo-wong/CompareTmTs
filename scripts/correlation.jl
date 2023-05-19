using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis

include(srcdir("correlation.jl"))

e5ds = ERA5Hourly(start=Date(1979),stop=Date(1979,5),path=datadir())
egeo = GeoRegion("LAT1")
correlationTmTs(e5ds,ERA5Region(egeo,gres=0.25))