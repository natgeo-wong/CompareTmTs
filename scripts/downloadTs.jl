using DrWatson
@quickactivate "CompareTmTs"

using ERA5Reanalysis

e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021),path=datadir())
egeo = ERA5Region(GeoRegion("GLB"),gres=0.25)
evar = SingleVariable("skt"); download(e5ds,evar,egeo)
evar = SingleVariable("t2m"); download(e5ds,evar,egeo)