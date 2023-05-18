using DrWatson
@quickactivate "CompareTmTs"
using ERA5Reanalysis

addGeoRegions(srcdir("addLatRegions.txt"))

e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021,12),path=datadir())
egeo = ERA5Region("GLB",gres=0.25)
evar = SingleVariable("Tm")
vgeo = [
    GeoRegion("LAT1"),GeoRegion("LAT2"),GeoRegion("LAT3"),GeoRegion("LAT4"),
    GeoRegion("LAT5"),GeoRegion("LAT6"),GeoRegion("LAT7"),GeoRegion("LAT8"),
    GeoRegion("LAT9"),GeoRegion("LAT10"),GeoRegion("LAT11"),GeoRegion("LAT12"),
    GeoRegion("LAT13"),GeoRegion("LAT14"),GeoRegion("LAT15"),GeoRegion("LAT16"),
    GeoRegion("LAT17"),GeoRegion("LAT18")
]

extract(vgeo,e5ds,evar,egeo)