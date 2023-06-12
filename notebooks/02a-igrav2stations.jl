### A Pluto.jl notebook ###
# v0.19.23

using Markdown
using InteractiveUtils

# ╔═╡ fc2fdbeb-f775-48af-9d93-82d8a79567a9
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ed931d06-d36a-4570-b4f3-c6eaea16ff59
begin
	@quickactivate "CompareTmTs"
	using DelimitedFiles
	using ParseIGRAv2
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ bcd1a6c2-ae3d-11ed-067b-6b7828037b63
md"
# 02a. IGRAv2 Stations
"

# ╔═╡ b1bc72b6-da83-4c30-b368-94eb9180b6c5
begin
	coasts = readdlm(datadir("GLB-l.txt"),comments=true,comment_char='#')
	x = coasts[:,1]
	y = coasts[:,2]
	md"Loading coastlines ..."
end

# ╔═╡ 6e0a5139-c17f-4d6c-9818-f55b07129442
begin
	stn_raw = stationinfodata(derived=false)[:,3:7]
	stn_drv = stationinfodata(derived=true)[:,3:7]
	md"Loading IGRAv2 station information ..."
end

# ╔═╡ a0dfdb6c-77dd-4464-9621-181c6fd294a4
begin
	lon_bin = -180 : 5 : 180
	lat_bin =  -90 : 5 : 90
	ind1 = (stn_raw[:,5] .>= 1979)
	ind2 = (stn_drv[:,5] .>= 1979)
	h_raw = fit(Histogram,(stn_raw[:,2],stn_raw[:,1]),(lon_bin,lat_bin))
	h_drv = fit(Histogram,(stn_drv[:,2],stn_drv[:,1]),(lon_bin,lat_bin))
	h_val = fit(Histogram,(stn_drv[ind2,2],stn_drv[ind2,1]),(lon_bin,lat_bin))
	md"Binning the density of stations by gridpoint ..."
end

# ╔═╡ ab7ed432-ab10-4938-ac94-905e2ced2b62
begin
	lat_den = cosd.(collect(-87.5 : 5 : 87.5))
	area = 2 * pi * 6.370 * lat_den * (6.370*2*pi/360*5)
end

# ╔═╡ 992f5acd-813d-4af8-b755-46244f3370ba
begin
	pplt.close(); f1,a1 = pplt.subplots(proj="moll",aspect=1.9,axwidth=2.5,ncols=3)

	a1[1].pcolormesh(lon_bin,lat_bin,h_raw.weights'./area,cmap="Blues",levels=vcat(0:10)*5e-2,extend="both")
	a1[2].pcolormesh(lon_bin,lat_bin,h_drv.weights'./area,cmap="Blues",levels=vcat(0:10)*5e-2,extend="both")
	c = a1[3].pcolormesh(lon_bin,lat_bin,h_val.weights'./area,cmap="Blues",levels=vcat(0:10)*5e-2,extend="both")
	
	# a1[1].scatter(stn_raw[:,2],stn_raw[:,1],s=2)
	# a1[1].scatter(stn_drv[:,2],stn_drv[:,1],s=2)
	# a1[2].scatter(stn_raw[ind1,2],stn_raw[ind1,1],s=2)
	# a1[2].scatter(stn_drv[ind2,2],stn_drv[ind2,1],s=2)
	a1[1].format(title="All Stations")
	a1[2].format(title="Stations with Derived Data")
	a1[3].format(title="Stations (Derived) ending after 1979")

	for ax in a1
		ax.format(coast=true)
	end

	f1.colorbar(c,label=L"Stations per 10$^6$ km$^2$")
	
	f1.savefig(plotsdir("02a-igrav2stations.png"),transparent=false,dpi=200)
	load(plotsdir("02a-igrav2stations.png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─b1bc72b6-da83-4c30-b368-94eb9180b6c5
# ╟─6e0a5139-c17f-4d6c-9818-f55b07129442
# ╟─a0dfdb6c-77dd-4464-9621-181c6fd294a4
# ╟─ab7ed432-ab10-4938-ac94-905e2ced2b62
# ╠═992f5acd-813d-4af8-b755-46244f3370ba
