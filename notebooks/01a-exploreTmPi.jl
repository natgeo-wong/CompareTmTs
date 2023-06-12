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
	using ERA5Reanalysis
	using TmPi
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ bcd1a6c2-ae3d-11ed-067b-6b7828037b63
md"
# 01a. Exploring the TmPi Dataset
"

# ╔═╡ 272b9a90-f85f-4817-91d9-b29e700c0cf1
evar_tm = SingleVariable("Tm")

# ╔═╡ ad28ba6e-1f4b-4bfc-8bcf-48c804936028
begin
	rds = readTm(Date(1981,9),path=datadir())
	lon = rds["longitude"][:]
	lat = rds["latitude"][:]
	Tm  = rds["Tm"][:]
	close(rds)
end

# ╔═╡ e10f8163-60b8-4039-8a9c-d43a0ada98a9
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=2,axwidth=3)
	
	c1 = a1[1].pcolormesh(lon,lat,Tm[:,:,313]',extend="both",levels=270:2:290)
	a1[1].scatter(-110.333+360,31.58,s=1)
	a1[1].format(coast=true,xlim=(-115,-105).+360,ylim=(27.5,32.5))
	f1.colorbar(c1)
	f1.savefig(plotsdir("01a-exploretmpi.png"),transparent=false,dpi=200)
	load(plotsdir("01a-exploretmpi.png"))
end

# ╔═╡ 5ab3a7af-a790-434a-a3cf-f6daf4de96a3
begin
	ads = readTm(Date(1995),path=datadir(),analysis=true)
	ana = ads["longitude"][:]
	ana = ads["latitude"][:]
	ana = ads["domain_yearly_mean_climatology"][:]
	close(ads)
end

# ╔═╡ c3f7b6a2-5326-4b26-902f-de33e4e46a73
begin
	pplt.close(); f2,a2 = pplt.subplots(proj="wintri",proj_kw=Dict("lon_0"=>180),aspect=1.7,axwidth=3)
	
	c2 = a2[1].pcolormesh(lon,lat,ana',extend="both",levels=270:2:290)
	a2[1].format(coast=true)
	f2.colorbar(c2)
	f2.savefig(plotsdir("01a-exploretmpi_anatest.png"),transparent=false,dpi=200)
	load(plotsdir("01a-exploretmpi_anatest.png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─272b9a90-f85f-4817-91d9-b29e700c0cf1
# ╠═ad28ba6e-1f4b-4bfc-8bcf-48c804936028
# ╠═e10f8163-60b8-4039-8a9c-d43a0ada98a9
# ╟─5ab3a7af-a790-434a-a3cf-f6daf4de96a3
# ╟─c3f7b6a2-5326-4b26-902f-de33e4e46a73
