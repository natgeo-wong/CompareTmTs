### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 30525f38-0614-47f6-81c8-1ccb36e91e43
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ c4aeba6c-22ad-11ee-1b25-b36ad2b696cb
begin
	@quickactivate "CompareTmTs"
	using ERA5Reanalysis
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ 8ed3908d-8e00-4b9a-a6fb-53b525425308
e5ds = ERA5Hourly(start=Date(2020),stop=Date(2020,12),path=datadir())

# ╔═╡ ed2faf78-6c4b-4660-9cdd-540da3084920
evar = SingleVariable("Tm")

# ╔═╡ b79ee1fa-b434-4257-8c34-2983daf27a85
egeo = ERA5Region("GLB",resolution=0.25)

# ╔═╡ 415af975-774d-4a56-ac34-4c4cf1e4adb6
begin
	eds = read(e5ds,evar,egeo,Date(2020))
	lon_e5 = eds["longitude"][:]
	lat_e5 = eds["latitude"][:]
	Tm_e5  = eds[evar.ID][:]
	close(eds)
end

# ╔═╡ f772ccab-89cb-426f-96d2-0170c3fa3c01
begin
	ds  = NCDataset(datadir("merra2Tm","merra2Tm-20200101.nc"))
	lon_m2 = ds["longitude"][:] .- 180
	lat_m2 = ds["latitude"][:]
	Tm_m2  = ds["Tm"][:]
	close(ds)
end

# ╔═╡ dc5d9398-2a75-41db-bb0b-c2425d900997
ggrd = RegionGrid(egeo.geo,lon_m2,lat_m2)

# ╔═╡ 2307d910-9010-4f15-96e7-8931b72e3a84
begin
	Tm = extractGrid(Tm_m2,ggrd)
	md"extractGrid"
end

# ╔═╡ 717e4659-af68-4e49-a452-664c284af8c6
begin
	pplt.close()
	fig,axs = pplt.subplots(axwidth=4,aspect=2,proj="robin",proj_kw=Dict("lon_0"=>180))
	
	# axs[1].pcolormesh(ggrd.lon,ggrd.lat,Tm[:,:,1]',levels=220:5:295)
	# axs[2].pcolormesh(lon_e5,lat_e5,dropdims(mean(Tm_e5[:,:,1:3],dims=3),dims=3)',levels=220:5:295)
	c = axs[1].pcolormesh(
		ggrd.lon[1:2:end],ggrd.lat,(Tm[1:2:end,:,1] .- dropdims(mean(Tm_e5[1:5:end,721:-2:1,1:3],dims=3),dims=3))',
		levels=vcat(-5:-1,-0.5,0.5,1:5),extend="both"
	)
	for ax in axs
		ax.format(coast=true)
	end

	fig.colorbar(c)
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─30525f38-0614-47f6-81c8-1ccb36e91e43
# ╟─c4aeba6c-22ad-11ee-1b25-b36ad2b696cb
# ╟─8ed3908d-8e00-4b9a-a6fb-53b525425308
# ╟─ed2faf78-6c4b-4660-9cdd-540da3084920
# ╟─b79ee1fa-b434-4257-8c34-2983daf27a85
# ╟─415af975-774d-4a56-ac34-4c4cf1e4adb6
# ╟─f772ccab-89cb-426f-96d2-0170c3fa3c01
# ╟─dc5d9398-2a75-41db-bb0b-c2425d900997
# ╟─2307d910-9010-4f15-96e7-8931b72e3a84
# ╟─717e4659-af68-4e49-a452-664c284af8c6
