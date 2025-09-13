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
begin
	ds  = NCDataset(datadir("e5gadiff-smoothed.nc"))
	Δ   = ds["Δ"][:]
	σ   = ds["σ"][:]
	close(ds)
end

# ╔═╡ 717e4659-af68-4e49-a452-664c284af8c6
begin
	pplt.close()
	fig,axs = pplt.subplots(nrows=2,axwidth=4,aspect=2,proj="robin",proj_kw=Dict("lon_0"=>180))

	lvls = vcat(-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5)
	# lvls = vcat(-5:-1,-0.5,0.5,1:5)*2
	
	c1 = axs[1].pcolormesh(0:2.5:357.5,90:-2:-90,-Δ',levels=lvls,extend="both")
	c2 = axs[2].pcolormesh(0:2.5:357.5,90:-2:-90,σ',levels=1:0.5:5,extend="both")
	
	for ax in axs
		ax.format(coast=true)
	end

	axs[1].colorbar(c1)
	axs[2].colorbar(c2)
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─30525f38-0614-47f6-81c8-1ccb36e91e43
# ╟─c4aeba6c-22ad-11ee-1b25-b36ad2b696cb
# ╟─8ed3908d-8e00-4b9a-a6fb-53b525425308
# ╟─717e4659-af68-4e49-a452-664c284af8c6
