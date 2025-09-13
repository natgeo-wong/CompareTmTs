### A Pluto.jl notebook ###
# v0.19.25

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
	using NCDatasets
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ bcd1a6c2-ae3d-11ed-067b-6b7828037b63
md"
# 01a. Exploring the TmPi Dataset
"

# ╔═╡ 4c55471e-1d8e-4fca-81fc-ad9dfee44c73
egeo = ERA5Region("GLB",gres=0.25)

# ╔═╡ 0f411cd8-30dc-4d00-a66f-b174eb0de1d0
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ 5cf63bdf-84e0-4b64-bd8f-c2efee783a88
begin
	coord = readdlm(datadir("GLB-l.txt"),comments=true,comment_char='#')
	clon = coord[:,1]; clat = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ 4fc01c2e-c1bc-46e0-b8cd-f56e51ec1e33
begin
	ds = NCDataset(datadir("mjophase.nc"))
	tm = ds["Tm"][:]
	close(ds)
end

# ╔═╡ f04272fa-fc91-450c-be17-84a89929607e
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=9,nrows=8,axwidth=6)

	lvls = -1:0.1:1
	
	c = axs[1].pcolormesh(lsd.lon,lsd.lat,tm[:,:,1]',levels=lvls,extend="both")
	axs[2].pcolormesh(lsd.lon,lsd.lat,tm[:,:,2]',levels=lvls,extend="both")
	axs[3].pcolormesh(lsd.lon,lsd.lat,tm[:,:,3]',levels=lvls,extend="both")
	axs[4].pcolormesh(lsd.lon,lsd.lat,tm[:,:,4]',levels=lvls,extend="both")
	axs[5].pcolormesh(lsd.lon,lsd.lat,tm[:,:,5]',levels=lvls,extend="both")
	axs[6].pcolormesh(lsd.lon,lsd.lat,tm[:,:,6]',levels=lvls,extend="both")
	axs[7].pcolormesh(lsd.lon,lsd.lat,tm[:,:,7]',levels=lvls,extend="both")
	axs[8].pcolormesh(lsd.lon,lsd.lat,tm[:,:,8]',levels=lvls,extend="both")

	for ax in axs
		ax.plot(clon,clat,c="k",lw=0.5)
		ax.format(ylim=(-20,20),xlim=(0,360),xlocator=0:60:360)
	end

	fig.colorbar(c,length=0.5)
	fig.savefig("test.png",transparent=false,dpi=300)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─4c55471e-1d8e-4fca-81fc-ad9dfee44c73
# ╟─0f411cd8-30dc-4d00-a66f-b174eb0de1d0
# ╟─5cf63bdf-84e0-4b64-bd8f-c2efee783a88
# ╟─4fc01c2e-c1bc-46e0-b8cd-f56e51ec1e33
# ╟─f04272fa-fc91-450c-be17-84a89929607e
