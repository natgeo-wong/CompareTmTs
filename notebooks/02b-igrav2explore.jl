### A Pluto.jl notebook ###
# v0.19.24

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
	using NCDatasets
	using ParseIGRAv2
	using Statistics
	
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

# ╔═╡ 37d0d82a-c735-429f-914b-993ff65cedeb
begin
	stndata = stationinfodata(derived=true)
	md"Loading station information ..."
end

# ╔═╡ bf8acaff-9f82-46ad-ae8d-17793a96a90e
begin
	ds = NCDataset(datadir("IGRAv2","Tm","compile-derived.nc"))
	stnTm = ds["Tm"][:]
	close(ds)
	ds = NCDataset(datadir("IGRAv2","era5Tm.nc"))
	eraTm = ds["Tm"][1:12:end,:]
	dt = ds["time"][1:12:end]
	close(ds)
	md"Loading the Tm data for ERA5 and IGRAv2 stations"
end

# ╔═╡ cf11d080-8ec3-4702-a971-cd3e05ffaf5b
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=2,axwidth=4,ncols=1)

	ii = 808
	idt = dt[.!isnan.(stnTm[:,ii])]
	a1[1].plot(dt,eraTm[:,ii])
	a1[1].plot(dt,stnTm[:,ii])
	a1[1].text(0.05,0.9,"lon = $(stndata[ii,4])",transform = a1[1].transAxes)
	a1[1].text(0.05,0.8,"lat = $(stndata[ii,3])",transform = a1[1].transAxes)
	a1[1].format(xlim=(minimum(idt),maximum(idt)))
	# a1[1].format(xlim=(Date(2002),Date(2003)))
	f1.savefig(plotsdir("02b-igrav2vsera5.png"),transparent=false,dpi=200)
	load(plotsdir("02b-igrav2vsera5.png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─b1bc72b6-da83-4c30-b368-94eb9180b6c5
# ╟─37d0d82a-c735-429f-914b-993ff65cedeb
# ╟─bf8acaff-9f82-46ad-ae8d-17793a96a90e
# ╠═cf11d080-8ec3-4702-a971-cd3e05ffaf5b
