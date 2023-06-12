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
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ bcd1a6c2-ae3d-11ed-067b-6b7828037b63
md"
# 01c. Basic Ts Statistics
"

# ╔═╡ 8ec33dcf-60ff-4a45-8e85-74601b8eb4bc
e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021),path=datadir())

# ╔═╡ faf7a70d-c07c-4519-b396-5fbc960416be
evar = SingleVariable("sst")

# ╔═╡ 420f34e4-3981-40c5-b2af-9fd031ba95e7
egeo = ERA5Region("TRP",gres=0.25)

# ╔═╡ 5a690017-7dde-4feb-86e6-78f68725f4aa
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ 172bf4b3-4038-4f11-b98a-621e1eb1e811
yrbeg = 1979

# ╔═╡ cacc5773-4891-4a34-ac66-662c623e1488
yrend = 2021

# ╔═╡ 5ab3a7af-a790-434a-a3cf-f6daf4de96a3
begin
	ayr = zeros(length(lsd.lon),length(lsd.lat))
	ast = zeros(length(lsd.lon),length(lsd.lat),yrend+1-yrbeg)
	amo = zeros(length(lsd.lon),length(lsd.lat),12)
	aitr= zeros(length(lsd.lon),length(lsd.lat),12)
	dhr = zeros(length(lsd.lon),length(lsd.lat),24)

	iyr = 0

	for yr = yrbeg : yrend

		global iyr += 1
		yds = read(e5ds,evar,egeo,Date(yr),analysis=true)
		global ayr += nomissing(yds["domain_yearly_mean_climatology"][:],NaN)
		global ast[:,:,iyr] = nomissing(yds["domain_yearly_mean_climatology"][:],NaN)
		global amo += nomissing(yds["domain_monthly_mean_climatology"][:],NaN)
		global dhr += nomissing(yds["domain_yearly_mean_hourly"][:],NaN)
		global aitr+= nomissing(yds["domain_monthly_maximum_climatology"][:],NaN) -
		              nomissing(yds["domain_monthly_minimum_climatology"][:],NaN)
		close(yds)
		
	end

	ayr /= (yrend+1-yrbeg)
	amo /= (yrend+1-yrbeg)
	aitr/= (yrend+1-yrbeg)
	dhr /= (yrend+1-yrbeg)
	
	aitr = dropdims(mean(aitr,dims=3),dims=3)
	amo  = dropdims(maximum(amo,dims=3)  .- minimum(amo,dims=3),dims=3)
	dhr  = dropdims(maximum(dhr,dims=3)  .- minimum(dhr,dims=3),dims=3)
	md"Compiling Statistics ..."
end

# ╔═╡ 4b089fcb-ede9-4eb1-b971-e79e42a5ed32
begin
	lrg = ones(yrend+1-yrbeg,2)
	lrg[:,2] .= collect(1979:2021)
	ydtrend = zeros(yrend+1-yrbeg)
	avar = zeros(length(lsd.lon),length(lsd.lat))
	for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
		iitm = @view ast[ilon,ilat,:]
		iitm = iitm .- lrg * (lrg \ iitm)
		avar[ilon,ilat] = maximum(iitm) .- minimum(iitm)
	end
end

# ╔═╡ c3f7b6a2-5326-4b26-902f-de33e4e46a73
begin
	pplt.close(); f2,a2 = pplt.subplots(
		aspect=1.9,axwidth=4,
		proj="robin",proj_kw=Dict("lon_0"=>180)
	)
	
	a2[1].contour(lsd.lon,lsd.lat,ayr'.-273.15,levels=[29,35,36],c="k",linestyle="--",lw=1)
	a2[1].format(coast=true)
	
	f2.savefig(plotsdir("01d-sstwarmpool.png"),transparent=false,dpi=200)
	load(plotsdir("01d-sstwarmpool.png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─8ec33dcf-60ff-4a45-8e85-74601b8eb4bc
# ╟─faf7a70d-c07c-4519-b396-5fbc960416be
# ╟─420f34e4-3981-40c5-b2af-9fd031ba95e7
# ╟─5a690017-7dde-4feb-86e6-78f68725f4aa
# ╠═172bf4b3-4038-4f11-b98a-621e1eb1e811
# ╠═cacc5773-4891-4a34-ac66-662c623e1488
# ╟─5ab3a7af-a790-434a-a3cf-f6daf4de96a3
# ╟─4b089fcb-ede9-4eb1-b971-e79e42a5ed32
# ╠═c3f7b6a2-5326-4b26-902f-de33e4e46a73
