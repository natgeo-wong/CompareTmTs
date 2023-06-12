### A Pluto.jl notebook ###
# v0.19.22

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
evar = SingleVariable("t2m")

# ╔═╡ 420f34e4-3981-40c5-b2af-9fd031ba95e7
egeo = ERA5Region("GLB",gres=0.25)

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
		global ayr += yds["domain_yearly_mean_climatology"][:]
		global ast[:,:,iyr] = yds["domain_yearly_mean_climatology"][:]
		global amo += yds["domain_monthly_mean_climatology"][:]
		global dhr += yds["domain_yearly_mean_hourly"][:]
		global aitr+= yds["domain_monthly_maximum_climatology"][:] -
		              yds["domain_monthly_minimum_climatology"][:]
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
		nrows=2,ncols=3,aspect=1.9,axwidth=2,wspace=[3,1],
		proj="robin",proj_kw=Dict("lon_0"=>180)
	)
	
	c2 = a2[1].pcolormesh(lsd.lon,lsd.lat,ayr',extend="both",levels=280:1:300)
	a2[1].format(coast=true,title=L"$\mu$")
	a2[1].colorbar(c2,loc="l",locator=280:10:300)
	
	c2 = a2[4].pcolormesh(lsd.lon,lsd.lat,(amo.+aitr.+avar.+dhr)',extend="both",levels=10:10:90)
	a2[4].format(coast=true,title=L"$\delta_t$")
	a2[4].colorbar(c2,loc="l")
	
	c2 = a2[2].pcolormesh(lsd.lon,lsd.lat,amo',extend="both",levels=5:5:45)
	a2[2].format(coast=true,title=L"$\delta_s$")
	# a2[2].colorbar(c2)
	
	c2 = a2[3].pcolormesh(lsd.lon,lsd.lat,aitr',extend="both",levels=5:5:45)
	a2[3].format(coast=true,title=L"$\delta_i$")
	a2[3].colorbar(c2)
	
	c2 = a2[5].pcolormesh(lsd.lon,lsd.lat,avar',extend="both",levels=1:9)
	a2[5].format(coast=true,title=L"$\delta_a$")
	
	c2 = a2[6].pcolormesh(lsd.lon,lsd.lat,dhr',extend="both",levels=1:9)
	a2[6].format(coast=true,title=L"$\delta_d$")
	a2[6].colorbar(c2)
	
	f2.savefig(plotsdir("01c-TsStats_$(evar.varID).png"),transparent=false,dpi=200)
	load(plotsdir("01c-TsStats_$(evar.varID).png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─8ec33dcf-60ff-4a45-8e85-74601b8eb4bc
# ╠═faf7a70d-c07c-4519-b396-5fbc960416be
# ╟─420f34e4-3981-40c5-b2af-9fd031ba95e7
# ╟─5a690017-7dde-4feb-86e6-78f68725f4aa
# ╠═172bf4b3-4038-4f11-b98a-621e1eb1e811
# ╠═cacc5773-4891-4a34-ac66-662c623e1488
# ╟─5ab3a7af-a790-434a-a3cf-f6daf4de96a3
# ╟─4b089fcb-ede9-4eb1-b971-e79e42a5ed32
# ╟─c3f7b6a2-5326-4b26-902f-de33e4e46a73
