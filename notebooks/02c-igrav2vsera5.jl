### A Pluto.jl notebook ###
# v0.19.26

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
	using GeoRegions
	using NCDatasets
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
	close(ds)
	Tmdiff = eraTm .- stnTm
	md"Calculating difference between ERA5 and IGRAv2 stations"
end

# ╔═╡ d7e9548e-cffd-47e8-b340-8a4db06e1118
begin
	stdTm = zeros(1535)
	for ii = 1 : 1535
		iiTm = @view Tmdiff[:,ii]
		iiTm = @view iiTm[iiTm.!==NaN]
		if length(iiTm) > (5*365)
			stdTm[ii] = std(iiTm)
		else
			stdTm[ii] = NaN
		end
	end
end

# ╔═╡ 484dacd4-7960-4e25-92dd-5452fe178d39
md"
### A. Calculating the Global Statistics for Comparison
"

# ╔═╡ cf11d080-8ec3-4702-a971-cd3e05ffaf5b
begin
	pplt.close(); f1,a1 = pplt.subplots(
		proj="moll",#proj_kw=Dict("lon_0"=>180),
		aspect=2,axwidth=4,ncols=1
	)

	ii = .!isnan.(stdTm)
	c = a1[1].scatter(
		stndata[ii,4],stndata[ii,3],c=stdTm[ii],
		s=2.5,levels=vcat(0.5:0.5:2.5,5:9),extend="both",
		cmap="balance",#cmap_kw=Dict("gamma"=>0.7)
	)

	for ax in a1
		ax.format(coast=true,ultitle="(a)")
	end

	f1.colorbar(c,label=L"$\sigma_{T_m}$ / K",length=0.7)
	
	f1.savefig(plotsdir("02c-igrav2vsera5.png"),transparent=false,dpi=400)
	load(plotsdir("02c-igrav2vsera5.png"))
end

# ╔═╡ 6004dd44-7379-4a5d-9b0f-7f66cf74d209
begin
	h = fit(Histogram, (stnTm[:],eraTm[:]),(200:325,200:325))
	a = Float64.(h.weights)
	a[iszero.(h.weights)] .= NaN
	md"Binning the data into a 2D histogram ..."
end

# ╔═╡ bee9301a-42d4-4382-b13f-ce72b9c2c87e
begin
	pplt.close(); f2,a2 = pplt.subplots()
	
	c2 = a2[1].pcolormesh(
		200:325,200:325,a',
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",
		levels=vcat(0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000)
	)
	a2[1].contour(
		200:325,200:325,a',levels=vcat(10,100,1000,10000,100000),
		c="k",lw=0.75,linestyle="--"
	)
	a2[1].plot([200,350],[200,350],c="k",lw=0.5)
	a2[1].format(
		xlim=(200,325),xlocator=200:25:325,xlabel=L"IGRAv2 Derived $T_m$",
		ylim=(200,325),ylocator=200:25:325,ylabel=L"ERA5 Derived $T_m$",
	)

	f2.colorbar(c2,label="Number of Data Points")
	f2.savefig(plotsdir("02c-igrav2era5scatter.png"),transparent=false,dpi=400)
	load(plotsdir("02c-igrav2era5scatter.png"))
end

# ╔═╡ 6ed17596-5be9-44aa-9e1f-5c77e9bde4c8
md"
### B. An Analysis by Geographic Region
"

# ╔═╡ c53bc1b8-3bd7-49ba-896a-1fc9356d1133
begin
	geoNA = [
		GeoRegion("AR6_NWN"),GeoRegion("AR6_NEN"),GeoRegion("AR6_WNA"),
		GeoRegion("AR6_CNA"),GeoRegion("AR6_ENA"),
		GeoRegion("AR6_NCA"),GeoRegion("AR6_SCA"),GeoRegion("AR6_CAR"),
	]
	geoSA = [
		GeoRegion("AR6_NWS"),GeoRegion("AR6_NSA"),GeoRegion("AR6_SAM"),
		GeoRegion("AR6_NES"),GeoRegion("AR6_SEA"),
		GeoRegion("AR6_SWS"),GeoRegion("AR6_SES"),GeoRegion("AR6_SSA"),
	]
	geoEU = [
		GeoRegion("AR6_NEU"),GeoRegion("AR6_WCE"),GeoRegion("AR6_MED"),
		GeoRegion("AR6_EEU"),
	]
	geoAS = [
		GeoRegion("AR6_RAR"),GeoRegion("AR6_WSB"),GeoRegion("AR6_ESB"),
		GeoRegion("AR6_RFE"),GeoRegion("AR6_ECA"),GeoRegion("AR6_TIB"),
		GeoRegion("AR6_SAS"),GeoRegion("AR6_EAS"),GeoRegion("AR6_SEA"),
	]
	geoAU = [
		GeoRegion("AR6_NAU"),GeoRegion("AR6_CAU"),GeoRegion("AR6_EAU"),
		GeoRegion("AR6_SAU"),GeoRegion("AR6_NZ")
	]
	geoAF = [
		GeoRegion("AR6_SAH"),GeoRegion("AR6_WAF"),GeoRegion("AR6_CAF"),
		GeoRegion("AR6_NEAF"),GeoRegion("AR6_SEAF"),GeoRegion("AR6_WSAF"),
		GeoRegion("AR6_ESAF"),GeoRegion("AR6_MDG"),GeoRegion("AR6_MED"),
	]
end

# ╔═╡ 0db7e996-5443-45b1-86cc-c0711c747bf6
function continent2stn(geovec,stndata)

	stnvec = zeros(length(stndata[:,1]))

	for igeo = 1 : length(geovec)

		geo = geovec[igeo]
		for istn = 1 : length(stnvec)
			stnvec[istn] += isinGeoRegion(Point2.(stndata[istn,4],stndata[istn,3]),geo,throw=false)
		end
		
	end

	return .!iszero.(stnvec)
	
end

# ╔═╡ afa9059c-f81d-49f8-b11d-62834a8982d8
begin
	stn_NA = continent2stn(geoNA,stndata)
	stn_SA = continent2stn(geoSA,stndata)
	stn_EU = continent2stn(geoEU,stndata)
	stn_AF = continent2stn(geoAF,stndata)
	stn_AS = continent2stn(geoAS,stndata)
	stn_AU = continent2stn(geoAU,stndata)
end

# ╔═╡ 5276c12e-3673-4a70-8842-4cdbde9f28bf
begin
	h_NA = fit(Histogram, (stnTm[:,stn_NA][:],eraTm[:,stn_NA][:]),(200:325,200:325))
	h_SA = fit(Histogram, (stnTm[:,stn_SA][:],eraTm[:,stn_SA][:]),(200:325,200:325))
	h_EU = fit(Histogram, (stnTm[:,stn_EU][:],eraTm[:,stn_EU][:]),(200:325,200:325))
	h_AF = fit(Histogram, (stnTm[:,stn_AF][:],eraTm[:,stn_AF][:]),(200:325,200:325))
	h_AS = fit(Histogram, (stnTm[:,stn_AS][:],eraTm[:,stn_AS][:]),(200:325,200:325))
	h_AU = fit(Histogram, (stnTm[:,stn_AU][:],eraTm[:,stn_AU][:]),(200:325,200:325))
end

# ╔═╡ 8d616b0a-3c80-47aa-bc8a-b8bd1215ec4b
function normalizeweights(raw)

	weights = raw / sum(raw)
	weights[iszero.(weights)] .= NaN

	return weights
	
end

# ╔═╡ 89f897d4-0b9b-4a18-a201-a9afb61772d9
begin
	weights_NA = normalizeweights(h_NA.weights)
	weights_SA = normalizeweights(h_SA.weights)
	weights_EU = normalizeweights(h_EU.weights)
	weights_AF = normalizeweights(h_AF.weights)
	weights_AS = normalizeweights(h_AS.weights)
	weights_AU = normalizeweights(h_AU.weights)
	md"Normalization of the respective weights ..."
end

# ╔═╡ decf7b15-cfe0-44b8-b906-9b190ba37568
begin
	pplt.close(); f3,a3 = pplt.subplots(ncols=3,nrows=2,axwidth=1.5,hspace=1,wspace=1)

	lvls = vcat(0.1,1,2,5,10,20,50,100,200)
	# lvls = 20:20:200
	
	c3 = a3[1].pcolormesh(
		200:325,200:325,weights_NA'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="both",levels=lvls
	)
	a3[1].contour(
		200:325,200:325,weights_NA'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[2].pcolormesh(
		200:325,200:325,weights_SA'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",levels=lvls
	)
	a3[2].contour(
		200:325,200:325,weights_SA'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[3].pcolormesh(
		200:325,200:325,weights_EU'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",levels=lvls
	)
	a3[3].contour(
		200:325,200:325,weights_EU'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[4].pcolormesh(
		200:325,200:325,weights_AS'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",levels=lvls
	)
	a3[4].contour(
		200:325,200:325,weights_AS'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[5].pcolormesh(
		200:325,200:325,weights_AF'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",levels=lvls
	)
	a3[5].contour(
		200:325,200:325,weights_AF'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[6].pcolormesh(
		200:325,200:325,weights_AU'*125^2,
		cmap="blues",cmap_kw=Dict("left"=>0.05),extend="max",levels=lvls
	)
	a3[6].contour(
		200:325,200:325,weights_AU'*125^2,levels=vcat(0.1,1,10,100),
		c="k",lw=0.75,linestyle="--"
	)

	a3[1].format(ultitle="(b) North America")
	a3[2].format(ultitle="(c) South America")
	a3[3].format(ultitle="(d) Europe")
	a3[4].format(ultitle="(e) Asia")
	a3[5].format(ultitle="(f) Africa")
	a3[6].format(ultitle="(g) Oceania")

	for ax in a3
		ax.plot([200,350],[200,350],c="k",a=0.2,lw=0.5,linestyle="--")
		ax.format(
			xlim=(200,325),xlocator=225:25:300,xlabel=L"IGRAv2 Derived $T_m$",
			ylim=(200,325),ylocator=225:25:300,ylabel=L"ERA5 Derived $T_m$",
		)
	end

	f3.colorbar(c3,label="Probability Density",length=0.7)
	f3.savefig(plotsdir("02c-igrav2era5scatter-region.png"),transparent=false,dpi=400)
	load(plotsdir("02c-igrav2era5scatter-region.png"))
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╠═ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─b1bc72b6-da83-4c30-b368-94eb9180b6c5
# ╟─37d0d82a-c735-429f-914b-993ff65cedeb
# ╟─bf8acaff-9f82-46ad-ae8d-17793a96a90e
# ╟─d7e9548e-cffd-47e8-b340-8a4db06e1118
# ╟─484dacd4-7960-4e25-92dd-5452fe178d39
# ╟─cf11d080-8ec3-4702-a971-cd3e05ffaf5b
# ╟─6004dd44-7379-4a5d-9b0f-7f66cf74d209
# ╟─bee9301a-42d4-4382-b13f-ce72b9c2c87e
# ╟─6ed17596-5be9-44aa-9e1f-5c77e9bde4c8
# ╟─c53bc1b8-3bd7-49ba-896a-1fc9356d1133
# ╠═0db7e996-5443-45b1-86cc-c0711c747bf6
# ╟─afa9059c-f81d-49f8-b11d-62834a8982d8
# ╠═5276c12e-3673-4a70-8842-4cdbde9f28bf
# ╠═8d616b0a-3c80-47aa-bc8a-b8bd1215ec4b
# ╟─89f897d4-0b9b-4a18-a201-a9afb61772d9
# ╟─decf7b15-cfe0-44b8-b906-9b190ba37568
