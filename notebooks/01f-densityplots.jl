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
	using PlutoUI
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

# ╔═╡ 63430452-9c10-4b83-8932-f47b7011fdb7
md"
### A. Loading LandSea masks for filtering
"

# ╔═╡ 5a690017-7dde-4feb-86e6-78f68725f4aa
lsd_GLB = getLandSea(ERA5Region("GLB",gres=0.25),path=datadir("emask"))

# ╔═╡ b1a55c36-834c-4b33-bb4a-d893ee86e1d4
lsd_TRP = getLandSea(ERA5Region("TRP"),path=datadir("emask"))

# ╔═╡ 5ab3a7af-a790-434a-a3cf-f6daf4de96a3
begin
	ds_Tm   = NCDataset(datadir("compile_Tm.nc"))
	ayr_Tm  = ds_Tm["mean"][:] .- 273.15
	dhr_Tm  = ds_Tm["variability_diurnal"][:]
	aitr_Tm = ds_Tm["variability_intraseasonal"][:]
	amo_Tm  = ds_Tm["variability_seasonal"][:]
	avar_Tm = ds_Tm["variability_interannual"][:]
	close(ds_Tm)
	
	ds_Ts   = NCDataset(datadir("compile_t2m.nc"))
	ayr_Ts  = ds_Ts["mean"][:] .- 273.15
	dhr_Ts  = ds_Ts["variability_diurnal"][:]
	aitr_Ts = ds_Ts["variability_intraseasonal"][:]
	amo_Ts  = ds_Ts["variability_seasonal"][:]
	avar_Ts = ds_Ts["variability_interannual"][:]
	close(ds_Ts)
	
	ds_sst  = NCDataset(datadir("compile_sst.nc"))
	sst = ds_sst["mean"][:]
	close(ds_sst)
	md"Loading compiled statistics ..."
end

# ╔═╡ e7187360-b1e7-4c14-b1ca-20fb1f2dd000
md"
### D. Correlation?
"

# ╔═╡ 0c1c859c-67e6-4a67-93b4-2f148df6fa8c
function binocnlnd(
	tmdata :: Array{<:Real,2},
	tsdata :: Array{<:Real,2},
	tmbin  :: AbstractRange,
	tsbin  :: AbstractRange,
	lsd    :: LandSea;
	dolatweight :: Bool = true
)

	wgtsocn = zeros(length(tmbin)-1,length(tsbin)-1)
	wgtslnd = zeros(length(tmbin)-1,length(tsbin)-1)
	wgtstot = zeros(length(tmbin)-1,length(tsbin)-1)
	for ilat = 1 : length(lsd.lat)
		if dolatweight
			wgt = cosd(lsd.lat[ilat])
		else
			wgt = 1
		end
		lsii = lsd.lsm[:,ilat]
		Tmii = tmdata[lsii.<0.4,ilat][:]
		Tsii = tsdata[lsii.<0.4,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtsocn[:,:] .+= hist.weights * wgt
		Tmii = tmdata[lsii.>0.4,ilat][:]
		Tsii = tsdata[lsii.>0.4,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtslnd[:,:] .+= hist.weights * wgt
		Tmii = tmdata[:,ilat][:]
		Tsii = tsdata[:,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtstot[:,:] .+= hist.weights * wgt
	end
	wgtslnd = wgtslnd ./ maximum(wgtslnd)
	wgtsocn = wgtsocn ./ maximum(wgtsocn)
	wgtstot = wgtstot ./ maximum(wgtstot)
	wgtsocn[wgtsocn.<0.1] .= NaN
	wgtslnd[wgtslnd.<0.1] .= NaN
	wgtstot[wgtstot.<0.1] .= NaN

	return wgtstot, wgtslnd, wgtsocn

end

# ╔═╡ c6c174ca-94b9-491d-bf4d-bafe59de3457
ggrd_TRP = RegionGrid(GeoRegion("TRP"),lsd_GLB.lon,lsd_GLB.lat)

# ╔═╡ c45e171d-b624-48e9-ac9d-e7bf106e06c3
ayr_Tm_TRP = extractGrid(ayr_Tm,ggrd_TRP)

# ╔═╡ 69b621ab-c4d8-4fb0-8949-1ba13420ba85
ayr_Ts_TRP = extractGrid(ayr_Ts,ggrd_TRP)

# ╔═╡ e2650c94-f37f-4ec7-b9eb-5caf3e7014d8
begin
	tmμbin = 7.5:0.2:17.5
	tsμbin = 20:0.2:30
	tmδdbin = 0 : 0.05 : 2.5
	tsδdbin = 0 : 0.2 : 15
	tmδibin = 0 : 0.5 : 20
	tsδibin = 0 : 0.5 : 25
	tmδsbin = 0 : 0.5 : 30
	tsδsbin = 0 : 0.5 : 40
	tmδabin = 0 : 0.05 : 2
	tsδabin = 0 : 0.05 : 2.5
end

# ╔═╡ 9c0c3af0-6c52-409f-9cd4-28a54b0cc270
begin
	GLBμ_tot,GLBμ_lnd,GLBμ_ocn = binocnlnd(ayr_Tm,ayr_Ts,tmμbin,tsμbin,lsd_GLB)
	GLBδd_tot,GLBδd_lnd,GLBδd_ocn = binocnlnd(dhr_Tm,dhr_Ts,tmδdbin,tsδdbin,lsd_GLB)
	GLBδi_tot,GLBδi_lnd,GLBδi_ocn = binocnlnd(aitr_Tm,aitr_Ts,tmδibin,tsδibin,lsd_GLB)
	GLBδs_tot,GLBδs_lnd,GLBδs_ocn = binocnlnd(amo_Tm,amo_Ts,tmδsbin,tsδsbin,lsd_GLB)
	GLBδa_tot,GLBδa_lnd,GLBδa_ocn = binocnlnd(avar_Tm,avar_Ts,tmδabin,tsδabin,lsd_GLB)
end

# ╔═╡ cca29898-12d4-444b-9470-707eebfe121f
TRP_tot,TRP_lnd,TRP_ocn = binocnlnd(
	ayr_Tm_TRP,ayr_Ts_TRP,tmμbin,tsμbin,lsd_TRP,dolatweight=true
)

# ╔═╡ 76f3ffc8-377d-4bac-86ee-a82cd26fa48f
begin
	pplt.close(); f4,a4 = pplt.subplots(axwidth=1.5,sharex=0,sharey=0)
	
	c = a4[1].pcolormesh(tmμbin,tsμbin,TRP_lnd',levels=(1:9)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	c = a4[1].pcolormesh(tmμbin,tsμbin,TRP_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))
	
	a4[1].colorbar(c)
	
	f4.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 3926b2be-8967-4dc5-b856-655d202fa12c
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=3,nrows=2,axwidth=1.5,sharex=0,sharey=0)
	
	c2_1 = a2[1].pcolormesh(tmμbin,tsμbin,GLBμ_lnd',levels=(1:9)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	c2_2 = a2[1].pcolormesh(tmμbin,tsμbin,GLBμ_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))
	
	a2[2].pcolormesh(tmδdbin,tsδdbin,GLBδd_lnd',levels=(1:9)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	a2[2].pcolormesh(tmδdbin,tsδdbin,GLBδd_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))
	
	a2[3].pcolormesh(tmδibin,tsδibin,GLBδi_lnd',levels=(1:9)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	a2[3].pcolormesh(tmδibin,tsδibin,GLBδi_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))
	
	a2[4].pcolormesh(tmδsbin,tsδsbin,GLBδs_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))
	a2[4].pcolormesh(tmδsbin,tsδsbin,GLBδs_lnd',levels=(1:9)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	
	a2[5].pcolormesh(tmδabin,tsδabin,GLBδa_lnd',levels=(1:5)/10,extend="max",cmap="greens",cmap_kw=Dict("left"=>0.1))
	a2[5].pcolormesh(tmδabin,tsδabin,GLBδa_ocn',levels=(1:9)/10,extend="max",cmap="blues",cmap_kw=Dict("left"=>0.1))

	# a2[4].format(xlim=(0,10),ylim=(0,10))
	
	f2.colorbar(c2_1,row=[1])
	f2.colorbar(c2_2,row=[2])
	f2.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─63430452-9c10-4b83-8932-f47b7011fdb7
# ╟─5a690017-7dde-4feb-86e6-78f68725f4aa
# ╟─b1a55c36-834c-4b33-bb4a-d893ee86e1d4
# ╟─5ab3a7af-a790-434a-a3cf-f6daf4de96a3
# ╟─e7187360-b1e7-4c14-b1ca-20fb1f2dd000
# ╠═0c1c859c-67e6-4a67-93b4-2f148df6fa8c
# ╟─c6c174ca-94b9-491d-bf4d-bafe59de3457
# ╠═c45e171d-b624-48e9-ac9d-e7bf106e06c3
# ╠═69b621ab-c4d8-4fb0-8949-1ba13420ba85
# ╠═e2650c94-f37f-4ec7-b9eb-5caf3e7014d8
# ╠═9c0c3af0-6c52-409f-9cd4-28a54b0cc270
# ╠═cca29898-12d4-444b-9470-707eebfe121f
# ╠═76f3ffc8-377d-4bac-86ee-a82cd26fa48f
# ╠═3926b2be-8967-4dc5-b856-655d202fa12c
