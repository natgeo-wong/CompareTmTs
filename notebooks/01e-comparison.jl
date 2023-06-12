### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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

# ╔═╡ 8ec33dcf-60ff-4a45-8e85-74601b8eb4bc
e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021),path=datadir())

# ╔═╡ 04aca35c-fff8-479e-92ba-2576b18d87fa
evar_Tm = SingleVariable("Tm")

# ╔═╡ faf7a70d-c07c-4519-b396-5fbc960416be
evar_t2m = SingleVariable("t2m")

# ╔═╡ fc0c2edd-2969-401f-86c8-6c9dbc0f78f9
evar_sst = SingleVariable("sst")

# ╔═╡ 420f34e4-3981-40c5-b2af-9fd031ba95e7
egeo = ERA5Region("GLB",gres=0.25)

# ╔═╡ 5719b0a5-c160-4400-8b13-fde3c34f9c2d
egeo_TRP = ERA5Region("TRP",gres=0.25)

# ╔═╡ 5a690017-7dde-4feb-86e6-78f68725f4aa
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ b1a55c36-834c-4b33-bb4a-d893ee86e1d4
lsd_TRP = getLandSea(egeo_TRP,path=datadir("emask"))

# ╔═╡ 5ab3a7af-a790-434a-a3cf-f6daf4de96a3
begin
	ds_Tm   = NCDataset(datadir("compile_Tm.nc"))
	ayr_Tm  = ds_Tm["mean"][:]
	dhr_Tm  = ds_Tm["variability_diurnal"][:]
	aitr_Tm = ds_Tm["variability_intraseasonal"][:]
	amo_Tm  = ds_Tm["variability_seasonal"][:]
	avar_Tm = ds_Tm["variability_interannual"][:]
	close(ds_Tm)
	
	ds_Ts   = NCDataset(datadir("compile_t2m.nc"))
	ayr_Ts  = ds_Ts["mean"][:]
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

# ╔═╡ 9e36cd2a-cc7c-43e0-9e5c-96aa146c7122
md"
### A. Basic side-by-side comparison of $T_m$ and $T_s$
"

# ╔═╡ 33ffc329-5dc7-4480-b989-8b88e63eb3b7
md"Make Image 1? $(@bind makeimage1 PlutoUI.Slider(0:1))"

# ╔═╡ c3f7b6a2-5326-4b26-902f-de33e4e46a73
begin
	if isone(makeimage1)
		pplt.close(); f1,a1 = pplt.subplots(
			nrows=6,ncols=2,aspect=1.9,axwidth=2,
			proj="robin",proj_kw=Dict("lon_0"=>180)
		)
	
		x = [60,75,95,175,160,175,115,70,60]
		y = [-5,5,10,5,0,-10,-10,-7,-5]
		
		c1_1 = a1[1].pcolormesh(lsd.lon,lsd.lat,ayr_Tm',extend="both",levels=250:2.5:300)
		a1[1].contour(lsd_TRP.lon,lsd_TRP.lat,sst'.-273.15,levels=[29,35,36],c="w",linestyle="--",lw=0.5)
		a1[1].format(title=L"$T_m$")
		a1[2].pcolormesh(lsd.lon,lsd.lat,ayr_Ts',extend="both",levels=250:2.5:300)
		a1[2].contour(lsd_TRP.lon,lsd_TRP.lat,sst'.-273.15,levels=[29,35,36],c="w",linestyle="--",lw=0.5)
		a1[2].format(title=L"$T_s$")
		a1[2].colorbar(c1_1,loc="r",locator=250:10:300,minorlocator=[])
	
		a1[1].format(
		    leftlabels=[
				L"(a) $\mu$",L"(b) $\delta_t$",L"(c) $\delta_s$",L"(d) $\delta_i$",L"(e) $\delta_d$",L"(f) $\delta_a$"
			],leftlabelrotation="horizontal",
		)
		
		c1_2 = a1[3].pcolormesh(
			lsd.lon,lsd.lat,(amo_Tm.+aitr_Tm.+avar_Tm.+dhr_Tm)',
			extend="both",levels=10:10:90
		)
		a1[4].pcolormesh(
			lsd.lon,lsd.lat,(amo_Ts.+aitr_Ts.+avar_Ts.+dhr_Ts)',
			extend="both",levels=10:10:90
		)
		a1[4].colorbar(c1_2,loc="r")
		
		c1_3 = a1[5].pcolormesh(lsd.lon,lsd.lat,amo_Tm',extend="both",levels=5:5:45)
		a1[6].pcolormesh(lsd.lon,lsd.lat,amo_Ts',extend="both",levels=5:5:45)
		a1[6].colorbar(c1_3,loc="r")
		
		c1_4 = a1[7].pcolormesh(lsd.lon,lsd.lat,aitr_Tm',extend="both",levels=5:5:45)
		a1[8].pcolormesh(lsd.lon,lsd.lat,aitr_Ts',extend="both",levels=5:5:45)
		a1[8].colorbar(c1_4,loc="r")
		
		c1_5 = a1[9].pcolormesh(lsd.lon,lsd.lat,dhr_Tm',extend="both",levels=1:9)
		a1[10].pcolormesh(lsd.lon,lsd.lat,dhr_Ts',extend="both",levels=1:9)
		a1[10].colorbar(c1_5,loc="r")
		
		c1_6 = a1[11].pcolormesh(lsd.lon,lsd.lat,avar_Tm',extend="both",levels=(1:9)/2)
		a1[12].pcolormesh(lsd.lon,lsd.lat,avar_Ts',extend="both",levels=(1:9)/2)
		a1[12].colorbar(c1_6,loc="r")
	
		for ax in a1
			ax.format(coast=true)
		end
		
		f1.savefig(plotsdir("01d-TsTm.png"),transparent=false,dpi=400)
		f1.savefig(plotsdir("01d-TsTm.pdf"),transparent=false,dpi=400)
		load(plotsdir("01d-TsTm.png"))
	else
		load(plotsdir("01d-TsTm.png"))
	end
end

# ╔═╡ f1abefad-f66b-41d9-91b2-24abdfbeceee
md"
### B. Comparison of the components of variability
"

# ╔═╡ 3cbadf62-69f5-4406-b775-107a84a3aa29
md"Make Image 2? $(@bind makeimage2 PlutoUI.Slider(0:1))"

# ╔═╡ af99619a-5b0f-435d-a61d-18876c83d54f
begin
	if isone(makeimage2)
		pplt.close(); f2,a2 = pplt.subplots(
			nrows=4,ncols=3,aspect=1.9,axwidth=2,wspace=[1,3],
			proj="robin",proj_kw=Dict("lon_0"=>180)
		)
		
		dTs_tot = amo_Ts.+aitr_Ts.+avar_Ts.+dhr_Ts
		dTm_tot = amo_Tm.+aitr_Tm.+avar_Tm.+dhr_Tm
	
		c2_1 = a2[1].pcolormesh(lsd.lon,lsd.lat,(amo_Tm./dTm_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		       a2[2].pcolormesh(lsd.lon,lsd.lat,(amo_Ts./dTs_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		c2_3 = a2[3].pcolormesh(lsd.lon,lsd.lat,(amo_Tm./dTm_tot.-amo_Ts./dTs_tot)'*100,levels=vcat(-50:10:-10,-5,5,10:10:50),extend="both")
	
		a2[4].pcolormesh(lsd.lon,lsd.lat,(aitr_Tm./dTm_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[5].pcolormesh(lsd.lon,lsd.lat,(aitr_Ts./dTs_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[6].pcolormesh(lsd.lon,lsd.lat,(aitr_Tm./dTm_tot.-aitr_Ts./dTs_tot)'*100,levels=vcat(-50:10:-10,-5,5,10:10:50),extend="both")
	
		a2[7].pcolormesh(lsd.lon,lsd.lat,(dhr_Tm./dTm_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[8].pcolormesh(lsd.lon,lsd.lat,(dhr_Ts./dTs_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[9].pcolormesh(lsd.lon,lsd.lat,(dhr_Tm./dTm_tot.-dhr_Ts./dTs_tot)'*100,levels=vcat(-50:10:-10,-5,5,10:10:50),extend="both")
	
		a2[10].pcolormesh(lsd.lon,lsd.lat,(avar_Tm./dTm_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[11].pcolormesh(lsd.lon,lsd.lat,(avar_Ts./dTs_tot)'*100,levels=10:10:90,extend="both",cmap="deep")
		a2[12].pcolormesh(lsd.lon,lsd.lat,(avar_Tm./dTm_tot.-avar_Ts./dTs_tot)'*100,levels=vcat(-50:10:-10,-5,5,10:10:50),extend="both")
	
		
		f2.colorbar(c2_1,loc="l",label="%",rows=[2,3],length=0.8,space=5)
		f2.colorbar(c2_3,loc="r",label="%",rows=[2,3],length=0.8)
	
		a2[1].format(leftlabels=[
			L"(a) $\delta_s/\delta_t$",L"(b) $\delta_i/\delta_t$",L"(c) $\delta_d/\delta_t$",L"(d) $\delta_a/\delta_t$"
		],leftlabelrotation="horizontal",title=L"T_m")
		a2[2].format(title=L"T_s")
		a2[3].format(title=L"T_m-T_s")
	
		for ax in a2
			ax.format(coast=true)
		end
		
		f2.savefig(plotsdir("01d-TsTmvariability.png"),transparent=false,dpi=400)
		f2.savefig(plotsdir("01d-TsTmvariability.pdf"),transparent=false,dpi=400)
		load(plotsdir("01d-TsTmvariability.png"))
	else
		load(plotsdir("01d-TsTmvariability.png"))
	end
end

# ╔═╡ b776a438-1ab5-4615-b0e5-e277b891b5da
md"
### C. Some double checking
"

# ╔═╡ 9d2eb0fa-be5f-4b9b-b2c2-195a005aad4d
md"Make Image 3? $(@bind makeimage3 PlutoUI.Slider(0:1))"

# ╔═╡ a28df267-51c8-4cfb-8c10-cf8a8ba38230
begin
	if isone(makeimage3)
		pplt.close(); f3,a3 = pplt.subplots(
			nrows=2,ncols=2,aspect=1.9,axwidth=2,
			proj="robin",proj_kw=Dict("lon_0"=>180)
		)

		lvls = vcat(0.1,0.2,1/3,0.5,2/3,0.8,0.9,10/9,10/8,1.5,2,3,5,10)
		
		c3 = a3[1].pcolormesh(lsd.lon,lsd.lat,(dhr_Tm./dhr_Ts)',levels=lvls,extend="both",cmap="RdBu_r")
		a3[2].pcolormesh(lsd.lon,lsd.lat,(aitr_Tm./aitr_Ts)',levels=lvls,extend="both",cmap="RdBu_r")
		a3[3].pcolormesh(lsd.lon,lsd.lat,(amo_Tm./amo_Ts)',levels=lvls,extend="both",cmap="RdBu_r")
		a3[4].pcolormesh(lsd.lon,lsd.lat,(avar_Tm./avar_Ts)',levels=lvls,extend="both",cmap="RdBu_r")
	
		for ax in a3
			ax.format(coast=true)
		end

		f3.colorbar(c3)
		f3.savefig(plotsdir("01d-TsTmRatio.png"),transparent=false,dpi=150)
		load(plotsdir("01d-TsTmRatio.png"))
	else
		load(plotsdir("01d-TsTmRatio.png"))
	end
end

# ╔═╡ e7187360-b1e7-4c14-b1ca-20fb1f2dd000
md"
### D. Correlation?
"

# ╔═╡ 7c7ca71a-40b8-4dd7-b0e6-aba1c6a3a573
begin
	tmbin = 265:0.5:295
	tsbin = 275:0.5:305
	wgtsocn = zeros(length(tmbin)-1,length(tsbin)-1)
	wgtslnd = zeros(length(tmbin)-1,length(tsbin)-1)
	wgtstot = zeros(length(tmbin)-1,length(tsbin)-1)
	for ilat = 241 : 481
		lsii = lsd.lsm[:,ilat]
		Tmii = ayr_Tm[lsii.<0.5,ilat][:]
		Tsii = ayr_Ts[lsii.<0.5,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtsocn[:,:] .+= hist.weights * cosd(lsd.lat[ilat])
		Tmii = ayr_Tm[lsii.>0.5,ilat][:]
		Tsii = ayr_Ts[lsii.>0.5,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtslnd[:,:] .+= hist.weights * cosd(lsd.lat[ilat])
		Tmii = ayr_Tm[:,ilat][:]
		Tsii = ayr_Ts[:,ilat][:]
		hist = fit(Histogram, (Tmii,Tsii),(tmbin,tsbin))
		wgtstot[:,:] .+= hist.weights * cosd(lsd.lat[ilat])
	end
	wgtslnd = wgtslnd ./ sum(wgtstot) * (length(tmbin)-1) * (length(tsbin)-1)
	wgtsocn = wgtsocn ./ sum(wgtstot) * (length(tmbin)-1) * (length(tsbin)-1)
	wgtstot = wgtstot ./ sum(wgtstot) * (length(tmbin)-1) * (length(tsbin)-1)
	wgtsocn[iszero.(wgtsocn)] .= NaN
	wgtslnd[iszero.(wgtslnd)] .= NaN
	wgtstot[iszero.(wgtstot)] .= NaN
end

# ╔═╡ 76f3ffc8-377d-4bac-86ee-a82cd26fa48f
begin
	pplt.close(); f4,a4 = pplt.subplots(axwidth=1.5,sharex=0,sharey=0)
	
	# c = a4[1].pcolormesh(210:310,215:315,wgtslnd',levels=vcat(1,2,5,10,20,50,100),extend="both",cmap="greens",cmap_kw=Dict("left"=>0.1))
	c = a4[1].pcolormesh(tmbin,tsbin,wgtsocn',levels=5:5:100,extend="both",cmap="blues",cmap_kw=Dict("left"=>0.1))
	# a4[2].scatter(dhr_Tm[:],dhr_Ts[:])
	# a4[3].scatter(aitr_Tm[:],aitr_Ts[:])
	# a4[4].scatter(amo_Tm[:],amo_Ts[:])
	# a4[5].scatter(avar_Tm[:],avar_Ts[:])

	a4[1].format(
		xlim=(minimum(tmbin),maximum(tmbin)),
		ylim=(minimum(tsbin),maximum(tsbin))
	)
	a4[1].colorbar(c)
	
	f4.savefig("test.png",transparent=false,dpi=300)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─8ec33dcf-60ff-4a45-8e85-74601b8eb4bc
# ╟─04aca35c-fff8-479e-92ba-2576b18d87fa
# ╠═faf7a70d-c07c-4519-b396-5fbc960416be
# ╟─fc0c2edd-2969-401f-86c8-6c9dbc0f78f9
# ╟─420f34e4-3981-40c5-b2af-9fd031ba95e7
# ╟─5719b0a5-c160-4400-8b13-fde3c34f9c2d
# ╟─5a690017-7dde-4feb-86e6-78f68725f4aa
# ╟─b1a55c36-834c-4b33-bb4a-d893ee86e1d4
# ╟─5ab3a7af-a790-434a-a3cf-f6daf4de96a3
# ╟─9e36cd2a-cc7c-43e0-9e5c-96aa146c7122
# ╟─33ffc329-5dc7-4480-b989-8b88e63eb3b7
# ╟─c3f7b6a2-5326-4b26-902f-de33e4e46a73
# ╟─f1abefad-f66b-41d9-91b2-24abdfbeceee
# ╟─3cbadf62-69f5-4406-b775-107a84a3aa29
# ╟─af99619a-5b0f-435d-a61d-18876c83d54f
# ╟─b776a438-1ab5-4615-b0e5-e277b891b5da
# ╟─9d2eb0fa-be5f-4b9b-b2c2-195a005aad4d
# ╟─a28df267-51c8-4cfb-8c10-cf8a8ba38230
# ╟─e7187360-b1e7-4c14-b1ca-20fb1f2dd000
# ╟─7c7ca71a-40b8-4dd7-b0e6-aba1c6a3a573
# ╟─76f3ffc8-377d-4bac-86ee-a82cd26fa48f
