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
	using ERA5Reanalysis
	using NCDatasets
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the CompareTmTs project..."
end

# ╔═╡ bcd1a6c2-ae3d-11ed-067b-6b7828037b63
md"
# 01a. Exploring the TmPi Dataset
"

# ╔═╡ 803ea2b9-1b0d-41de-82bd-bdd6a7444de9
e5ds = ERA5Hourly(start=Date(1979,1),stop=Date(2021,12),path=datadir())

# ╔═╡ f35cc155-ad5b-47ba-9702-094ff4df4982
ereg = ERA5Region("GLB",resolution=1)

# ╔═╡ 78cf32dd-7c9b-4990-83c2-1dbfa22dc8f3
egeo = RectRegion("TRP_TMP","GLB","Temporary Tropics",[20,-25,360,0],savegeo=false)

# ╔═╡ 0f411cd8-30dc-4d00-a66f-b174eb0de1d0
lsd = getLandSea(e5ds,ereg)

# ╔═╡ 2a49afd5-f41e-4ac9-b7f5-b40bdd072bf3
lsd_TMP = getLandSea(e5ds,ERA5Region(egeo,resolution=1))

# ╔═╡ 95bf6b0e-7395-4c4b-8ef7-4e59ca37c21a
evar_tm = SingleVariable("Tm")

# ╔═╡ e320d474-1e70-488b-963a-744764bf6e87
evar_ts = SingleVariable("t2m")

# ╔═╡ 5cf63bdf-84e0-4b64-bd8f-c2efee783a88
begin
	coord = readdlm(datadir("GLB-l.txt"),comments=true,comment_char='#')
	clon = coord[:,1]; clat = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ 4fc01c2e-c1bc-46e0-b8cd-f56e51ec1e33
begin
	ds = NCDataset(datadir("tmtscorr.nc"))
	ρ  = ds["rho"][:]
	close(ds)
	ds = NCDataset(datadir("tmtscorr_GLB.nc"))
	ρd = ds["rho"][:]
	close(ds)
end

# ╔═╡ 1baec3c8-687a-45e3-8b8b-524738a0df9b
ggrd = RegionGrid(egeo,lsd.lon,lsd.lat)

# ╔═╡ ac94d629-5ad5-49e6-8829-a3f854463063
ρ_TRP = extractGrid(ρ,ggrd)

# ╔═╡ 762b6473-d18c-49f2-9265-3e9c0ec977c6
ρd_TRP = extractGrid(ρd,ggrd)

# ╔═╡ f04272fa-fc91-450c-be17-84a89929607e
begin
	pplt.close(); fig,axs = pplt.subplots(
		aspect=2,nrows=3,axwidth=2,#wspace=[1,2.5],
		# proj="robin",proj_kw=Dict("lon_0"=>180),
	)

	lvls = vcat(.05,0.1,0.15,0.2:0.1:0.8,0.85,0.9,0.95)
	lvls2 = -0.2:0.05:0.2
	
	c1 = axs[1].pcolormesh(lsd.lon,lsd.lat,abs.(ρ'),levels=lvls,cmap="RdBu_r",extend="both")
	axs[2].pcolormesh(lsd.lon,lsd.lat,abs.(ρd'),levels=lvls,cmap="RdBu_r",extend="both")
	c2 = axs[3].pcolormesh(lsd.lon,lsd.lat,abs.(ρ').-abs.(ρd'),levels=lvls2,cmap="RdBu_r",extend="both")
	axs[3].contour(lsd.lon,lsd.lat,abs.(ρ'),levels=[0.4],linestyle=":",c="k")

	for ax in axs
		ax.plot(clon,clat,c="k",lw=0.25)
		ax.format(xlim=(90,165),ylim=(-20,20))
	end

	axs[1].format(title="(a) Hourly",suptitle=L"Correlation between $T_m$ and $T_s$")
	axs[2].format(title="(b) Daily")
	axs[3].format(title="(c) Hourly - Daily")

	# fig.colorbar(c1,loc="b",length=0.8,locator=lvls,col=(1,2))
	# fig.colorbar(c2,loc="b",length=0.8,locator=-0.2:0.1:0.2,col=3)
	fig.savefig(plotsdir("01g-tstmcorr.png"),transparent=false,dpi=300)
	load(plotsdir("01g-tstmcorr.png"))
end

# ╔═╡ 6e6b521f-43ea-4491-80ae-f61a34b3cd41
begin
	pplt.close(); f2,a2 = pplt.subplots()
	
	a2[1].scatter(ρ[lsd.lsm.<0.5],ρd[lsd.lsm.<0.5],s=2)
	a2[1].scatter(ρ[lsd.lsm.>0.5],ρd[lsd.lsm.>0.5],s=2)
	a2[1].scatter(ρ_TRP[lsd_TMP.lsm.>0.5],ρd_TRP[lsd_TMP.lsm.>0.5],s=2)
	a2[1].format(xlabel=L"$\rho_h$",ylabel=L"$\rho_d$")
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─bcd1a6c2-ae3d-11ed-067b-6b7828037b63
# ╟─fc2fdbeb-f775-48af-9d93-82d8a79567a9
# ╟─ed931d06-d36a-4570-b4f3-c6eaea16ff59
# ╟─803ea2b9-1b0d-41de-82bd-bdd6a7444de9
# ╟─f35cc155-ad5b-47ba-9702-094ff4df4982
# ╠═78cf32dd-7c9b-4990-83c2-1dbfa22dc8f3
# ╟─0f411cd8-30dc-4d00-a66f-b174eb0de1d0
# ╟─2a49afd5-f41e-4ac9-b7f5-b40bdd072bf3
# ╟─95bf6b0e-7395-4c4b-8ef7-4e59ca37c21a
# ╟─e320d474-1e70-488b-963a-744764bf6e87
# ╟─5cf63bdf-84e0-4b64-bd8f-c2efee783a88
# ╟─4fc01c2e-c1bc-46e0-b8cd-f56e51ec1e33
# ╟─1baec3c8-687a-45e3-8b8b-524738a0df9b
# ╠═ac94d629-5ad5-49e6-8829-a3f854463063
# ╠═762b6473-d18c-49f2-9265-3e9c0ec977c6
# ╠═f04272fa-fc91-450c-be17-84a89929607e
# ╠═6e6b521f-43ea-4491-80ae-f61a34b3cd41
