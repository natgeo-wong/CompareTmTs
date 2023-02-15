### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ dc2574f0-d0e3-11ec-3365-5d97655fccec
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Loading DrWatson to ensure project reproducibility"
end

# ╔═╡ 00d7653f-d878-4c80-a952-3240382a179f
begin
	@quickactivate "CompareTmTs"
	using NCDatasets
	using Printf
	using Statistics
	
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow
	
	pplt = pyimport("proplot")
	md"Loading modules for CompareTmTs project ..."
end

# ╔═╡ c543aab7-84d1-4314-9a9d-49b7554c2e7a
md"
### A. Let's Load some RCE data ...
"

# ╔═╡ 58994a96-16be-4e36-a389-c3fc03bfd4cd
begin
	fol = datadir("RCEProfile")
	p   = zeros(64)
	q   = zeros(64)
	t   = zeros(64)
	for inc = 1 : 10
		istr  = @sprintf("%02d",inc)
		ids   = NCDataset(joinpath(fol,"RCE_ExploreWTGSpace-Control-member$(istr).nc"))
		p[:] += ids["p"][:]
		q[:] += dropdims(mean(ids["QV"][:],dims=2),dims=2)
		t[:] += dropdims(mean(ids["TABS"][:],dims=2),dims=2)
	end
	p = p / 10
	q = q / 10
	t = t / 10
	md"Loading RCE data ..."
end

# ╔═╡ b7526a77-a1b6-48ed-884e-c558bc54c930
md"
### B. Here goes plotting ...
"

# ╔═╡ df0dbb14-b42b-495f-bc89-046929fcfbe9
begin
	arr = [
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[1,1,8,8,2,3,9,9,10,10,4,5,11,11],
		[0,0,6,6,6,6,6,6,7,7,7,7,7,7]
	]
	md"Set up plotting arrays ..."
end

# ╔═╡ a7ed378d-077a-4fa0-9720-7115bb171699
function hump(radius)
	x  = -1 : 0.01 : 1
	nx = length(x)
	y  = zeros(nx)
	for ix = 1 : nx
		if abs(x[ix]) < radius
			y[ix] = cos(x[ix] / radius * pi/2) * 0.8
		end
	end
	return x,y
end

# ╔═╡ e368cd1d-2a88-496c-bc1a-8810ab13693e
begin
	pplt.close(); fig,axs = pplt.subplots(
		arr,wspace=[0,1,0,0,0.1,0,0,1,0,0,0.1,0,0],
		aspect=0.5,axwidth=1,sharey=0
	)

	axs[1].plot(t,p)
	axs[1].fill_between([180,320],125,70,color="gray",alpha=0.1)
	axs[1].text(215,100,"Tropopause")
	axs[1].format(ylim=(1000,25),yscale="log",xlim=(180,320),ylabel="Pressure / hPa",xlabel=L"$T_v$ / K",ltitle=L"(a) Large-Scale $T_v$")
	
	axs[2].fill_between([-2,0],2,1.75,color="gray",alpha=0.1)
	axs[2].arrow(-0.2,2.95,0,-0.95,width=0.1,head_width=0.25,head_length=0.1,edgecolor="none",facecolor="k",length_includes_head=true)
	axs[2].text(-0.8,2.5,L"w_{wtg}")
	axs[2].format(xlim=(-1,0),ltitle="(b) Warm Perturbation",ultitle="(i)")
	
	axs[3].fill_between([0,2],2,1.75,color="gray",alpha=0.1)
	axs[3].plot(sin.(((1.9:0.01:3).-1.9)*pi/1.1)/1.2,1.9:0.01:3,c="red3",lw=1,linestyle="--")
	axs[3].plot(sin.(((1.9:0.01:3).-1.9)*pi/1.1)/5,1.9:0.01:3,c="red9")
	axs[3].arrow(0.8,2.45,-0.6,0,width=0.02,head_width=0.05,head_length=0.15,edgecolor="none",facecolor="gray4",length_includes_head=true)
	axs[3].format(grid=false,xlabel=L"$\delta T_v$")
	
	axs[4].fill_between([-2,0],2.1,1.85,color="gray",alpha=0.1)
	axs[4].plot(-sin.(((2.05:0.01:3).-2.05)*pi/0.95)/1.2,2.05:0.01:3,c="blue3",lw=1,linestyle="--")
	axs[4].plot(-sin.(((2.05:0.01:3).-2.05)*pi/0.95)/5,2.05:0.01:3,c="blue9")
	axs[4].arrow(-0.8,2.525,0.6,0,width=0.02,head_width=0.05,head_length=0.15,edgecolor="none",facecolor="gray4",length_includes_head=true)
	axs[4].format(grid=false,ltitle="(c) Cool Perturbation",xlabel=L"$\delta T_v$",ultitle="(i)")
	
	axs[5].fill_between([0,2],2.1,1.85,color="gray",alpha=0.1)
	axs[5].arrow(0.2,2.15,0,0.8,width=0.1,head_width=0.25,head_length=0.1,edgecolor="none",facecolor="k",length_includes_head=true)
	axs[5].text(0.35,2.55,L"w_{wtg}")
	axs[5].format(xlim=(0,1))
	
	axs[6].plot(hump(0.1)[1],hump(0.1)[2],c="red3",lw=1,linestyle="--")
	axs[6].plot(hump(0.8)[1],hump(0.8)[2]*0.25,c="red9")
	axs[6].format(ylabel=L"$\delta T_v$",ltitle="(ii)")
	
	axs[7].plot(hump(0.1)[1],-hump(0.1)[2],c="blue3",lw=1,linestyle="--")
	axs[7].plot(hump(0.8)[1],-hump(0.8)[2]*0.25,c="blue9")
	axs[7].format(ltitle="(ii)")
	
	for ii in 6 : 7
		    axs[ii].arrow(0.2,0.5,0.6,0,width=0.1,head_width=0.25,head_length=0.05,edgecolor="none",facecolor="gray4",length_includes_head=true)
		    axs[ii].arrow(0.2,-0.5,0.6,0,width=0.1,head_width=0.25,head_length=0.05,edgecolor="none",facecolor="gray4",length_includes_head=true)
		    axs[ii].arrow(-0.2,0.5,-0.6,0,width=0.1,head_width=0.25,head_length=0.05,edgecolor="none",facecolor="gray4",length_includes_head=true)
		    axs[ii].arrow(-0.2,-0.5,-0.6,0,width=0.1,head_width=0.25,head_length=0.05,edgecolor="none",facecolor="gray4",length_includes_head=true)
		axs[ii].text(-0.95,0.35,"+")
		axs[ii].text(-0.95,-.75,L"-")
		axs[ii].text(0.9,0.35,"+")
		axs[ii].text(0.9,-.75,L"-")
		axs[ii].format(
			xloc="neither",xticks=[],ytickloc="neither",
			ylim=(-1,1),xlim=(-1,1),yticks=[0],yticklabels=[]
		)
	end
	
	for ii in [8 10]
	    for addii = 1.55 : 0.15 : 2.95
	        axs[ii].plot(-1:0.01:0.8,sin.((-1:0.01:0.8)*pi*5.5)*0.02 .+addii,c="gray4")
	    end
	end
	
	for ii in [9 11]
	    for addii = 1.55 : 0.15 : 2.95
	        axs[ii].plot(-0.8:0.01:1,sin.((-0.8:0.01:1)*pi*5.5)*0.02 .+addii,c="gray4")
	    end
	end
	
	for ii in 8 : 11
	    axs[ii].format(xloc="neither",yloc="neither",xticks=[],yticks=[],xlim=(-1,1))
	end
	
	for ii in [2,5]
	    axs[ii].format(xloc="neither",yloc="neither",xticks=[],yticks=[])
	end
	
	axs[3].format(yticks=[],yloc="left",xloc="bottom",xlim=(0,1.1),xticklabels=["0","+"])
	axs[4].format(yticks=[],yloc="right",xloc="bottom",xlim=(-1.1,0),xticklabels=["0","-","0"])
	
	for ii = vcat(2:5,8:11)
	    axs[ii].format(ylim=[3,1.4],yticks=[],yticklabels=[])
	end
	
	fig.savefig(plotsdir("wtgfigure.png"),transparent=false,dpi=300)
	load(plotsdir("wtgfigure.png"))
end

# ╔═╡ Cell order:
# ╟─dc2574f0-d0e3-11ec-3365-5d97655fccec
# ╟─00d7653f-d878-4c80-a952-3240382a179f
# ╟─c543aab7-84d1-4314-9a9d-49b7554c2e7a
# ╟─58994a96-16be-4e36-a389-c3fc03bfd4cd
# ╟─b7526a77-a1b6-48ed-884e-c558bc54c930
# ╟─df0dbb14-b42b-495f-bc89-046929fcfbe9
# ╟─a7ed378d-077a-4fa0-9720-7115bb171699
# ╟─e368cd1d-2a88-496c-bc1a-8810ab13693e
