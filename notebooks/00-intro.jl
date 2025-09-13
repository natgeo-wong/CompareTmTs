### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ d20942a4-51f2-11eb-0977-33f2082de974
md"
# CompareTmTs: Comparing Tm and Ts from ERA5-based Products

This project is twofold: (1) it creates a dataset for Tm and Pi based on the formulae from Davis et al. [1985], and (2) determines if the linear $T_m$-$T_s$ relationship first canonized by Bevis et al. [1992] is applicable outside of the United States.

Notebook List Breakdown:
1. Basic analysis on the calculated $T_m$ and $T_s$ based off ERA5 reanalysis products
2. A look at IGRAv2 radiosonde products and validation of the $T_m$ products using radiosonde data
3. Comparison of the $T_m$ and $T_s$ relationship
4. WTG figure
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╟─d20942a4-51f2-11eb-0977-33f2082de974
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
