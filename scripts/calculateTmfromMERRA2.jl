using DrWatson
@quickactivate "CompareTmTs"

include(srcdir("merra2.jl"))

merraTm(start=Date(2020),stop=Date(2020,12,31))