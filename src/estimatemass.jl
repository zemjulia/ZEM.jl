import CrPlots
import FEHM
import Interpolations
import Mads
import PyCall
import PyPlot

@PyCall.pyimport scipy.interpolate as scint

case = "cond05a"
madsfilename = "w01-v04.mads"

newdir = string(case, "_massestimate")
if !isdir(newdir)
	cp(case, newdir)
else
	warn("directory already exists -- not copying a fresh one")
end
startingdir = pwd()
cd(newdir)
if !isdefined(:md)
	md = Mads.loadmadsfile(madsfilename)
else
	warn("reusing existing md rather than reloading")
end
Mads.writeparameters(md)
include(joinpath(newdir, "runmodelcode.jl"))
if !isfile("w01.hyco")
	setupmaindir()
	preprocess(:wl)
else
	warn("skipping setupmaindir, because w01.hyco already exists")
end
#update the data file to have avs and the times we want
cd("wl")
rm("w01.data")
datalines = readlines("../wl.data")
f = open("w01.data", "w")
write(f, datalines[1])
write(f, "cont\navs 1 1e20\nc\nendavs\n")
for i = 2:length(datalines)
	write(f, datalines[i])
end
close(f)
if !isfile("w01.avs_log")
	runmodelindir(:wl)
else
	warn("skipping runmodel because w01.avs_log already exists")
end
if !isfile("concentrations.jld")
	FEHM.avs2jld("../../../smoothgrid2/tet.geo", "w01", "concentrations.jld"; timefilter=t->t > 366 && t < 20000)
else
	warn("skipping avs2jld because concentrations.jld already exists")
end
crs, times, xs, ys, zs = JLD.load("concentrations.jld", "Cr", "times", "xs", "ys", "zs")
cd("..")

function fehmday2year(d)
	return 1964 + (d / 365.25)
end
const topnodes = map(Int, readdlm("../../smoothgrid2/out_top.nodes_xyz")[:, 1])
const xy2topnode = Dict{Tuple{Float64, Float64}, Int}()
for m in topnodes
	xy2topnode[(xs[m], ys[m])] = m
end
const topnode2i = Dict(zip(topnodes, 1:length(topnodes)))

function integratealongz(cr)
	local integrals = Array(Float64, length(topnodes))
	local zcoords = Array(Array{Float64, 1}, length(topnodes))
	local crvals = Array(Array{Float64, 1}, length(topnodes))
	for i = 1:length(topnodes)
		zcoords[i] = Array{Float64}(0)
		crvals[i] = Array{Float64}(0)
	end
	for i = 1:length(cr)
		x, y, z, crval = xs[i], ys[i], zs[i], cr[i]
		tni = topnode2i[xy2topnode[(x, y)]]
		push!(zcoords[tni], z)
		push!(crvals[tni], max(0, crval - 5))
	end
	for i = 1:length(integrals)
		integrals[i] = 0.0
		integrals[i] = 0.5 * crvals[i][1] * (zcoords[i][2] - zcoords[i][1])
		for j = 2:length(zcoords[i]) - 1
			integrals[i] += 0.5 * crvals[i][j] * (zcoords[i][j + 1] - zcoords[i][j - 1])
		end
		integrals[i] += 0.5 * crvals[i][end] * (zcoords[i][end] - zcoords[i][end - 1])
	end
	return integrals
end

const n = 1000
const xmin = CrPlots.getwelllocation("R-01")[1] - 250
const xmax = CrPlots.getwelllocation("R-36")[1] + 250
const ymin = CrPlots.getwelllocation("R-34")[2] - 250
const ymax = CrPlots.getwelllocation("R-43")[2] + 1000
const gridxs = [x for x = range(xmin, xmax; length=n), y = range(ymin, ymax; length=n)]
const gridys = [y for x = range(xmin, xmax; length=n), y = range(ymin, ymax; length=n)]

function integratealongxy(zintegrals)
	gd1 = scint.griddata(hcat(xs[topnodes], ys[topnodes]), zintegrals, (gridxs, gridys), fill_value=0.0, method="linear")
	return mean(gd1) * (maximum(gridxs) - minimum(gridxs)) * (maximum(gridys) - minimum(gridys))
end

#if true
if !isdefined(:integral3ds)
	integral3ds = Array(Float64, length(times))
	for (i, t) in enumerate(times)
		zintegrals = integratealongz(crs[i])
		integral3ds[i] = integratealongxy(zintegrals)
	end
else
	warn("skipping integral3ds")
end

porosity = 0.3
#porosity = 0.07686779642653208
#integral3ds have units of (μg/L-pores)*m^3-space
literspermetercubed = 1000
mass = integral3ds * literspermetercubed * porosity * 1e-9#convert to kg

for (i, t) in enumerate(times)
	@show round(Int, fehmday2year(t)), round(Int, mass[i])
end
fig, ax = PyPlot.subplots()
ax[:plot](map(x->fehmday2year(x), times), mass)
ax[:set_ylabel]("Cr mass [kg]")
ax[:set_xlabel]("Year")
ax[:set_xlim](1970, 2017)
fig[:savefig]("mass.pdf")
display(fig); println()
PyPlot.close(fig)

cd(startingdir)
