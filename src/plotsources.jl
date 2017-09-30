import FEHM
import CrPlots
import DataStructures

if isdefined(:workingdir)#get this from sanitycheck if it is there
	basedir = workingdir
else
	basedir = "cond05"
end
figdir = "$basedir/figs"
if !isfile("$basedir/w01.jld")
	timefilter = t->t > 0 && t < (2029 - 1964) * 365.25 && abs(round(t / 365.25) - t / 365.25) < 1e-2
	FEHM.avs2jld("$basedir/wl/w01.00001_geo", "$basedir/wl/w01", "$basedir/w01.jld"; timefilter=timefilter)
end
xs, ys, zs, times, crdata = JLD.load("$basedir/w01.jld", "xs", "ys", "zs", "times", "Cr")
topnodes = map(Int, readdlm("../smoothgrid2/out_top.nodes_xyz")[:, 1])
if !isdir(figdir)
	mkdir(figdir)
end

x0 = CrPlots.getwelllocation("R-62")[1] - 250
x1 = CrPlots.getwelllocation("R-36")[1] + 250
y0 = CrPlots.getwelllocation("R-50")[2] - 500
y1 = CrPlots.getwelllocation("R-43")[2] + 250
boundingbox = (x0, y0, x1, y1)
plotwells = [["CrPZ-$i" for i in 1:5]; ["CrIN-$i" for i in 1:5]; ["CrEX-$i" for i in [1, 3]]; ["R-61", "R-42", "R-28", "R-50", "R-44", "R-45", "R-13", "SIMR-2", "R-43", "R-62", "R-15", "PM-03", "O-04", "R-11", "R-36"]]
lowerlimit = 1000
upperlimit = 100000

function loadsourcedata(filename)
	data = readdlm(filename; skipstart=15)
	sourcedata = fill(NaN, length(xs))
	sourcex0sdict = DataStructures.DefaultDict(()->Float64[])
	sourcey0sdict = DataStructures.DefaultDict(()->Float64[])
	for i = 1:size(data, 1)
		if data[i, 1] > 0#if the data is being defined at a node (not the zones defined by the boundaries/R-36)
			concentration = data[i, 4]
			t0 = data[i, 5]
			node = Int(data[i, 1])
			sourcedata[node] = concentration
			push!(sourcex0sdict[t0], xs[node])
			push!(sourcey0sdict[t0], ys[node])
		end
	end
	times = Float64[]
	sourcex0s = Float64[]
	sourcey0s = Float64[]
	for t0 in keys(sourcex0sdict)
		push!(times, t0 / 365.25 + 1964)
		@show mean(sourcex0sdict[t0])
		@show length(sourcex0sdict[t0])
		push!(sourcex0s, mean(sourcex0sdict[t0]))
		push!(sourcey0s, mean(sourcey0sdict[t0]))
		@show t0, sourcex0sdict[t0][1]
	end
	return sourcedata, times, sourcex0s, sourcey0s
end

sourcedata, sourcetimes, sourcex0s, sourcey0s = loadsourcedata("$basedir/wl.trac")

fig, ax, img = CrPlots.crplot(boundingbox, xs[topnodes], ys[topnodes], log10.(sourcedata[topnodes]); upperlimit=log10(upperlimit), lowerlimit=log10(lowerlimit))
for i = 1:length(sourcetimes)
	@show i
	@show sourcetimes[i]
	@show round(Int, 1964 + sourcetimes[i] / 365.25)
	ax[:text](sourcex0s[i], sourcey0s[i], round(Int, sourcetimes[i]); weight="bold", fontsize=18)
end
CrPlots.addwells(ax, plotwells; alpha=0.25)
CrPlots.addcbar(fig, img, "log_10\nCr [ppb]", CrPlots.getticks([log10(lowerlimit), log10(upperlimit)]))
CrPlots.addmeter(ax, boundingbox[3] - 1100, boundingbox[2] + 0, [250, 500, 1000], ["250m", "500m", "1km"])
display(fig); println()
fig[:savefig]("$figdir/sources.png", dpi=120)
PyPlot.close(fig)
