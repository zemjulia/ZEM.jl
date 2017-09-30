import Interpolations
import FEHM
import JLD
import CrPlots

if isdefined(:workingdir)#get this from sanitycheck if it is there
	basedir = workingdir
else
	basedir = "cond11_scenarios"
end
case = split(splitdir(basedir)[end], '_')[1]

if !isdefined(:scenarios)
	scenarios = JLD.load(basedir * ".jld")
end

moviedir = "$basedir/scenario_movies"

xs = nothing
ys = nothing
zs = nothing
times = Any[]
crdatas = Dict()
times = collect(365.25:365.25/8:23376.0)[36:end]
#times = [times[1], times[div(end, 2)], times[end]]
for scenarioname in keys(scenarios)
	if !isfile("$basedir/w01_$(scenarioname).jld")
		timefilter = t->t > 0 && t < (2029 - 1964) * 365.25 && abs(round(t / 365.25) - t / 365.25) < 1e-2
		FEHM.avs2jld("$basedir/wl_$(scenarioname)/w01.00001_geo", "$basedir/wl_$(scenarioname)/w01", "$basedir/w01_$(scenarioname).jld"; timefilter=timefilter)
	end
	lasttimes = times
	xs, ys, zs, thesetimes, thiscrdata = JLD.load("$basedir/w01_$(scenarioname).jld", "xs", "ys", "zs", "times", "Cr")
	itp = Interpolations.interpolate((thesetimes,), thiscrdata, Interpolations.Gridded(Interpolations.Linear()))
	crdatas[scenarioname] = map(t->itp[t], times)
end
topnodes = map(Int, readdlm("../smoothgrid2/out_top.nodes_xyz")[:, 1])
if !isdir(moviedir)
	mkdir(moviedir)
end

x0 = CrPlots.getwelllocation("R-62")[1] - 250
x1 = CrPlots.getwelllocation("R-36")[1] + 250
y0 = CrPlots.getwelllocation("R-50")[2] - 500
y1 = CrPlots.getwelllocation("R-43")[2] + 250
boundingbox = (x0, y0, x1, y1)
plotwells = [["CrPZ-$i" for i in 1:5]; ["CrIN-$i" for i in 1:6]; ["CrEX-$i" for i in 1:3]; ["R-61", "R-42", "R-28", "R-50", "R-44", "R-45", "R-13", "SIMR-2", "R-43", "R-62", "R-15", "PM-03", "O-04", "R-11", "R-36"]]
lowerlimit = 50
upperlimit = 1000

for scenarioname in keys(scenarios)
	crdata = crdatas[scenarioname]
	#set up the data for the top of the water table
	topdata = Array{Any}(length(times))
	for i = 1:length(times)
		topdata[i] = crdata[i][topnodes]
	end
	#set up the data for the average over the top 10 meters of the water table
	topdict = Dict()
	thirteenmeternodes = Array{Array{Int, 1}}(length(xs))
	for node in topnodes
		topdict[(xs[node], ys[node])] = node
		thirteenmeternodes[node] = Int[]
	end
	for (i, x, y, z) in zip(1:length(xs), xs, ys, zs)
		topnode = topdict[(x, y)]
		if z >= zs[topnode] - 13
			push!(thirteenmeternodes[topnode], i)
		end
	end
	thirteenmeterdata = deepcopy(topdata)
	for i in 1:length(times)
		for (j, topnode) in enumerate(topnodes)
			thirteenmeterdata[i][j] = mean(crdata[i][thirteenmeternodes[topnode]])
		end
	end
	for (prefix, data) in zip(["top_$(scenarioname)", "thirteenmeters_$(scenarioname)"], [topdata, thirteenmeterdata])
		@sync @parallel for (i, t) in collect(enumerate(times))
			@show 1964 + (t / 365.25)
			fig, ax, img = CrPlots.crplot(boundingbox, xs[topnodes], ys[topnodes], data[i]; upperlimit=upperlimit, lowerlimit=lowerlimit)
			CrPlots.addwells(ax, plotwells)
			CrPlots.addcbar(fig, img, "Cr [ppb]", CrPlots.getticks([lowerlimit, upperlimit]))
			CrPlots.addmeter(ax, boundingbox[3] - 1100, boundingbox[2] + 0, [250, 500, 1000], ["250m", "500m", "1km"])
			CrPlots.addpbar(fig, ax, (t - minimum(times)) / (maximum(times) - minimum(times)), "Year $(@sprintf("%.2f",1964 + t / 365.25))")
			#display(fig); println()
			fig[:savefig]("$moviedir/$(case)_$prefix$(lpad(i, 4, 0)).png", dpi=120)
			PyPlot.close(fig)
		end
		run(`ffmpeg -i $moviedir/$(case)_$prefix%04d.png -vcodec libx264 -pix_fmt yuv420p -f mp4 -y $moviedir/$(case)_$(prefix).mp4`) # good for latex
	end
end

for scenarioname1 in keys(scenarios), scenarioname2 in ["no_action"]
	if scenarioname1 == scenarioname2
		continue
	end
	crdata = crdatas[scenarioname1] - crdatas[scenarioname2]
	#set up the data for the top of the water table
	topdata = Array{Any}(length(times))
	for i = 1:length(times)
		topdata[i] = crdata[i][topnodes]
	end
	#set up the data for the average over the top 10 meters of the water table
	topdict = Dict()
	thirteenmeternodes = Array{Array{Int, 1}}(length(xs))
	for node in topnodes
		topdict[(xs[node], ys[node])] = node
		thirteenmeternodes[node] = Int[]
	end
	for (i, x, y, z) in zip(1:length(xs), xs, ys, zs)
		topnode = topdict[(x, y)]
		if z >= zs[topnode] - 13
			push!(thirteenmeternodes[topnode], i)
		end
	end
	thirteenmeterdata = deepcopy(topdata)
	for i in 1:length(times)
		for (j, topnode) in enumerate(topnodes)
			thirteenmeterdata[i][j] = mean(crdata[i][thirteenmeternodes[topnode]])
		end
	end
	for (prefix, data) in zip(["top_$(scenarioname1)_minus_$(scenarioname2)", "thirteenmeters_$(scenarioname1)_minus_$(scenarioname2)"], [topdata, thirteenmeterdata])
		diffupperlimit = maximum(map(maximum, data))
		difflowerlimit = minimum(map(minimum, data))
		@sync @parallel for (i, t) in collect(enumerate(times))
			@show 1964 + (t / 365.25)
			fig, ax, img = CrPlots.crplot(boundingbox, xs[topnodes], ys[topnodes], data[i]; upperlimit=diffupperlimit, lowerlimit=difflowerlimit)
			CrPlots.addwells(ax, plotwells)
			CrPlots.addcbar(fig, img, "Cr [ppb]", collect(linspace(difflowerlimit, diffupperlimit, 5)))
			CrPlots.addmeter(ax, boundingbox[3] - 1100, boundingbox[2] + 0, [250, 500, 1000], ["250m", "500m", "1km"])
			CrPlots.addpbar(fig, ax, (t - minimum(times)) / (maximum(times) - minimum(times)), "Year $(@sprintf("%.2f",1964 + t / 365.25))")
			#display(fig); println()
			fig[:savefig]("$moviedir/$(case)_$prefix$(lpad(i, 4, 0)).png", dpi=120)
			PyPlot.close(fig)
		end
		run(`ffmpeg -i $moviedir/$(case)_$prefix%04d.png -vcodec libx264 -pix_fmt yuv420p -f mp4 -y $moviedir/$(case)_$(prefix).mp4`) # good for latex
	end
end
