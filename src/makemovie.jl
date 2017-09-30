import FEHM
import CrPlots

if isdefined(:workingdir)#get this from sanitycheck if it is there
	basedir = workingdir
else
	basedir = "cond05"
end
moviedir = "$basedir/movie"
if !isfile("$basedir/w01.jld")
	timefilter = t->t > 0 && t < (2029 - 1964) * 365.25 && abs(round(t / 365.25) - t / 365.25) < 1e-2
	FEHM.avs2jld("$basedir/wl/w01.00001_geo", "$basedir/wl/w01", "$basedir/w01.jld"; timefilter=timefilter)
end
xs, ys, zs, times, crdata = JLD.load("$basedir/w01.jld", "xs", "ys", "zs", "times", "Cr")
topnodes = map(Int, readdlm("../smoothgrid2/out_top.nodes_xyz")[:, 1])
if !isdir(moviedir)
	mkdir(moviedir)
end

x0 = CrPlots.getwelllocation("R-62")[1] - 250
x1 = CrPlots.getwelllocation("R-36")[1] + 250
y0 = CrPlots.getwelllocation("R-50")[2] - 500
y1 = CrPlots.getwelllocation("R-43")[2] + 250
boundingbox = (x0, y0, x1, y1)
plotwells = [["CrPZ-$i" for i in 1:5]; ["CrIN-$i" for i in 1:5]; ["CrEX-$i" for i in [1, 3]]; ["R-61", "R-42", "R-28", "R-50", "R-44", "R-45", "R-13", "SIMR-2", "R-43", "R-62", "R-15", "PM-03", "O-04", "R-11", "R-36"]]
lowerlimit = 50
upperlimit = 1000

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
for (prefix, data) in zip(["top", "thirteenmeters"], [topdata, thirteenmeterdata])
	framenum = 1
	for (i, t) in enumerate(times)
		if i == 1 || (i > 1 && times[i] - times[i - 1] > 365.25 / 2)
			fig, ax, img = CrPlots.crplot(boundingbox, xs[topnodes], ys[topnodes], data[i]; upperlimit=upperlimit, lowerlimit=lowerlimit)
			CrPlots.addwells(ax, plotwells)
			CrPlots.addcbar(fig, img, "Cr [ppb]", CrPlots.getticks([lowerlimit, upperlimit]))
			CrPlots.addmeter(ax, boundingbox[3] - 1100, boundingbox[2] + 0, [250, 500, 1000], ["250m", "500m", "1km"])
			CrPlots.addpbar(fig, ax, (t - minimum(times)) / (maximum(times) - minimum(times)), "Year $(round(Int, 1964 + t / 365.25))")
			#display(fig); println()
			fig[:savefig]("$moviedir/$prefix$(lpad(framenum, 4, 0)).png", dpi=120)
			framenum += 1
			PyPlot.close(fig)
			@show 1964 + (t / 365.25)
		end
		#@show i
	end
	run(`ffmpeg -i $moviedir/$prefix%04d.png -vcodec libx264 -pix_fmt yuv420p -f mp4 -y $moviedir/$(prefix)_$(basedir).mp4`)
end
