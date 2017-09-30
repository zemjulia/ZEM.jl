import DataStructures
import Distributions
import Interpolations
import Kriging
@everywhere begin
	include("checknode.jl")
	include("gettrends.jl")
	include("myreadlines.jl")
	dirs = [:wl, "crex1", "crex3", "crin1", "crin2", "crin3", "crin4", "crin5", "o04", "pm01", "pm02", "pm03", "pm04", "pm05", "r28", "r42"]
	#dirs = [:wl]
	calibwells = ["crex1","crex3","crpz1","crpz2_1","crpz3","crpz4","crpz5","r1","r11","r13","r15","r28","r33_1","r33_2","r35b","r36","r42","r43_1","r43_2","r44_1","r44_2","r45_1","r45_2","r50_1","r50_2","r61_1","r61_2","r62","simr2","crin1","crin2","crin3","crin4","crin5","crin6","crex2","crpz2_2"]
	calibzones = ["620101","620301","630101","630201","630301","630401","630501","4010001","4110001","4130001","4150001","4280001","4330001","4330002","4350201","4360001","4420001","4430001","4430002","4440001","4440002","4450001","4450002","4500001","4500002","4610001","4610002","4620001","7050201","612101","611201","610301","618401","610501","610601","620201","630202"]
	starttimedict = DataStructures.DefaultDict(366)
	crtimedict = Dict("4010001"=>[], "4110001"=>collect(2005:2020), "4130001"=>collect(2002:2020), "4150001"=>[2000; 2001; collect(2003:2020)], "4280001"=>collect(2005:2020), "4330001"=>collect(2005:2020), "4330002"=>collect(2005:2020), "4350201"=>collect(2007:2020), "4360001"=>collect(2008:2020), "4420001"=>collect(2008:2020), "4430001"=>collect(2008:2020), "4430002"=>collect(2008:2020), "4440001"=>collect(2009:2020), "4440002"=>collect(2009:2020), "4450001"=>collect(2009:2020), "4450002"=>collect(2009:2020), "4500001"=>collect(2010:2020), "4500002"=>collect(2010:2020), "4610001"=>collect(2011:2020), "4610002"=>collect(2011:2020), "4620001"=>collect(2012:2020), "612101"=>collect(2010:2020), "611201"=>collect(2010:2020), "610301"=>collect(2010:2020), "618401"=>collect(2010:2020), "610501"=>collect(2010:2020), "610601"=>collect(2010:2020), "620101"=>collect(2010:2020), "620201"=>collect(2010:2020), "620301"=>collect(2010:2020), "630101"=>collect(2010:2020), "630201"=>collect(2010:2020), "630202"=>collect(2010:2020), "630301"=>collect(2010:2020), "630401"=>collect(2010:2020), "630501"=>collect(2010:2020), "7050201"=>collect(2010:2020))
	function year2fehmday(year)
		return (year - 1964) * 365.25
	end
	function parsezones()
		info("open ../../smoothgrid2/well_screens.zonn")
		f = open("../../smoothgrid2/well_screens.zonn")
		lines = myreadlines(f)
		nnumlines = Int64[]
		for i = 1:length(lines)
			if contains(lines[i], "nnum")
				push!(nnumlines, i)
			end
		end
		zonenumbers = Array{Int64}(length(nnumlines))
		for i = 1:length(nnumlines)
			zonenumbers[i] = parse(Int, split(lines[nnumlines[i] - 1])[1])
		end
		nodenumbers = Array{Array{Int64, 1}}(length(nnumlines))
		for i = 1:length(nnumlines)
			nodenumbers[i] = Int64[]
			numnodes = parse(Int, lines[nnumlines[i] + 1])
			j = 2
			while(length(nodenumbers[i]) < numnodes)
				newnodes = map(x->parse(Int, x), split(lines[nnumlines[i] + j]))
				nodenumbers[i] = [nodenumbers[i]; newnodes]
				j += 1
			end
		end
		close(f)
		d1 = Dict()
		d2 = Dict()
		for i = 1:length(nnumlines)
			d2[zonenumbers[i]] = Int64[]
			for j = 1:length(nodenumbers[i])
				d1[nodenumbers[i][j]] = zonenumbers[i]
				push!(d2[zonenumbers[i]], nodenumbers[i][j])
			end
		end
		return d1, d2
	end
	function postprocess(dir::AbstractString)
		info("open w01_head_his.dat")
		f = open(string(dir, "/w01_head_his.dat"))
		lines = myreadlines(f)
		close(f)
		zoneheaders = map(x -> string(parse(Int, x[2:8])), split(lines[2], "Zone")[2:end])
		#splitlines = map(line -> split(line, " ", 0, false), lines[5:end])
		floats = map(line -> map(float, split(line, " "; limit=0, keep=false)[1:end-1]), lines[5:end])
		#find the columns that corresponds to the zones we want to calibrate against
		#columns = Array(Int64, length(calibzones))
		columns = fill(-1, length(calibzones))
		initialheads = Array{Float64}(length(calibzones))
		if floats[end][1] == 1830
			for i = 1:length(calibzones)
				for j = 1:length(zoneheaders)
					if zoneheaders[j] == calibzones[i]
						columns[i] = j + 1
						if length(dir) >= 3 && (dir[1:3] == "r44" || dir[1:3] == "r45" || dir[1:3] == "cre")
							obsstarttime = 0
						else
							obsstarttime = starttimedict[zoneheaders[j]]
						end
						#get the initial head
						k = 1
						while floats[k][1] < obsstarttime
							k += 1
						end
						initialheads[i] = floats[k][j + 1]
					end
				end
			end
			if minimum(columns) == -1
				minind = indmin(columns)
				error("Didn't find a column corresponding to calibration zone: $(calibzones[minind]) for $dir")
			end
			ddinterps = Dict()
			for i = 1:length(calibzones)
				times = map(j->floats[j][1], 1:length(floats))
				drawdowns = map(j->initialheads[i] - floats[j][columns[i]], 1:length(floats))
				ddinterps[calibzones[i]] = Interpolations.interpolate((times,), drawdowns, Interpolations.Gridded(Interpolations.Linear()))
			end
			info("open obs")
			f = open(string(dir, ".obs"), "w")
			for t = 0:1830
				write(f, @sprintf("%.3f time ", Float64(t)))
				for j = 1:length(calibzones)
					write(f, string(ddinterps[calibzones[j]][t], " "))
				end
				write(f, "\n")
			end
		else#this if/else is needed for when FEHM fails to run -- if it fails to run, we set all the drawdowns to one million, so the optimization doesn't go there
			info("open obs bad")
			f = open(string(dir, ".obs"), "w")
			for time = 0:1830
				write(f, @sprintf("%.3f time ", Float64(time)))
				for j = 1:length(calibzones)
					write(f, "1.0e6 ")
				end
				write(f, "\n")
			end
		end
		close(f)
	end

	function postprocess(dir::Symbol)
		node2zone, zone2nodes = parsezones()
		dir = string(dir)
		info("postproc")
		f = open(string(dir, "/w01_head_his.dat"))
		lineswl = myreadlines(f)
		close(f)
		f = open(string(dir, "/w01_Cr.dat"))
		linescr = myreadlines(f)
		close(f)
		nodeheaders = map(x -> parse(Int, x[2:8]), split(lineswl[2], "Node")[2:end])
		timesandwls = map(line -> map(float, split(line, " "; limit=0, keep=false)[1:end-1]), lineswl[5:end])
		timesandcrs = map(line -> map(float, split(line, " "; limit=0, keep=false)[1:end-1]), linescr[5:end])
		#find the columns that corresponds to the zones we want to calibrate against
		columns = Array{Int64}(length(calibzones))
		calibzonewls = Array{Float64}(length(calibzones))
		calibzonecrs = Array{Array{Float64}}(length(calibzones))
		for i = 1:length(calibzones)
			calibzonecrs[i] = zeros(length(crtimedict[calibzones[i]]))
			goodnodes = zone2nodes[parse(Int, calibzones[i])]
			runningtotalwl = 0.
			numterms = 0
			for j = 1:length(nodeheaders)
				if nodeheaders[j] in goodnodes
					runningtotalwl += timesandwls[end][j + 1]
					k = 1
					l = 1
					for year in crtimedict[calibzones[i]]
						fehmday = year2fehmday(year)
						#go through the rows looking for the right time
						while abs(fehmday - timesandcrs[k][1]) > 1
							k += 1
						end
						calibzonecrs[i][l] += timesandcrs[k][j + 1]
						l += 1
					end
					numterms += 1
				end
			end
			if length(goodnodes) != numterms
				warn("$(numterms) nodes were found in zone $(calibzones[i]), but there are $(length(goodnodes)) nodes in that zone. Do you mean to include only the top one?")
			end
			calibzonewls[i] = runningtotalwl / numterms
			calibzonecrs[i] = calibzonecrs[i] / numterms
		end
		f = open("allobs.yaml", "w")
		concobs = Dict()
		for j = 1:length(calibzones)
			for i = 1:length(calibzonecrs[j])
				key = string("conc_", calibwells[j], "_", crtimedict[calibzones[j]][i])
				val = calibzonecrs[j][i]
				concobs[key] = val
				write(f, string(key, ": ", val, "\n"))
			end
		end
		slopes = gettrends(concobs)
		for key in keys(slopes)
			write(f, string(key, ": ", slopes[key], "\n"))
		end
		for j = 1:length(calibzones)
			write(f, string("wl_", calibwells[j], ": ", calibzonewls[j], "\n"))
		end
		close(f)
	end

	function preprocess(dir::AbstractString)
		info("preprocessor setup $dir")
		try
			rm(dir, force=true, recursive=true)
			mkdir(dir)
			symlinkdir("w01.data", dir)
			symlinkdir("w01.files", dir)
			symlinkdir("w01.node", dir)
			symlinkdir("w01.node1", dir)
			symlinkdir("w01.nodes", dir)
			symlinkdir("w01.ppor", dir)
			symlinkdir("w01.hyco", dir)
			symlinkdir("w01.flow", dir)
			symlinkdir("nop.temp", dir)
			symlinkdir("infiltration.zonn", dir)
			symlink("../" * dir * ".boun", dir * "/w01.boun")
		catch
			warn("preprocessor setup $dir failed")
		end
	end

	function preprocess(dir::Symbol)
		info("preprocessor setup $dir symbol")
		dir = string(dir)
		preprocess(dir)
		try
			rm(dir * "/w01.data", force=true)
			symlink("../" * dir * ".data", dir * "/w01.data")
			rm(dir * "/w01.files", force=true)
			symlink("../" * dir * ".files", dir * "/w01.files")
			symlinkdir("w01.fin-i", dir)
			symlinkdir("wl.trac", dir)
			rm(dir * "/w01.flow", force=true)
			symlinkdir("wl.flow", dir)
		catch
			warn("preprocessor setup $dir symbol failed")
		end
	end

	function runmodelindir(dir::Symbol)
		runmodelindir(string(dir))
	end

	function runmodelindir(dir::AbstractString)
		info("run $dir")
		run(pipeline(`bash -c "cd $dir; xfehm w01.files"`, stdout=DevNull, stderr=DevNull))
		# run(`bash -c "cd $dir; xfehm w01.files"`)
		# Mads.runcmd(`bash -c "cd $dir; xfehm w01.files"`)
	end
 *
	function symlinkdir(filename::String, dir::String)
		symlink(abspath(filename), dir * "/" * filename)
	end
end

function setupmaindir()
	try
		rm("w01.hyco", force=true)
		rm("w01.ppor", force=true)
		rm("w01.nodes", force=true)
	catch
		None
	end
	info("setupmaindir")
	#set up the .trac file
	tracfile = open("wl.trac", "w")
	tracheader = open("wl.trac.header")
	lines = myreadlines(tracheader)
	close(tracheader)
	for line in lines
		write(tracfile, line)
	end
	wlflowfile = open("wl.flow", "w")
	wlflowheaderlines = myreadlines("wl.flow.header")
	for line in wlflowheaderlines
		write(wlflowfile, line)
	end
	nodexyz = readdlm("../../smoothgrid2/out_top.nodes_xyz")
	nodes = map(x->round(Int, x), nodexyz[:, 1])
	xs = nodexyz[:, 2]
	ys = nodexyz[:, 3]
	sourceinfo = readdlm("sourceinfo")
	sourcedists = Array{Distributions.Distribution}(size(sourceinfo, 1))
	infils = Array{Float64}(size(sourceinfo, 1))
	strengths = Array{Float64}(size(sourceinfo, 1))
	t0s = Array{Float64}(size(sourceinfo, 1))
	tfs = Array{Float64}(size(sourceinfo, 1))
	x0s = Array{Float64}(size(sourceinfo, 1))
	y0s = Array{Float64}(size(sourceinfo, 1))
	for i = 1:size(sourceinfo, 1)
		strength, x0, y0, xr, yr, cov, t0, tf, infil = sourceinfo[i, :]
		strengths[i] = strength
		t0s[i] = t0
		tfs[i] = tf
		x0s[i] = x0
		y0s[i] = y0
		sourcedists[i] = Distributions.MvNormal([x0, y0], [xr * xr cov; cov yr * yr])
		infils[i] = infil
	end
	for j = 1:length(nodes)
		x, y = xs[j], ys[j]
		writeit = false
		t0 = Inf
		tf = -Inf
		cmax = 0.
		sumctimesinfil = 0.
		totalinfil = 0.
		for i = 1:length(sourcedists)
			if Distributions.pdf(sourcedists[i], [x, y]) / Distributions.pdf(sourcedists[i], [x0s[i], y0s[i]]) >= 0.6065306597126334 # go out 1 standard deviation
				c = strengths[i] * Distributions.pdf(sourcedists[i], [x, y]) / Distributions.pdf(sourcedists[i], [x0s[i], y0s[i]])
				totalinfil += infils[i]
				sumctimesinfil += c * infils[i]
				t0 = min(t0, t0s[i])
				tf = max(tf, tfs[i])
				writeit = true
			end
		end
		if writeit
			write(tracfile, string(nodes[j], " ", nodes[j], " 1 ", sumctimesinfil / totalinfil, " ", t0, " ", tf, "\n"))
			write(wlflowfile, string(nodes[j], " ", nodes[j], " 1 ", totalinfil, " 1 0\n"))
		end
	end
	write(tracfile, "-3 0 0 5 0 23376\n")
	write(tracfile, "-5 0 0 5 0 23376\n")
	write(tracfile, "-336000 0 0 5 0 23376\n\n")
	close(tracfile)
	write(wlflowfile, "\n")
	close(wlflowfile)

	#set up the nodefile
	_, zone2nodes = parsezones()
	nodes = Int[]
	for k in keys(zone2nodes)
		for node in zone2nodes[k]
			push!(nodes, node)
		end
	end
	f = open("w01.nodes", "w")
	write(f, "nodes\n")
	write(f, "$(length(nodes))\n")
	for i = 1:length(nodes)
		write(f, "$(nodes[i])\n")
	end
	write(f, "\n\n")
	close(f)

	info("hyco")
	#read the parameter file that is output by mads
	f = open("w01.hyco_hetero")
	line = chomp(readline(f))
	vals = map(x -> float(x), split(line, " "))#load the hyco vals
	kx = vals[1:div(end, 3)]
	ky = vals[div(end, 3) + 1:2 * div(end, 3)]
	kz = vals[2 * div(end, 3) + 1:end]
	close(f)
	f = open("pilotpoints.txt")
	pilotpoints = transpose(readdlm(f))
	close(f)
	xyz = readdlm("../../smoothgrid2/tet.xyz")'
	krigingparams = readdlm("krigingparams")
	mycov = h->Kriging.sphericalcov(h, krigingparams[1], krigingparams[2])
	allkxs = Kriging.krige(xyz, pilotpoints, kx, mycov)
	allkys = Kriging.krige(xyz, pilotpoints, ky, mycov)
	allkzs = Kriging.krige(xyz, pilotpoints, kz, mycov)
	#write the file that is read by fehm
	f = open("w01.hyco", "w")
	write(f, "hyco\n")
	for i = 1:length(allkxs)
		write(f, "$i $i 1 $(10 ^ allkxs[i]) $(10 ^ allkys[i]) $(10 ^ allkzs[i])\n")
	end
	write(f, "\n")
	close(f)
	#now do the ppor
	info("ppor")
	f = open("w01.stor_hetero")
	line = chomp(readline(f))
	ppor = map(x -> float(x), split(line, " "))#vals should contain the value of the xy hycos at each of the pilot points follow by a single value giving the factor for the vertical hyco
	close(f)
	allppors = Kriging.krige(xyz, pilotpoints, ppor, mycov)
	f = open("w01.ppor", "w")
	write(f, "ppor\n")
	write(f, "-1\n")
	for i = 1:length(allppors)
		write(f, "$i $i 1 $(10 ^ allppors[i])\n")
	end
	write(f, "\n")
	close(f)
end
