import PyPlot
import YAML
import RobustPmap
import Mads
@everywhere include("setscenarios.jl")
@everywhere include("myreadlines.jl")

@everywhere begin
	names2zones = Dict()
	names2zones["CrIN-1"] = 612101
	names2zones["CrIN-2"] = 611201
	names2zones["CrIN-3"] = 610301
	names2zones["CrIN-4"] = 618401
	names2zones["CrIN-5"] = 610501
	names2zones["CrIN-6"] = 610601
	names2zones["CrEX-1"] = 620101
	names2zones["CrEX-2"] = 620206
	names2zones["CrEX-3"] = 620301
	names2zones["R-61#1"] = 4610001
	names2zones["R-61#2"] = 4610002
	modelstartyear = 1964
	function runscenario(scenarioname, wellsinuse, pntstartyear, gpmrates, duration)
		gpm2kgpersec = 0.06309
		kgpersecrates = gpmrates * gpm2kgpersec

		println("$scenarioname")
		if !isdir("wl_$scenarioname")# || true
			initialdir = pwd()
			cp("wl", "wl_$scenarioname"; remove_destination=true)
			cd("wl_$scenarioname")
			info("Set up FEHM")
			rm("w01.data")
			datalines = myreadlines("../wl.data")
			f = open("w01.data", "w")
			write(f, datalines[1])
			write(f, "cont\navs 3 1e20\nc\nendavs\n")
			i = 2
			while datalines[i] != "time\n"
				write(f, datalines[i])
				i += 1
			end
			write(f, datalines[i])#write "time\n"
			i += 1
			write(f, datalines[i])#write group 1 of time macro
			i += 1
			pntstarttime = (pntstartyear - modelstartyear) * 365.25
			pntendtimes = pntstarttime + scenarios[scenarioname]["duration"] * 365.25
			dittimes = unique([[pntstarttime]; pntendtimes])
			ditsteps = map(x->365.25, dittimes)
			while datalines[i] != "\n"#write dit groups of time macro
				splitline = split(datalines[i])
				currentday = parse(Float64, splitline[1])
				if (currentday / 365.25) + modelstartyear + 1e-2 >= pntstartyear
					for j = 0:19
						if !(currentday + .5 * 36.525 * j in dittimes)
							push!(dittimes, currentday + .5 * 36.525 * j)
							push!(ditsteps, .5 * 36.525)
						end
					end
				else
					if !(parse(Float64, split(datalines[i])[1]) in dittimes)
						push!(dittimes, parse(Float64, split(datalines[i])[1]))
						push!(ditsteps, parse(Float64, split(datalines[i])[end]))
					end
				end
				i += 1
			end
			ditsortedindices = sort(collect(1:length(dittimes)), by=i->dittimes[i])
			for j = 1:length(dittimes)
				write(f, "$(dittimes[ditsortedindices[j]]) -1 1 1 $(ditsteps[ditsortedindices[j]])\n")
			end
			for j = i:length(datalines)
				write(f, datalines[j])
			end
			close(f)
			rm("w01.boun")
			bounlines = myreadlines("../wl.boun")
			bounendline = filter(i->bounlines[i] == "end\n", 1:length(bounlines))[1]
			f = open("w01.boun", "w")
			for i = 1:bounendline - 1
				if !contains(bounlines[i], "year")
					write(f, bounlines[i])
				end
			end
			for i = 1:length(wellsinuse)
				write(f, "model $(i + 2) $(wellsinuse[i])\n")
				write(f, "tran\ntime\n3\n")
				write(f, "0\n$((pntstartyear - modelstartyear) * 365.25)\n$((pntstartyear - modelstartyear + duration[i]) * 365.25)\n")
				write(f, "dsw\n0\n$(kgpersecrates[i])\n0\n")
			end
			write(f, "end\n")
			for i = 1:length(wellsinuse)
				write(f, "-$(names2zones[wellsinuse[i]]) 0 0 $(i + 2)\n")
			end
			for i = bounendline+1:length(bounlines)
				write(f, bounlines[i])
			end
			close(f)
			rm("wl.trac")
			traclines = myreadlines("../wl.trac")
			f = open("wl.trac", "w")
			write(f, traclines[1])
			write(f, "0.0 1.1 1.0E-3 1.\n")
			for i = 3:length(traclines) - 1
				write(f, traclines[i])
			end
			for i = 1:length(wellsinuse)
				if gpmrates[i] < 0
					write(f, "-$(names2zones[wellsinuse[i]]) 0 0 5 $((pntstartyear - modelstartyear) * 365.25) $((pntstartyear - modelstartyear + duration[i]) * 365.25)\n")
				end
			end
			write(f, "\n")
			close(f)
			info("Run FEHM ...")
			@time run(pipeline(`xfehm w01.files`, stdout=DevNull, stderr=DevNull))
			run(`ln -s $geofilename w01.00001_geo`)
			cd(initialdir)
		else
			println("skip running $scenarioname")
		end
	end

	function postprocess_full(dir::Symbol)
		node2zone, zone2nodes = parsezones()
		dir = string(dir)
		info("postproc_drawdowns")
		f = open(string(dir, "/w01_head_his.dat"))
		lineswl = myreadlines(f)
		close(f)
		f = open(string(dir, "/w01_Cr.dat"))
		linescr = myreadlines(f)
		close(f)
		nodeheaders = map(x -> parse(Int, x[2:8]), split(lineswl[2], "Node")[2:end])
		timesandwls = map(line -> map(float, split(line, " "; limit=0, keep=false)[1:end-1]), lineswl[11:end])
		timesandcrs = map(line -> map(float, split(line, " "; limit=0, keep=false)[1:end-1]), linescr[11:end])
		wltimes = map(x->x[1], timesandwls) / 365.25 + modelstartyear
		crtimes = map(x->x[1], timesandcrs) / 365.25 + modelstartyear
		wls = Dict()
		crs = Dict()
		for cz in calibzones
			wls[cz] = zeros(length(wltimes))
			crs[cz] = zeros(length(crtimes))
		end
		#find the columns that corresponds to the zones we want to calibrate against
		for i = 1:length(calibzones)
			goodnodes = zone2nodes[parse(Int, calibzones[i])]
			numfoundnodes = 0
			for j = 1:length(nodeheaders)
				if nodeheaders[j] in goodnodes
					numfoundnodes += 1
					wls[calibzones[i]] += map(x->x[j + 1], timesandwls)
					crs[calibzones[i]] += map(x->x[j + 1], timesandcrs)
				end
			end
			@assert numfoundnodes == length(goodnodes)
			wls[calibzones[i]] /= length(goodnodes)
			crs[calibzones[i]] /= length(goodnodes)
		end
		return wltimes, wls, crtimes, crs
	end

	function ppscenarios(scenarios)
		cd(workingdir)
		plotdir = "scenario_plots"
		if !isdir(plotdir)
			mkdir(plotdir)
		end
		crs = Dict()
		wls = Dict()
		scenarionames = sort(collect(keys(scenarios)))
		wltimes = Dict()
		crtimes = Dict()
		for scenarioname in scenarionames
			#postprocess(Symbol("wl_" * scenarioname))
			#mv("allobs.yaml", scenarioname * ".yaml"; remove_destination=true)
			#crdata[scenarioname] = YAML.load(open(scenarioname * ".yaml"))
			wltimes[scenarioname], wls[scenarioname], crtimes[scenarioname], crs[scenarioname] = postprocess_full(Symbol("wl_" * scenarioname))
		end
		wells = [["crin$i" for i in 1:6]; ["crin6e", "crin6w", "crex1","crex2","crex3","crpz1","crpz2_1","crpz3","crpz4","crpz5","r11","r13","r15","r28","r33_1","r33_2","r35b","r36","r42","r43_1","r43_2","r44_1","r44_2","r45_1","r45_2","r50_1","r50_2","r61_1","r61_2","r62","simr2"]]
		shortname2zone = Dict(zip(calibwells, calibzones))
		ts = collect(2005:2027)
		for i = 1:length(wells)
			wlfig, wlax = PyPlot.subplots()
			crfig, crax = PyPlot.subplots()
			crfig[:set_size_inches](8, 5)
			maxobs = 0.
			for scenarioname in scenarionames
				#=
				obs = Float64[]
				goodts = Float64[]
				for t in ts
					if haskey(crdata[scenarioname], string("conc_", wells[i], "_", t))
						push!(obs, crdata[scenarioname][string("conc_", wells[i], "_", t)])
						push!(goodts, t)
					end
				end
				=#
				crax[:plot](crtimes[scenarioname], crs[scenarioname][shortname2zone[wells[i]]], scenarios[scenarioname]["color"], linewidth=3, alpha=0.5)
				wlax[:plot](wltimes[scenarioname], wls[scenarioname][shortname2zone[wells[i]]], scenarios[scenarioname]["color"], linewidth=3, alpha=0.5)
				maxobs = max(maxobs, maximum(crs[scenarioname][shortname2zone[wells[i]]]))
			end
			crax[:plot](ts, map(x->50, ts), "k--", linewidth=3)
			crax[:set_title](wells[i])
			crax[:set_xlabel]("Year")
			crax[:set_ylabel]("Concentration [ppb]")
			if maxobs < 50
				crax[:set_ylim](0, 60)
			end
			crax[:set_xlim](ts[1], ts[end])
			crax[:legend](map(x->scenarios[x]["legend"], scenarionames), loc=6)
			crfig[:tight_layout]()
			crfig[:savefig](string("$plotdir/$(case)_btcs_", wells[i], ".pdf"))
			PyPlot.close(crfig)
			wlax[:set_title](wells[i])
			wlax[:set_xlabel]("Year")
			wlax[:set_ylabel]("Water level [m]")
			wlax[:set_xlim](2016.9, 2018)
			wlax[:legend](map(x->scenarios[x]["legend"], scenarionames))
			wlfig[:tight_layout]()
			wlfig[:savefig](string("$plotdir/$(case)_wls_", wells[i], ".pdf"))
			PyPlot.close(wlfig)
		end
	end
end

if length(ARGS) == 0
	defaultmadsfilename = "cond11_scenarios/w01-v10.mads"
	push!(ARGS, defaultmadsfilename)
	warn("Using default Mads file: $defaultmadsfilename")
end
if length(ARGS) < 2
	defaultgeofilename = "../../../smoothgrid2/tet.geo"
	push!(ARGS, defaultgeofilename)
	warn("Using default grid geometry file: $defaultgeofilename")
end
if !isfile(ARGS[1])
	error("Mads file $(ARGS[1]) is missing!")
end
info("Mads file: $(ARGS[1])")
info("Grid file: $(ARGS[2])")

sp = splitdir(ARGS[1])
case = split(sp[end-1], '_')[1]
info("Case: $(case)")
workingdir = abspath(joinpath(sp[1:end - 1]...))
madsfilename = sp[end]
geofilename = ARGS[2]
scenarios = setscenarios(workingdir * ".jld")
@eval @everywhere scenarios = $scenarios
@eval @everywhere geofilename = $geofilename
@eval @everywhere workingdir = $workingdir
@eval @everywhere case = $case

@everywhere startdir = pwd()
@everywhere cd(workingdir)

if !isdefined(:md)
	info("Loading Mads input file ..")
	md = Mads.loadmadsfile(madsfilename)
end
Mads.writeparameters(md)

include(abspath("gettrends.jl"))
include(abspath("checknode.jl"))
include(abspath("runmodelcode.jl"))

push!(calibwells, "simr2")
push!(calibzones, "7050201")
starttimedict["7050201"] = 365
crtimedict["7050201"] = collect(2005:2027)
push!(calibwells, "crin6e")
push!(calibzones, "610603")
starttimedict["610603"] = 365
crtimedict["610603"] = collect(2005:2027)
push!(calibwells, "crin6w")
push!(calibzones, "610602")
starttimedict["610602"] = 365
crtimedict["610602"] = collect(2005:2027)
for k in keys(crtimedict)
	crtimedict[k] = collect(2005:2027)
end

if isfile("w01.hyco")
	warn("skipping main dir setup")
else
	setupmaindir()
	preprocess(:wl)
end

RobustPmap.rpmap(x->runscenario(x, scenarios[x]["wellsinuse"], scenarios[x]["pntstartyear"], scenarios[x]["gpmrates"], scenarios[x]["duration"]), keys(scenarios))
#map(x->runscenario(x, scenarios[x]["wellsinuse"], scenarios[x]["pntstartyear"], scenarios[x]["gpmrates"], scenarios[x]["duration"]), keys(scenarios))
ppscenarios(scenarios)
@everywhere cd(startdir)
