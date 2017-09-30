import Mads
using LaTeXStrings
import PyPlot
import PyCall
@PyCall.pyimport aquiferdb as db
rundirs = ["anew01"]
filenames = ["w01-v04.mads"]
if !isdefined(:mds)
	doruns = true
	info("Reading mads files ...")
	mds = Array(Any, length(rundirs))
	for i in 1:length(rundirs)
		info("Reading problem $(rundirs[i]) ...")
		mds[i] = Mads.loadmadsfile("$(rundirs[i])/$(filenames[i])")
		mds[i]["Filename"] = filenames[i]
	end
end
wellnames = ["R-11", "R-13", "R-15", "R-28", "R-33#1", "R-33#2", "R-35b", "R-36", "R-42", "R-43#1", "R-43#2", "R-44#1", "R-44#2", "R-45#1", "R-45#2", "R-50#1", "R-50#2", "R-61#1", "R-61#2", "R-62", "CrEX-1", "CrEX-3", "CrIN-1", "CrIN-2", "CrIN-3", "CrIN-4", "CrIN-5"]
plotgroups = Any[["r11"], ["r13"], ["r15"], ["r28"], ["r33_1"], ["r33_2"], ["r35b"], ["r36"], ["r42"], ["r43_1"], ["r43_2"], ["r44_1"], ["r44_2"], ["r45_1"], ["r45_2"], ["r50_1"], ["r50_2"], ["r61_1"], ["r61_2"], ["r62"], ["crex1"], ["crex3"], ["crin1"], ["crin2"], ["crin3"], ["crin4"], ["crin5"]]

if isdefined(:doruns) && doruns
	numsamples = 2000
	seed = 2016
	samples = Array(Any, length(rundirs))
	llhoods = Array(Any, length(rundirs))
	fs = Array(Any, length(rundirs))
	for i = 1:length(rundirs)
		cd(rundirs[i])
		mds[i]["Restart"] = true
		fs[i] = Mads.forward(mds[i])
		cd("..")
	end
	doruns = false
end

if !isdefined(:organizedata)
	organizedata = true
end
if organizedata
	conc_keys = Mads.filterkeys(mds[1]["Observations"], "conc")
	fgroupobs = map(rd->Dict(), rundirs)
	grouptimes = map(wellname->Int[], wellnames)
	obs2time(obsname) = parse(Int, split(obsname, "_")[end])
	for i = 1:length(conc_keys)
		for j = 1:length(wellnames)
			for k = 1:length(plotgroups[j])
				if ismatch(Regex(string("conc_", plotgroups[j][k], "_\\d+")), conc_keys[i])
					push!(grouptimes[j], obs2time(conc_keys[i]))
				end
			end
		end
	end
	grouptimes = map(x->sort(unique(filter(t->t<=2017, x))), grouptimes)
	for i = 1:length(wellnames)
		for j = 1:length(rundirs)
			fgroupobs[j][wellnames[i]] = Dict()
		end
		for t in grouptimes[i]
			for j = 1:length(rundirs)
				fgroupobs[j][wellnames[i]][t] = median(map(k->fs[j][string("conc_", plotgroups[i][k], "_", t)], 1:length(plotgroups[i])))
			end
		end
	end
	organizedata = false
end

if !isdefined(:dodatabase)
	dodatabase = true
end
if dodatabase
	dbobs = Dict()
	dbobstimes = Dict()
	wellname(fullname) = split(fullname, "#")[1]
	portdesc(fullname) = length(split(fullname, "#")) == 1 ? "SINGLE COMPLETION" : split(fullname, "#")[2]
	db.connecttodb()
	for i = 1:length(wellnames)
		if wellnames[i][1] == 'R'
			dbobs[wellnames[i]], dbobstimes[wellnames[i]] = db.getchromiumconcentrations(wellname(wellnames[i]), portdesc(wellnames[i]))
			ind = find(dbobstimes[wellnames[i]].>Date(2013, 7))
		end
	end
	db.disconnectfromdb()
	dodatabase = false
end

plotrows = 5
plotcols = 6
facecolors = ["k", "r", "g", "b", "y", "c"]
fig, axs = PyPlot.subplots(plotrows, plotcols, sharey=false, sharex=true, figsize=(18, 10))
for j = 1:length(grouptimes)
	for i = 1:length(rundirs)
		axs[j][:plot](map(year->Date(year, 7), grouptimes[j]), map(t->fgroupobs[i][wellnames[j]][t], grouptimes[j]), facecolors[i])
	end
	if wellnames[j][1] == 'R'
		axs[j][:plot](dbobstimes[wellnames[j]], dbobs[wellnames[j]], "r.")
	end
	axs[j][:tick_params](axis="both", which="major", labelsize=8)
	axs[j][:set_title](wellnames[j])
	axs[j][:set_ylabel](L"Cr [ppb]")
end
preds = [ [1, 2, 3, 4] [2, 3, 4, 5] [2, 2.1, 2.2, 2.3] [0.6, 0.9, 1.2, 1.3] [0.8, 1.2, 1.4, 1.6] ]
axs[end-2][:set_title]("Legend: predictions")
for i = 1:length(rundirs)
	axs[end-2][:plot](map(year->Date(year, 7), [2001, 2006, 2011, 2016]), preds[:,i], facecolors[i], label=rundirs[i])
end
axs[end-2][:legend](fontsize=4)
axs[end-2][:tick_params](axis="both", which="major", labelsize=8)
axs[end][:set_title]("Legend: observations")
axs[end][:plot](map(year->Date(year, 7), collect(2001:2016)), log10(.5 * (9.99 * collect(linspace(1, 4, 16))) + .5 * (10.01 + 10 * collect(linspace(1, 4, 16))) + 5 * randn(16)), "r.")
axs[end][:tick_params](axis="both", which="major", labelsize=8)
PyPlot.setp(axs, xticks=map(year->Date(year, 7), [2001, 2006, 2011, 2016]))
PyPlot.tight_layout()
display(fig); println()
PyPlot.savefig(rundirs[1] * "/waffle2014.pdf")
PyPlot.close(fig)
