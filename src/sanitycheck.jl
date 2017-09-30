import FEHM

if length(ARGS) == 0
	defaultmadsfilename = "cond11/w01-v10.mads"
	push!(ARGS, defaultmadsfilename)
	warn("you didn't give me a filename, so I'm going with the default -- $defaultmadsfilename")
end
if length(ARGS) < 2
	defaultgeofilename = "../../../smoothgrid2/tet.geo"
	push!(ARGS, defaultgeofilename)
	warn("you didn't give me a geo filename, so I'm going with the default -- $defaultgeofilename")
end
if !isfile(ARGS[1])
	error("file $(ARGS[1]) is missing!")
end
sp = splitdir(ARGS[1])
workingdir = joinpath(sp[1:end - 1]...)
madsfilename = sp[end]
geofilename = ARGS[2]
import Mads

startdir = pwd()

cd(workingdir)
if !isdefined(:md)
	info("Loading Mads input file ..")
	md = Mads.loadmadsfile(madsfilename)
end

info("Plot the model predictions ...")
if !isdefined(:results)
	results = Mads.forward(md)
end
if !isdir("figs")
	mkdir("figs")
end
include("$workingdir/plotresults.jl")

info("Generate the AVS files ...")
Mads.writeparameters(md)
rmslines = readlines("runmodel.jl")
goodrmslines = filter(x->!contains(x, "pmap("), rmslines)
if length(rmslines) != length(goodrmslines) + 3
	error("there are more than 3 pmaps in your run-model-split.jl, but this script assumes there are only 3")
end
textcode = join(goodrmslines, "")
#eval(parse(textcode))
f = open("temp_run-model-split.jl", "w")
write(f, textcode)
close(f)
include("$workingdir/temp_run-model-split.jl")
rm("temp_run-model-split.jl")
map(preprocess, dirs)
rm("wl/w01.data")
datalines = readlines("wl.data")
f = open("wl/w01.data", "w")
write(f, datalines[1])
write(f, "cont\navs 3 1e20\nc\nendavs\n")
for i = 2:length(datalines)
	write(f, datalines[i])
end
close(f)
cd("wl")

info("Run FEHM ...")
if !isfile("w01.00001_con_node.avs")
	run(pipeline(`xfehm w01.files`, stdout=DevNull, stderr=DevNull))
end
run(`ln -s $geofilename w01.00001_geo`)
cd("../..")

include("makemovie.jl")

cd(startdir)

include("plotsources.jl")
