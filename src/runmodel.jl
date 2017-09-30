import RobustPmap
proc = spawn(`ls`)
timedwait(() -> process_exited(proc), 10.) # 10 seconds
if process_running(proc)
	kill(proc)
	exit()
end
info("start runmodel.jl")
include("runmodelcode.jl")
info("setup")
setupmaindir()
info("pre")
checknode()
RobustPmap.rpmap(preprocess, dirs)
info("run")
checknode()
RobustPmap.rpmap(runmodelindir, dirs)
info("post")
checknode()
RobustPmap.rpmap(postprocess, dirs)
info("done")
:done
