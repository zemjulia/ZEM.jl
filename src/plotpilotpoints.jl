import CrPlots

basedir = "cond01"
numcondpoints = 37
pps = readdlm("$basedir/pilotpoints.txt")
x0 = CrPlots.bgx0
x1 = CrPlots.bgx1
y0 = CrPlots.bgy0
y1 = CrPlots.bgy1
boundingbox = (x0, y0, x1, y1)
fig, ax = CrPlots.crplot(boundingbox)
plotwells = ["R-42", "R-11","R-44","R-35b","R-36","R-62","R-33","SIMR-2","R-50","R-28","R-45","R-13","CrEX-1","CrEX-3","R-01","R-15","R-43","R-61", "PM-03", "PM-04", "PM-05", "O-04"]
CrPlots.addwells(ax, plotwells)
#CrPlots.addpoints(ax, pps[1:numcondpoints, 1:2]'; colorstring="r.", markersize=30)
CrPlots.addpoints(ax, pps[numcondpoints+1:end, 1:2]'; colorstring="b.", markersize=15)
display(fig); println()
fig[:savefig]("$basedir/pilotpoints.png")
PyPlot.close(fig)
