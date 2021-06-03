import JLD

function setscenarios(filename::AbstractString="scenarios.jld")
	colors = ["b", "g", "r", "c", "m", "y", "k", "b--", "g--", "r--", "c--", "m--", "y--", "k--"]
	scenarios = Dict()
	scenarioname = "no_action"
	scenarios[scenarioname] = Dict()
	scenarios[scenarioname]["legend"] = "No Action"
	scenarios[scenarioname]["wellsinuse"] = []
	scenarios[scenarioname]["pntstartyear"] = 2017
	scenarios[scenarioname]["gpmrates"] = []
	scenarios[scenarioname]["duration"] = []
	scenarios[scenarioname]["color"] = "b"

	colorindex = 2
	scenarioname = "scenario_1_20170916"
	scenarios[scenarioname] = Dict()
	scenarios[scenarioname]["legend"] = "Test 1: CrIN 1,3,5 injection"
	scenarios[scenarioname]["wellsinuse"] = ["CrEX-1", "CrEX-2", "CrEX-3", "CrIN-1", "CrIN-3", "CrIN-5"]
	scenarios[scenarioname]["pntstartyear"] = 2017
	scenarios[scenarioname]["gpmrates"] = [75, 75, 75, -75, -75, -75]
	scenarios[scenarioname]["duration"] = [7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25]
	scenarios[scenarioname]["color"] = colors[colorindex]; colorindex += 1
	scenarioname = "scenario_2_20170916"
	scenarios[scenarioname] = Dict()
	scenarios[scenarioname]["legend"] = "Test 2: CrIN 2,4,6 injection"
	scenarios[scenarioname]["wellsinuse"] = ["CrEX-1", "CrEX-2", "CrEX-3", "CrIN-2", "CrIN-4", "CrIN-6"]
	scenarios[scenarioname]["pntstartyear"] = 2017
	scenarios[scenarioname]["gpmrates"] = [75, 75, 75, -75, -75, -75]
	scenarios[scenarioname]["duration"] = [7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25, 7 / 365.25]
	scenarios[scenarioname]["color"] = colors[colorindex]; colorindex += 1
	#=
	scenarioname = "scenario_3_20170916"
	scenarios[scenarioname] = Dict()
	scenarios[scenarioname]["legend"] = "CrIN-6 injection"
	scenarios[scenarioname]["wellsinuse"] = ["CrIN-6"]
	scenarios[scenarioname]["pntstartyear"] = 2017
	scenarios[scenarioname]["gpmrates"] = [75]
	scenarios[scenarioname]["duration"] = [1]
	scenarios[scenarioname]["color"] = colors[colorindex]; colorindex += 1
	scenarioname = "scenario_4_20170916"
	scenarios[scenarioname] = Dict()
	scenarios[scenarioname]["legend"] = "CrIN-6 extraction"
	scenarios[scenarioname]["wellsinuse"] = ["CrIN-6"]
	scenarios[scenarioname]["pntstartyear"] = 2017
	scenarios[scenarioname]["gpmrates"] = [-75]
	scenarios[scenarioname]["duration"] = [1]
	scenarios[scenarioname]["color"] = colors[colorindex]; colorindex += 1
	=#
	JLD.save(filename, scenarios)
	return scenarios
end
