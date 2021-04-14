* FIGURE 2 - Monte Carlo Simulations: Impact of Economy Size on the Equilibrium Distribution of Cutoffs (Constrained / Truncated DA)

*** Execute Matlab script figure_2.m before running this program ***

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/figures/figure_2a.xlsx"
if _rc { 
	di "{red: ERROR: Cannot execute figure_2.do. Execute Matlab script figure_2.m before running this program}"
	exit
	}

do "./do/school_cutoffs.do"

* Figure 2(a): 100 students, 6 schools

import excel using "${folder}\figures\figure_2a.xlsx", firstrow clear
Cutoffs_graph_new, ceiling(20) step(5) s2(4.25 0.18 6.5 0.18) samp(100)
graph export "${folder}\figures\figure_2a.pdf", replace

* Figure 2(b): 500 students, 6 schools

import excel using "${folder}\figures\figure_2b.xlsx", firstrow clear
Cutoffs_graph_new, ceiling(20) step(5) s2(7.4 0.16 15.2 0.16) samp(500)
graph export "${folder}\figures\figure_2b.pdf", replace

graph close _all
