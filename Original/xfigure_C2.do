* FIGURE C2 - Monte Carlo Simulations: Equilibrium Distribution of School Cutoffs (6 schools, 500 students)

*** Execute Matlab script figure_C2.m before running this program ***

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/figures/xfigure_C2a.xlsx"
if _rc { 
	di "{red: ERROR: Cannot execute xfigure_C2.do. Execute Matlab script xfigure_C2.m before running this program}"
	exit
	}
	
do "./do/school_cutoffs.do"

* Figure C2(a): Constrained/truncated DA

import excel using "${folder}\figures\xfigure_C2a.xlsx", firstrow clear
Cutoffs_graph_new, ceiling(20) step(5) s2(7.4 0.16 15.2 0.16) samp(500)
graph export "${folder}\figures\xfigure_C2a.pdf", replace


* Figure C2(b): Unconstrained DA with cost

import excel using "${folder}\figures\xfigure_C2b.xlsx", firstrow clear
Cutoffs_graph_new, ceiling(20) step(5) s2(7.4 0.16 15.2 0.16) samp(500)
graph export "${folder}\figures\xfigure_C2b.pdf", replace
