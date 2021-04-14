* FIGURE C3 - Monte Carlo Simulations: Impact of Economy Size on the Equilibrium Distribution of Cutoffs (Constrained/Truncated DA)

*** Execute Matlab script figure_C4.m before running this program ***

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/figures/xfigure_C3a.xlsx"
if _rc { 
	di "{red: ERROR: Cannot execute xfigure_C3.do. Execute Matlab script xfigure_C3.m before running this program}"
	exit
	}
	
do "./do/school_cutoffs.do"

* Figure C3(a): 100 students

import excel using "${folder}\figures\xfigure_C3a.xlsx", firstrow clear
Cutoffs_graph, ceiling(50) step(10) s2(3.7 0.13 10.5 0.13) samp(100)
graph export "${folder}\figures\xfigure_C3a.pdf", replace

* Figure C3(b): 200 students

import excel using "${folder}\figures\xfigure_C3b.xlsx", firstrow clear
Cutoffs_graph, ceiling(50) step(10) s2(3.7 0.13 12.5 0.13) samp(200)
graph export "${folder}\figures\xfigure_C3b.pdf", replace

* Figure C3(c): 500 students

import excel using "${folder}\figures\xfigure_C3c.xlsx", firstrow clear
Cutoffs_graph, ceiling(50) step(10) s2(6.5 0.15 20.5 0.15) samp(500)
graph export "${folder}\figures\xfigure_C3c.pdf", replace

* Figure C3(d): 5,000 students

import excel using "${folder}\figures\xfigure_C3d.xlsx", firstrow clear
Cutoffs_graph, ceiling(50) step(10) s2(30 0.185 40.5 0.185)  samp(5000)
graph export "${folder}\figures\xfigure_C3d.pdf", replace
