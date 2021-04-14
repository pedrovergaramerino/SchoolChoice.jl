* FIGURE 5 - Goodness of Fit: Observed vs. Simulated Cutoffs

*** Execute Matlab script figure_5.m before running this program ***

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute figure_5.do. This program requires proprietary data - check ReadMe file}"
	exit
	}
else {
capture confirm file "${folder}/figures/figure_5.xlsx"
if _rc { 
	di "{red: ERROR: Cannot execute figure_5.do. Execute Matlab script figure_5.m before running this program}"
	exit
	}
}

import excel using "${folder}\figures\figure_5.xlsx", firstrow

twoway ///
(scatter observed id, connect(l) lcolor(gs9) msymbol(circle_hollow) mcolor(gs9) msize(medlarge) ) ///
(scatter TT id,  lp(solid) msymbol(dh) mcolor(navy) lcolor(navy) msize(medlarge) ) ///
(scatter ST id,  msymbol(X)  mcolor(red) msize(medlarge) ) ///
(scatter MEI id,  msymbol(+)  mcolor(green) msize(medlarge)), ///
xlabel(1 `"s1"' 2 `"s2"' 3 `"s3"' 4 `"s4"' 5 `"s5"' 6 `"s6"' 7 `"s7"' 8 `"s8"' 9 `"s9"' 10 `"s10"' 11 `"s11"') ///
xtitle("School", size(medlarge) height(5)) ///
ytitle("Observed / Simulated Cutoff", size(medlarge) height(5)) ///
graphregion( fcolor(white) lcolor(white)) ///
ylabel(0(0.1)0.8, format(%2.1fc)) ///
legend ( pos(11) ring(0) col(1) rowgap(2) size(medsmall) ///
	   order(1 "Observed cutoffs" - "Simulated cutoffs from estimates based on:" 2 "(1) Weak truth-telling (WTT)" 3 "(2) Stability" 4 "(3) Stability and undominated strategies") )	  

graph export "${folder}\figures\figure_5.pdf", replace

graph close _all
