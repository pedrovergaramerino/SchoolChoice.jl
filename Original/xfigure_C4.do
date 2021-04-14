* FIGURE C4 - Monte Carlo Simulations: Impact of the Marginal Cost of Applying to Schools on Equilibrium Outcomes (500 Students, 6 Schools)

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/figures/xfigure_C4.xlsx"
if _rc { 
	di "{red: ERROR: Cannot execute xfigure_C4.do. Execute Matlab script xfigure_C4.m before running this program}"
	exit
	}

import excel using "${folder}\figures\xfigure_C4.xlsx", firstrow

twoway ///
(scatter pct_stable id, connect(l) lcolor(navy) msymbol(d) mcolor(navy) msize(medlarge) yaxis(1)) ///
(scatter pct_WTT id, connect(l) lcolor(red) msymbol(circle_hollow) mcolor(red) msize(medlarge) yaxis(1)) ///
(scatter m_ROL_length id, connect(l) lcolor(green) msymbol(S) mcolor(green) msize(medlarge) yaxis(2)), ///
xlabel(1 `"0"' 2 `"1e-6"' 3 `"1e-3"' 4 `"0.01"' 5 `"0.1"' 6 `"1.0"') ///
xtitle("Marginal application cost", size(medlarge) height(5)) ///
graphregion( fcolor(white) lcolor(white)) ///
ylabel(0(0.2)1, format(%2.1fc)) ///
ylabel(0(1)6 , axis(2)) ///
ytitle("", axis(2)) ///
legend ( pos(6) ring(1) bmargin(small) col(1) rowgap(3) size(medsmall)  ///
	   order(3 "Average length of ROL (right axis)" 1 "Share of students assigned to favorite feasible school (left axis)" 2 "Share of weakly truth-telling students (left axis)" ///
	   ) )
	   
graph export "${folder}\figures\xfigure_C4.pdf", replace
