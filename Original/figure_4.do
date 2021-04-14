* FIGURE 4 - School Cutoffs in 2012 and 2013

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute figure_4.do. This program requires proprietary data - check ReadMe file}"
	exit
	}

set more off

use "${folder}\data\district_sud.dta", clear

** Keeping only cutoff variables
keep sch_id sch_cutoff sch_cutoff_12
bysort sch_id : keep if _n==1
replace sch_cutoff=round(sch_cutoff, 0.001)

** Creating labels for graphic
gen sch_label=""
forvalues i=1(1)11 {
replace sch_label="s `i'" if sch_id==`i'
}

** Controling the position of label
capture drop pos
gen pos=1
replace pos=12 if sch_i==2
replace pos=3 if sch_id==3
replace pos=4 if sch_id==1
twoway (scatter  sch_cutoff sch_cutoff_12,  sort msymbol(o) mcolor(black) mlabel(sch_label)  mlabvposition(pos) mlabsize(vsmall) mlabcolor(black)) ///
 (function y=x, range(0 0.8) msize(small) lcolor(gs10)  lpattern(dash)), legend(off) ///
ytitle("School cutoffs in 2013",size(small)) ylabel(,labsize(small)) xtitle("School cutoffs in 2012",size(small)) xlabel(,labsize(small)) ///
 graphregion(color(white))  bgcolor(white)  aspectratio(1)
 
graph export "${folder}\figures\figure_4.pdf", replace

graph close _all


 
 
