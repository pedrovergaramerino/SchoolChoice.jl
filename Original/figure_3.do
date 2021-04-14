* FIGURE 3 - Fraction of Students Ranking Each of the Four Most Selective Schools in the Southern District of Paris, by Distance to School Cutoff

* N.B.: requires package egenmore
* To install, type "ssc install egenmore" in command window

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute figure_3.do. This program requires proprietary data - check ReadMe file}"
	exit
	}

set more off

use "${folder}\data\district_sud.dta", clear

* Define the best school in terms of cutoffs 
gsort - sch_cut_pts sch_id
egen r_school=group(sch_cut_pts sch_id), missing
egen max_rank=max(r_school)
gen sch_rank=max_rank-r_school+1
sort sch_rank
drop r_school 
label var sch_rank "rank of school by cutoff"

* Drop low income students (always above cutoff due to low income bonus)
drop if stu_lowinc==1


* Generate variables for graphical analysis 

local nb=max_rank 
forvalues i=1(1) `nb'{

	* Rescale distance to the score so that 0 corresponds to the school cutoff in points
	gen dist_cut_`i'=stu_priority_pts-sch_cut_pts if sch_rank==`i'
	label var dist_cut_`i' "distance to cutoff `i'"
	
			
	* Define whether a student has the school of a rank `i' in his list*
	gen sch_rk`i'_ok=(sch_rank==`i' & choice_rk~=0)
	egen sch_rk`i'_in_rol=sum(sch_rk`i'_ok),by(stu_id)
	drop sch_rk`i'_ok
	label var sch_rk`i'_in_rol "school in sudent's rol"
	*su dist_cut_`i', d
		
	*Generate the dummy for being above the school cutoff
	gen d`i'=(dist_cut_`i'>=0 & dist_cut_`i'~=.)
	replace d`i'=. if dist_cut_`i'==.
	tab d`i', miss
	
	rename dist_cut_`i' dist`i'
	gen d`i'_dist`i'=dist`i'*d`i'

}

* Export Figure

capture program drop figure
program define figure 
		set more off
		preserve 

		local mylabel : variable label `1'		
		
		*grouping the data into cells on each side of the cutoff 

		gen bin = floor(dist`2'/10)*10 + 10/2 
  
		sort bin
		egen bin_mean = mean(`1'), by(bin) 
		egen stderror = semean(`1'), by(bin) 


		*compute the 95% CI (binomial) 
		egen max_bin=max(bin)
		egen min_bin=min(bin)
		gen ub=.
		gen lb=.
		gen nb=.
		local nbmax=max_bin
		local nbmin=min_bin
		forvalues x=`nbmin'(10)`nbmax' {
			quietly {
			ci proportions `1' if bin==`x'
			replace ub=r(ub) if bin==`x' 
			replace lb=r(lb) if bin==`x' 
			replace nb=r(N)  if bin==`x' 
			}
			}

		
		* 50/50 BW
		gen BW50`2'=(bin>-50 & bin<50 & sch_rank==`2')
		
		* Drop bins with less than 10 observations: 
		replace bin_mean=. if nb<=10
		replace ub=. if nb<=10
		replace lb=. if nb<=10
		
		*generate Figure 
				 		
		graph twoway (scatter bin_mean bin if BW50`2' == 1, sort msymbol(o) mcolor(black)) /// 
		(line ub bin if BW50`2' == 1,sort lcolor(gs10) lpattern(dash))  /// 
		(line lb bin if  BW50`2' == 1,sort lcolor(gs10) lpattern(dash)),  ///  
		xline(0) legend(off) xtitle("Distance between student priority index and school cutoff" "(original scale in points)", size(small) ) ytitle("Fraction of students ranking school `4'", size(small) ) ///
	    graphregion(color(white))  bgcolor(white)	title(School with `3', size(small)) /// /*nodraw*/ 
		name(F_SCH`2'_NB50, replace) 
	
		restore
end

figure sch_rk1_in_rol 1 "the highest cutoff (School 11)" 11

figure sch_rk2_in_rol 2 "the 2nd highest cutoff (School 9)" 9

figure sch_rk3_in_rol 3 "the 3rd highest cutoff (School 10)" 10

figure sch_rk4_in_rol 4 "the 4th highest cutoff (School 7)" 7

graph combine F_SCH1_NB50 F_SCH2_NB50 F_SCH3_NB50 F_SCH4_NB50 ,col(2) row(4) xcommon ycommon  graphregion(color(white)) 

graph export "${folder}\figures\figure_3.pdf", replace

graph close _all
