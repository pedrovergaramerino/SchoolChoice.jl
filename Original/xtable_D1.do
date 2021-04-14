* TABLE D1 - Fraction of Students Ranking Each of the Four Most Selective Schools in the Southern District of Paris, by Distance to School Cutoff

* N.B.: requires package outreg2
* To install, type "ssc install outreg2" in command window

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute table_D1.do. This program requires proprietary data - check ReadMe file}"
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

* Generate variables for distance to cutoff

local nb=max_rank 
forvalues i=1(1) `nb'{

	* Rescale distance to the score so that 0 corresponds to the school cutoff (in points, not proportion) 
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
	
	gen dist_`i'_above=dist`i'/100*d`i'
	gen dist_`i'_below=dist`i'/100*(1-d`i')
	
	*Generate a variable for school_priority pts and interaction above and below cutoff
	gen priority_`i'=stu_priority_pts if sch_rank==`i'
	label var priority_`i' "student priority points for school `i'"
	gen priority_`i'_below=priority_`i'*(1-d`i')
	gen priority_`i'_above=priority_`i'*d`i'
	
	rename d`i' above_cutoff_`i'

}

rename stu_dnb_pct stu_dnb
rename stu_dnbm_pct stu_dnbm
rename stu_dnbf_pct stu_dnbf

replace stu_dnbm=stu_dnbm
replace stu_dnbf=stu_dnbf
replace stu_dnb=stu_dnb

gen stu_dnb2=stu_dnb^2
gen stu_dnbf2=stu_dnbf^2
gen stu_dnbm2=stu_dnbm^2

* Program for regressions
forvalues i=1(1) 1{


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)')  cttop(without controls) nocons  excel replace


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' stu_dnbf  stu_dnbm if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)') cttop(dnb controls) nocons  excel


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' stu_dnbf stu_dnbm stu_dnbf2  stu_dnbm2  distance closest_sch collocated_sch stu_highses if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)')  cttop(all controls) nocons  excel

}


forvalues i=2(1) 2{


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)')  cttop(without controls)  nocons excel 


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' stu_dnbf  stu_dnbm if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)') cttop(dnb controls) nocons   excel


regress sch_rk`i'_in_rol   dist_`i'_below dist_`i'_above  above_cutoff_`i' stu_dnbf stu_dnbm stu_dnbf2  stu_dnbm2  distance closest_sch collocated_sch stu_highses if stu_lowinc==0 , robust
test  dist_`i'_below dist_`i'_above  above_cutoff_`i'
outreg2 using "${folder}\tables\xtable_D1.out", dec(3) adds(R2_adjusted, `e(r2_a)', F-test, r(F), Prob > F, `r(p)')  cttop(all controls) nocons  excel

}

erase "${folder}\tables\xtable_D1.out"

* Note : output table D1.xml to be opened in Excel

