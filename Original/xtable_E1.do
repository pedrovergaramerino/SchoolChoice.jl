* TABLE E1 -  Assigned and Unassigned Students in the Southern District of Paris

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute table_E1.do. This program requires proprietary data - check ReadMe file}"
	exit
	}

set more off

use "${folder}\data\district_sud.dta", clear

bysort stu_id : keep if _n==1
keep stu_age stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc stu_enrolment stu_assigned stu_id
gen nb_inparis=1 if stu_enrolment~=4
gen in_assigned=(stu_enrolment==1)
gen in_otherpub=(stu_enrolment==2)
gen in_private=(stu_enrolment==3)

*stats 
tabstat  /*carc student*/ stu_age stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc , by(stu_assigned)   s(mean count) f(%5.2f)    varwidth(20)
*** Stats on students's assignment (we only take into  students who remained in the Parisian LEA) 
bysort stu_assigned : ta stu_enrolment if stu_enrolment~=4, freq

*table 
gen nb_tot=1
collapse  stu_age stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc (sum) nb_inparis in_assigned in_otherpub in_private nb_tot , by(stu_assigned)
replace in_assigned=in_assigned/nb_inparis
replace in_otherpub=in_otherpub/nb_inparis
replace in_private=in_private/nb_inparis
drop nb_inparis
*save "${folder}\tables\xtable_E1.dta", replace
outsheet using "${folder}\tables\xtable_E1.xls", replace
