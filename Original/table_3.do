* TABLE 3 - High School Applicants in the Southern District of Paris: Summary Statistics

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute table_3.do. This program requires proprietary data - check ReadMe file}"
	exit
	}

set more off

use "${folder}\data\district_sud.dta", clear

/* Definition of caracteristics of students'choice and assignment */

*Number of district choices
gen district_ch=(choice_rk>0)
egen nb_district_ch=sum(district_ch), by(stu_id)

* Assigned to first-choice school
gen first_a=(choice_rk==1 & sch_assignment==1)
egen first_ass=sum(first_a), by(stu_id)

/* Info on assigned school*/ 
* Distance to school
gen dist_a=distance if sch_assignment==1
egen distance_ass=mode(dist_a), by(stu_id)

* Mean DNB Math score in high school (2012)
gen sch_dnbm_12_a=sch_dnbm if sch_assignment==1
egen sch_dnbm_ass=mode(sch_dnbm_12_a), by(stu_id)

* Mean DNB French score in high school (2012)
gen sch_dnbf_12_a=sch_dnbf if sch_assignment==1
egen sch_dnbf_ass=mode(sch_dnbf_12_a), by(stu_id)

* Mean composite score in high school (2012)
gen sch_dnb_12_a=sch_dnb if sch_assignment==1
egen sch_dnb_ass=mode(sch_dnb_12_a), by(stu_id)

* Mean proportion of high SES students (2012)

gen sch_high_ses_a=sch_highses if sch_assignment==1
egen sch_highses_ass=mode(sch_high_ses_a), by(stu_id) 


/** we keep one obs per student (first choice) */
bysort stu_id : keep if choice_rk==1 
save "${folder}\tables\stats_sud.dta",replace

***** Stats displayed in stata results log *****

tabstat  /*carc student*/ stu_age stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc nb_district_ch stu_assigned  first_ass /*  
*/   /*carc 1st choice*/ distance sch_dnbf sch_dnbm  sch_dnb sch_highses  /*
 */ /*carc adm*/ distance_ass sch_dnbf_ass sch_dnbm_ass  sch_dnb_ass sch_highses_ass, /*
 */   s(mean sd min max count) f(%5.2f)  columns(statistics)  varwidth(20)
 
**** OPTIONAL Export file *****
foreach x in  mean sd min max count {
use "${folder}\tables\stats_sud.dta", clear
collapse (`x') /*carc student*/ stu_age stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc  nb_district_ch stu_assigned first_ass /*  
*/   /*carc 1st choice*/ distance sch_dnbf sch_dnbm  sch_dnb sch_highses  /*
 */ /*carc adm*/ distance_ass sch_dnbf_ass sch_dnbm_ass  sch_dnb_ass sch_highses_ass
replace stu_age=round(stu_age,.1)
replace nb_district_ch=round(nb_district_ch,.1)
foreach j in stu_female stu_dnbf_pct stu_dnbm_pct  stu_dnb_pct stu_highses stu_lowinc stu_assigned first_ass distance sch_dnbf sch_dnbm  sch_dnb sch_highses distance_ass sch_dnbf_ass sch_dnbm_ass  sch_dnb_ass sch_highses_ass  {
replace `j'=round(`j',0.01)
}
gen stat="`x'"
order stat, first
save "${folder}\tables\stats_stu_SOUTH_`x'.dta",replace
}

use "${folder}\tables\stats_stu_SOUTH_mean.dta", clear
append using "${folder}\tables\stats_stu_SOUTH_sd.dta"
append using "${folder}\tables\stats_stu_SOUTH_min.dta"
append using "${folder}\tables\stats_stu_SOUTH_max.dta"
append using "${folder}\tables\stats_stu_SOUTH_count.dta"

replace stat="N" if stat=="count"

* label student's carac*
label var stu_age "Age"
label var stu_female "Female"
label var stu_dnbf_pct "French Score"
label var stu_dnbm_pct  "Math Score"
label var stu_dnb_pct "high SES"
label var stu_lowinc "With low-income bonus"

*label choice carac*
label var nb_district_ch "Number of choices within district"
label var stu_assigned "Assigned to a within-district school"
label var  first_ass "Assigned to a first choice school"

*Attributes of first choice *

label var distance "Distance (km)"
label var sch_dnbf "Mean student French score"
label var sch_dnbm "Mean student math score"
label var sch_dnb "Mean student composite score"
label var sch_highses "Fraction high SES in school"

*Attributes of assigned school *
label var distance_ass "Distance (km) for assigned school"
label var sch_dnbf_ass "Mean student French score (assigned school)"  
label var sch_dnbm_ass "Mean student math score (assigned school)"
label var sch_dnb_ass "Mean student composite score (assigned school)"
label var sch_highses_ass "Fraction high SES in school (assigned school)"

*save "${folder}\tables\table_3.dta", replace
foreach x in  mean sd min max count {
erase "${folder}\tables\stats_stu_SOUTH_`x'.dta"
}

erase "${folder}\tables\stats_sud.dta"

outsheet using "${folder}\tables\table_3.xls", replace

