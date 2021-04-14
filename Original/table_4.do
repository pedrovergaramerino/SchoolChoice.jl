* TABLE 4 -  High Schools in the Southern District of Paris: Summary Statistics

clear

if "${folder}" == "" {
	display "{red: ERROR: Specify path to folder containing replication files in change_directory.do and execute script}"
}
else{
	cd "$folder"
}
capture confirm file "${folder}/data/nul"
if _rc { 
	di "{red: ERROR: Cannot execute table_4.do. This program requires proprietary data - check ReadMe file}"
	exit
	}


set more off

use "${folder}\data\district_sud.dta", clear

* count number of assigned students
egen sch_count=sum(sch_assignment), by(sch_id)

* calculate share of students ranking the school
gen inrol=(choice_rk~=0)
egen nb_inrol=sum(inrol), by(sch_id)
gen count_stu=1
egen nb_tot=sum(count_stu), by(sch_id)
gen sch_inrol=nb_inrol/nb_tot
bysort sch_id : keep if _n==1

* create stats table
keep sch_id sch_dnbf sch_dnbm sch_dnb sch_highses sch_capacity sch_count sch_cutoff sch_inrol
foreach j in sch_id sch_dnbf sch_dnbm sch_dnb sch_highses sch_inrol {
replace `j'=round(`j',0.01)
}
replace sch_cutoff=round(sch_cutoff,0.001)
order sch_id sch_dnbf sch_dnbm sch_dnb sch_highses sch_capacity sch_count  sch_cutoff  sch_inrol, first

** Add admission cutoffs, fraction of rol ranking it, count 

* save "${folder}\tables\table_4.dta", replace
outsheet using "${folder}\tables\table_4.xls", replace

