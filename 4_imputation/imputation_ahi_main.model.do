// Project Sleep apnea and gut microbiota // 

// The aim of this script is to imput AHI values for the participants that have
// information on T90 and ODI, but not on AHI

// Imputation is performed using multiple imputation and the pmm method 

// The data was initially prepared and exported from R 

cd "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"

capture log close

log using imputation_ahi_mainmodel.txt, replace

di c(current_date)
di c(current_time)

set more off
clear 

cap postclose myfile 

postfile myfile str20 MGS double rho p_value N using "cor_ahi_imput_mgs.dta", replace

use "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/pheno.dta", clear

foreach mgs of varlist HG3A*{

	 
	di "`mgs'"

	use "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/pheno.dta", clear 

	drop if Alkohol == . | smokestatus == . // drop obs with missing information on basic model vars


	** gen Miss_ahi = missing(ahi) // Flag variable for missing AHI
	
	// Imputation 


*local var_to_keep Miss_ahi t90 age Sex Alkohol smokestatus* plate* shannon BMI Fibrer Energi_kcal leisurePA* educat* placebirth* visit_month* metformin ppi hypermed dyslipmed `mgs' rank*

	mi set mlong 

	mi register imputed ahi

	*mi register regular `var_to_keep'

	mi impute pmm ahi odi age Sex i.smokestatus Alkohol shannon BMI t90 i.visit_month i.plate `mgs', add(10) knn(5)

// Diagnostics for imputation 

** mi xeq 1/2: summarize ahi if Miss_ahi==0; summarize ahi if Miss_ahi==1; summarize ahi 

** qui mi xeq 1: twoway (kdensity ahi if Miss_ahi==0) || ///
**	(kdensity ahi if Miss_ahi==1) || ///
**	(kdensity ahi if Miss_ahi), legend(label(1 "observed") label(2 "imputed") label(3 "completed"))

// Partial rank correlation 


	mi passive: egen X_rank_imp = rank(ahi)
	mi passive: egen Y_rank_imp = rank(`mgs')
	
	foreach var in age Alkohol shannon BMI {
		qui mi passive: egen rank_`var' = rank(`var')
	}
	
	local main_model Sex rank_age rank_Alkohol smokestatus_* plate_* rank_shannon rank_BMI

	qui mi estimate,saving(model1,replace): regress X_rank_imp `main_model'
	mi predict Xxb_imp using model1
	mi passive: gen Xres_imp=X_rank_imp-Xxb_imp
	qui mi estimate,saving(model2,replace): regress Y_rank_imp `main_model'
	mi predict Yxb_imp using model2
	mi passive: gen Yres_imp=Y_rank_imp-Yxb_imp
	
	mi passive: egen Xres_imp_std=std(Xres_imp)
	mi passive: egen Yres_imp_std=std(Yres_imp)

	mi estimate: regress Yres_imp_std Xres_imp_std,dof(3301)
	
	
	matrix A = r(table)


post myfile ("`mgs'") (A[1,1]) (A[4,1]) (e(N)) 

}

postclose myfile

di c(current_date)
di c(current_time)


use cor_ahi_imput_mgs.dta, clear 

gen exposure = "ahi" 

export delimited using "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/cor_ahi_imput_mgs.tsv", delim(tab) replace datafmt

di c(current_date)
di c(current_time)

capture log close

exit






