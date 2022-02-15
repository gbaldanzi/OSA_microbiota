// Project Sleep apnea and gut microbiota // 

// The aim of this script is to imput AHI values for the participants that have
// information on T90 and ODI, but not on AHI

// Imputation is performed using multiple imputation and the pmm method 

// The data was initially prepared and exported from R 

// After imputation, the script will investigate the correlation between AHI and the 
// gut microbiota species in the full model 

capture log close

log using imputation_ahi_fullmodel.txt, replace
set more off
clear 

use "/home/baldanzi/Sleep_apnea/Results/sig_mgs_basicmodel.dta"

levelsof sigmgs, clean local(sigMGS)

cd "/home/baldanzi/Sleep_apnea/Results/"

cap postclose myfile 

postfile myfile str20 MGS double rho p_value N using "/home/baldanzi/Sleep_apnea/Results/cor2_ahi_imput_mgs.dta", replace

use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta", clear

foreach mgs of varlist `sigMGS' {


	use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta",  clear 

	// Variables with names too long for Stata

	rename leisurePA_regularandmodrtactvty leisuraPA_regmod

	rename leisurePA_regularexercisortrnng leisurePA_regtrning

	rename educat_lowersecondaryeducation educat_lowsec

	rename  educat_uncmpltdprmryorlwrscndry educat_unc

	rename educat_uppersecondaryeducation educat_uppersec

	local var_full_model Alkohol smokestatus Fibrer Energi_kcal leisurePA educat placebirth visit_month metformin hypermed dyslipmed ppi
	
	foreach var of varlist `var_full_model'{
	drop if `var' == . 
	}

	foreach var in age Alkohol shannon BMI Fibrer Energi_kcal {
		egen rank_`var' = rank(`var')
	}

	gen Miss_ahi = missing(ahi) // Flag variable for missing AHI


	// Imputation 


*local var_to_keep Miss_ahi t90 age Sex Alkohol smokestatus* plate* shannon BMI Fibrer Energi_kcal leisurePA* educat* placebirth* visit_month* metformin ppi hypermed dyslipmed `mgs' rank*

	mi set mlong 

	mi register imputed ahi

	*mi register regular `var_to_keep'

	mi impute pmm ahi odi age i.Sex i.smokestatus Alkohol shannon BMI t90 i.visit_month i.plate i.leisurePA i.educat i.placebirth ///
	metformin ppi hypermed dyslipmed Fibrer Energi_kcal `mgs', add(10) knn(5)

// Diagnostics for imputation 

** mi xeq 1/2: summarize ahi if Miss_ahi==0; summarize ahi if Miss_ahi==1; summarize ahi 

** qui mi xeq 1: twoway (kdensity ahi if Miss_ahi==0) || ///
**	(kdensity ahi if Miss_ahi==1) || ///
**	(kdensity ahi if Miss_ahi), legend(label(1 "observed") label(2 "imputed") label(3 "completed"))

// Partial rank correlation 


	local full_model Sex rank_* smokestatus_* plate_* leisurePA_* educat_* placebirth_* visit_month_* metformin hypermed dyslipmed ppi
	
	mi passive: egen X_rank = rank(ahi)
	mi passive: egen Y_rank = rank(`mgs')

	mi xeq: regress Y_rank `full_model' ;  predict res1, r
	mi xeq: regress X_rank `full_model' ;  predict res2, res

	mi xeq: egen res1_std = std(res1)
	mi xeq: egen res2_std = std(res2)

	mi estimate: regress res1_std res2_std

	matrix A = r(table)
	di A[4,1]
	di e(N)


post myfile ("`mgs'") (A[1,1]) (A[4,1]) (e(N)) 

}

postclose myfile


use "/home/baldanzi/Sleep_apnea/Results/cor2_ahi_imput_mgs.dta", clear 

gen exposure = "ahi" 

export delimited using "/home/baldanzi/Sleep_apnea/Results/cor2_ahi_imput_mgs.tsv", delim(tab) replace datafmt



capture log close 

exit 
