// Project Sleep apnea and gut microbiota // 

// The aim of this script is to imput AHI values for the participants that have
// information on T90 and ODI, but not on AHI

// Imputation is performed using multiple imputation and the pmm method 

// The data was initially prepared and exported from R 

cd "/home/baldanzi/Sleep_apnea/Results/imp3"

capture log close

log using imputation_ahi_basicmodel_3.txt , replace
set more off
clear 

cap postclose myfile 

postfile myfile str20 MGS double rho p_value N using "cor_ahi_imput_mgs_3.dta", replace

use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta", clear 

foreach mgs of varlist HG3A_0801-HG3A_1200{

	di `mgs'
	
	use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta", clear 

	// Variables with names too long for Stata

	rename leisurePA_regularandmodrtactvty leisuraPA_regmod

	rename leisurePA_regularexercisortrnng leisurePA_regtrning

	rename educat_lowersecondaryeducation educat_lowsec

	rename  educat_uncmpltdprmryorlwrscndry educat_unc

	rename educat_uppersecondaryeducation educat_uppersec


	drop if Alkohol == . | smokestatus == . // drop obs with missing information on basic model vars


	foreach var in age Alkohol shannon BMI Fibrer Energi_kcal {
		egen rank_`var' = rank(`var')
	}

	gen Miss_ahi = missing(ahi) // Flag variable for missing AHI
	
	// Imputation 

	mi set mlong 

	mi register imputed ahi

	mi impute pmm ahi odi age i.Sex i.smokestatus Alkohol shannon BMI t90 i.visit_month i.plate `mgs', add(10) knn(5)

// Diagnostics for imputation 

** mi xeq 1/2: summarize ahi if Miss_ahi==0; summarize ahi if Miss_ahi==1; summarize ahi 

** qui mi xeq 1: twoway (kdensity ahi if Miss_ahi==0) || ///
**	(kdensity ahi if Miss_ahi==1) || ///
**	(kdensity ahi if Miss_ahi), legend(label(1 "observed") label(2 "imputed") label(3 "completed"))

// Partial rank correlation 


	local basic_model Sex rank_age rank_Alkohol smokestatus_* plate_* rank_shannon 
	
	mi passive: egen X_rank = rank(ahi)
	mi passive: egen Y_rank = rank(`mgs')

	mi xeq: regress Y_rank `basic_model' ;  predict res1, r


	mi xeq: regress X_rank `basic_model' ;  predict res2, res

	mi xeq: egen res1_std = std(res1)
	mi xeq: egen res2_std = std(res2)

	mi estimate: regress res1_std res2_std

	matrix A = r(table)
	di A[4,1]
	di e(N)


post myfile ("`mgs'") (A[1,1]) (A[4,1]) (e(N)) 

}

postclose myfile


use "cor_ahi_imput_mgs_3.dta", clear 

gen exposure = "ahi" 

export delimited using "/home/baldanzi/Sleep_apnea/Results/cor_ahi_imput_mgs_3.tsv", delim(tab) replace datafmt


capture log close

exit

