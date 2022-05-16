// Project Sleep apnea and gut microbiota // 

// The aim of this script is to impute AHI values for the participants that have
// information on T90 and ODI, but not on AHI

// Imputation is performed using multiple imputation and the pmm method 

// The data was initially prepared and exported from R 

cd "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"

capture log close

log using imputation_ahi_4.txt, replace

di c(current_date)
di c(current_time)

set seed 7

set more off
clear 

cap postclose myfile_4 

postfile myfile_4 str20 MGS double rho p_value N using "cor_ahi_imput_mgs_4.dta", replace

use "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/pheno_4.dta", clear

// Imputation and analysis were peformed in a loop so that we did not have to include 
// all species in the same imputation equatation. Therefore, for every species 
// we performed an imputation step followed by the partial Spearman correlation


foreach mgs of varlist HG3A*{

	
	di "`mgs'"

	use "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/pheno_4.dta", clear 
	
	** Flag Miss ahi observations 
	** gen Miss_ahi = 0 if ahi !=.
	** replace Miss_ahi =1 if ahi==.
	
	// Imputation 


	mi set mlong 

	mi register imputed ahi

	*mi register regular `var_to_keep'

	mi impute pmm ahi odi t90 age Sex smokestatus* Alkohol shannon BMI  month* plate* Fibrer Energi_kcal leisurePA* educat* placebirth* WaistHip `mgs', add(10) knn(5)

// Diagnostics for imputation 

 ** mi xeq 1/2: summarize ahi if Miss_ahi==0; summarize ahi if Miss_ahi==1; summarize ahi 

 ** qui mi xeq 1: twoway (kdensity ahi if Miss_ahi==0) || ///
**	(kdensity ahi if Miss_ahi==1) || ///
**	(kdensity ahi if Miss_ahi), legend(label(1 "observed") label(2 "imputed") label(3 "completed"))

// Partial rank correlation 


	// Create ranks for exposure and outcome
	mi passive: egen X_rank_imp = rank(ahi)
	mi passive: egen Y_rank_imp = rank(`mgs')
	
	// Create ranks for continuous covariates
	foreach var in age Alkohol BMI Fibrer Energi_kcal {
		qui mi passive: egen rank_`var' = rank(`var')
	}
	
	
	// * is used to capture all dummy variables for that covariate
	local ext_model Sex rank_* smokestatus* plate* educat* leisurePA* month* placebirth*

	qui mi estimate,saving(model1_4,replace): regress X_rank_imp `ext_model'
	mi predict Xxb_imp using model1_4
	mi passive: gen Xres_imp=X_rank_imp-Xxb_imp
	qui mi estimate,saving(model2_4,replace): regress Y_rank_imp `ext_model'
	mi predict Yxb_imp using model2_4
	mi passive: gen Yres_imp=Y_rank_imp-Yxb_imp
	
	mi passive: egen Xres_imp_std=std(Xres_imp)
	mi passive: egen Yres_imp_std=std(Yres_imp)

	mi estimate: regress Yres_imp_std Xres_imp_std,dof(3170)
	
	
	matrix A = r(table)


post myfile_4 ("`mgs'") (A[1,1]) (A[4,1]) (e(N)) 

}

postclose myfile_4

di c(current_date)
di c(current_time)


use cor_ahi_imput_mgs_4.dta, clear 

gen exposure = "ahi" 

// Export results back to R

export delimited using "cor_ahi_imput_mgs_4.tsv", delim(tab) replace datafmt

di c(current_date)
di c(current_time)

capture log close

exit






