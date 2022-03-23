// Project Sleep apnea and gut microbiota // 

// The aim of this script is to imput AHI values for the participants that have
// information on T90 and ODI, but not on AHI

// Imputation is performed using multiple imputation and the pmm method 

// The data was initially prepared and exported from R 

// After imputation, the script will investigate the correlation between AHI and the 
// gut microbiota species in the full model 

clear

capture log close

cd "/home/baldanzi/Sleep_apnea/Results/"

log using imputation_ahi_fullmodel.txt, replace

di c(current_date)
di c(current_time)

set more off
clear 

use "/home/baldanzi/Sleep_apnea/Results/sig_mgs_basicmodel.dta"    // import species that were assoacited with AHI (imputed), T90 or ODI

levelsof sigmgs, clean local(sigMGS)  // pass species names to a local named sigMGS

cap postclose myfile 

postfile myfile str20 MGS double rho p_value N using "cor2_ahi_imput_mgs.dta", replace

use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta", clear

foreach mgs of varlist `sigMGS' {

	di "`mgs'"

	use "/home/baldanzi/Datasets/sleep_SCAPIS/pheno.dta",  clear 
	
	// Drop obs with missing information on full model vars
	
	local var_full_model Alkohol smokestatus Fibrer Energi_kcal leisurePA educat placebirth visit_month metformin hypermed dyslipmed ppi
		
	foreach var of varlist `var_full_model'{
	di "`var'"
	drop if `var' == . 
	}

	** gen Miss_ahi = missing(ahi) // Flag variable for missing AHI


	// Imputation 

	mi set mlong 

	mi register imputed ahi

	mi impute pmm ahi odi age i.Sex i.smokestatus Alkohol shannon BMI t90 i.visit_month i.plate i.leisurePA i.educat i.placebirth ///
	metformin ppi hypermed dyslipmed Fibrer Energi_kcal `mgs', add(10) knn(5)


// Partial rank correlation 
	
	mi passive: egen X_rank_imp = rank(ahi)
	mi passive: egen Y_rank_imp = rank(`mgs')
	
	foreach var in age Alkohol shannon BMI Fibrer Energi_kcal {
		qui mi passive: egen rank_`var' = rank(`var') 
	}
	
	local full_model Sex rank_* smokestatus_* plate_* leisurePA_* educat_* placebirth_* visit_month_* metformin hypermed dyslipmed ppi
	
	qui mi estimate, saving(model1, replace): regress X_rank_imp `full_model'
	mi predict Xxb_imp using model1 
	mi passive: gen Xres_imp=X_rank-Xxb_imp
	
	qui mi estimate, saving(model2, replace): regress Y_rank_imp `full_model'
	mi predict Yxb_imp using model2
	mi passive: gen Yres_imp=Y_rank_imp-Yxb_imp
	
	mi passive: egen Xres_imp_std=std(Xres_imp)
	mi passive: egen Yres_imp_std=std(Yres_imp)
	
	mi estimate: regress Yres_imp_std Xres_imp_std, dof(3106)

	matrix A = r(table)
	di A[4,1]
	di e(N)


post myfile ("`mgs'") (A[1,1]) (A[4,1]) (e(N)) 

}

postclose myfile

di c(current_date)
di c(current_time)


use "/home/baldanzi/Sleep_apnea/Results/cor2_ahi_imput_mgs.dta", clear 

gen exposure = "ahi" 

export delimited using "/home/baldanzi/Sleep_apnea/Results/cor2_ahi_imput_mgs.tsv", delim(tab) replace datafmt

di c(current_date)
di c(current_time)

capture log close 

exit 
