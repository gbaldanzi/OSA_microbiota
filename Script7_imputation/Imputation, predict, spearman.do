cls
clear
set obs 250
gen X=rnormal()
gen Z=rnormal()
gen Y=rnormal()+0.8*X+0.5*Z
gen U=runiform()
replace Y=. if runiform()<0.2

egen X_rank=rank(X)
egen Y_rank=rank(Y)
egen Z_rank=rank(Z)

noi pcorr Y_rank X_rank Z_rank

mi set wide
mi register imputed Y
mi impute chained (regress) Y=X Z,add(10) //unnecessarily many imputations in practice, but good for simulation

mi passive: egen X_rank_imp=rank(X) //mi xeq can be used, but only works as it should if the variable is already created and registred as passive...a bit tricky for the egen-variables.
mi passive: egen Y_rank_imp=rank(Y)
mi passive: egen Z_rank_imp=rank(Z)

mi estimate,saving(model1,replace): regress X_rank_imp Z_rank_imp
mi predict Xxb_imp using model1
mi passive: gen Xres_imp=X_rank_imp-Xxb_imp
mi estimate,saving(model2,replace): regress Y_rank_imp Z_rank_imp
mi predict Yxb_imp using model2

mi passive: gen Yres_imp=Y_rank_imp-Yxb_imp
mi passive: egen Xres_imp_std=std(Xres_imp)
mi passive: egen Yres_imp_std=std(Yres_imp)

mi estimate: regress Yres_imp Xres_imp
noi mi estimate: regress Yres_imp_std Xres_imp_std,dof(247)  
//calculates partial spearman with MI. Note the DOF option! It's because when we run regress on residuals, the observations are no longer independent.
//(among other things, the residuals have to sum to 0). The DOF will be the number of observation minus 2 (slope, intercept) minus the variables you've "partialled out" when creating the residuals - in this case only Z_rank_imp
//That is 250-2-1=247 (instead of the automatic 248). In your case, it should make a somewhat bigger difference since there's a lot of adjustment variables (batches, for example).

/*
*findit mibeta
*noi mibeta Yres_imp Xres_imp  //slightly different since Rubin's rules are performed on transformed correlation coefficients, which is probably preferable.
//Unfortunately, this function doesn't produce confidence intervals - and I think combining imputation and bootstrap is really too much.
//But you can loook at the mibeta function if you have the time (it's not part of official Stata, but created by the StataCorp team), and make sure that the estimates
//don't deviate too much from the mi estimate: regress-estimate of the correlation.
*/

/*
//completed data (the Stata suggestion: same same as just running partial spearman in the first place [there's some minor differences compared to pcorr due to the missingness)
mi xeq: regress X_rank_imp Z_rank_imp; predict resX,r
mi xeq: regress Y_rank_imp Z_rank_imp; predict resY,r 
mi passive: egen resXstd=std(resX)
mi passive: egen resYstd=std(resY)

noi mi estimate: regress resYstd resXstd
*noi pcorr Y_rank X_rank Z_rank

regress Y_rank Z_rank
predict resY_compl,r
regress X_rank Z_rank
predict resX_compl,r 
egen resXstd_compl=std(resX_compl)
egen resYstd_compl=std(resY_compl)

noi regress resYstd_compl resXstd_compl
*/
