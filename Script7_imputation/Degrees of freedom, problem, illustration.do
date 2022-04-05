//extreme example of the degrees of freedom-problem
cls
clear
set obs 250
gen X=rnormal()
forv i=1/240{
	gen Z`i'=rnormal()
	egen Z`i'_rank=rank(Z`i')
}
gen Y=rnormal() //no real association with anything

egen X_rank=rank(X)
egen Y_rank=rank(Y)

noi pcorr Y_rank X_rank Z*_rank

regress Y_rank Z*_rank //we partial away 240 variables
predict resY_compl,r
regress X_rank Z*_rank
predict resX_compl,r 
egen resXstd_compl=std(resX_compl)
egen resYstd_compl=std(resY_compl)

noi regress resYstd_compl resXstd_compl //same point estimate as pcorr, but much too low p-value
noi regress resYstd_compl resXstd_compl,dof(8) //...much better
