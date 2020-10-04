*!balance version 0.3.0
*!Written 05Oct2020
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program balance_patients, eclass byable(recall)
	version 15.1
	syntax [if] [in], ///
		range_r1(numlist integer >1 min=2 max=2 sort) ///
		range_r2(numlist integer >1 min=2 max=2 sort) ///
		MAXNsubpops(numlist integer >1 min=1 max=1) COVar(varname numeric) ///
		[ Plot noSHOWResults noCLeanup ]
	
	/* Options:
	   --------
		 range_r1(numlist integer >1 min=2 max=2 sort)
		                          --> range of r1 values
		 range_r2(numlist integer >1 min=2 max=2 sort)
		                          --> range of r2 values
		 maxnsubpops(numlist integer >1 min=1 max=1)
		                          --> maximum number of subpopulations to consider
		 covar(varname numeric)   --> covariate to use for generating the
                                  subpopulations
		 plot								      --> produce a plot of the results
		 noshowresults				    --> do not show the result summary
		 nocleanup						    --> Mata temporary objects are not removed
															    (undocumented)
	 */

	local cmdline : list clean 0

	/* Mark sample to use */
	marksample __touse__
	/* End of marking sample to use */
	
	/* Parse range_r1 numlist */
	tokenize "`range_r1'"
	local r1_min "`1'"
	local r1_max "`2'"
	/* End of parsing range_r1 numlist */

	/* Parse range_r2 numlist */
	tokenize "`range_r2'"
	local r2_min "`1'"
	local r2_max "`2'"
	/* End of parsing range_r2 numlist */

	/* Perform calculations */
	tempname cvrt freqdist allfreqs k i res_calc
	mata: `cvrt' = st_data(., "covar")
	mata: `freqdist' = uniqrows(`cvrt', 1)
	mata: `freqdist' = (`freqdist', runningsum(`freqdist'[., 2]))
	mata: `k' = rows(`freqdist')
	mata: `allfreqs' = J(`k', `k', .)
	mata: `allfreqs'[1, .] = `freqdist'[., 3]'
	mata: for (`i' = 2; `i' <= `k'; `i'++) ///
		  `allfreqs'[`i', `i'..`k'] = `freqdist'[`i'..`k', 3]' - ///
			  J(1, `k' - `i' + 1, `freqdist'[`i' - 1, 3])
	mata: `res_calc' = balance_p_calc(strtoreal(st_local("r1_min")), ///
	  strtoreal(st_local("r1_max")), strtoreal(st_local("r2_min")), ///
		strtoreal(st_local("r2_max")), strtoreal(st_local("maxnsubpops")), ///
    `freqdist', `allfreqs')
	/* End of performinf calculations */

	/* Display results */
	if ("`showresults'" == "") {
		mata: balance_p_print(`res_calc', strtoreal(st_local("r1_min")), ///
	  strtoreal(st_local("r1_max")), strtoreal(st_local("r2_min")), ///
		strtoreal(st_local("r2_max")), strtoreal(st_local("maxnsubpops")), ///
		  `freqdist', `allfreqs', st_data(., "`covar'", "`__touse__'"))
	}
	/* End of displaying results */
	
	/* Plot results */
	if ("`plot'" != "") {
		display
		display as error "warning: the plot option is not yet implemented."
	}
	/* End of plotting results */
	
	/* Return values */
	quietly count if `__touse__'
	local nobs = r(N)
	ereturn post, obs(`nobs') esample(`__touse__')

	ereturn local covariate "`covar'"	
	ereturn local maxnsubpops "`maxnsubpops'"
	ereturn local r2_range "`range_r2'"
	ereturn local r1_range "`range_r1'"
	ereturn local cmd "balance_patients"
	ereturn local cmdline "balance_patients `cmdline'"
	ereturn local title "Balanced subpopulations determination for STEPP (patients)"
	
	mata: st_numscalar("e(r1_best)", `res_calc'.r1best)
	mata: st_numscalar("e(r2_best)", `res_calc'.r2best)
	mata: st_numscalar("e(var_best)", `res_calc'.varbest)
	mata: st_numscalar("e(nsubpops_best)", `res_calc'.indbest)
	mata: st_matrix("e(all_res)", `res_calc'.resmat)
	/* End of returning values */

	/* Clean up */
	if ("`cleanup'" == "") {
// 		capture mata: cleanup()   // all objectes are deleted apart from __steppes__
//		mata: st_rclear()
	}
	/* End of cleaning up */

end
