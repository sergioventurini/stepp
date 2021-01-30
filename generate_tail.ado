*!generate_tail version 0.1
*!Written 30Jan2021
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program generate_tail, rclass
	version 15.1
	syntax [if] [in], ///
	  COVariate(varname numeric) Nsub(numlist integer >0 max=1) ///
		[ Dir(string) noSHOWResults noCLeanup ]
	
	/* Options:
	   --------
		 covariate(varname numeric)
													--> covariate to use for generating the subpopulations
		 dir(string)			    --> direction of tail-oriented windows; either 'LE'
															(default) or 'GE'; subpopulations with covariate
															values less than or equal/greater than or equal
															to the generated values
		 noshowresults				--> do not show the result summary
		 nocleanup						--> Mata temporary objects are not removed
															(undocumented)
	 */

	local cmdline : list clean 0
	local varlist : list clean varlist
	
	/* Mark sample to use */
	marksample __touse__
	markout `__touse__' `covariate'
	/* End of marking sample to use */
	
	/* Check direction */
	local dir = strupper("`dir'")
	if !("`dir'" == "LE" | "`dir'" == "GE") {
		display as error "dir can either be 'LE' or 'GE'"
		exit
	}
	/* End of checking model type */
	
	/* Generate values */
	tempname stwin r1_mat np_mat
	mata: `stwin' = gen_tailwin(st_data(., "`covariate'"), `nsub', "`dir'")
	mata: st_matrix("`r1_mat'", `stwin'.v)
	mata: st_matrix("`np_mat'", `stwin'.np)
	local nsubpop: rowsof `r1_mat'
	local i = 1
	quietly summarize `covariate'
	if ("`dir'" == "LE") {
		local r2_val = r(min)
	}
	else if ("`dir'" == "GE") {
		local r2_val = r(max)
	}
	while (`i' <= `nsubpop') {
		local v = strofreal(`r1_mat'[`i', 1])
		if ("`dir'" == "LE") {
			local r1_tmp `r1_tmp' `v'
			local r2_tmp `r2_tmp' `r2_val'
		}
		else if ("`dir'" == "GE") {
			local r1_tmp `r1_tmp' `r2_val'
			local r2_tmp `r2_tmp' `v'
		}
		local np = strofreal(`np_mat'[`i', 1])
		local np_tmp `np_tmp' `np'
		local ++i
	}
	/* End of generating values */
	
	/* Display results */
	if ("`showresults'" == "") {
		display as result "[Note: covariate values for " _continue
		display as result "tail-oriented windows are available in " _continue
		display as result "stored results (type 'return list')]"
	}
	/* End of displaying results */
	
	/* Return values */
	quietly count if `__touse__'
	return scalar nobs = r(N)
	return scalar nsubpop = `nsubpop' + 1
	return local minpatspop `"`r1_tmp'"'
	return local patspop `"`r2_tmp'"'
	return local npats `"`np_tmp'"'
	return local covariate "`covariate'"
	return local dir "`dir'"
	/* End of returning values */
	
	/* Clean up */
	if ("`cleanup'" == "") {
// 		capture mata: cleanup()   // all objectes are deleted apart from __steppes__
	}
	/* End of cleaning up */
end
