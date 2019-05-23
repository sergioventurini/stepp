*!stepp version 0.1.0
*!Written 22May2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program stepp, byable(onecall)
	version 15.1
	syntax [anything] [if] [in] [, * ]
	
	if replay() {
		if ("`e(cmd)'" != "stepp") {
			error 301
		}
		if (_by()) {
			error 190
		}
		
		local loptions = length(`"`options'"')
		local type = e(type)
		if (`"`options'"' == substr("summary", 1, max(1, `loptions'))) {
			mata: __steppes__.summary(strupper("`type'"))
		}
		else if (`"`options'"' == "") {
			mata: __steppes__.print(strupper("`type'"), 1, 1, 1)
		}
		else {
			display as error "unrecognized option `options'"
		}
		exit
	}
	
	if _by() {
		local BY `"by `_byvars' `_byrc0':"'
	}
	if (_caller() < 8) {
		local version : display "version " string(_caller()) ", missing :"
	}
	else {
		local version : display "version " string(_caller()) " :"
	}
	
	`version' `BY' Estimate `0'  // version is not necessary
end

program Estimate, eclass byable(recall)
	version 15.1
	syntax varlist(min=2 max=2 numeric) [if] [in], ///
		TYpe(string) PATspop(numlist integer >0 max=1) ///
		MINPatspop(numlist integer >0 max=1) TRts(numlist integer ascending min=2) ///
		COVSubpop(varname numeric) [ WINType(string) WINBasedon(string) ///
		Failure(varname numeric) COMPrisk(varname numeric) COVariates(varlist) ///
		TImepoint(numlist >0 max=1) FAmily(string) Link(string) noCONStant ///
		noTEst NPerm(numlist integer >0 max=1) Seed(numlist integer >0 max=1) ///
		Eps(real 0.00001) noCLeanup ]
	
	/* Options:
	   --------
		 varlist							--> variables list to analyze providing the response
															and the treatment indicator
		 type(string)					--> model type (either 'km', 'ci' or 'glm')
		 patspop(numlist integer >0 max=1)
													--> number of patients in each subpopulation (r2)
		 minpatspop(numlist integer >0 max=1)
													--> largest number of patients in common among
															consecutive subpopulations (r1)
		 trts(numlist integer min=2)
													--> treatments list
		 covsubpop(varname numeric)
													--> covariate to use for generating the subpopulations
		 wintype(string)			--> window type (either 'sliding' or 'tail-oriented')
		 winbasedon(string)		--> what window is based on (either 'all' or 'event')
		 failure(varname numeric)
													--> censoring variable (to use with type = "km")
		 comprisk(varname numeric)
													--> competing risk variable (to use with type = "ci")
		 covariates(varlist)
													--> list of additional covariates (to use with
															type = "glm")
		 timepoint(numlist >0 max=1)
													--> timepoint at which to make the analysis
		 family(string)				--> glm family (either 'gaussian', 'binomial' or
															'poisson')
		 link(string)					--> glm link (either 'identity', 'logit', 'probit' or
															'log')
		 noconstant						--> do not include the constant (to use with
															type = "glm")
		 notest								--> do not perform the permutation test
		 nperm(numlist integer >0 max=1)
													--> number of replications in the permutaton test
		 seed(numlist integer >0 max=1)
													--> random seed
     eps(real 0.00001)    --> small value to add to zero times when (type == "km")
		 nocleanup						--> Mata temporary objects are not removed
															(undocumented)
	 */

	local cmdline : list clean 0
	local varlist : list clean varlist
	
	/* Mark sample to use */
	marksample __touse__
	/* End of marking sample to use */
	
	/* Check model type */
	local type = strlower("`type'")
	if !("`type'" == "km" | "`type'" == "ci" | "`type'" == "glm") {
		display as error "type can either be 'km', 'ci' or 'glm'"
		exit
	}
	if ("`type'" == "glm") {
		if ("`family'" == "") {
			display as error "the family argument must be provided"
			exit
		}
		if !("`family'" == "gaussian" | "`family'" == "binomial" | ///
			"`family'" == "poisson") {
			display as error "family can either be 'gaussian', 'binomial' or 'poisson'"
			exit
		}
		if ("`link'" == "") {
			display as error "the link argument must be provided"
			exit
		}
		if !("`link'" == "identity" | "`link'" == "logit" | ///
			"`link'" == "probit" | "`link'" == "log") {
			display as error "the link argument can either be 'identity', 'logit', 'probit' or 'log'"
			exit
		}
		if ("`family'" == "gaussian") {
			if ("`link'" != "identity") {
				display as error "gaussian models allow only to specify an identity link"
				exit
			}
		}
		if ("`family'" == "binomial") {
			if !("`link'" == "logit" | "`link'" == "probit") {
				display as error "binomial models allow only to specify either a logit or a probit link"
				exit
			}
		}
		if ("`family'" == "poisson") {
			if ("`link'" != "log") {
				display as error "poisson models only allow to specify a log link"
				exit
			}
		}
	}
	/* End of checking model type */
	
	/* Parse window attributes */
	if ("`wintype'" == "") local wintype = "sliding"
	if ("`winbasedon'" == "") local winbasedon = "all"
	if !("`wintype'" == "sliding" | "`wintype'" == "tail-oriented" ) {
		display as error "window type can only be set to 'sliding' or 'tail-oriented'"
		exit
	}
	if !("`winbasedon'" == "all" | "`winbasedon'" == "event" ) {
		display as error "window based-on attribute can only be set to 'all' or 'event'"
		exit
	}
	/* End parsing window attributes */
	
	/* Parse list of variables provided */
	tokenize `"`varlist'"'
	local response `1'
	local trt `2'
	tempname ntrts_data
	mata: `ntrts_data' = length(uniqrows(st_data(., "`trt'", "`__touse__'")))
	mata: st_local("ntrts_data", strofreal(`ntrts_data'))
	if (`ntrts_data' > 10) {
		display as error "number of treatments too large; make sure you " _continue
		display as error "didn't swap the response variable and the " _continue
		display as error "treatment indicator"
		exit
	}
	/* End of parsing list of variables provided */

	/* Parse trts numlist */
	tempname is_trts_uncorrect
	mata: `is_trts_uncorrect' = ///
		any(strtoreal(tokens(st_local("trts")))' :!= uniqrows(st_data(., "`trt'", "`__touse__'")))
	mata: st_local("is_trts_uncorrect", strofreal(`is_trts_uncorrect'))
	if (`is_trts_uncorrect') {
		display as error "treatment labels provided in the trts() option " _continue
		display as error "do not correspond to those included in the " _continue
		display as error "treatment column"
		exit
	}
	/* End of parsing trts numlist */

	/* Parse noconstant option */
	if (("`constant'" == "") & ("`type'" == "glm")) {
		if ("`covariates'" != "") {
			local covariates "_cons_ `covariates'"
		}
		else {
			local covariates "_cons_"
		}
	}
	/* End of parsing noconstant option */
	
	/* Parse trts option */
	tempname trts_vec ntrts
	mata: `trts_vec' = strtoreal(tokens("`trts'"))
	mata: `ntrts' = length(`trts_vec')
	mata: st_local("ntrts", strofreal(`ntrts'))
	if (`ntrts' != `ntrts_data') {
		display as error "number of treatments provided in trts option " _continue
		display as error "(`ntrts') is different than the number of treatments included in the"
		display as error "treatment indicator provided (`ntrts_data')"
		exit
	}
	/* End of parsing trts option */
	
  /* Parse eps option */
  if ("`type'" == "km") {
    if (`eps' < 0) {
      display as error "eps option must be strictly positive"
      exit
    }
    
    // generate a new (fake) response with eps added (used when there are times equal to zero)
    local eps_added = 0
    tempvar response_eps
    quietly count if (`response' == 0) & `__touse__'
    if (r(N) > 0) {
      local eps_added = 1
      quietly generate double `response_eps' = `response' if `__touse__'
			quietly replace  `response_eps' = `response' + `eps' if ///
				(`response_eps' == 0) & `__touse__'
      local response_orig = "`response'"
      local response = "`response_eps'"
    }
  }
  /* End of parsing eps option */
  
	/* Generate model's subpopulations */
	tempname stwin subp
	mata: `stwin' = stwin_wrap("`wintype'", strtoreal("`minpatspop'"), ///
		strtoreal("`patspop'"), "`winbasedon'")
	mata: `subp' = stsubpop_wrap(&`stwin', st_data(., "`covsubpop'", "`__touse__'"))
	/* End of generating model's subpopulations */

	/* Get model's estimates */
	tempname failvec compriskvec
	if ("`type'" == "km") {
		mata: `failvec' = st_data(., "`failure'", "`__touse__'")
		mata: `compriskvec' = J(1, 0, .)
	}
	else if ("`type'" == "ci") {
		mata: `failvec' = J(1, 0, .)
		mata: `compriskvec' = st_data(., "`comprisk'", "`__touse__'")
	}
	else if ("`type'" == "glm") {
		mata: `failvec' = J(1, 0, .)
		mata: `compriskvec' = J(1, 0, .)
	}
	mata: __steppes__ = steppes_wrap(&`subp', strupper("`type'"), ///
		st_data(., "`response'", "`__touse__'"), ///
		st_data(., "`trt'", "`__touse__'"), ///
		`trts_vec', `failvec', `compriskvec', strtoreal("`timepoint'"), ///
		tokens("`covariates'"), "`family'", "`link'", "`__touse__'")
	/* End of getting model's estimates */
	
	/* Perform permutation test */
	if ("`test'" == "") {
		if ("`type'" == "km") {
			mata: steppes_test_wrap(&__steppes__, strupper("`type'"), ///
				strtoreal("`nperm'"), 1, strtoreal("`seed'"), "`response'", ///
				"`failure'", "`trt'", 1, "`covariates'", "`__touse__'")
		}
		else if (("`type'" == "ci") | ("`type'" == "glm")) {
			mata: steppes_test_wrap(&__steppes__, strupper("`type'"), ///
				strtoreal("`nperm'"), 1, strtoreal("`seed'"), "", "", "", ., "", ///
				"`__touse__'")
		}
	}
	/* End of performing permutation test */
	
	/* Display results */
	mata: __steppes__.print(strupper("`type'"), 1, 1, 1)
	if ("`type'" == "km") {
		if (`eps_added') {
			display as text "Note: the 'eps' option value has been added to the `response_orig' variable"
			display as text "      to avoid the exclusion of times equal to zero"
		}
	}
	/* End of displaying results */
	
	/* Return values */
	local props "`test'"
	quietly count if `__touse__'
	local nobs = r(N)
	ereturn post, obs(`nobs') esample(`__touse__') properties(`props')
	if ("`nperm'" != "") {
		ereturn scalar nperm = `nperm'
	}
	if ("`type'" != "glm") {
		ereturn scalar timepoint = `timepoint'
	}
	ereturn scalar r1 = `minpatspop'
	ereturn scalar r2 = `patspop'
	mata: st_numscalar("e(nsubpop)", __steppes__.subpop->nsubpop)
	mata: st_numscalar("e(ntrts)", `ntrts')

	ereturn local title "Subpopulation treatment effect pattern plot (STEPP) analysis"
	ereturn local estat_cmd "stepp_estat"
	ereturn local cmdline "stepp `cmdline'"
	ereturn local cmd "stepp"
	if ("`type'" == "km") {
		if (`eps_added') {
			ereturn local responsevar "`response_orig'"
			ereturn scalar eps = `eps'
		}
		else {
			ereturn local responsevar "`response'"
		}
  }
  else {
    ereturn local responsevar "`response'"
  }
	ereturn local trtvar "`trt'"
	if ("`type'" == "km") {
		ereturn local censorvar "`failure'"
	}
	else if ("`type'" == "ci") {
		ereturn local compriskvar "`comprisk'"
	}
	else if ("`type'" == "glm") {
		ereturn local family "`family'"
		ereturn local link "`link'"
		ereturn local covariates "`covariates'"
		ereturn local noconstant "`noconstant'"
	}
	mata: st_global("e(wintype)", __steppes__.subpop->win->type)
	mata: st_global("e(winbasedon)", __steppes__.subpop->win->basedon)
	ereturn local covsubpop "`covsubpop'"
	ereturn local type "`type'"
	
	mata: return_results(&__steppes__, "`test'")
	mata: st_matrix("e(subpop)", __steppes__.subpop->subpop)
	mata: st_matrix("e(npatsub)", __steppes__.subpop->npatsub)
	mata: st_matrix("e(medianz)", __steppes__.subpop->medianz)
	mata: st_matrix("e(minc)", __steppes__.subpop->minc)
	mata: st_matrix("e(maxc)", __steppes__.subpop->maxc)
	tempname trts_vec_ret
	mata: st_matrix("`trts_vec_ret'", `trts_vec')
	local trts_nm ""
	local j = 1
	foreach trt in `trts' {
		local trts_nm "`trts_nm' trt`j'"
		local ++j
	}
	matrix colnames `trts_vec_ret' = `trts_nm'
	ereturn matrix trts = `trts_vec_ret'
	/* End of returning values */
	
	/* Clean up */
	if ("`cleanup'" == "") {
// 		capture mata: cleanup()   // all objectes are deleted apart from __steppes__
		mata: st_rclear()
	}
	/* End of cleaning up */
end

program disp_corr
	version 15.1
	syntax , matrix(string) [ Title(string) CUToff(real 0) ]

	/* Options:
	   --------
		 matrix(string)							--> matrix containing the numbers to display
		 title(string)							--> table's main title
		 cutoff(real 0)							--> do not show correlation smaller than cutoff
	  */
	
	local skip0 = 0

	if ("`title'" != "") {
		display
		display as text _skip(`skip0') "`title'"
	}

	local allvars : rownames `matrix'
	tokenize `allvars'
	local nvar : word count `allvars'
	local j0 = 1
	while (`j0' <= `nvar') {
		local j1 = min(`j0' + 9, `nvar')
		local j = `j0'
		local l = 9*(`j1' - `j0' + 1)
		display as text "{hline 13}{c TT}{hline `l'}"
		display as text _skip(13) "{c |}" _continue
		while (`j' <= `j1') {
			display as text %9s abbrev("``j''", 8) _continue
			local j = `j' + 1
		}
		display as text _newline "{hline 13}{c +}{hline `l'}"

		local i = `j0'
		while (`i' <= `nvar') {
			display as text %12s abbrev("``i''", 12) " {c |} " _continue
			local j = `j0'
			while (`j' <= min(`j1', `i')) {
				if (!missing(`matrix'[`i', `j']) & abs(`matrix'[`i', `j']) < `cutoff') {
					local c`j' = .
				}
				else {
					local c`j' = `matrix'[`i', `j']
				}
				local j = `j' + 1
			}
			local j = `j0'
			while (`j' <= min(`j1', `i')) {
				if (`c`j'' < .) {
					display as result " " %7.4f `c`j'' " " _continue
				}
				else {
					display as result _skip(9) _continue
				}
				local j = `j' + 1
			}
			display
			local i = `i' + 1
		}
		local j0 = `j0' + 10
		display as text "{hline 13}{c BT}{hline `l'}"
	}
end
