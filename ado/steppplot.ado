*!steppplot version 0.3.0
*!Written 08Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program steppplot
	version 14.2

	if ("`e(cmd)'" != "stepp") {
		error 301
	}
	
	tempname results eresults
	_return hold `results'
	_return restore `results' , hold
	_estimates hold `eresults', copy restore
	capture nobreak noisily {
		_steppplot `0'
	}
	local rc = _rc
	_return restore `results'
	_estimates unhold `eresults'
	if (`rc') {
		exit `rc'
	}
end

program _steppplot
	syntax [ , Subpop Trteff Diff Ratio All Conf(numlist >0 <1 min=1 max=1) ///
		Pointwise noPOPsize ]
	
	gettoken subcmd rest : 0 , parse(", ")

	local 0 `", `options'"'
	
	tempvar __touse__
	quietly generate `__touse__' = e(sample)
	
	local pstyle "solid dash dash_dot shortdash shortdash_dot longdash longdash_dot dot"
	local lcol "gs0 gs7 gs13 gs3 gs11 gs5 gs9 gs15 gs1"
	local mstyle "O D T S Oh Dh Th Sh"
	local props = e(properties)
	local notest "notest"
	local isnotest : list notest in props
	local dq `"""'
	
	if ("`conf'" == "") local conf = 0.95
	
	tempname alpha zcrit
	mata: `alpha' = 1 - strtoreal(st_local("conf"))
	if ("`pointwise'" != "") {
		mata: st_numscalar("`zcrit'", invnormal(1 - `alpha'/2))
	}
	else {
		mata: st_numscalar("`zcrit'", ///
			invnormal(1 - `alpha'/(2*st_numscalar("e(nsubpop)"))))
	}
	
	if (("`all'" != "") | ("`rest'" == "")) {
		local subpop "subpop"
		local trteff "trteff"
		local diff "diff"
		local ratio "ratio"
	}
	
	local type = e(type)
	local family = e(family)
	local link = e(link)
	
	if ("`subpop'" != "") {
		tempname allvals nsubpop medians trts npatsub
		local ntrts = e(ntrts)
		local nsubpop = e(nsubpop)
		matrix `medians' = e(medianz)
		matrix `trts' = e(trts)
		matrix `npatsub' = e(npatsub)

		mata: `allvals' = J(`nsubpop', 4, .)
		mata: `allvals'[., 1] = (1::`nsubpop')
		mata: `allvals'[., 2] = st_matrix("e(medianz)")
		mata: `allvals'[., 3] = st_matrix("e(minc)")
		mata: `allvals'[., 4] = st_matrix("e(maxc)")
		
		forvalues i = 1/`nsubpop' {
			mata: st_local("mini", strofreal(`allvals'[`i', 3]))
			mata: st_local("maxi", strofreal(`allvals'[`i', 4]))
			local allname "`allname' function `i', range(`mini' `maxi') lpattern(solid) lcolor(gs0) lwidth(medthick) ||"
			local xlabi = `medians'[`i', 1]
			local xlabs "`xlabs' `xlabi'"
			local sizei = strofreal(`npatsub'[`i', 1])
			local npatsubi `"(n=`sizei')"'
			local xannot `"`xannot' `i' `xlabi' `"`npatsubi'"'"'
			local ylabs `"`ylabs' `i' `"Subpopulation `i'"'"'
		}
		local ymax = `nsubpop'*1.05
		local maxi = `mini' + 1e-10
		local allname "`allname' function `ymax', range(`mini' `maxi') ||"
		
		preserve
		
		quietly drop _all
		tempname allvals_d
		mata: st_matrix("`allvals_d'", `allvals')
		quietly svmat double `allvals_d', names("value")
		local allname "`allname' scatter value1 value2, msymbol(O) mcolor(gs0) mfcolor(gs0) mlcolor(gs0)"
		
		tempname subpop_plot
		local xtitle "Subpopulations by median covariate"
		twoway `allname' ||, ylabel(`ylabs', angle(0) labsize(small)) ///
			xlabel(`xlabs', labsize(small) angle(90)) ///
			text(`xannot', size(small) placement(north) justification(center) ///
			alignment(top) margin(small)) xtitle(`xtitle', margin(medsmall)) ///
			legend(off) scheme(sj) name(`subpop_plot', replace)
		
		restore
	}
	
	if ("`trteff'" != "") {
		tempname medians trts npatsub
		local ntrts = e(ntrts)
		local nsubpop = e(nsubpop)
		matrix `medians' = e(medianz)
		matrix `trts' = e(trts)
		matrix `npatsub' = e(npatsub)
		local st_pvalue "supremum test p-value ="
		local allnames ""
		local xlabs ""
		local allconnect ""
		local xannot ""
		local allcolor ""
		local allmarker ""
		local allwidth ""
		local allpattern ""
		local alllabels ""
		local allorder ""

		preserve

		quietly drop _all
		tempname xvalues
		mata: st_addobs(`nsubpop')
		mata: (void) st_addvar("double", `xvalues' = st_tempname())
		mata: st_store(., `xvalues', st_matrix("`medians'"))
		mata: st_local("xvalues", `xvalues')
		
		tempname skmObs
		matrix `skmObs' = J(`nsubpop', `ntrts', .)
		local trt1 = `trts'[1, 1]
		forvalues j = 1/`ntrts' {
			matrix `skmObs'[1, `j'] = e(sObs_`j')*100
			local allnames "`allnames' skmObs`j'"
			local allconnect "`allconnect' l"
			local trtj = `trts'[1, `j']
			local legend "`legend' label(`j' Treatment `trtj')"
			local lwidth "`lwidth' medthick"
			if ((!`isnotest') & (`j' > 1)) {
				local jm1 = `j' - 1
				local pvalue = e(pvalue_`jm1')
				if (`pvalue' > 0) {
					mata: st_local("pvalue", strofreal(`pvalue', "%10.2e"))
				}
				if (`j' < `ntrts') {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1') -"
				}
				else {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1')"
				}
			}
		}
		quietly svmat double `skmObs', names("skmObs")
		forvalues i = 1/`nsubpop' {
			local xlabi = `medians'[`i', 1]
			local xlabs "`xlabs' `xlabi'"
			local sizei = strofreal(`npatsub'[`i', 1])
			local npatsubi `"(n=`sizei')"'
			local xannot `"`xannot' 0 `xlabi' `"`npatsubi'"'"'
		}
		
		tempname trteff_plot
		local ytitle "Treatment effect estimates"
		local xtitle "Subpopulations by median covariate"
		if ("`popsize'" != "") local xannot ""
		twoway scatter `allnames' `xvalues', ylabel(0(20)100) xlabel(`xlabs') ///
			ytitle(`ytitle', margin(small)) xtitle(`xtitle', margin(medsmall)) ///
			connect(`allconnect') legend(`legend') lpattern(`pstyle') ///
			lwidth(`lwidth') lcolor(`lcol') mcolor(`lcol') mfcolor(`lcol') ///
			mlcolor(`lcol') scheme(sj) ///
			text(`xannot', orientation(vertical) size(small) placement(north) ///
			justification(center) alignment(top)) note("`st_pvalue'") ///
			name(`trteff_plot', replace)
		
		restore
	}
	
	if ("`diff'" != "") {
		tempname medians trts npatsub
		local ntrts = e(ntrts)
		local nsubpop = e(nsubpop)
		matrix `medians' = e(medianz)
		matrix `trts' = e(trts)
		matrix `npatsub' = e(npatsub)
		local st_pvalue "supremum test p-value ="
		local allnames ""
		local xlabs ""
		local allconnect ""
		local xannot ""
		local allcolor ""
		local allmarker ""
		local allwidth ""
		local allpattern ""
		local alllabels ""
		local allorder ""

		preserve

		quietly drop _all
		tempname xvalues
		mata: st_addobs(`nsubpop')
		mata: (void) st_addvar("double", `xvalues' = st_tempname())
		mata: st_store(., `xvalues', st_matrix("`medians'"))
		mata: st_local("xvalues", `xvalues')
		
		tempname skmObs skmObsV dskmObs dskmObsSE dskmObsSE_tmp ///
			dskmObsL dskmObsU dskmObsL_tmp dskmObsU_tmp
		matrix `skmObs' = J(`nsubpop', `ntrts', .)
		mata: `skmObsV' = J(`nsubpop', `ntrts', .)
		matrix `dskmObs' = J(`nsubpop', `ntrts' - 1, .)
		matrix `dskmObsSE' = J(`nsubpop', `ntrts' - 1, .)
		mata: `dskmObsSE_tmp' = J(`nsubpop', `ntrts' - 1, .)
		matrix `dskmObsL' = J(`nsubpop', `ntrts' - 1, .)
		matrix `dskmObsU' = J(`nsubpop', `ntrts' - 1, .)
		mata: `dskmObsL_tmp' = J(`nsubpop', `ntrts' - 1, .)
		mata: `dskmObsU_tmp' = J(`nsubpop', `ntrts' - 1, .)
		local trt1 = `trts'[1, 1]
		forvalues j = 1/`ntrts' {
			matrix `skmObs'[1, `j'] = e(sObs_`j')*100
			mata: `skmObsV'[., `j'] = (st_matrix("e(sSE_`j')")*100):^2
			local trtj = `trts'[1, `j']
			if (`j' > 1) {
				matrix `dskmObs'[1, `j' - 1] = `skmObs'[1..., 1] - `skmObs'[1..., `j']
				mata: `dskmObsSE_tmp'[., `j' - 1] = sqrt(`skmObsV'[., 1] + `skmObsV'[., `j'])
				mata: `dskmObsL_tmp'[., `j' - 1] = ///
					st_matrix("`dskmObs'")[., `j' - 1] - st_numscalar("`zcrit'")*`dskmObsSE_tmp'[., `j' - 1]
				mata: `dskmObsU_tmp'[., `j' - 1] = ///
					st_matrix("`dskmObs'")[., `j' - 1] + st_numscalar("`zcrit'")*`dskmObsSE_tmp'[., `j' - 1]
				
				local jm1 = `j' - 1
				local skmObs_nm "dskmObs`jm1' dskmObsL`jm1' dskmObsU`jm1'"
				local allnames "`allnames' `skmObs_nm'"
				local allconnect "`allconnect' l l l"
				local markertouse : word `jm1' of `mstyle'
				local allmarker "`allmarker' `markertouse' `markertouse' `markertouse'"
				local allwidth "`allwidth' medthick medthick medthick"
				local colortouse : word `jm1' of `lcol'
				local allcolor "`allcolor' `colortouse' `colortouse' `colortouse'"
				local patterntouse : word `jm1' of `pstyle'
				local allpattern "`allpattern' `patterntouse' `patterntouse' `patterntouse'"
				local j1 = 1 + 3*(`jm1' - 1)
				local j2 = 2 + 3*(`jm1' - 1)
				local j3 = 3 + 3*(`jm1' - 1)
				local alllabels `"`alllabels' label(`j1' Treatment `trtj' vs. `trt1') label(`j2' `dq'`dq') label(`j3' `dq'`dq')"'
				local allorder "`allorder' `j1'"
			}
			if ((!`isnotest') & (`j' > 1)) {
				local jm1 = `j' - 1
				local pvalue = e(pvalue_`jm1')
				if (`pvalue' > 0) {
					mata: st_local("pvalue", strofreal(`pvalue', "%10.2e"))
				}
				if (`j' < `ntrts') {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1') -"
				}
				else {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1')"
				}
			}
		}
		
		mata: st_matrix("`dskmObsSE'", `dskmObsSE_tmp')
		mata: st_matrix("`dskmObsL'", `dskmObsL_tmp')
		mata: st_matrix("`dskmObsU'", `dskmObsU_tmp')
		quietly svmat double `skmObs', names("skmObs")
		quietly svmat double `dskmObs', names("dskmObs")
		quietly svmat double `dskmObsSE', names("dskmObsSE")
		quietly svmat double `dskmObsL', names("dskmObsL")
		quietly svmat double `dskmObsU', names("dskmObsU")
		forvalues i = 1/`nsubpop' {
			local xlabi = `medians'[`i', 1]
			local xlabs "`xlabs' `xlabi'"
			local sizei = strofreal(`npatsub'[`i', 1])
			local npatsubi `"(n=`sizei')"'
			local xannot `"`xannot' -100 `xlabi' `"`npatsubi'"'"'
		}
		
		tempname diff_plot
		local ytitle "Treatment effect estimate differences"
		local xtitle "Subpopulations by median covariate"
		if ("`popsize'" != "") local xannot ""
		twoway scatter `allnames' `xvalues', ylabel(-100(20)100) xlabel(`xlabs') ///
			ytitle(`ytitle', margin(small)) xtitle(`xtitle', margin(medsmall)) ///
			connect(`allconnect') legend(`alllabels' order(`allorder')) ///
			lpattern(`allpattern') lwidth(`allwidth') lcolor(`allcolor') ///
			mcolor(`allcolor') mfcolor(`allcolor') mlcolor(`allcolor') ///
			msymbol(`allmarker') scheme(sj) ///
			text(`xannot', orientation(vertical) size(small) placement(north) ///
			justification(center) alignment(top)) note("`st_pvalue'") ///
			yline(0, lwidth(thin)) name(`diff_plot', replace)
		
		restore
	}
	
	if ("`ratio'" != "") {
		tempname medians trts npatsub
		local ntrts = e(ntrts)
		local nsubpop = e(nsubpop)
		matrix `medians' = e(medianz)
		matrix `trts' = e(trts)
		matrix `npatsub' = e(npatsub)
		local st_pvalue "supremum test p-value ="
		local allnames ""
		local xlabs ""
		local allconnect ""
		local xannot ""
		local allcolor ""
		local allmarker ""
		local allwidth ""
		local allpattern ""
		local alllabels ""
		local allorder ""

		preserve

		quietly drop _all
		tempname xvalues
		mata: st_addobs(`nsubpop')
		mata: (void) st_addvar("double", `xvalues' = st_tempname())
		mata: st_store(., `xvalues', st_matrix("`medians'"))
		mata: st_local("xvalues", `xvalues')
		
		tempname logHR logHR_tmp logHRSE logHRSE_tmp HR HRL HRU HRL_tmp HRU_tmp
		matrix `logHR' = J(`nsubpop', `ntrts' - 1, .)
		mata: `logHR_tmp' = J(`nsubpop', `ntrts' - 1, .)
		matrix `logHRSE' = J(`nsubpop', `ntrts' - 1, .)
		mata: `logHRSE_tmp' = J(`nsubpop', `ntrts' - 1, .)
		matrix `HRL' = J(`nsubpop', `ntrts' - 1, .)
		matrix `HRU' = J(`nsubpop', `ntrts' - 1, .)
		mata: `HRL_tmp' = J(`nsubpop', `ntrts' - 1, .)
		mata: `HRU_tmp' = J(`nsubpop', `ntrts' - 1, .)
		local ntrtsm1 = `ntrts' - 1
		local trt1 = `trts'[1, 1]
		forvalues j = 1/`ntrtsm1' {
			mata: `logHR_tmp'[., `j'] = st_matrix("e(logHR_`j')")
			mata: `logHRSE_tmp'[., `j'] = st_matrix("e(logHRSE_`j')")
			local trtj = `trts'[1, `j' + 1]
			mata: `HRL_tmp'[., `j'] = exp(`logHR_tmp'[., `j'] - ///
				st_numscalar("`zcrit'")*`logHRSE_tmp'[., `j'])
			mata: `HRU_tmp'[., `j'] = exp(`logHR_tmp'[., `j'] + ///
				st_numscalar("`zcrit'")*`logHRSE_tmp'[., `j'])
			
			local Ratios_nm "HR`j' HRL`j' HRU`j'"
			local allnames "`allnames' `Ratios_nm'"
			local allconnect "`allconnect' l l l"
			local markertouse : word `j' of `mstyle'
			local allmarker "`allmarker' `markertouse' `markertouse' `markertouse'"
			local allwidth "`allwidth' medthick medthick medthick"
			local colortouse : word `j' of `lcol'
			local allcolor "`allcolor' `colortouse' `colortouse' `colortouse'"
			local patterntouse : word `j' of `pstyle'
			local allpattern "`allpattern' `patterntouse' `patterntouse' `patterntouse'"
			local j1 = 1 + 3*(`j' - 1)
			local j2 = 2 + 3*(`j' - 1)
			local j3 = 3 + 3*(`j' - 1)
			local alllabels `"`alllabels' label(`j1' Treatment `trtj' vs. `trt1') label(`j2' `dq'`dq') label(`j3' `dq'`dq')"'
			local allorder "`allorder' `j1'"
			local jp1 = `j' + 1
			if (!`isnotest') {
				local pvalue = e(pvalue_`j')
				if (`pvalue' > 0) {
					mata: st_local("pvalue", strofreal(`pvalue', "%10.2e"))
				}
				if (`j' < `ntrtsm1') {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1') -"
				}
				else {
					local st_pvalue "`st_pvalue' `pvalue' (`trtj' vs. `trt1')"
				}
			}
		}
		
		mata: st_matrix("`logHR'", `logHR_tmp')
		mata: st_matrix("`logHRSE'", `logHRSE_tmp')
		mata: st_matrix("`HR'", exp(`logHR_tmp'))
		mata: st_matrix("`HRL'", `HRL_tmp')
		mata: st_matrix("`HRU'", `HRU_tmp')
		quietly svmat double `logHR', names("logHR")
		quietly svmat double `logHRSE', names("logHRSE")
		quietly svmat double `HR', names("HR")
		quietly svmat double `HRL', names("HRL")
		quietly svmat double `HRU', names("HRU")
		forvalues i = 1/`nsubpop' {
			local xlabi = `medians'[`i', 1]
			local xlabs "`xlabs' `xlabi'"
			local sizei = strofreal(`npatsub'[`i', 1])
			local npatsubi `"(n=`sizei')"'
			local xannot `"`xannot' 0 `xlabi' `"`npatsubi'"'"'
		}
		
		tempname ratio_plot maxHR
		if ("`type'" != "glm") {
			if ("`type'" == "km") {
				local ytitle "Hazard ratios"
			}
			else if ("`type'" == "ci") {
				local ytitle "Subdistribution hazard ratios"
			}
		}
		else {
			if ("`family'" == "gaussian") {
				local ytitle "Ratios"
			}
			else if ("`family'" == "binomial") {
				local ytitle "Odds ratios"
			}
			else if ("`family'" == "poisson") {
				local ytitle "Risk ratios"
			}
		}
		local xtitle "Subpopulations by median covariate"
		mata: st_local("yrange", invtokens(strofreal(pretty((0 \ vec(`HRU_tmp')), 7))'))
		if ("`popsize'" != "") local xannot ""
		twoway scatter `allnames' `xvalues', ylabel(`yrange') xlabel(`xlabs') ///
			ytitle(`ytitle', margin(small)) xtitle(`xtitle', margin(medsmall)) ///
			connect(`allconnect') legend(`alllabels' order(`allorder')) ///
			lpattern(`allpattern') lwidth(`allwidth') lcolor(`allcolor') ///
			mcolor(`allcolor') mfcolor(`allcolor') mlcolor(`allcolor') ///
			msymbol(`allmarker') scheme(sj) ///
			text(`xannot', orientation(vertical) size(small) placement(north) ///
			justification(center) alignment(top)) note("`st_pvalue'") ///
			yline(1, lwidth(thin)) name(`ratio_plot', replace)
		
		restore
	}
	
	/* Clean up */
	if ("`cleanup'" == "") {
		capture mata: cleanup()   // all objectes are deleted apart from __steppes__
		mata: st_rclear()
	}
	/* End of cleaning up */
end
