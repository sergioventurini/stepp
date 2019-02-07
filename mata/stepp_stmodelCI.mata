*!stepp_stmodelCI version 0.1.0
*!Written 08Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void stmodelCI_class::new()
{
	class_name = "stmodelCI"
}

string scalar stmodelCI_class::getclass()
{
	return(class_name)
}

struct steppes_effect scalar stmodelCI_class::estimate(
	pointer(class stsubpop_class scalar) scalar sp, string scalar touse)
{
	struct steppes_effect scalar estimate
	real scalar nsubpop, timePoint, ntrts, i, j, index, temp_oObs_1, temp_oObs_j
	real vector survTime, coltrt, trts, type, txassign, sel, sel_ind, temp_Obs,
		temp_HRs, temp_Obs_1, temp_Obs_j, temp_SE_1, temp_SE_j
	real matrix subpop
	class AssociativeArray scalar Obs, TrtEff, HRs, Ratios, result, temp_TrtEff,
		temp_Ratios, temp_result
	
	nsubpop = sp->nsubpop
	subpop = sp->subpop
	coltrt = this.coltrt
	survTime = this.coltime
	trts = this.trts
	timePoint = this.timePoint
	type = this.coltype

	ntrts = length(trts)

	// set up the return structure
	Obs.reinit("string", 1)
	Obs.put("sObs", J(nsubpop, 1, 0))
	Obs.put("sSE", J(nsubpop, 1, 0))
	Obs.put("oObs", 0)
	Obs.put("oSE", 0)
	TrtEff.reinit("real", 1)
	for (i = 1; i <= ntrts; i++) {
		TrtEff.put(i, Obs)
	}

	HRs.reinit("string", 1)
	HRs.put("skmw", 0)
	HRs.put("logHR", J(nsubpop, 1, 0))
	HRs.put("logHRSE", J(nsubpop, 1, 0))
	HRs.put("ologHR", 0)
	HRs.put("ologHRSE", 0)
	HRs.put("logHRw", 0)
	Ratios.reinit("real", 1)
	for (i = 1; i <= (ntrts - 1); i++) {
		Ratios.put(i, HRs)
	}

	// 	estimate the abs and rel trt effects comparing trt 1 with trt j
	for (j = 2; j <= ntrts; j++) {
		txassign = J(length(coltrt), 1, -1)
		txassign[selectindex(coltrt :== trts[1])] =
			J(length(selectindex(coltrt :== trts[1])), 1, 1)
		txassign[selectindex(coltrt :== trts[j])] =
			J(length(selectindex(coltrt :== trts[j])), 1, 0)

		sel = (txassign :== 1 :| txassign :== 0)

		// for each subgroup i
		for (i = 1; i <= nsubpop; i++) {
			sel_ind = selectindex((subpop[., i] :== 1) :& (txassign :!= -1))
			result = cuminc_HR(survTime[sel_ind], type[sel_ind], txassign[sel_ind])
			
			temp_result = result.get("1 1")
			if (max(temp_result.get("time")) >= timePoint) {
				index = sum(temp_result.get("time") :<= timePoint)
				temp_TrtEff = TrtEff.get(1)
				temp_Obs = temp_TrtEff.get("sObs")
				temp_Obs[i] = temp_result.get("est")[index]
				temp_TrtEff.put("sObs", temp_Obs)
				temp_Obs = temp_TrtEff.get("sSE")
				temp_Obs[i] = sqrt(temp_result.get("var")[index])
				temp_TrtEff.put("sSE", temp_Obs)
				TrtEff.put(1, temp_TrtEff)
			}
			temp_result = result.get("0 1")
			if (max(temp_result.get("time")) >= timePoint) {
				index = sum(temp_result.get("time") :<= timePoint)
				temp_TrtEff = TrtEff.get(j)
				temp_Obs = temp_TrtEff.get("sObs")
				temp_Obs[i] = temp_result.get("est")[index]
				temp_TrtEff.put("sObs", temp_Obs)
				temp_Obs = temp_TrtEff.get("sSE")
				temp_Obs[i] = sqrt(temp_result.get("var")[index])
				temp_TrtEff.put("sSE", temp_Obs)
				TrtEff.put(j, temp_TrtEff)
			}
			
			temp_Ratios = Ratios.get(j - 1)
			temp_HRs = temp_Ratios.get("logHR")
			temp_HRs[i] = result.get("ome")/result.get("omevar")
			temp_Ratios.put("logHR", temp_HRs)
			temp_HRs = temp_Ratios.get("logHRSE")
			temp_HRs[i] = sqrt(1/result.get("omevar"))
			temp_Ratios.put("logHRSE", temp_HRs)
			Ratios.put(j - 1, temp_Ratios)
		}
		temp_Ratios = Ratios.get(j - 1)
		temp_Ratios.put("logHRw",
			sum(temp_Ratios.get("logHR") :/ temp_Ratios.get("logHRSE")))
		temp_TrtEff = TrtEff.get(1)
		temp_Obs_1 = temp_TrtEff.get("sObs")
		temp_SE_1 = temp_TrtEff.get("sSE")
		temp_TrtEff = TrtEff.get(j)
		temp_Obs_j = temp_TrtEff.get("sObs")
		temp_SE_j = temp_TrtEff.get("sSE")
		temp_Ratios.put("skmw",
			sum((temp_Obs_1 - temp_Obs_j) :/ sqrt(temp_SE_1:^2 + temp_SE_j:^2)))
		Ratios.put(j - 1, temp_Ratios)

		// estimate the overall effect
		sel_ind = selectindex(sel)
		result = cuminc_HR(survTime[sel_ind], type[sel_ind], txassign[sel_ind])
		
			temp_result = result.get("1 1")
			if (max(temp_result.get("time")) >= timePoint) {
				index = length(temp_result.get("time")[selectindex(temp_result.get("time") :<= timePoint)])
				temp_TrtEff = TrtEff.get(1)
				temp_TrtEff.put("oObs", temp_result.get("est")[index])
				temp_TrtEff.put("oSE", sqrt(temp_result.get("var")[index]))
				TrtEff.put(1, temp_TrtEff)
		}
			temp_result = result.get("0 1")
			if (max(temp_result.get("time")) >= timePoint) {
				index = sum(temp_result.get("time") :<= timePoint)
				temp_TrtEff = TrtEff.get(j)
				temp_TrtEff.put("oObs", temp_result.get("est")[index])
				temp_TrtEff.put("oSE", sqrt(temp_result.get("var")[index]))
				TrtEff.put(j, temp_TrtEff)
		}
		temp_Ratios = Ratios.get(j - 1)
		temp_Ratios.put("ologHR", result.get("ome")/result.get("omevar"))
		temp_Ratios.put("ologHRSE", sqrt(1/result.get("omevar")))
		Ratios.put(j - 1, temp_Ratios)

		temp_TrtEff = TrtEff.get(1)
		temp_Obs_1 = temp_TrtEff.get("sObs")
		temp_oObs_1 = temp_TrtEff.get("oObs")
		temp_TrtEff = TrtEff.get(j)
		temp_Obs_j = temp_TrtEff.get("sObs")
		temp_oObs_j = temp_TrtEff.get("oObs")
		if (hasmissing(temp_Obs_1) | hasmissing(temp_Obs_j) |
			hasmissing(temp_oObs_1) | hasmissing(temp_oObs_j)) {
			printf("\n")
      printf("{err}Unable to estimate survival time at " +
        strofreal(timePoint) + " time-unit(s) for trt " + strofreal(j) + "\n")
      printf("{err}because there are too few events within one or more subpopulation(s).\n")
      printf("{err}The problem may be avoided by constructing larger ")
      printf("{err}subpopulations \n")
      printf("{err}and/or by selecting a different timepoint for estimation.\n")
			_error(3000)
		}
	}

	// estimate (cumulative incidence) object to be returned
	estimate.model = "CIe"
	estimate.ntrts = ntrts
	estimate.TrtEff = &TrtEff
	estimate.Ratios = &Ratios

	return(estimate)
}

struct steppes_result scalar stmodelCI_class::test(
	pointer(class stsubpop_class scalar) scalar sp, real scalar nperm,
	struct steppes_effect scalar effect, | real scalar showstatus,
	real scalar seed, string scalar touse)
{
	struct steppes_result scalar test
	class AssociativeArray scalar Res, PermTest, result, temp_result, temp_TrtEff,
		temp_Ratios, temp_PermTest
	class stwin_class scalar win
	real scalar nsubpop, timePoint, ntrts, j, no, p, terminate, Ntemp,
		skip, i, overallsCI1, overallsCI2, index, overallSlogHR, s, rm,
		n_el
	real vector coltrt, survTime, trts, type, txassign, sel_ind, IndexSet1,
		IndexSet2, sel, sCI1, sCI2, sCISE1, sCISE2, slogHRs, slogHRSE, sObs, oObs,
		logHR, ologHR
	real matrix osubpop, differences, logHRs, Subpop, permuteSubpop, subpop
	string scalar todisp, spaces
	struct ssigma_struct scalar sigmas
	
	if (showstatus == .) showstatus = 1
	if (seed != .) rseed(seed)
	
	if (nperm > 0) {
		win = *sp->win
		nsubpop = sp->nsubpop
		osubpop = sp->subpop
		coltrt = this.coltrt
		survTime = this.coltime
		trts = this.trts
		timePoint = this.timePoint
		type = this.coltype

		ntrts = length(trts)

	  PermTest.reinit("string", 1)
		PermTest.put("sigma", J(nperm, nsubpop, 0))
		PermTest.put("HRsigma", J(nperm, nsubpop, 0))
		PermTest.put("pvalue", 0)
		PermTest.put("chi2pvalue", 0)
		PermTest.put("HRpvalue", 0)
		
		Res.reinit("real", 1)
		for (j = 1; j <= (ntrts - 1); j++) {
			Res.put(j, PermTest)
		}

		for (j = 2; j <= ntrts; j++) {
			// do the permutations --> compare trt1 with the rest
			
			txassign = J(length(coltrt), 1, -1)
		  sel_ind = selectindex(coltrt :== trts[1])
			txassign[sel_ind] = J(length(sel_ind), 1, 1)
		  sel_ind = selectindex(coltrt :== trts[j])
		  txassign[sel_ind] = J(length(sel_ind), 1, 0)
			
		  // set up the intial values
			differences = J(nperm, nsubpop, 0)
			logHRs = J(nperm, nsubpop, 0)
			no = 0
			p = 0
			terminate = 0
			Ntemp = rows(osubpop)
			IndexSet1 = selectindex(txassign :== 1)
			IndexSet2 = selectindex(txassign :== 0)
			
			if (showstatus) {
				printf("{txt}\n")
				printf("{txt}Computing the pvalue with ")
				printf("{res}%f", nperm)
				printf("{txt} number of permutations comparing trt %f with trt %f\n\n",
					trts[j], trts[1])
				printf("{hline 4}{c +}{hline 3} 1 ")
				printf("{hline 3}{c +}{hline 3} 2 ")
				printf("{hline 3}{c +}{hline 3} 3 ")
				printf("{hline 3}{c +}{hline 3} 4 ")
				printf("{hline 3}{c +}{hline 3} 5\n")
			}
			skip = 0
			while (no < nperm) {
				if (showstatus) {
					if (mod(no + 1, 50) == 0) {
						todisp = strtrim(strofreal(no + 1, "%9.0f"))
						if (strlen(todisp) < strlen(strofreal(nperm))) {
							skip = strlen(strofreal(nperm)) - strlen(todisp)
						}
						spaces = ""
						for (i = 1; i <= (skip + 3); i++) {
							spaces = spaces + char(32)
						}
						printf("{txt}." + spaces + todisp + "\n")
						skip = 0
					}
					else if (no + 1 == nperm) {
						printf("{txt}.\n")
					}
					else {
						printf("{txt}.")
					}
					displayflush()
				}
				
				sel = selectindex((txassign :== 1) :| (txassign :== 0))
				
				Subpop = osubpop
				permuteSubpop = J(Ntemp, nsubpop, 0)
				permuteSubpop[IndexSet1, .] = Subpop[jumble(IndexSet1), .]
				permuteSubpop[IndexSet2, .] = Subpop[jumble(IndexSet2), .]
				subpop = permuteSubpop
				sCI1 = J(nsubpop, 1, .)
				sCI2 = J(nsubpop, 1, .)
				sCISE1 = J(nsubpop, 1, .)
				sCISE2 = J(nsubpop, 1, .)
				overallsCI1 = .
				overallsCI2 = .
				slogHRs = J(nsubpop, 1, 0)
				slogHRSE = J(nsubpop, 1, 0)

				for (i = 1; i <= nsubpop; i++) {
					sel_ind = selectindex((subpop[., i] :== 1) :& (txassign :!= -1))
					result = cuminc_HR(survTime[sel_ind], type[sel_ind], txassign[sel_ind])
					
					temp_result = result.get("1 1")
					if (max(temp_result.get("time")) >= timePoint) {
						index = sum(temp_result.get("time") :<= timePoint)
						sCI1[i] = temp_result.get("est")[index]
						sCISE1[i] = sqrt(temp_result.get("var")[index])
					}
					temp_result = result.get("0 1")
					if (max(temp_result.get("time")) >= timePoint) {
						index = sum(temp_result.get("time") :<= timePoint)
						sCI2[i] = temp_result.get("est")[index]
						sCISE2[i] = sqrt(temp_result.get("var")[index])
					}
					slogHRs[i] = result.get("ome")/result.get("omevar")
					slogHRSE[i] = sqrt(1/result.get("omevar"))
				}

				result = cuminc_HR(survTime[sel], type[sel], txassign[sel])
				temp_result = result.get("1 1")
				if (max(temp_result.get("time")) >= timePoint) {
					index = sum(temp_result.get("time") :<= timePoint)
					overallsCI1 = temp_result.get("est")[index]
				}
				temp_result = result.get("0 1")
				if (max(temp_result.get("time")) >= timePoint) {
					index = sum(temp_result.get("time") :<= timePoint)
					overallsCI2 = temp_result.get("est")[index]
				}
				overallSlogHR = result.get("ome")/result.get("omevar")

				if (!hasmissing(sCI1) & !hasmissing(sCI2) &
					!hasmissing(overallsCI1) & !hasmissing(overallsCI2)) {
					no++
					p++
					for (s = 1; s <= nsubpop; s++) {
						differences[p, s] = (sCI1[s] - sCI2[s]) - (overallsCI1 - overallsCI2)
						logHRs[p, s] = slogHRs[s] - overallSlogHR
					}
				}
				terminate++
				if (terminate >= nperm + 10000) {
					printf("{err}After permuting %f plus 10000 times, ", nperm)
					printf("{err}the program was unable to generate the permutation ")
					printf("{err}distribution based on %f permutations of the data.\n", nperm)
					printf("{err}Consider creating larger subpopulations or selecting a ")
					printf("{err}different timepoint for estimation.\n")
					_error(3000)
				}
			}
			
			// generating the sigmas and pvalues
			temp_TrtEff = (*effect.TrtEff).get(1)
			sObs = temp_TrtEff.get("sObs")
			temp_TrtEff = (*effect.TrtEff).get(j)
			sObs = sObs - temp_TrtEff.get("sObs")
			temp_TrtEff = (*effect.TrtEff).get(1)
			oObs = temp_TrtEff.get("oObs")
			temp_TrtEff = (*effect.TrtEff).get(j)
			oObs = oObs - temp_TrtEff.get("oObs")

			temp_Ratios = (*effect.Ratios).get(j - 1)
			logHR = temp_Ratios.get("logHR")
			ologHR = temp_Ratios.get("ologHR")

			if (win.type == "tail-oriented") {
				// remove the trivial case of the full cohort
				rm = length(win.r1) + 1
				n_el = cols(differences)
				sel_ind = selectindex((1..n_el) :!= rm)
				differences = differences[., sel_ind]
				n_el = length(sObs)
				sel_ind = selectindex((1::n_el) :!= rm)
				sObs = sObs[sel_ind]
				n_el = length(oObs)
				sel_ind = selectindex((1::n_el) :!= rm)
				oObs = oObs[sel_ind]
				n_el = cols(logHRs)
				sel_ind = selectindex((1..n_el) :!= rm)
				logHRs = logHRs[., sel_ind]
				n_el = length(logHR)
				sel_ind = selectindex((1::n_el) :!= rm)
				logHR = logHR[sel_ind]
				n_el = length(ologHR)
				sel_ind = selectindex((1::n_el) :!= rm)
				ologHR = ologHR[sel_ind]
			}

			temp_PermTest = Res.get(j - 1)
			sigmas = ssigma(differences)
			temp_PermTest.put("sigma", sigmas.sigma)
			temp_PermTest.put("chi2pvalue", ppv(differences, sigmas.sigmainv, sObs,
				oObs, nperm))
			temp_PermTest.put("pvalue", ppv2(differences, sObs, oObs, nperm))

			sigmas = ssigma(logHRs)
			temp_PermTest.put("HRsigma", sigmas.sigma)
			temp_PermTest.put("HRpvalue", ppv2(logHRs, logHR, ologHR, nperm))
			Res.put(j - 1, temp_PermTest)
		}
		
		test.model = "CIt"
		test.ntrts = ntrts
		test.Res = &Res
	}
	
	return(test)
}

void print_estimate_CI(class steppes_class scalar x, real scalar timePoint,
	real vector trts)
{
	real scalar j, nsubpop, diff_oObs, diff_oSE
	real matrix temp
	string scalar title, fcname, subtitle
	string vector rnames
	string matrix cnames
	class AssociativeArray scalar TrtEff, temp_TrtEff, Ratios, temp_Ratios
	
	TrtEff = (*x.effect.TrtEff)
	Ratios = (*x.effect.Ratios)
	nsubpop = x.subpop->nsubpop
	
	for (j = 1; j <= x.effect.ntrts; j++) {
		// print out the cumulative incidence results
		title = "Cumulative incidence estimates for treatment group " +
			strofreal(trts[j]) + " at time point " + strofreal(timePoint)
		temp_TrtEff = TrtEff.get(j)
		temp = J(nsubpop, 3, .)
		temp[., 1] = (1::nsubpop)
		temp[., 2] = temp_TrtEff.get("sObs")
		temp[., 3] = temp_TrtEff.get("sSE")
		rnames = strofreal(temp[., 1])
		fcname = ("" \ "Subpopulation")
		cnames = ("Cumulative", "Standard" \ "incidence", "error")
		
		if (x.subpop->win->type != "tail-oriented") {
 			temp = (temp \ (., temp_TrtEff.get("oObs"), temp_TrtEff.get("oSE")))
			rnames = (rnames \ "Overall")
		}

		if (x.subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 12, title, "", 
				nsubpop, 0, 4, 0, 0,
				selectindex(temp[., 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 12, title, "",
				(nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
	}

	for (j = 2; j <= x.effect.ntrts; j++) {
		if (j == 2) {
			title = "Cumulative incidence differences at time point " +
				strofreal(timePoint)
		}
		else {
			title = ""
		}
		subtitle = "{txt}trt " + strofreal(trts[j]) + " vs. trt " +
			strofreal(trts[1])

		temp = J(nsubpop, 3, .)
		temp[., 1] = (1::nsubpop)
		temp_TrtEff = TrtEff.get(1)
		temp[., 2] = temp_TrtEff.get("sObs")
		temp_TrtEff = TrtEff.get(j)
		temp[., 2] = round(temp[., 2] - temp_TrtEff.get("sObs"), .0001)
		temp_TrtEff = TrtEff.get(1)
		temp[., 3] = temp_TrtEff.get("sSE"):^2
		temp_TrtEff = TrtEff.get(j)
		temp[., 3] = round(sqrt(temp[., 3] + temp_TrtEff.get("sSE"):^2), .0001)
		rnames = strofreal(temp[., 1])
		fcname = ("" \ "" \ "Subpopulation")
		cnames = ("Cumulative", "" \ "incidence", "Standard" \ "difference", "error")

		if (x.subpop->win->type != "tail-oriented") {
 			temp_TrtEff = TrtEff.get(1)
			diff_oObs = temp_TrtEff.get("oObs")
			diff_oSE = temp_TrtEff.get("oSE")
			temp_TrtEff = TrtEff.get(j)
			diff_oObs = diff_oObs - temp_TrtEff.get("oObs")
			diff_oSE = sqrt(diff_oSE^2 + temp_TrtEff.get("oSE")^2)
			temp = (temp \ (., diff_oObs, diff_oSE))
			rnames = (rnames \ "Overall")
		}
			
		if (x.subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 12, title,
				subtitle, nsubpop, 0, 4, 0, 0,
				selectindex(temp[., 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 12, title,
				subtitle, (nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
		
		if (j == 2) {
			title = "Hazard ratio estimates"
		}
		else {
			title = ""
		}
		subtitle = ""
		
		temp = J(nsubpop, 4, .)
		temp[., 1] = (1::nsubpop)
		temp_Ratios = Ratios.get(j - 1)
		temp[., 2] = temp_Ratios.get("logHR")
		temp[., 3] = temp_Ratios.get("logHRSE")
		temp[., 4] = exp(temp[., 2])
		rnames = strofreal(temp[., 1])
		fcname = "Subpopulation"
		cnames = ("Log HR", "Std. err.", "Hazard ratio")

		if (x.subpop->win->type != "tail-oriented") {
			temp = (temp \ (., temp_Ratios.get("ologHR"),
				temp_Ratios.get("ologHRSE"), exp(temp_Ratios.get("ologHR"))))
			rnames = (rnames \ "Overall")
		}
			
		if (x.subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..4], rnames, cnames, fcname, 15, 14, title,
				subtitle, nsubpop, 0, 4, 0, 0,
				selectindex(temp[., 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..4], rnames, cnames, fcname, 15, 14, title,
				subtitle, (nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
	}
}

void print_cov_CI(class steppes_class scalar stobj, real scalar timePoint,
	real vector trts)
{
	real scalar nsubpop, j, ns, longest
	class AssociativeArray scalar temp_t, t
	string vector names
	struct steppes_result scalar result_empty
	
	nsubpop = stobj.subpop->nsubpop
	names = J(nsubpop, 1, "")
	for (j = 1; j <= nsubpop; j++) {
		names[j] = "SP" + strofreal(j) + "-Overall"
	}

	longest = min((max(strlen(names)), 11))
	if (stobj.result != result_empty) {
		temp_t = (*stobj.result.Res)
		for (j = 1; j <= (stobj.result.ntrts - 1); j++) {
			ns = stobj.subpop->nsubpop
			if (stobj.subpop->win->type == "tail-oriented") ns--
			printf("\n")
			printf("{txt}The covariance matrix of the cumulative incidence ")
			printf("{txt}differences at {res}%f {txt}time units for the ", timePoint)
			printf("{res}%f {txt}subpopulations is:\n", ns)
			printf("{txt}trt %f vs. trt %f\n", trts[j + 1], trts[1])
			t = temp_t.get(j)
			print_matrix(t.get("sigma"), names, names, "", ., longest)

			printf("\n")
			printf("{txt}The covariance matrix of the log hazard ratios for the ")
			printf("{res}%f {txt}subpopulations is:\n", ns)
			print_matrix(t.get("HRsigma"), names, names, "", ., longest)
		}
	}
}

void print_stat_CI(class steppes_class scalar stobj, real vector trts)
{
	real scalar j
	class AssociativeArray scalar temp_t, t
	struct steppes_result scalar result_empty
	
	if (stobj.result != result_empty) {
		temp_t = (*stobj.result.Res)
		for (j = 1; j <= (stobj.result.ntrts - 1); j++) {
			t = temp_t.get(j)
			printf("\n")
			printf("{txt}Supremum test results\n")
			printf("{txt}trt %f vs. trt %f\n", trts[j + 1], trts[1])
			printf("{txt}Interaction p-value based on cumulative incidence estimates: {res}%f\n", t.get("pvalue"))
			printf("{txt}Interaction p-value based on hazard ratio estimates: {res}%f\n", t.get("HRpvalue"))

			printf("\n")
			printf("{txt}Chi-square test results\n")
			printf("{txt}Interaction p-value based on cumulative incidence estimates: {res}%f\n", t.get("chi2pvalue"))
		}
	}
}

void stmodelCI_class::print(class steppes_class scalar stobj,
	real scalar estimate, real scalar cov, real scalar test)
{
	// 1. estimates
	if (estimate) {
		print_estimate_CI(stobj, timePoint, trts)
	}

	// 2. covariance matrices
	if (cov) {
		print_cov_CI(stobj, timePoint, trts)
	}

	// 3. supremum test and chi-square test results
	if (test) {
		print_stat_CI(stobj, trts)
	}
	
	displayflush()
}

end
