*!stepp_stmodelKM version 0.1.0
*!Written 08Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void stmodelKM_class::new()
{
	class_name = "stmodelKM"
}

string scalar stmodelKM_class::getclass()
{
	return(class_name)
}

struct steppes_effect scalar stmodelKM_class::estimate(
	pointer(class stsubpop_class scalar) scalar sp, string scalar touse)
{
	struct steppes_effect scalar estimate
	real scalar nsubpop, timePoint, ntrts, i, j, overallSkmObs, overallSkmSE,
		overallLogHR, overallLogHRSE, logHRw, skmw
	real vector survTime, coltrt, trts, censor, txassign, skmObs, skmSE,sel,
		logHR, logHRSE, sel_ind
	real matrix subpop
	class AssociativeArray scalar trteff, overalltrteff, Obs, TrtEff, HRs, Ratios,
		temp_TrtEff, temp_Ratios, LogRank, temp_TrtEff_1, temp_TrtEff_j
	
	nsubpop = sp->nsubpop
	subpop = sp->subpop
	coltrt = this.coltrt
	survTime = this.survTime
	trts = this.trts
	timePoint = this.timePoint
	censor = this.censor

	ntrts = length(trts)
	txassign = J(length(coltrt), 1, .)
	for (j = 1; j <= ntrts; j++) {
		txassign[selectindex(coltrt :== trts[j])] = J(sum(coltrt :== trts[j]), 1, j)
	}
	
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

	temp_TrtEff.reinit("string", 1)
	
	// estimate treatment effect for each treatment j
	for (j = 1; j <= ntrts; j++) {
		skmObs = J(nsubpop, 1, 0)
		skmSE  = J(nsubpop, 1, 0)
		for (i = 1; i <= nsubpop; i++) {
			trteff = tpest(kmest(
				survTime[selectindex(txassign :== j :& subpop[., i] :== 1)],
				censor[selectindex(txassign :== j :& subpop[., i] :== 1)]), timePoint)
			skmObs[i] = minmax((trteff.get("est"), 0), 1)[2]
			skmSE[i]  = sqrt(trteff.get("var"))
		}
		overalltrteff = tpest(kmest(
			survTime[selectindex(txassign :== j)],
			censor[selectindex(txassign :== j)]), timePoint)
		overallSkmObs = minmax((overalltrteff.get("est"), 0), 1)[2]
		overallSkmSE  = sqrt(overalltrteff.get("var"))
		
		// check to make sure that subpopulation estimates are OK
		if ((sum(skmObs :== .) > 0) | (overallSkmObs == .)) {
			printf("\n")
			printf("{err}Unable to estimate survival time at " +
				strofreal(timePoint) + " time-unit(s) for trt " + strofreal(j) + "\n")
			printf("{err}because there are too few events within one or more subpopulation(s).\n")
			printf("{err}The problem may be avoided by constructing larger ")
			printf("{err}subpopulations \n")
			printf("{err}and/or by selecting a different timepoint for estimation.\n")
			_error(3000)
		}
		
		temp_TrtEff.put("sObs", skmObs)
		temp_TrtEff.put("sSE", skmSE)
		temp_TrtEff.put("oObs", overallSkmObs)
		temp_TrtEff.put("oSE", overallSkmSE)
		TrtEff.put(j, temp_TrtEff)
	}

	// 	estimate the relative treatment effect comparing each one with trt 1
	for (j = 2; j <= ntrts; j++) {
		txassign = J(length(coltrt), 1, -1)
		txassign[selectindex(coltrt :== trts[1])] =
			J(length(selectindex(coltrt :== trts[1])), 1, 1)
		txassign[selectindex(coltrt :== trts[j])] =
			J(length(selectindex(coltrt :== trts[j])), 1, 0)

		logHR = J(nsubpop, 1, 0)
		logHRSE = J(nsubpop, 1, 0)
		sel = (txassign :== 1 :| txassign :== 0)

		// for each subgroup i
		for (i = 1; i <= nsubpop; i++) {
			sel_ind = selectindex((subpop[., i] :== 1) :& (sel :== 1))
			LogRank = survdiff_fit(survTime[sel_ind], censor[sel_ind],
				as_factor(txassign[sel_ind]))
			logHR[i] = -(LogRank.get("observed")[1] - LogRank.get("expected")[1])/
				LogRank.get("var")[1, 1]
			logHRSE[i] = sqrt(1/LogRank.get("var")[1, 1])
		}

		// estimate the overall effect
		LogRank = survdiff_fit(survTime[selectindex(sel)],
			censor[selectindex(sel)], as_factor(txassign[selectindex(sel)]))
		overallLogHR = -(LogRank.get("observed")[1] - LogRank.get("expected")[1]) /
			LogRank.get("var")[1, 1]
		overallLogHRSE = sqrt(1/LogRank.get("var")[1, 1])
		logHRw = sum(logHR :/ logHRSE)
		temp_TrtEff_1 = TrtEff.get(1)
		temp_TrtEff_j = TrtEff.get(j)
		skmw = sum((temp_TrtEff_1.get("sObs") - temp_TrtEff_j.get("sObs")) :/
			sqrt(temp_TrtEff_1.get("sSE"):^2 + temp_TrtEff_j.get("sSE"):^2))
		
		temp_Ratios = Ratios.get(j - 1)
		temp_Ratios.put("skmw", skmw)
		temp_Ratios.put("logHR", logHR)
		temp_Ratios.put("logHRSE", logHRSE)
		temp_Ratios.put("ologHR", overallLogHR)
		temp_Ratios.put("ologHRSE", overallLogHRSE)
		temp_Ratios.put("logHRw", logHRw)
		Ratios.put(j - 1, temp_Ratios)
	}

	// return the estimate as a structure
	estimate.model = "KMe"
	estimate.ntrts = ntrts
	estimate.TrtEff = &TrtEff
	estimate.Ratios = &Ratios

	return(estimate)
}

struct steppes_result scalar stmodelKM_class::test(
	pointer(class stsubpop_class scalar) scalar sp, real scalar nperm,
	struct steppes_effect scalar effect, string scalar timevar,
	string scalar failvar, string scalar trtvar, | real scalar showstatus,
	real scalar Cox, string vector MM, real scalar seed, string scalar touse)
{
	struct steppes_result scalar test
	class AssociativeArray scalar Res, PermTest, seff1, seff2, LogRank,
		trteff, temp_TrtEff, temp_Ratios, temp_PermTest
	class stwin_class scalar win
	real scalar nsubpop, timePoint, ntrts, j, no, p, sel_idx, selo_idx, terminate,
		Ntemp, skip, i, overallSkm1, overallSkm2, slogHRw, overallSLogHR,
		overallSLogHRSE, s, rm, n_el
	real vector coltrt, survTime, trts, censor, txassign, sel, sel_ind, IndexSet1,
		IndexSet2, skm1, skm2, skmSE1, skmSE2, slogHRs, slogHRSEs, selo, selo_ind,
		sObs, oObs, logHR, ologHR
	real matrix osubpop, differences, logHRs, Subpop, permuteSubpop, subpop
	string scalar fmla, xnam, sel_name, selo_name, todisp, spaces
	struct ssigma_struct scalar sigmas
	
	if (showstatus == .) showstatus = 1
	if (Cox == .) Cox = 0
	if (seed != .) rseed(seed)
	
	if (nperm > 0) {
		win = *sp->win
		nsubpop = sp->nsubpop
		osubpop = sp->subpop
		coltrt = this.coltrt
		survTime = this.survTime
		trts = this.trts
		timePoint = this.timePoint
		censor = this.censor

		ntrts = length(trts)
		txassign = J(length(coltrt), 1, -1)
		
		for (j = 1; j <= ntrts; j++) {
			txassign[selectindex(coltrt :== trts[j])] =
				J(sum(coltrt :== trts[j]), 1, j)
		}

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

		stata("preserve", 1)
		for (j = 2; j <= ntrts; j++) {
			// do the permutations --> compare trt1 with the rest
			differences = J(nperm, nsubpop, 0)
			logHRs = J(nperm, nsubpop, 0)
			no = 0
			p = 0
			terminate = 0
			Ntemp = rows(osubpop)
			IndexSet1 = (1::Ntemp)[selectindex(txassign :== 1)]
			IndexSet2 = (1::Ntemp)[selectindex(txassign :== j)]
			
			if (Cox) {
				if (MM == J(1, 0, "")) {
				 fmla = trtvar
				}
				else {
					xnam = invtokens(MM)
					fmla = trtvar + " " + xnam
				}
				stata("quietly stset " + timevar + ", failure(" + failvar + ")", 1)
			}
			
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
				
				Subpop = osubpop
				permuteSubpop = J(Ntemp, nsubpop, 0)
				permuteSubpop[selectindex(txassign :== 1), .] =
					Subpop[jumble(IndexSet1), .]
				permuteSubpop[selectindex(txassign :== j), .] =
					Subpop[jumble(IndexSet2), .]
				subpop = permuteSubpop
				skm1 = J(nsubpop, 1, 0)
				skm2 = J(nsubpop, 1, 0)
				skmSE1 = J(nsubpop, 1, 0)
				skmSE2 = J(nsubpop, 1, 0)
				slogHRs = J(nsubpop, 1, 0)
				slogHRSEs = J(nsubpop, 1, 0)

				for (i = 1; i <= nsubpop; i++) {
					seff1 = tpest(kmest(
						survTime[selectindex(txassign :== 1 :& subpop[., i] :== 1)],
						censor[selectindex(txassign :== 1 :& subpop[., i] :== 1)]),
						timePoint)
					seff2 = tpest(kmest(
						survTime[selectindex(txassign :== j :& subpop[., i] :== 1)],
						censor[selectindex(txassign :== j :& subpop[., i] :== 1)]),
						timePoint)
					skm1[i] = max((seff1.get("est"), 0))
					skm2[i] = max((seff2.get("est"), 0))
					skmSE1[i] = sqrt(seff1.get("var"))
					skmSE2[i] = sqrt(seff2.get("var"))
					sel = (txassign :== 1 :| txassign :== j) :& (subpop[., i] :== 1)
					sel_ind = selectindex(sel)
				
					if (Cox) {
						sel_idx = st_addvar("int", sel_name = st_tempname())
						st_store(selectindex(st_data(., touse)), sel_idx, sel)
						stata("quietly stcox " + fmla + " if " + sel_name, 1)
						if (MM == J(1, 0, "")) {
							slogHRs[i] = st_matrix("e(b)")
							slogHRSEs[i] = sqrt(st_matrix("e(V)"))
						} else {
							slogHRs[i] = st_matrix("e(b)")[1]
							slogHRSEs[i] = sqrt(st_matrix("e(V)")[1, 1])
						}
						st_dropvar(sel_name)
					} else {
						LogRank = survdiff_fit(survTime[sel_ind], censor[sel_ind],
							as_factor(txassign[sel_ind]))
						slogHRs[i] = -(LogRank.get("observed")[1] - LogRank.get("expected")[1])/
							LogRank.get("var")[1, 1]
						slogHRSEs[i] = sqrt(1/LogRank.get("var")[1, 1])
					}
				}
				
				selo = (txassign :== 1 :| txassign :== j)
				selo_ind = selectindex(selo)

				slogHRw	= sum(slogHRs :/ slogHRSEs)

				trteff = tpest(kmest(survTime[selectindex(txassign :== 1)],
					censor[selectindex(txassign :== 1)]), timePoint)
				overallSkm1 = max((trteff.get("est"), 0))
				trteff = tpest(kmest(survTime[selectindex(txassign :== j)],
					censor[selectindex(txassign :== j)]), timePoint)
				overallSkm2 = max((trteff.get("est"), 0))

				if (Cox) {
					selo_idx = st_addvar("int", selo_name = st_tempname())
					st_store(selectindex(st_data(., touse)), selo_idx, selo)
					stata("quietly stcox " + fmla + " if " + selo_name, 1)
					if (MM == J(1, 0, "")) {
						overallSLogHR = st_matrix("e(b)")
						overallSLogHRSE = sqrt(st_matrix("e(V)"))
					} else {
						overallSLogHR = st_matrix("e(b)")[1]
						overallSLogHRSE = sqrt(st_matrix("e(V)")[1, 1])
					}
					st_dropvar(selo_name)
				} else {
					LogRank = survdiff_fit(survTime[selo_ind], censor[selo_ind],
						as_factor(txassign[selo_ind]))
					overallSLogHR = -(LogRank.get("observed")[1] - LogRank.get("expected")[1])/
						LogRank.get("var")[1, 1]
					overallSLogHRSE = sqrt(1/LogRank.get("var")[1, 1])
				}

				if (!hasmissing(skm1) & !hasmissing(skm2) &
					!hasmissing(overallSkm1) & !hasmissing(overallSkm2)) {
					no++
					p++
					for (s = 1; s <= nsubpop; s++) {
						differences[p, s] = (skm1[s] - skm2[s]) - (overallSkm1 - overallSkm2)
						logHRs[p, s] = slogHRs[s] - overallSLogHR
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
			if (Cox) {
				stata("quietly stset, clear", 1)
				st_rclear()
				st_eclear()
			}
			stata("restore", 1)
			
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
		
		test.model = "KMt"
		test.ntrts = ntrts
		test.Res = &Res
	}
	
	pragma unused slogHRw
	pragma unused overallSLogHRSE
	
	return(test)
}

void print_estimate_KM(class steppes_class scalar x, real scalar timePoint,
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
		// print out the Kaplan-Meier results
		title = "Survival estimates for treatment group " +
			strofreal(trts[j]) + " at time point " + strofreal(timePoint)
		temp_TrtEff = TrtEff.get(j)
		temp = J(nsubpop, 3, .)
		temp[., 1] = (1::nsubpop)
		temp[., 2] = temp_TrtEff.get("sObs")
		temp[., 3] = temp_TrtEff.get("sSE")
		rnames = strofreal(temp[., 1])
		fcname = ("" \ "Subpopulation")
		cnames = ("Survival", "Standard" \ "probability", "error")
		
		if (x.subpop->win->type != "tail-oriented") {
 			temp = (temp \ (., temp_TrtEff.get("oObs"), temp_TrtEff.get("oSE")))
			rnames = (rnames \ "Overall")
		}

		if (x.subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 14, title, "", 
				nsubpop, 0, 4, 0, 0,
				selectindex(temp[, 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 14, title, "",
				(nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
	}

	printf("\n")
	printf("Survival differences and hazard ratio estimates\n")
	for (j = 2; j <= x.effect.ntrts; j++) {
		if (j == 2) {
			title = "Survival differences at time point " + strofreal(timePoint)
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
		fcname = ("" \ "Subpopulation")
		cnames = ("Survival", "Standard" \ "difference", "error")

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
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 14, title,
				subtitle, nsubpop, 0, 4, 0, 0,
				selectindex(temp[, 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 14, title,
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
				selectindex(temp[, 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..4], rnames, cnames, fcname, 15, 14, title,
				subtitle, (nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
	}
}

void print_cov_KM(class steppes_class scalar stobj, real scalar timePoint,
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
			printf("{txt}The covariance matrix of the Kaplan-Meier ")
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

void print_stat_KM(class steppes_class scalar stobj, real vector trts)
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
			printf("{txt}Interaction p-value based on Kaplan-Meier estimates: {res}%f\n", t.get("pvalue"))
			printf("{txt}Interaction p-value based on hazard ratio estimates: {res}%f\n", t.get("HRpvalue"))

			printf("\n")
			printf("{txt}Chi-square test results\n")
			printf("{txt}Interaction p-value based on Kaplan-Meier estimates: {res}%f\n", t.get("chi2pvalue"))
		}
	}
}

void stmodelKM_class::print(class steppes_class scalar stobj,
	real scalar estimate, real scalar cov, real scalar test)
{
	// 1. estimates
	if (estimate) {
		print_estimate_KM(stobj, timePoint, trts)
	}

	// 2. covariance matrices
	if (cov) {
		print_cov_KM(stobj, timePoint, trts)
	}

	// 3. supremum test and chi-square test results
	if (test) {
		print_stat_KM(stobj, trts)
	}
	
	displayflush()
}

end
