*!stepp_stmodelGLM version 0.1.0
*!Written 14Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

void stmodelGLM_class::new()
{
	class_name = "stmodelGLM"
// 	if ((glm != "binomial") & (glm != "poisson") & (glm != "gaussian")) {
// 		printf("{err}the GLM model can either be 'binomial', 'poisson' or 'gaussian'\n")
// 		_error(3001)
// 	}
}

string scalar stmodelGLM_class::getclass()
{
	return(class_name)
}

real scalar stmodelGLM_class::validity()
{
	real valid
	
	valid = 1
	
	if ((glm == "binomial") & (link == "")) {
		this.link = "logit"
	}
	if ((glm == "poisson") & (link == "")) {
		this.link = "log"
	}
	if ((glm == "gaussian") & (link == "")) {
		this.link = "identity"
	}
	
	if ((glm != "binomial") & (glm != "poisson") & (glm != "gaussian")) {
		printf("{err}the GLM model can either be 'binomial', 'poisson' or 'gaussian'\n")
		valid = 0
	}
	if ((glm == "binomial") & ((link != "logit") & (link != "probit"))) {
		printf("{err}the binomial model only allows 'logit' or 'probit' as the link function\n")
		valid = 0
	}
	if ((glm == "poisson") & (link != "log")) {
		printf("{err}the poisson model only allows 'log' as the link function\n")
		valid = 0
	}
	if ((glm == "gaussian") & (link != "identity")) {
		printf("{err}the gaussian model only allows 'identity' as the link function\n")
		valid = 0
	}
	
	if (glm == "") {
		printf("the glm model type must be provided\n")
		valid = 0
	}
	
	return(valid)
}

struct steppes_effect scalar stmodelGLM_class::estimate(
	pointer(class stsubpop_class scalar) scalar sp, string scalar touse)
{
	struct steppes_effect scalar estimate
	real scalar nsubpop, ntrts, i, j, OuTcOmE_idx, txassign_idx, has_intercept,
		cstart, seli_idx, logHRw, skmw, h, p1_fit, p1_se, p2_fit, p2_se, op1_fit,
		op1_se, op2_fit, op2_se, oexp1, oexp2
	real vector OuTcOmE, coltrt, trts, txassign, DIFF, DIFFSE, seli, seli_ind, cm,
		tmp_sObs_1, tmp_sObs_j, tmp_sSE_1, tmp_sSE_j, tmp_logHR, tmp_logHRSE, onv1
	real matrix subpop, mm1, om1
	class AssociativeArray scalar Obs, TrtEff, HRs, Ratios, temp_Ratios,
		temp_TrtEff_1, temp_TrtEff_j
	string scalar glm, link, model,OuTcOmE_nm, txassign_nm, fmla, xnam, seli_name,
		model_glm, cm_txt, onv1_txt
	string vector covariate, cm_nm, onv1_nm

	if (this.validity()) {
		nsubpop = sp->nsubpop
		subpop = sp->subpop
		coltrt = this.coltrt
		OuTcOmE = this.colY
		trts = this.trts
		covariate = this.MM
		glm = this.glm
		link = this.link
		debug = this.debug

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

		// create a formula for the model
		OuTcOmE_idx = st_addvar("double", OuTcOmE_nm = st_tempname())
		st_store(selectindex(st_data(., touse)), OuTcOmE_idx, OuTcOmE)
		txassign_idx = st_addvar("double", txassign_nm = st_tempname())
		st_store(selectindex(st_data(., touse)), txassign_idx, txassign)
		if (covariate != "") {
			has_intercept = (covariate[1] == "_cons_")
			if (has_intercept) {
				if  (length(covariate) > 1) {
					cstart = 2
					xnam = (txassign_nm, covariate[cstart::length(covariate)])
					fmla = OuTcOmE_nm + " " + invtokens(xnam)
				}
				else {
					cstart = 1
					xnam = txassign_nm
					fmla = OuTcOmE_nm + " " + invtokens(xnam)
				}
			}
			else {
				cstart = 1
				xnam = (txassign_nm, covariate[cstart::length(covariate)])
				fmla = OuTcOmE_nm + " " + invtokens(xnam)
			}
		}
		else {
			cstart = 2   // the default model has an intercept
			fmla = OuTcOmE_nm + " " + txassign_nm
		}
		
		// estimate treatment effect for each treatment j
		for (j = 2; j <= ntrts; j++) {
			DIFF = J(nsubpop, 1, 0)
			DIFFSE = J(nsubpop, 1, 0)
			tmp_sObs_1 = J(nsubpop, 1, 0)
			tmp_sObs_j = J(nsubpop, 1, 0)
			tmp_sSE_1 = J(nsubpop, 1, 0)
			tmp_sSE_j = J(nsubpop, 1, 0)
			temp_TrtEff_1 = TrtEff.get(1)
			temp_TrtEff_j = TrtEff.get(j)
			tmp_logHR = J(nsubpop, 1, 0)
			tmp_logHRSE = J(nsubpop, 1, 0)
			temp_Ratios = Ratios.get(j - 1)
		
			seli = (txassign :== 1 :| txassign :== j)

			for (i = 1; i <= nsubpop; i++) {
				seli_ind = (subpop[., i] :== 1) :& (seli :== 1)
				seli_idx = st_addvar("int", seli_name = st_tempname())
				st_store(selectindex(st_data(., touse)), seli_idx, seli_ind)

				if (glm == "gaussian") {
					if (link == "identity") {
						model_glm = "regress " + fmla + " if " + seli_name + ","
					}
				}
				else if (glm == "binomial") {
					if (link == "logit") {
						model_glm = "logit " + fmla + " if " + seli_name + ","
					}
				}
				else if (glm == "poisson") {
					if (link == "logarithm") {
						model_glm = "poisson " + fmla + " if " + seli_name + ","
					}
				}
				else {
					model_glm = "glm " + fmla + " if " + seli_name + ", family(" + glm +
						") link(" + link + ")"
				}
				if (!has_intercept) model_glm = model_glm + " noconstant"
				stata("quietly " + model_glm, 1)

		    if (covariate != "") {
					if (has_intercept) {
						if (length(covariate) > 1) {
							mm1 = st_data(., covariate[selectindex(covariate :!= "_cons_")],
								seli_name)
							mm1 = (J(rows(mm1), 1, 1), mm1)
						}
						else {
							mm1 = J(sum(seli_ind), 1, 1)
						}
					}
					else {
						mm1 = st_data(., covariate, seli_name)
					}
		      mm1 = (txassign[selectindex(seli_ind)], mm1)
		      cm = mean(mm1)
					cm_nm = txassign_nm
		      if (has_intercept) {
						if (cols(mm1) > 2) {
							cm = (cm[1], cm[3..cols(mm1)])
							cm_nm = (cm_nm, covariate[2..length(covariate)])
						}
						else {
							cm = cm[1]
						}
					}
					else {
						cm_nm = (cm_nm, covariate)
					}
		      cm[1] = 1  	// treatment indicator is 1
					cm_txt = ""
					for (h = 1; h < length(cm); h++) {
						cm_txt = cm_txt + cm_nm[h] + " = " + strofreal(cm[h]) + " "
					}
					cm_txt = cm_txt + cm_nm[length(cm)] + " = " + strofreal(cm[length(cm)])
					stata("quietly margins , at(" + cm_txt + ") predict(xb)", 1)
					p1_fit = st_matrix("r(b)")
					p1_se = sqrt(st_matrix("r(V)"))
		      cm[1] = j	  // treatment indicator is j
					cm_txt = ""
					for (h = 1; h < length(cm); h++) {
						cm_txt = cm_txt + cm_nm[h] + " = " + strofreal(cm[h]) + " "
					}
					cm_txt = cm_txt + cm_nm[length(cm)] + " = " + strofreal(cm[length(cm)])
					stata("quietly margins , at(" + cm_txt + ") predict(xb)", 1)
					p2_fit = st_matrix("r(b)")
					p2_se = sqrt(st_matrix("r(V)"))
		    } else {
					stata("quietly margins , at(" + txassign_nm + " = 1) predict(xb)", 1)
					p1_fit = st_matrix("r(b)")
					p1_se = sqrt(st_matrix("r(V)"))
					stata("quietly margins , at(" + txassign_nm + " = " + strofreal(j) + ") predict(xb)", 1)
					p2_fit = st_matrix("r(b)")
					p2_se = sqrt(st_matrix("r(V)"))
		    }
				st_dropvar(seli_name)
				
				if (glm == "gaussian") {
					tmp_sObs_1[i] = p1_fit
					tmp_sSE_1[i] = p1_se

					tmp_sObs_j[i] = p2_fit
					tmp_sSE_j[i] = p2_se

					DIFF[i] = st_matrix("e(b)")[1]
					DIFFSE[i] = sqrt(st_matrix("e(V)")[1, 1])
					
					tmp_logHR[i] = log(tmp_sObs_1[i])-log(tmp_sObs_j[i])
					// use delta method
					tmp_logHRSE[i] = sqrt((tmp_sSE_1[i]^2)/(tmp_sObs_1[i]^2) + (tmp_sSE_j[i]^2)/(tmp_sObs_j[i]^2))
				} else if (glm == "binomial") {
					tmp_sObs_1[i] = 1/(1 + exp(-p1_fit))
					tmp_sSE_1[i] = p1_se*(tmp_sObs_1[i]^2)/exp(p1_fit)   // delta method

					tmp_sObs_j[i] = 1/(1 + exp(-p2_fit))
					tmp_sSE_j[i] = p2_se*(tmp_sObs_j[i]^2)/exp(p2_fit)   // delta method

					DIFF[i] = tmp_sObs_1[i] - tmp_sObs_j[i]
					DIFFSE[i] =  sqrt(tmp_sSE_1[i]^2 + tmp_sSE_j[i]^2)
					
					tmp_logHR[i] = st_matrix("e(b)")[1]   // log OR from the model coefficient
					tmp_logHRSE[i] = sqrt(st_matrix("e(V)")[1, 1])
				} else if (glm == "poisson") {
					tmp_sObs_1[i] = exp(p1_fit)
					tmp_sSE_1[i] = p1_se*tmp_sObs_1[i]   // delta method

					tmp_sObs_j[i] = exp(p2_fit)
					tmp_sSE_j[i] = p2_se*tmp_sObs_j[i]   // delta method

					DIFF[i] = tmp_sObs_1[i] - tmp_sObs_j[i]
					DIFFSE[i] =  sqrt(tmp_sSE_1[i]^2 + tmp_sSE_j[i]^2)

					tmp_logHR[i] = st_matrix("e(b)")[1]   // log RR from the model coefficient
					tmp_logHRSE[i] = sqrt(st_matrix("e(V)")[1, 1])
				}
			}
			temp_TrtEff_1.put("sObs", tmp_sObs_1)
			temp_TrtEff_1.put("sSE", tmp_sSE_1)
			temp_TrtEff_j.put("sObs", tmp_sObs_j)
			temp_TrtEff_j.put("sSE", tmp_sSE_j)
			
			temp_Ratios.put("logHR", tmp_logHR)
			temp_Ratios.put("logHRSE", tmp_logHRSE)			
			skmw = sum(DIFF :/ DIFFSE)
			temp_Ratios.put("skmw", skmw)
		  logHRw = sum(tmp_logHR :/ tmp_logHRSE)
			temp_Ratios.put("logHRw", logHRw)
			
			// estimate the overall effect
			if (glm == "guassian") {
				if (link == "identity") {
					model_glm = "regress " + fmla + ","
				}
			}
			if (glm == "binomial") {
				if (link == "logit") {
					model_glm = "logit " + fmla + ","
				}
			}
			else if (glm == "poisson") {
				if (link == "logarithm") {
					model_glm = "poisson " + fmla + ","
				}
			}
			else {
				model_glm = "glm " + fmla + ", family(" + glm + ") link(" + link + ")"
			}
			if (!has_intercept) model_glm = model_glm + " noconstant"
			stata("quietly " + model_glm, 1)

			if (covariate != "") {
				if (has_intercept) {
					if (length(covariate) > 1) {
						om1 = st_data(., covariate[selectindex(covariate :!= "_cons_")])
						om1 = (J(rows(om1), 1, 1), om1)
					}
					else {
						om1 = J(length(txassign), 1, 1)
					}
				}
				else {
					om1 = st_data(., covariate)
				}
				om1 = (txassign, om1)
				onv1 = mean(om1)
				onv1_nm = txassign_nm
				if (has_intercept) {
					if (cols(om1) > 2) {
						onv1 = (onv1[1], onv1[3..cols(om1)])
						onv1_nm = (onv1_nm, covariate[2..length(covariate)])
					}
					else {
						onv1 = onv1[1]
					}
				}
				else {
					onv1_nm = (onv1_nm, covariate)
				}
				onv1[1] = 1  	// treatment indicator is 1
				onv1_txt = ""
				for (h = 1; h < length(onv1); h++) {
					onv1_txt = onv1_txt + onv1_nm[h] + " = " + strofreal(onv1[h]) + " "
				}
				onv1_txt = onv1_txt + onv1_nm[length(onv1)] + " = " + strofreal(onv1[length(onv1)])
				stata("quietly margins , at(" + onv1_txt + ") predict(xb)", 1)
				op1_fit = st_matrix("r(b)")
				op1_se = sqrt(st_matrix("r(V)"))
				onv1[1] = j	  // treatment indicator is j
				onv1_txt = ""
				for (h = 1; h < length(onv1); h++) {
					onv1_txt = onv1_txt + onv1_nm[h] + " = " + strofreal(onv1[h]) + " "
				}
				onv1_txt = onv1_txt + onv1_nm[length(onv1)] + " = " + strofreal(onv1[length(onv1)])
				stata("quietly margins , at(" + onv1_txt + ") predict(xb)", 1)
				op2_fit = st_matrix("r(b)")
				op2_se = sqrt(st_matrix("r(V)"))
			} else {
				stata("quietly margins , at(" + txassign_nm + " = 1) predict(xb)", 1)
				op1_fit = st_matrix("r(b)")
				op1_se = sqrt(st_matrix("r(V)"))
				stata("quietly margins , at(" + txassign_nm + " = " + strofreal(j) + ") predict(xb)", 1)
				op2_fit = st_matrix("r(b)")
				op2_se = sqrt(st_matrix("r(V)"))
			}
			
			if (glm == "gaussian") {
				temp_TrtEff_1.put("oObs", op1_fit)
				temp_TrtEff_1.put("oSE", op1_se)

				temp_TrtEff_j.put("oObs", op2_fit)
				temp_TrtEff_j.put("oSE", op2_se)

				temp_Ratios.put("ologHR", log(temp_TrtEff_1.get("oObs")) -
					temp_TrtEff_j.get("oObs"))
				// delta method
				temp_Ratios.put("ologHRSE", sqrt(
					(temp_TrtEff_1.get("oSE")^2)/(temp_TrtEff_1.get("oObs")^2) +
					(temp_TrtEff_j.get("oSE")^2)/(temp_TrtEff_j.get("oObs")^2)))
			} else if (glm == "binomial") {
				oexp1 = exp(op1_fit)
				temp_TrtEff_1.put("oObs", oexp1/(1 + oexp1))
				temp_TrtEff_1.put("oSE", op1_se*(oexp1/((1 + oexp1)^2)))   // delta method

				oexp2 = exp(op2_fit)
				temp_TrtEff_j.put("oObs", oexp2/(1 + oexp2))
				temp_TrtEff_j.put("oSE", op2_se*(oexp2/((1 + oexp2)^2)))   // delta method

				temp_Ratios.put("ologHR", st_matrix("e(b)")[1])
				temp_Ratios.put("ologHRSE", sqrt(st_matrix("e(V)")[1, 1]))
			} else if (glm == "poisson") {
				temp_TrtEff_1.put("oObs", exp(op1_fit))
				temp_TrtEff_1.put("oSE", op1_se*exp(op1_fit))   // delta method

				temp_TrtEff_j.put("oObs", exp(op2_fit))
				temp_TrtEff_j.put("oSE", op2_se*exp(op2_fit))   // delta method

				temp_Ratios.put("ologHR", st_matrix("e(b)")[1])
				temp_Ratios.put("ologHRSE", sqrt(st_matrix("e(V)")[1, 1]))
			}

			if (any(tmp_sObs_1 :== .) | any(tmp_sObs_j :== .) |
				(temp_TrtEff_1.get("oObs") == .) | (temp_TrtEff_j.get("oObs") == .)) {
				printf("\n")
				printf("{err}Unable to estimate the effect because there are too " + 
					"few events within one or more subpopulation(s).\n")
				printf("{err}The problem may be avoided by constructing larger " +
					"subpopulations.\n")
				_error(3000)
			}

			TrtEff.put(j, temp_TrtEff_j)
			Ratios.put(j - 1, temp_Ratios)
		}
		TrtEff.put(1, temp_TrtEff_1)
		
		st_dropvar(OuTcOmE_nm)
		st_dropvar(txassign_nm)
		st_rclear()
		st_eclear()

		if (glm == "gaussian") {
			model = "GLMGe"
		}
		else if (glm == "binomial") {
			model = "GLMBe"
		}
		else if (glm == "poisson") {
			model = "GLMPe"
		}

		// return the estimate as a structure
		estimate.model = model
		estimate.ntrts = ntrts
		estimate.TrtEff = &TrtEff
		estimate.Ratios = &Ratios

		return(estimate)
	}
	else {
		printf("{err}unknown model\n")
		_error(3001)
	}
}

struct steppes_result scalar stmodelGLM_class::test(
	pointer(class stsubpop_class scalar) scalar sp, real scalar nperm,
	struct steppes_effect scalar effect, | real scalar showstatus,
	real scalar seed, string scalar touse)
{
	struct steppes_result scalar test
	class AssociativeArray scalar Res, PermTest, temp_TrtEff, temp_Ratios,
		temp_PermTest
	class stwin_class scalar win
	real scalar nsubpop, ntrts, j, no, p, terminate, Ntemp, skip, i, overallSglm1,
		overallSglm2, overalllogRatio, s, rm, n_el, has_intercept, cstart,
		OuTcOmE_idx, txassign_idx, seli_idx, h, p1_fit, p1_se, p2_fit, p2_se,
		selj_idx
	real vector coltrt, trts, OuTcOmE, txassign, selj, seli, IndexSet1, IndexSet2,
		sglm1, sglm2, sglmSE1, sglmSE2, slogratios, slogratioSEs, sObs, oObs, logHR,
		ologHR, nv1, nv2, sel_ind
	real matrix osubpop, differences, logratios, Subpop, permuteSubpop, subpop,
		mm1
	string scalar fmla, xnam, todisp, spaces, OuTcOmE_nm, txassign_nm, seli_name,
		model_glm, nv_txt, selj_name
	string vector covariate, nv_nm
	struct ssigma_struct scalar sigmas
	
	if (showstatus == .) showstatus = 1
	if (seed != .) rseed(seed)
	
	if (nperm > 0) {
		win = *sp->win
		nsubpop = sp->nsubpop
		osubpop = sp->subpop
		coltrt = this.coltrt
		OuTcOmE = this.colY
		trts = this.trts
		covariate = this.MM
		glm = this.glm
		link = this.link
		debug = (this.debug != 0)

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

		// create a formula for the model
		OuTcOmE_idx = st_addvar("double", OuTcOmE_nm = st_tempname())
		st_store(selectindex(st_data(., touse)), OuTcOmE_idx, OuTcOmE)
		txassign_idx = st_addvar("double", txassign_nm = st_tempname())
		st_store(selectindex(st_data(., touse)), txassign_idx, txassign)
		if (covariate != "") {
			has_intercept = (covariate[1] == "_cons_")
			if (has_intercept) {
				if  (length(covariate) > 1) {
					cstart = 2
					xnam = (txassign_nm, covariate[cstart::length(covariate)])
					fmla = OuTcOmE_nm + " " + invtokens(xnam)
				}
				else {
					cstart = 1
					xnam = txassign_nm
					fmla = OuTcOmE_nm + " " + invtokens(xnam)
				}
			}
			else {
				cstart = 1
				xnam = (txassign_nm, covariate[cstart::length(covariate)])
				fmla = OuTcOmE_nm + " " + invtokens(xnam)
			}
		}
		else {
			cstart = 2   // the default model has an intercept
			fmla = OuTcOmE_nm + " " + txassign_nm
		}
		
		for (j = 2; j <= ntrts; j++) {
			// do the permutations --> compare trt1 with the rest
			differences = J(nperm, nsubpop, 0)
			logratios = J(nperm, nsubpop, 0)
			no = 0
			p = 0
			terminate = 0
			Ntemp = rows(osubpop)
			IndexSet1 = (1::Ntemp)[selectindex(txassign :== 1)]
			IndexSet2 = (1::Ntemp)[selectindex(txassign :== j)]
					
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
				sglm1 = J(nsubpop, 1, 0)
				sglm2 = J(nsubpop, 1, 0)
				sglmSE1 = J(nsubpop, 1, 0)
				sglmSE2 = J(nsubpop, 1, 0)
				slogratios = J(nsubpop, 1, 0)
				slogratioSEs = J(nsubpop, 1, 0)

				selj = (txassign :== 1 :| txassign :== j)
				
				for (i = 1; i <= nsubpop; i++) {
					seli = selj :& (subpop[., i] :== 1)
					seli_idx = st_addvar("int", seli_name = st_tempname())
					st_store(selectindex(st_data(., touse)), seli_idx, seli)

					if (glm == "gaussian") {
						if (link == "identity") {
							model_glm = "regress " + fmla + " if " + seli_name + ","
						}
					}
					else if (glm == "binomial") {
						if (link == "logit") {
							model_glm = "logit " + fmla + " if " + seli_name + ","
						}
					}
					else if (glm == "poisson") {
						if (link == "logarithm") {
							model_glm = "poisson " + fmla + " if " + seli_name + ","
						}
					}
					else {
						model_glm = "glm " + fmla + " if " + seli_name + ", family(" + glm +
							") link(" + link + ")"
					}
					if (!has_intercept) model_glm = model_glm + " noconstant"
					stata("quietly " + model_glm, 1)
				
					if (covariate != "") {
						if (has_intercept) {
							if (length(covariate) > 1) {
								mm1 = st_data(., covariate[selectindex(covariate :!= "_cons_")],
									seli_name)
								mm1 = (J(rows(mm1), 1, 1), mm1)
							}
							else {
								mm1 = J(sum(seli), 1, 1)
							}
						}
						else {
							mm1 = st_data(., covariate, seli_name)
						}
						mm1 = (txassign[selectindex(seli)], mm1)
						nv1 = mean(mm1)
						nv_nm = txassign_nm
						if (has_intercept) {
							if (cols(mm1) > 2) {
								nv1 = (nv1[1], nv1[3..cols(mm1)])
								nv_nm = (nv_nm, covariate[2..length(covariate)])
							}
							else {
								nv1 = nv1[1]
							}
						}
						else {
							nv_nm = (nv_nm, covariate)
						}
						nv1[1] = 1  	// treatment indicator is 1
						nv_txt = ""
						for (h = 1; h < length(nv1); h++) {
							nv_txt = nv_txt + nv_nm[h] + " = " + strofreal(nv1[h]) + " "
						}
						nv_txt = nv_txt + nv_nm[length(nv1)] + " = " + strofreal(nv1[length(nv1)])
						stata("quietly margins , at(" + nv_txt + ") predict(xb)", 1)
						p1_fit = st_matrix("r(b)")
						p1_se = sqrt(st_matrix("r(V)"))
						nv2 = nv1
						nv2[1] = j	  // treatment indicator is j
						nv_txt = ""
						for (h = 1; h < length(nv2); h++) {
							nv_txt = nv_txt + nv_nm[h] + " = " + strofreal(nv2[h]) + " "
						}
						nv_txt = nv_txt + nv_nm[length(nv2)] + " = " + strofreal(nv2[length(nv2)])
						stata("quietly margins , at(" + nv_txt + ") predict(xb)", 1)
						p2_fit = st_matrix("r(b)")
						p2_se = sqrt(st_matrix("r(V)"))
					} else {
						stata("quietly margins , at(" + txassign_nm + " = 1) predict(xb)", 1)
						p1_fit = st_matrix("r(b)")
						p1_se = sqrt(st_matrix("r(V)"))
						stata("quietly margins , at(" + txassign_nm + " = " + strofreal(j) + ") predict(xb)", 1)
						p2_fit = st_matrix("r(b)")
						p2_se = sqrt(st_matrix("r(V)"))
					}
					st_dropvar(seli_name)

					if (glm == "gaussian") {
						sglm1[i] = p1_fit
						sglmSE1[i] = p1_se
						sglm2[i] = p2_fit
						sglmSE2[i] = p2_se
						slogratios[i] = log(sglm1[i]) - log(sglm2[i])
						slogratioSEs[i] = sqrt((sglmSE1[i]^2)/(sglm1[i]^2) +
							(sglmSE2[i]^2)/(sglm2[i]^2))
					}
					else if (glm == "binomial") {
						sglm1[i] = 1/(1 + exp(-p1_fit))
						sglmSE1[i] = p1_se*((sglm1[i]^2)/exp(p1_fit))   // delta method
						sglm2[i] = 1/(1 + exp(-p2_fit))
						sglmSE2[i] = p2_se*((sglm2[i]^2)/exp(p2_fit))   // delta method
						slogratios[i] = st_matrix("e(b)")[1]
						slogratioSEs[i] = sqrt(st_matrix("e(V)")[1, 1])
					}
					else if (glm == "poisson") {
						sglm1[i] = exp(p1_fit)
						sglmSE1[i] = p1_se*sglm1[i]   // delta method
						sglm2[i] = exp(p2_fit)
						sglmSE2[i] = p2_se*sglm2[i]   // delta method
						slogratios[i] = st_matrix("e(b)")[1]
						slogratioSEs[i] = sqrt(st_matrix("e(V)")[1, 1])
					}
					else {
						printf("{err}unsupported model\n")
						_error(3000)
					}
				}
				
				selj_idx = st_addvar("int", selj_name = st_tempname())
				st_store(selectindex(st_data(., touse)), selj_idx, selj)

				if (glm == "gaussian") {
					if (link == "identity") {
						model_glm = "regress " + fmla + " if " + selj_name + ","
					}
				}
				else if (glm == "binomial") {
					if (link == "logit") {
						model_glm = "logit " + fmla + " if " + selj_name + ","
					}
				}
				else if (glm == "poisson") {
					if (link == "logarithm") {
						model_glm = "poisson " + fmla + " if " + selj_name + ","
					}
				}
				else {
					model_glm = "glm " + fmla + " if " + selj_name + ", family(" + glm +
						") link(" + link + ")"
				}
				if (!has_intercept) model_glm = model_glm + " noconstant"
				stata("quietly " + model_glm, 1)
			
				if (covariate != "") {
					if (has_intercept) {
						if (length(covariate) > 1) {
							mm1 = st_data(., covariate[selectindex(covariate :!= "_cons_")],
								selj_name)
							mm1 = (J(rows(mm1), 1, 1), mm1)
						}
						else {
							mm1 = J(sum(selj), 1, 1)
						}
					}
					else {
						mm1 = st_data(., covariate, selj_name)
					}
					mm1 = (txassign[selectindex(selj)], mm1)
					nv1 = mean(mm1)
					nv_nm = txassign_nm
					if (has_intercept) {
						if (cols(mm1) > 2) {
							nv1 = (nv1[1], nv1[3..cols(mm1)])
							nv_nm = (nv_nm, covariate[2..length(covariate)])
						}
						else {
							nv1 = nv1[1]
						}
					}
					else {
						nv_nm = (nv_nm, covariate)
					}
					nv1[1] = 1  	// treatment indicator is 1
					nv_txt = ""
					for (h = 1; h < length(nv1); h++) {
						nv_txt = nv_txt + nv_nm[h] + " = " + strofreal(nv1[h]) + " "
					}
					nv_txt = nv_txt + nv_nm[length(nv1)] + " = " + strofreal(nv1[length(nv1)])
					stata("quietly margins , at(" + nv_txt + ") predict(xb)", 1)
					p1_fit = st_matrix("r(b)")
					p1_se = sqrt(st_matrix("r(V)"))
					nv2 = nv1
					nv2[1] = j	  // treatment indicator is j
					nv_txt = ""
					for (h = 1; h < length(nv2); h++) {
						nv_txt = nv_txt + nv_nm[h] + " = " + strofreal(nv2[h]) + " "
					}
					nv_txt = nv_txt + nv_nm[length(nv2)] + " = " + strofreal(nv2[length(nv2)])
					stata("quietly margins , at(" + nv_txt + ") predict(xb)", 1)
					p2_fit = st_matrix("r(b)")
					p2_se = sqrt(st_matrix("r(V)"))
				} else {
					stata("quietly margins , at(" + txassign_nm + " = 1) predict(xb)", 1)
					p1_fit = st_matrix("r(b)")
					p1_se = sqrt(st_matrix("r(V)"))
					stata("quietly margins , at(" + txassign_nm + " = " + strofreal(j) + ") predict(xb)", 1)
					p2_fit = st_matrix("r(b)")
					p2_se = sqrt(st_matrix("r(V)"))
				}
				st_dropvar(selj_name)

				if (glm == "gaussian") {
					overallSglm1 = p1_fit
					overallSglm2 = p2_fit
					overalllogRatio = log(overallSglm1[i]) - log(overallSglm2[i])
				}
				else if (glm == "binomial") {
					overallSglm1 = 1/(1 + exp(-p1_fit))
					overallSglm2 = 1/(1 + exp(-p2_fit))
					overalllogRatio = st_matrix("e(b)")[1]
				}
				else if (glm == "poisson") {
					overallSglm1 = exp(p1_fit)
					overallSglm2 = exp(p2_fit)
					overalllogRatio = st_matrix("e(b)")[1]
				}

				if (!hasmissing(sglm1) & !hasmissing(sglm2) &
					!hasmissing(overallSglm1) & !hasmissing(overallSglm2)) {
					no++
					p++
					for (s = 1; s <= nsubpop; s++) {
						differences[p, s] = (sglm1[s] - sglm2[s]) - (overallSglm1 - overallSglm2)
						logratios[p, s] = slogratios[s] - overalllogRatio
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
				n_el = cols(logratios)
				sel_ind = selectindex((1..n_el) :!= rm)
				logratios = logratios[., sel_ind]
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

			sigmas = ssigma(logratios)
			temp_PermTest.put("HRsigma", sigmas.sigma)
			temp_PermTest.put("HRpvalue", ppv2(logratios, logHR, ologHR, nperm))
			Res.put(j - 1, temp_PermTest)
		}
		
		st_dropvar(OuTcOmE_nm)
		st_dropvar(txassign_nm)
		st_rclear()
		st_eclear()
		
		test.model = "GLMt"
		test.ntrts = ntrts
		test.Res = &Res
	}
	
	return(test)
}

void stmodelGLM_class::subgroup(real vector subsample)
{
	this.coltrt  = this.coltrt[selectindex(subsample)]
	this.colY    = this.colY[selectindex(subsample)]
}

void print_estimate_GLM(class steppes_class scalar x, string scalar family,
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
		// print out the GLM results
		if (family == "gaussian") {
			title = "Effect estimates for treatment group " + strofreal(trts[j])
			cnames = ("", "Standard" \ "Effect", "error")
		}
		else if (family == "binomial") {
			title = "Risk estimates for treatment group " + strofreal(trts[j])
			cnames = ("", "Standard" \ "Risk", "error")
		}
		else if (family == "poisson") {
			title = "Effect estimates for treatment group " + strofreal(trts[j])
			cnames = ("", "Standard" \ "Effect", "error")
		}
		temp_TrtEff = TrtEff.get(j)
		temp = J(nsubpop, 3, .)
		temp[., 1] = (1::nsubpop)
		temp[., 2] = temp_TrtEff.get("sObs")
		temp[., 3] = temp_TrtEff.get("sSE")
		rnames = strofreal(temp[., 1])
		fcname = ("" \ "Subpopulation")
		
		if (x.subpop->win->type != "tail-oriented") {
 			temp = (temp \ (., temp_TrtEff.get("oObs"), temp_TrtEff.get("oSE")))
			rnames = (rnames \ "Overall")
		}

		if (x.subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 10, title, "", 
				nsubpop, 0, 4, 0, 0,
				selectindex(temp[, 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 10, title, "",
				(nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
	}

	printf("\n")
	printf("Effect differences and ratio estimates\n")
	for (j = 2; j <= x.effect.ntrts; j++) {
		printf("\n")
		printf("{txt}trt " + strofreal(trts[j]) + " vs. trt " +
			strofreal(trts[1]) + "\n")
		
		if (family == "gaussian") {
			title = "Effect differences"
			cnames = ("Effect", "Standard" \ "difference", "error")
		}
		else if (family == "binomial") {
			title = "Risk differences"
			cnames = ("Risk", "Standard" \ "difference", "error")
		}
		else if (family == "poisson") {
			title = "Effect differences"
			cnames = ("Effect", "Standard" \ "difference", "error")
		}
		if (j > 2) title = ""
		subtitle = ""

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
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 13, title,
				subtitle, nsubpop, 0, 4, 0, 0,
				selectindex(temp[, 1] :== (length(x.subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..3], rnames, cnames, fcname, 15, 13, title,
				subtitle, (nsubpop, nsubpop + 1), 0, 4, 0, 0, 0)
		}
		
		if (family == "gaussian") {
			title = "Effect ratios"
			cnames = ("Log effect", "Standard", "Efect" \ "ratio", "error", "ratio")
		}
		else if (family == "binomial") {
			title = "Odds ratios"
			cnames = ("Log odds", "Standard", "" \ "ratio", "error", "Odds ratio")
		}
		else if (family == "poisson") {
			title = "Relative effect"
			cnames = ("Log", "Standard", "Relative" \ " effect", "error", "effect")
		}
		if (j > 2) title = ""
		subtitle = ""
		
		temp = J(nsubpop, 4, .)
		temp[., 1] = (1::nsubpop)
		temp_Ratios = Ratios.get(j - 1)
		temp[., 2] = temp_Ratios.get("logHR")
		temp[., 3] = temp_Ratios.get("logHRSE")
		temp[., 4] = exp(temp[., 2])
		rnames = strofreal(temp[., 1])
		fcname = ("" \ "Subpopulation")

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

void print_cov_GLM(class steppes_class scalar stobj, real vector trts)
{
	real scalar nsubpop, j, ns
	class AssociativeArray scalar temp_t, t
	string vector names
	struct steppes_result scalar result_empty
	
	nsubpop = stobj.subpop->nsubpop
	names = J(nsubpop, 1, "")
	for (j = 1; j <= nsubpop; j++) {
// 		names[j] = "SP" + strofreal(j) + "-Overall"
		names[j] = "SP" + strofreal(j)
	}

	if (stobj.result != result_empty) {
		temp_t = (*stobj.result.Res)
		for (j = 1; j <= (stobj.result.ntrts - 1); j++) {
			ns = stobj.subpop->nsubpop
			if (stobj.subpop->win->type == "tail-oriented") ns--
			printf("\n")
			printf("{txt}The covariance matrix of the effect ")
			printf("{txt}difference estimates for the {res}%f ", ns)
			printf("{txt}subpopulations is:\n")
			printf("{txt}trt %f vs. trt %f\n", trts[j + 1], trts[1])
			t = temp_t.get(j)
			st_matrix("__tmp_matrix__", t.get("sigma"))
			st_matrixrowstripe("__tmp_matrix__", (J(nsubpop, 1, ""), names))
			st_matrixcolstripe("__tmp_matrix__", (J(nsubpop, 1, ""), names))
			stata("disp_corr, matrix(__tmp_matrix__)", 0)
			stata("matrix drop __tmp_matrix__", 1)

			printf("\n")
			printf("{txt}The covariance matrix of the log effect ratios for the ")
			printf("{res}%f {txt}subpopulations is:\n", ns)
			st_matrix("__tmp_matrix__", t.get("HRsigma"))
			st_matrixrowstripe("__tmp_matrix__", (J(nsubpop, 1, ""), names))
			st_matrixcolstripe("__tmp_matrix__", (J(nsubpop, 1, ""), names))
			stata("disp_corr, matrix(__tmp_matrix__)", 0)
			stata("matrix drop __tmp_matrix__", 1)
		}
	}
}

void print_stat_GLM(class steppes_class scalar stobj, real vector trts)
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
			printf("{txt}Interaction p-value based on effect difference estimates: {res}%f\n", t.get("pvalue"))

			printf("\n")
			printf("{txt}Chi-square test results\n")
			printf("{txt}Interaction p-value based on effect difference estimates: {res}%f\n", t.get("chi2pvalue"))
		}
	}
}

void stmodelGLM_class::print(class steppes_class scalar stobj,
	real scalar estimate, real scalar cov, real scalar test)
{
	// 1. estimates
	if (estimate) {
		print_estimate_GLM(stobj, glm, trts)
	}

	// 2. covariance matrices
	if (cov) {
		print_cov_GLM(stobj, trts)
	}

	// 3. supremum test and chi-square test results
	if (test) {
		print_stat_GLM(stobj, trts)
	}
	
	displayflush()
}

end
