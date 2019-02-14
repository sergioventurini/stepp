/* -------------------------------------------------------------------------- */
// NOTE: examples from T1 to T5 were used for testing in the very early
//			 development of the package

// cd "/Users/Sergio/Dropbox (Personal)/stepp"

// /* Example T1 */
// /* ---------- */
// // The following example checks if the cinc() function works properly
// use ./data/simdataCI, clear

// sort time

// mata:
// 	y = st_data(., "time")
// 	ic = st_data(., "type")
// 	ic_idx = selectindex(ic :> 0)
// 	ic[ic_idx] = J(length(ic_idx), 1, 1)
// 	icc = st_data(., "trt")
// 	n2 = length(uniqrows(y[selectindex(ic :!= 0)]))
// 	n2 = 2*n2 + 2
// 	x = f = v = J(n2, 1, .)
	
// 	cinc(y, ic, icc, x, f, v)	
// end

// capture confirm new variable x
// if (_rc == 0) {
// 	quietly generate x = .
// }
// else {
// 	quietly replace x = .
// }
// capture confirm new variable f
// if (_rc == 0) {
// 	quietly generate f = .
// }
// else {
// 	quietly replace f = .
// }
// capture confirm new variable v
// if (_rc == 0) {
// 	quietly generate v = .
// }
// else {
// 	quietly replace v = .
// }
// mata:
// 	n = length(y)
// 	st_store(., "x", x[|1\n|])
// 	st_store(., "f", f[|1\n|])
// 	st_store(., "v", v[|1\n|])
// end

// twoway scatter f x, scheme(sj)

// /* Example T2 */
// /* ---------- */
// // The following example checks if the crstm() function works properly
// use ./data/simdataCI, clear

// sort time

// mata:
// 	time = st_data(., "time")
// 	type = st_data(., "type")
// 	cencode = 0
// 	censind = (type :!= cencode)
// 	causeind = (type :== 1)
// 	causeind = 2*censind - causeind
// 	trt = st_data(., "trt")
// 	strata = J(length(time), 1, 1)
// 	rho = 0
// 	ng = length(uniqrows(trt))
// 	ng1 = ng - 1
// 	s = J(ng1, 1, 0)
// 	vs = J(ng1, ng1, 0)
	
// 	crstm(time, causeind, trt, strata, rho, s, vs)
	
// 	s
// 	vs
// end

// /* Example T3 */
// /* ---------- */
// use ./data/simdataCI, clear

// mata:
// 	stwin = stwin_class()
// 	stwin.r1 = 200
// 	stwin.r2 = 300
// 	stwin.validity()

// 	subp = stsubpop_class()
// 	subp.win = &stwin
// 	subp.win->r1
// 	subp.colvar = st_data(., "covariate")
// 	p = &subp
// 	subp.generate(p)

// 	subp.validity()
	
// // 	subp.win->type
// // 	subp.win->r1
// // 	subp.win->r2
// // 	subp.win->basedon
// // 	subp.colvar
// // 	subp.nsubpop
// // 	subp.subpop
// // 	subp.npatsub
// // 	subp.medianz
// // 	subp.minz
// // 	subp.maxz
// // 	subp.minc
// // 	subp.maxc
// // 	subp.init

// 	covariate = st_data(., "covariate")
// 	mergevec = J(8, 1, 0)
// 	mergevec[(1, 3, 6)] = J(3, 1, 1)
// // 	subp_new = subp->merge(p, mergevec)   // doesn't work in interactive mode

// 	class stsubpop_class scalar merge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real vector mergevec)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->merge(sp, mergevec)
		
// 		return(subp_new)
// 	}
	
// 	subp_new = merge_interactive(p, mergevec)
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// 	class stsubpop_class scalar edge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real scalar j,
// 		string scalar side)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->edge(sp, j, side)
		
// 		return(subp_new)
// 	}
	
// //	subp_new = edge_interactive(p, 2, "L")
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// 	subp_new.summary()

// // 	tmp = gen_tow_struct()
// // 	tmp = gen_tailwin(covariate, 8, "LE")
// // 	tmp.v
// // 	tmp.np
// // 	tmp = gen_tailwin(covariate, 8, "RE")
// // 	tmp.v
// // 	tmp.np

// 	stmodelCI_ex = stmodelCI_class()
// 	stmodelCI_ex.coltrt = st_data(., "trt")
// 	stmodelCI_ex.coltime = st_data(., "time")
// 	stmodelCI_ex.coltype = st_data(., "type")
// 	stmodelCI_ex.trts = uniqrows(st_data(., "trt"))
// 	stmodelCI_ex.timePoint = 1.0
	
// 	steppes_ex = steppes_class()
// 	steppes_ex.estimate(&subp, &stmodelCI_ex)
// 	steppes_ex.nperm = 25
//  	steppes_ex.test("CI", &subp, 1, 101)
	
// 	/*
// 	// subpop
// 	// - stwin
// 		 steppes_ex.subpop->win->type
// 		 steppes_ex.subpop->win->r1
// 		 steppes_ex.subpop->win->r2
// 		 steppes_ex.subpop->win->basedon
// 	// - colvar
// 		 steppes_ex.subpop->colvar
// 	// - nsubpop
// 		 steppes_ex.subpop->nsubpop
// 	// - subpop
// 		 steppes_ex.subpop->subpop
// 	// - npatsub
// 		 steppes_ex.subpop->npatsub
// 	// - medianz
// 		 steppes_ex.subpop->medianz
// 	// - minz
// 		 steppes_ex.subpop->minz
// 	// - maxz
// 		 steppes_ex.subpop->maxz
// 	// - minc
// 		 steppes_ex.subpop->minc
// 	// - maxc
// 		 steppes_ex.subpop->maxc
		 
// 	// model
// 	steppes_ex_mod = stmodelCI_class()
// 	steppes_ex_mod = steppes_ex.model 	// needed to recover the correct class of
// 																			// the model item in a steppes object
// 	// steppes_ex_mod->class_name
// 	// - coltrt
// 		 steppes_ex_mod->coltrt
// 	// - coltime
// 		 steppes_ex_mod->coltime
// 	// - coltype
// 		 steppes_ex_mod->coltype
// 	// - trts
// 		 steppes_ex_mod->trts
// 	// - timePoint
// 		 steppes_ex_mod->timePoint

// 	// effect
// 	effect_stepp = &steppes_ex.effect
// 	// - model
// 		 effect_stepp->model
// 	// - ntrts
// 		 effect_stepp->ntrts
// 	// - TrtEff
// 		 TrtEff_1 = effect_stepp->TrtEff->get(1)
// 		 TrtEff_j = effect_stepp->TrtEff->get(2)
// 	//   * sObs
// 			 TrtEff_1.get("sObs")
// 	//   * sSE
// 			 TrtEff_1.get("sSE")
// 	//   * oObs
// 			 TrtEff_j.get("oObs")
// 	//   * oSE
// 			 TrtEff_j.get("oSE")
// 	// - Ratios
// 		 Ratios = effect_stepp->Ratios->get(1)
// 	//   * skmw
// 			 Ratios.get("skmw")
// 	//   * logHR
// 			 Ratios.get("logHR")
// 	//   * logHRSE
// 			 Ratios.get("logHRSE")
// 	//   * ologHR
// 			 Ratios.get("ologHR")
// 	//   * ologHRSE
// 			 Ratios.get("ologHRSE")
// 	//   * logHRw
// 			 Ratios.get("logHRw")
	
// 	// result
// 	result_stepp = &steppes_ex.result
// 	// - model
// 		 result_stepp->model
// 	// - ntrts
// 		 result_stepp->ntrts
// 	// - Res
// 		 Res = result_stepp->Res->get(1)
// 	//   * sigma
// 			 Res.get("sigma")
// 	//   * HRsigma
// 			 Res.get("HRsigma")
// 	//   * pvalue
// 			 Res.get("pvalue")
// 	//   * chi2pvalue
// 			 Res.get("chi2pvalue")
// 	//   * HRpvalue
// 			 Res.get("HRpvalue")

// 	// nperm
// 	steppes_ex.nperm
// 	*/

// 	steppes_ex.print("CI", 1, 1, 1)
// end

// /* Example T4 */
// /* ---------- */
// use ./data/simdataKM, clear

// // stset time, failure(censor)
// // sts generate y = s
// // // sts graph, ci
// // // sts test trt, mat(u V)
// // stcox trt

// // mata:
// // 	time = st_data(., "_t")
// // 	censor = st_data(., "_d")
// // 	trt = st_data(., "trt")

// // 	res = survdiff_fit(time, censor, trt, J(length(time), 1, 1), 0)
// // 	res.get("observed")
// // 	res.get("expected")
// // 	res.get("var")
// // end

// // mata:
// // 	res = kmest(st_data(., "time"), st_data(., "censor"))
// // 	res[1].time
// // 	res[1].est
// // 	res[1].var
// // 	res[1].tint
// // 	res[1].ndd
// // 	trteff = tpest(res, 4)
// // 	trteff.get("est")
// // 	trteff.get("var")
// // end

// mata:
// 	stwin = stwin_class()
// 	stwin.r1 = 200
// 	stwin.r2 = 300
// 	stwin.validity()

// 	subp = stsubpop_class()
// 	subp.win = &stwin
// 	subp.win->r1
// 	subp.colvar = st_data(., "covariate")
// 	p = &subp
// 	subp.generate(p)

// 	subp.validity()
	
// // 	subp.win->type
// // 	subp.win->r1
// // 	subp.win->r2
// // 	subp.win->basedon
// // 	subp.colvar
// // 	subp.nsubpop
// // 	subp.subpop
// // 	subp.npatsub
// // 	subp.medianz
// // 	subp.minz
// // 	subp.maxz
// // 	subp.minc
// // 	subp.maxc
// // 	subp.init

// 	covariate = st_data(., "covariate")
// 	mergevec = J(8, 1, 0)
// 	mergevec[(1, 3, 6)] = J(3, 1, 1)
// // 	subp_new = subp->merge(p, mergevec)   // doesn't work in interactive mode

// 	class stsubpop_class scalar merge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real vector mergevec)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->merge(sp, mergevec)
		
// 		return(subp_new)
// 	}
	
// 	subp_new = merge_interactive(p, mergevec)
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// 	class stsubpop_class scalar edge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real scalar j,
// 		string scalar side)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->edge(sp, j, side)
		
// 		return(subp_new)
// 	}
	
// //	subp_new = edge_interactive(p, 2, "L")
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// // 	subp_new.summary()

// // 	tmp = gen_tow_struct()
// // 	tmp = gen_tailwin(covariate, 8, "LE")
// // 	tmp.v
// // 	tmp.np
// // 	tmp = gen_tailwin(covariate, 8, "RE")
// // 	tmp.v
// // 	tmp.np

// 	stmodelKM_ex = stmodelKM_class()
// 	stmodelKM_ex.coltrt = st_data(., "trt")
// 	stmodelKM_ex.survTime = st_data(., "time")
// 	stmodelKM_ex.censor = st_data(., "censor")
// 	stmodelKM_ex.trts = uniqrows(st_data(., "trt"))
// 	stmodelKM_ex.timePoint = 4.0
	
// 	steppes_ex = steppes_class()
// 	steppes_ex.estimate(&subp, &stmodelKM_ex)
// 	steppes_ex.nperm = 25
//  	steppes_ex.test("KM", &subp, 1, 101, "time", "censor", "trt", 1, J(1, 0, ""))
// //  	steppes_ex.test("KM", &subp, 1, 101, "time", "censor", "trt", 1, "covariate")
	
// // 	// subpop
// // 	// - stwin
// // 		 steppes_ex.subpop->win->type
// // 		 steppes_ex.subpop->win->r1
// // 		 steppes_ex.subpop->win->r2
// // 		 steppes_ex.subpop->win->basedon
// // 	// - colvar
// // 		 steppes_ex.subpop->colvar
// // 	// - nsubpop
// // 		 steppes_ex.subpop->nsubpop
// // 	// - subpop
// // 		 steppes_ex.subpop->subpop
// // 	// - npatsub
// // 		 steppes_ex.subpop->npatsub
// // 	// - medianz
// // 		 steppes_ex.subpop->medianz
// // 	// - minz
// // 		 steppes_ex.subpop->minz
// // 	// - maxz
// // 		 steppes_ex.subpop->maxz
// // 	// - minc
// // 		 steppes_ex.subpop->minc
// // 	// - maxc
// // 		 steppes_ex.subpop->maxc
		 
// // 	// model
// // 	steppes_ex_mod = stmodelKM_class()
// // 	steppes_ex_mod = steppes_ex.model 	// needed to recover the correct class of
// // 																			// the model item in a steppes object
// // 	// steppes_ex_mod->class_name
// // 	// - coltrt
// // 		 steppes_ex_mod->coltrt
// // 	// - coltime
// // 		 steppes_ex_mod->coltime
// // 	// - colcensor
// // 		 steppes_ex_mod->colcensor
// // 	// - trts
// // 		 steppes_ex_mod->trts
// // 	// - timePoint
// // 		 steppes_ex_mod->timePoint

// // 	// effect
// // 	effect_stepp = &steppes_ex.effect
// // 	// - model
// // 		 effect_stepp->model
// // 	// - ntrts
// // 		 effect_stepp->ntrts
// // 	// - TrtEff
// // 		 TrtEff_1 = effect_stepp->TrtEff->get(1)
// // 		 TrtEff_j = effect_stepp->TrtEff->get(2)
// // 	//   * sObs
// // 			 TrtEff_1.get("sObs")
// // 	//   * sSE
// // 			 TrtEff_1.get("sSE")
// // 	//   * oObs
// // 			 TrtEff_j.get("oObs")
// // 	//   * oSE
// // 			 TrtEff_j.get("oSE")
// // 	// - Ratios
// // 		 Ratios = effect_stepp->Ratios->get(1)
// // 	//   * skmw
// // 			 Ratios.get("skmw")
// // 	//   * logHR
// // 			 Ratios.get("logHR")
// // 	//   * logHRSE
// // 			 Ratios.get("logHRSE")
// // 	//   * ologHR
// // 			 Ratios.get("ologHR")
// // 	//   * ologHRSE
// // 			 Ratios.get("ologHRSE")
// // 	//   * logHRw
// // 			 Ratios.get("logHRw")
	
// // 	// result
// // 	result_stepp = &steppes_ex.result
// // 	// - model
// // 		 result_stepp->model
// // 	// - ntrts
// // 		 result_stepp->ntrts
// // 	// - Res
// // 		 Res = result_stepp->Res->get(1)
// // 	//   * sigma
// // 			 Res.get("sigma")
// // 	//   * HRsigma
// // 			 Res.get("HRsigma")
// // 	//   * pvalue
// // 			 Res.get("pvalue")
// // 	//   * chi2pvalue
// // 			 Res.get("chi2pvalue")
// // 	//   * HRpvalue
// // 			 Res.get("HRpvalue")

// // 	// nperm
// // 	steppes_ex.nperm

// 	steppes_ex.print("KM", 1, 1, 1)
// end

// /* Example T5 */
// /* ---------- */
// use ./data/aspirin, clear

// misstable pattern
// drop if missing(AD) | missing(AL)

// keep if DOSE == 0 | DOSE == 81

// generate trtA = cond(DOSE == 81, 1, 0)
// generate ADorLE = cond(AD == 1 | AL == 1, 1, 0)

// mata:
// 	stwin = stwin_class()
// 	stwin.r1 = 30
// 	stwin.r2 = 100
// 	stwin.validity()

// 	subp = stsubpop_class()
// 	subp.win = &stwin
// 	subp.win->r1
// 	subp.colvar = st_data(., "AGE")
// 	p = &subp
// 	subp.generate(p)

// 	subp.validity()
	
// // 	subp.win->type
// // 	subp.win->r1
// // 	subp.win->r2
// // 	subp.win->basedon
// // 	subp.colvar
// // 	subp.nsubpop
// // 	subp.subpop
// // 	subp.npatsub
// // 	subp.medianz
// // 	subp.minz
// // 	subp.maxz
// // 	subp.minc
// // 	subp.maxc
// // 	subp.init

// 	covariate = st_data(., "AGE")
// 	mergevec = J(8, 1, 0)
// 	mergevec[(1, 3, 6)] = J(3, 1, 1)
// // 	subp_new = subp->merge(p, mergevec)   // doesn't work in interactive mode

// 	class stsubpop_class scalar merge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real vector mergevec)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->merge(sp, mergevec)
		
// 		return(subp_new)
// 	}
	
// 	subp_new = merge_interactive(p, mergevec)
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// 	class stsubpop_class scalar edge_interactive(
// 		pointer(class stsubpop_class scalar) scalar sp, real scalar j,
// 		string scalar side)
// 	{
// 		class stsubpop_class scalar subp_new
		
// 		subp_new = sp->edge(sp, j, side)
		
// 		return(subp_new)
// 	}
	
// //	subp_new = edge_interactive(p, 2, "L")
// // 	subp_new.win->type
// // 	subp_new.win->r1
// // 	subp_new.win->r2
// // 	subp_new.win->basedon
// // 	subp_new.colvar
// // 	subp_new.nsubpop
// // 	subp_new.subpop
// // 	subp_new.npatsub
// // 	subp_new.medianz
// // 	subp_new.minz
// // 	subp_new.maxz
// // 	subp_new.minc
// // 	subp_new.maxc
// // 	subp_new.init

// // 	subp_new.summary()

// // 	tmp = gen_tow_struct()
// // 	tmp = gen_tailwin(covariate, 8, "LE")
// // 	tmp.v
// // 	tmp.np
// // 	tmp = gen_tailwin(covariate, 8, "RE")
// // 	tmp.v
// // 	tmp.np

// 	stmodelGLM_ex = stmodelGLM_class()
// 	stmodelGLM_ex.coltrt = st_data(., "trtA")
// 	stmodelGLM_ex.colY = st_data(., "ADorLE")
// 	stmodelGLM_ex.trts = uniqrows(st_data(., "trt"))
// 	stmodelGLM_ex.MM = ""
// // 	stmodelGLM_ex.MM = ("_cons_")
// // 	stmodelGLM_ex.MM = ("_cons_", "G")
// // 	stmodelGLM_ex.MM = "G"
// 	stmodelGLM_ex.glm = "binomial"
// 	stmodelGLM_ex.link = "logit"
// 	stmodelGLM_ex.debug = 0
// 	stmodelGLM_ex.validity()
	
// 	steppes_ex = steppes_class()
// 	steppes_ex.estimate(&subp, &stmodelGLM_ex)
// 	steppes_ex.nperm = 25
//  	steppes_ex.test("GLM", &subp, 1, 101)

// // 	// subpop
// // 	// - stwin
// // 		 steppes_ex.subpop->win->type
// // 		 steppes_ex.subpop->win->r1
// // 		 steppes_ex.subpop->win->r2
// // 		 steppes_ex.subpop->win->basedon
// // 	// - colvar
// // 		 steppes_ex.subpop->colvar
// // 	// - nsubpop
// // 		 steppes_ex.subpop->nsubpop
// // 	// - subpop
// // 		 steppes_ex.subpop->subpop
// // 	// - npatsub
// // 		 steppes_ex.subpop->npatsub
// // 	// - medianz
// // 		 steppes_ex.subpop->medianz
// // 	// - minz
// // 		 steppes_ex.subpop->minz
// // 	// - maxz
// // 		 steppes_ex.subpop->maxz
// // 	// - minc
// // 		 steppes_ex.subpop->minc
// // 	// - maxc
// // 		 steppes_ex.subpop->maxc
		 
// // 	// model
// // 	steppes_ex_mod = stmodelGLM_class()
// // 	steppes_ex_mod = steppes_ex.model 	// needed to recover the correct class of
// // 																			// the model item in a steppes object
// // 	// steppes_ex_mod->class_name
// // 	// - coltrt
// // 		 steppes_ex_mod->coltrt
// // 	// - colY
// // 		 steppes_ex_mod->colY
// // 	// - trts
// // 		 steppes_ex_mod->trts
// // 	// - MM
// // 		 steppes_ex_mod->MM
// // 	// - glm
// // 		 steppes_ex_mod->glm
// // 	// - link
// // 		 steppes_ex_mod->link
// // 	// - debug
// // 		 steppes_ex_mod->debug

// // 	// effect
// // 	effect_stepp = &steppes_ex.effect
// // 	// - model
// // 		 effect_stepp->model
// // 	// - ntrts
// // 		 effect_stepp->ntrts
// // 	// - TrtEff
// // 		 TrtEff_1 = effect_stepp->TrtEff->get(1)
// // 		 TrtEff_j = effect_stepp->TrtEff->get(2)
// // 	//   * sObs
// // 			 TrtEff_1.get("sObs")
// // 	//   * sSE
// // 			 TrtEff_1.get("sSE")
// // 	//   * oObs
// // 			 TrtEff_1.get("oObs")
// // 	//   * oSE
// // 			 TrtEff_1.get("oSE")
// // 	//   * sObs
// // 			 TrtEff_j.get("sObs")
// // 	//   * sSE
// // 			 TrtEff_j.get("sSE")
// // 	//   * oObs
// // 			 TrtEff_j.get("oObs")
// // 	//   * oSE
// // 			 TrtEff_j.get("oSE")
// // 	// - Ratios
// // 		 Ratios = effect_stepp->Ratios->get(1)
// // 	//   * skmw
// // 			 Ratios.get("skmw")
// // 	//   * logHR
// // 			 Ratios.get("logHR")
// // 	//   * logHRSE
// // 			 Ratios.get("logHRSE")
// // 	//   * ologHR
// // 			 Ratios.get("ologHR")
// // 	//   * ologHRSE
// // 			 Ratios.get("ologHRSE")
// // 	//   * logHRw
// // 			 Ratios.get("logHRw")
	
// // 	// result
// // 	result_stepp = &steppes_ex.result
// // 	// - model
// // 		 result_stepp->model
// // 	// - ntrts
// // 		 result_stepp->ntrts
// // 	// - Res
// // 		 Res = result_stepp->Res->get(1)
// // 	//   * sigma
// // 			 Res.get("sigma")
// // 	//   * HRsigma
// // 			 Res.get("HRsigma")
// // 	//   * pvalue
// // 			 Res.get("pvalue")
// // 	//   * chi2pvalue
// // 			 Res.get("chi2pvalue")
// // 	//   * HRpvalue
// // 			 Res.get("HRpvalue")

// // 	// nperm
// // 	steppes_ex.nperm

// 	steppes_ex.print("GLM", 1, 1, 1)
// end
/* -------------------------------------------------------------------------- */


/* Example 1 */
/* --------- */
sysuse simdataKM, clear

stepp time trt, covsubpop(covariate) failure(censor) type(km) patspop(300) ///
	minpatspop(200) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, trteff conf(95) //nopop

/* Example 2 */
/* --------- */
sysuse simdataCI, clear

stepp time trt, covsubpop(covariate) comprisk(type) type(ci) patspop(300) ///
	minpatspop(200) trts(0 1) timepoint(1.0) nperm(25) //notest

steppplot, all conf(95) //nopop

/* Example 3 */
/* --------- */
sysuse aspirin, clear

drop if missing(AD) | missing(AL)

keep if DOSE == 0 | DOSE == 81

generate trtA = cond(DOSE == 81, 1, 0)
generate ADorLE = cond(AD == 1 | AL == 1, 1, 0)

stepp ADorLE trtA, covsubpop(AGE) type(glm) patspop(100) ///
	minpatspop(30) trts(0 1) nperm(10) family(binomial) link(logit)

/*
stepp ADorLE trtA, covsubpop(AGE) type(glm) patspop(100) ///
	minpatspop(30) trts(0 1) nperm(10) covariates(G) family(binomial) ///
	link(logit)
*/

/*
stepp ADorLE trtA, covsubpop(AGE) type(glm) patspop(100) ///
	minpatspop(30) trts(0 1) nperm(10) covariates(G) noconstant ///
	family(binomial) link(logit)
*/

steppplot, all conf(95) trtlabs(1 "Placebo" 2 "81 mg aspirin") ///
	xtitle("Subpopulations by median age") ytitle(Risk) //nopop

/* Example 4 */
/* --------- */
sysuse bigKM, clear

stepp time trt, covsubpop(ki67) failure(event) type(km) patspop(150) ///
	minpatspop(50) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) ///
	xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival") //nopop

/* Example 5 */
/* ---------- */
sysuse bigCI, clear

stepp time trt, covsubpop(ki67) comprisk(event) type(ci) patspop(150) ///
	minpatspop(50) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) ///
	xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival") //nopop
