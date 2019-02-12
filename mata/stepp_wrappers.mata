*!stepp_wrappers version 0.1.0
*!Written 12Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

class stwin_class scalar stwin_wrap(string scalar type, real scalar r1,
	real scalar r2, string scalar basedon)
{
	class stwin_class scalar stwin
	
	stwin.type = type
	stwin.r1 = r1
	stwin.r2 = r2
	stwin.basedon = basedon
	
	if (!stwin.validity()) {
		printf("{err}problems encountered in the generation of the stepp window ")
		printf("{err}object\n")
		_error(3000)
	}
	
	return(stwin)
}

class stsubpop_class scalar stsubpop_wrap(
	pointer(class stwin_class scalar) scalar wp, real vector colvar)
{
	class stsubpop_class scalar subp
	
	subp.win = wp
	subp.colvar = colvar
	subp.generate(&subp)
	
	if (!subp.validity()) {
		printf("{err}problems encountered in the generation of the stepp ")
		printf("{err}subpopulations\n")
		_error(3000)
	}
	
	return(subp)
}

class steppes_class scalar steppes_wrap(
	pointer(class stsubpop_class scalar) scalar sp, string scalar type,
	real vector response, real vector trt, real vector trts, | ///
	real vector failure, real vector comprisk, real scalar timepoint,
	string vector covariates, string scalar family, string scalar link,
	string scalar touse)
{
	class steppes_class scalar steppes
	
	if (type == "KM") {
		class stmodelKM_class scalar stmodelKM
		stmodelKM.survTime = response
		stmodelKM.coltrt = trt
		stmodelKM.censor = failure
		stmodelKM.trts = trts
		stmodelKM.timePoint = timepoint
		steppes.estimate(sp, &stmodelKM, touse)
	}
	else if (type == "CI") {
		class stmodelCI_class scalar stmodelCI
		stmodelCI.coltime = response
		stmodelCI.coltrt = trt
		stmodelCI.coltype = comprisk
		stmodelCI.trts = trts
		stmodelCI.timePoint = timepoint
		steppes.estimate(sp, &stmodelCI, touse)
	}
	else if (type == "GLM") {
		class stmodelGLM_class scalar stmodelGLM
		stmodelGLM.colY = response
		stmodelGLM.coltrt = trt
		stmodelGLM.trts = trts
		stmodelGLM.glm = family
		stmodelGLM.link = link
		stmodelGLM.MM = covariates
		stmodelGLM.debug = 0
		steppes.estimate(sp, &stmodelGLM, touse)
	}
	
	return(steppes)
}

void steppes_test_wrap(
	pointer(class steppes_class scalar) scalar stobj, string scalar type,
	real scalar nperm, real scalar showstatus, real scalar seed, | ///
	string scalar timevar, string scalar failvar, string scalar trtvar,
	real scalar Cox, string vector MM, string scalar touse)
{
	stobj->nperm = nperm
	if (type == "KM") {
		stobj->test(type, stobj->subpop, showstatus, seed, timevar, failvar, trtvar,
			Cox, MM, touse)
	}
	else if ((type == "CI") | (type == "GLM")) {
		stobj->test(type, stobj->subpop, showstatus, seed, "", "", "", ., "", touse)
	}
}

end
