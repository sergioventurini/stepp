*!stepp_steppes version 0.1.0
*!Written 14Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

void steppes_class::new()
{
	class_name = "steppes"
}

string scalar steppes_class::getclass()
{
	return(class_name)
}

void steppes_class::estimate(
	pointer(class stsubpop_class scalar) scalar sp,
	pointer(class stmodel_class scalar) scalar model, string scalar touse)
{
	struct steppes_effect scalar effect
	
	effect = model->estimate(sp, touse)

	this.subpop = sp
	this.model = model
	this.effect = effect
}

void steppes_class::test(string scalar type,
	pointer(class stsubpop_class scalar) scalar sp, | real scalar showstatus,
	real scalar seed, string scalar timevar, string scalar failvar,
	string scalar trtvar, real scalar Cox, string vector MM, string scalar touse)
{
	if ((subpop == NULL) | (this.model == NULL)) {
		printf("{err}warning: you have to estimate the effect first before testing for interaction!\n")
	}
	else {
		struct steppes_result scalar result
	if (type == "CI") {
		pointer(class stmodelCI_class scalar) scalar modelCI
		modelCI = model
		result = model->test(sp, this.nperm, this.effect, showstatus, seed, touse)
	}
	else if (type == "KM") {
		pointer(class stmodelKM_class scalar) scalar modelKM
		modelKM = model
		result = model->test(sp, this.nperm, this.effect, timevar, failvar,
			trtvar, showstatus, Cox, MM, seed, touse)
	}
	else if (type == "GLM") {
		pointer(class stmodelGLM_class scalar) scalar modelGLM
		modelGLM = model
		result = model->test(sp, this.nperm, this.effect, showstatus, seed, touse)
	}

		this.nperm = nperm
		this.result = result
	}
}

void steppes_class::print(string scalar type, real scalar estimate,
	real scalar cov, real scalar test)
{
	real scalar ntrts, n, j, nj
	real vector coltrt, trts
	struct steppes_result scalar result_empty
	
	ntrts = effect.ntrts
	n = 0

	if (type == "CI") {
		pointer(class stmodelCI_class scalar) scalar modelCI
		modelCI = model
		coltrt = modelCI->coltrt
		trts = modelCI->trts
	}
	else if (type == "KM") {
		pointer(class stmodelKM_class scalar) scalar modelKM
		modelKM = model
		coltrt = modelKM->coltrt
		trts = modelKM->trts
	}
	else if (type == "GLM") {
		pointer(class stmodelGLM_class scalar) scalar modelGLM
		modelGLM = model
		coltrt = modelGLM->coltrt
		trts = modelGLM->trts
	}
		
	printf("{txt}\n")
	for (j = 1; j <= ntrts; j++) {
 		nj = sum(coltrt :== trts[j])
		n  = n + nj
		printf("{txt}Sample size in treatment %f: {res}%f\n", trts[j], nj)
	}
	printf("{txt}Total sample size (excluding missing data): {res}%f\n", n)
	model->print(this, estimate, cov, test)

	if (this.result != result_empty) {
		printf("\n")
		printf("{txt}Note: p-values are not adjusted for multiple testing\n")
	}
}

void steppes_class::summary(string scalar type)
{
	real scalar nsubpop, ntrts, j, i
	real vector colvar, minc, maxc, trts, coltrt, txassign, trtj
	real matrix temp
	pointer(class stsubpop_class scalar) scalar subpop
	pointer(class stmodel_class scalar) scalar model
	string scalar title
	string vector fcname, rnames, cnames
	
	this.subpop->summary()

	// print number of patients in each subpopulation for each treatment
	if (this.subpop->init) {
		subpop = this.subpop
		nsubpop = subpop->nsubpop
		colvar = subpop->colvar
		minc = subpop->minc
		maxc = subpop->maxc

		model = this.model
		if (type == "CI") {
			pointer(class stmodelCI_class scalar) scalar modelCI
			modelCI = model
			coltrt = modelCI->coltrt
			trts = modelCI->trts
		}
		else if (type == "KM") {
			pointer(class stmodelKM_class scalar) scalar modelKM
			modelKM = model
			coltrt = modelKM->coltrt
			trts = modelKM->trts
		}
		else if (type == "GLM") {
			pointer(class stmodelGLM_class scalar) scalar modelGLM
			modelGLM = model
			coltrt = modelGLM->coltrt
			trts = modelGLM->trts
		}
		ntrts = length(trts)
		txassign = J(length(coltrt), 1, -1)

		for (j = 1; j <= ntrts; j++) {
			txassign[selectindex(coltrt :== trts[j])] = J(sum(coltrt :== trts[j]), 1, j)
		}

		title = "Treatments sample size information (with only specified treatments):"
		temp = J(nsubpop, ntrts + 2, .)
		temp[., 1] = (1::nsubpop)
		for (i = 1; i <= nsubpop; i++) {
			trtj = J(ntrts, 1, 0)
			for (j = 1; j <= ntrts; j++) {
				trtj[j] = trtj[j] +
					sum((colvar :>= minc[i]) :& (colvar :<= maxc[i]) :& (txassign :== j))
			}
			temp[|i, 2 \ i, (1 + ntrts)|] = trtj'
			temp[i, ntrts + 2] = sum(trtj)
		}
		rnames = strofreal(temp[., 1])
		fcname = ("Subpopulation")
		for (j = 1; j <= ntrts; j++) {
			cnames = (cnames, "trt " + strofreal(j))
		}
		cnames = (cnames, "Total")
		
		if (subpop->win->type == "tail-oriented") {
			printf("\n")
			print_table(temp[., 2..(ntrts + 2)], rnames, cnames, fcname, 15, 10, title, "", 
				nsubpop, 0, 0, 0, 0,
				selectindex(temp[., 1] :== (length(subpop->win->r1) + 1)))
		}
		else {
			printf("\n")
			print_table(temp[., 2..(ntrts + 2)], rnames, cnames, fcname, 15, 10, title, "", 
				nsubpop, 0, 0, 0, 0, 0)
		}
	}
}

void steppes_class::plot()
{

}

end
