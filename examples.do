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
