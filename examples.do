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

steppplot, all conf(95) ryscale(-1 2.5) //nopop

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
	xtitle("Subpopulations by median age") ytitle(Risk) ///
	ryscale(-2 5) //nopop

/* Example 4 */
/* --------- */
sysuse bigKM, clear

stepp time trt, covsubpop(ki67) failure(event) type(km) patspop(150) ///
	minpatspop(50) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) ///
	xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival") ryscale(-3 7.5) //nopop

/* Example 5 */
/* ---------- */
sysuse bigCI, clear

stepp time trt, covsubpop(ki67) comprisk(event) type(ci) patspop(150) ///
	minpatspop(50) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) ///
	xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival") tyscale(-20 50) ///
	dyscale(-100 40) ryscale(-3 11) //nopop

/* Example 6 */
/* --------- */
sysuse simdataKM, clear

quietly replace time = 0 in 3
quietly replace time = 0 in 30
count if time == 0

stepp time trt, covsubpop(covariate) failure(censor) type(km) patspop(300) ///
	minpatspop(200) trts(1 2) timepoint(4.0) nperm(25) eps(.0001) //notest

/* Example 7 */
/* --------- */
use ./data/bcdm, clear

encode arm, generate(arm_new)

count if (dfs == 0) & (trial == "IBCSG-9")

stepp dfs arm_new if trial == "IBCSG-9", covsubpop(ervalue) ///
	failure(dfsevent) type(km) patspop(300) minpatspop(200) trts(11 12) ///
	timepoint(20) nperm(25) eps(.0001) //notest

steppplot

/* Example 8 (event-based) */
/* ----------------------- */
sysuse bigKM, clear

stepp time trt, type(ci) wintype(sliding_events) mineventspop(10) eventspop(20) minsubpops(5) trts(1 2) covsubp(ki67) timepoint(4) comprisk(event) notest

stepp time trt, covsubpop(ki67) comprisk(event) type(ci) eventspop(20) ///
	mineventspop(10) trts(1 2) timepoint(4) nperm(25) ///
	wintype("sliding_events") minsubpops(5) //noshowsubpops //noshowresults //notest

steppplot, all conf(95) tyscale(-13 29) dyscale(-35 35) ryscale(-1 3.2) //nopop

/* Example 9 (event-based) */
/* ----------------------- */
sysuse bigKM, clear

stepp time trt, covsubpop(ki67) failure(event) type(km) eventspop(20) ///
	mineventspop(10) trts(1 2) timepoint(4) nperm(25) ///
	wintype("sliding_events") minsubpops(5) //notest

steppplot, all conf(95) //nopop
