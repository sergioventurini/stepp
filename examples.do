/* Example 1 */
/* --------- */
sysuse simdataKM, clear

stepp time trt, covsubpop(covariate) failure(censor) type(km) patspop(300) ///
	minpatspop(200) trts(1 2) timepoint(4.0) nperm(25) //notest

steppplot, trteff //nopop

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

steppplot, diff conf(15) trtlabs(1 Taxmoxifen 2 Letrozole) ///
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

stepp time trt, type(ci) wintype(sliding_events) mineventspop(10) ///
  eventspop(20) minsubpops(5) trts(1 2) covsubp(ki67) timepoint(4) comprisk(event) notest

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

/* Example 10 (balanced subpopulations - patients) */
/* ----------------------------------------------- */
sysuse balance_example, clear

// mata: covar = st_data(., "covar")
// mata: frq = uniqrows(covar, 1)
// mata: frq = (frq, runningsum(frq[., 2]))
// mata: k = rows(frq)
// mata:
//   allfrq = J(k, k, .)
//   allfrq[1, .] = frq[., 3]'
//   for (i = 2; i <= k; i++) {
//     allfrq[i, i..k] = frq[i..k, 3]' - J(1, k - i + 1, frq[i - 1, 3])
//   }
// end
// mata: r1 = 300
// mata: r2 = 950
// mata: maxnsubpops = 50
// mata: resunb = unbalance(r1, r2, maxnsubpops, frq, allfrq)
// mata: resunb.var
// mata: resunb.nsubpops
// mata: resunb.subpops

balance_patients, range_r1(300 500) range_r2(950 1050) ///
  maxnsubpops(50) covar(covar) //plot

// matrix resmat = e(all_res)
// svmat resmat
// twoway contour resmat3 resmat2 resmat1 if resmat4 == 7, heatmap

/* Example 11 (single group) */
/* ------------------------- */
sysuse simdataKM, clear

stepp time if trt == 2, covsubpop(covariate) failure(censor) type(km) ///
  patspop(300) minpatspop(200) timepoint(4.0)

steppplot, subpop
steppplot, trteff //nopop
steppplot

/* Example 12 (single group) */
/* ------------------------- */
sysuse simdataCI, clear

stepp time if trt == 0, covsubpop(covariate) comprisk(type) type(ci) ///
  patspop(300) minpatspop(200) timepoint(1.0)

steppplot, subpop

/* Example 13 (single group) */
/* ------------------------- */
sysuse aspirin, clear

drop if missing(AD) | missing(AL)

keep if DOSE == 0 | DOSE == 81

generate trtA = cond(DOSE == 81, 1, 0)
generate ADorLE = cond(AD == 1 | AL == 1, 1, 0)

stepp ADorLE if trtA == 0, covsubpop(AGE) type(glm) patspop(100) ///
	minpatspop(30) nperm(10) family(binomial) link(logit)

steppplot, all xtitle("Subpopulations by median age") ytitle(Risk)

/* Example 14 (single group) */
/* ------------------------- */
sysuse bigKM, clear

stepp time, covsubpop(ki67) failure(event) type(km) patspop(150) ///
	minpatspop(50) timepoint(4.0)

steppplot, all xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival")

/* Example 15 (single group) */
/* ------------------------- */
sysuse bigCI, clear

stepp time, covsubpop(ki67) comprisk(event) type(ci) patspop(150) ///
	minpatspop(50) timepoint(4.0)

steppplot, all xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ///
	ytitle("4-year disease free survival") tyscale(-20 50)

/* Example 16 (single group) */
/* ------------------------- */
sysuse simdataKM, clear

quietly replace time = 0 in 3
quietly replace time = 0 in 30
count if time == 0

stepp time, covsubpop(covariate) failure(censor) type(km) patspop(300) ///
	minpatspop(200) timepoint(4.0)

/* Example 17 (single group) */
/* ------------------------- */
use ./data/bcdm, clear

count if (dfs == 0) & (trial == "IBCSG-9")

stepp dfs if trial == "IBCSG-9", covsubpop(ervalue) ///
	failure(dfsevent) type(km) patspop(300) minpatspop(200) ///
	timepoint(20) nperm(25) eps(.0001)

steppplot

/* Example 18 (missing values) */
/* --------------------------- */
sysuse bigKM, clear

replace ki67 = . in 1/80
 
stepp time, covsubpop(ki67) failure(event) type(km) patspop(150) ///
  minpatspop(50) timepoint(4) nperm(50)

steppplot

/* Example 19 (tail-oriented windows) */
/* ---------------------------------- */
sysuse bigKM, clear

local nsub_tmp = 10
// generate_tail, covariate(ki67) nsub(`nsub_tmp') dir("LE")
generate_tail, covariate(ki67) nsub(`nsub_tmp') dir("GE")

stepp time trt, covsubpop(ki67) failure(event) type(km) ///
	trts(1 2) timepoint(4.0) nperm(10) wintype("tail-oriented") ///
	patspop("`r(patspop)'") minpatspop("`r(minpatspop)'")
