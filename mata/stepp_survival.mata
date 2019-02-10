*!stepp_survdiff version 0.1.0
*!Written 11Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void survdiff2(real scalar ntot, real scalar ngroup, real scalar nstrat,
	real scalar rho, real vector time, real vector status, real vector group,
	real vector strata, real vector obs, real vector exp, real vector var,
	real vector risk, real vector kaplan)
{
	/* Description:
		 ------------
		 Function that computes quantities related to the G-rho familiy of tests
	*/
	
	/* Source:
		 -------
		 Adapted from the survdiff2() C code in the survival R package
	*/
	
	/* Arguments:
		 ----------
		 INPUTS
		 - ntot				--> real scalar providing the number of observations
		 - ngroup			-->	real scalar providing the number of groups
		 - nstrat			-->	real scalar providing the number of strata
		 - rho				-->	real scalar providing the rho power to use in performing
											the G-test
		 - time				--> real vector of follow up times
		 - status			--> real vector of status indicator (see the Surv() R function
											help)
		 - group			--> real vector providing the group each observation belongs
											to
		 - strata			--> real vector providing the stratum each observation belongs
											to
		 
		 OUTPUTS
		 - obs				--> real vector of weighted observed number of events in each
											group; if there are strata, this will be a matrix with
											one column per stratum
		 - exp				--> real vector of weighted expected number of events in each
											group; if there are strata, this will be a matrix with
											one column per stratum
		 - var				--> real vector of (stacked) variances of the test
		 - risk				--> real vector of times of the cumulative incidence function
		 - kaplan			--> real vector of times of the cumulative incidence function
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/
	
	real scalar i, j, k, kk, n, istart, koff, km, nrisk, wt, tmp, deaths

	istart = 1
	koff = 0
	for (i = 1; i <= (ngroup*ngroup); i++) {
		var[i] = 0
	}
	for (i = 1; i <= (nstrat*ngroup); i++) {
		obs[i] = 0
		exp[i] = 0
	}

	while (istart < ntot) {  /* loop over the strata */
		for (i = 1; i <= ngroup; i++) {
			risk[i] = 0
		}

		/* last obs of this strata */
		for (i = istart; i <= ntot; i++) {
			if (strata[i] == 1) {
				break
			}
		}
		n = i

		/*
		Compute the k-m, which is only needed if rho!=0
		We want it set up as a left-continuous function (unusual)
		*/
		if (rho != 0) {
			km = 1
			for (i = istart; i <= n; ) {
				kaplan[i] = km
				nrisk = n - i + 1
				deaths = status[i]
				for (j = (i + 1); j <= n; j++) {
					if (time[j] == time[i]) {
						kaplan[j] = km
						deaths = deaths + status[j]
					}
					else {
						break
					}
				}
				km = km*(nrisk - deaths)/nrisk
				i = j
			}
		}

		/*
		Now for the actual test
		*/
		for (i = n; i >= istart; i--) {
			if (rho == 0) {
				wt = 1
			}
			else {
				wt = kaplan[i]^rho
			}

			deaths = 0
			for (j = i; j >= istart; j--) {
				if (time[j] == time[i]) {
					k = group[j]
					deaths = deaths + status[j]
					risk[k] = risk[k] + 1
					obs[k + koff] = obs[k + koff] + status[j]*wt
				}
				else {
					break
				}
			}
			i = j + 1
			nrisk = n - i + 1

			if (deaths > 0) {  /* a death time */
				for (k = 1; k <= ngroup; k++) {
					exp[k + koff] = exp[k + koff] + wt*deaths*risk[k]/nrisk
				}

				if (nrisk == 1) {
					continue  /* only 1 subject, so no variance */
				}
				kk = 0
				for (j = 1; j <= ngroup; j++) {
					tmp = (wt^2)*deaths*risk[j]*(nrisk - deaths)/(nrisk*(nrisk - 1))
					var[kk + j] = var[kk + j] + tmp
					for (k = 1; k <= ngroup; k++) {
						var[kk + 1] = var[kk + 1] - tmp*risk[k]/nrisk
						kk++
					}
				}
			}
		}
		istart = n
		koff = koff + ngroup
	}
}

class AssociativeArray scalar survdiff_fit(real vector time, real vector censor,
	real vector x, | real vector strat, real scalar rho)
{
	/* Description:
		 ------------
		 Function that computes quantities related to the G-rho familiy of tests
	*/
	
	/* Source:
		 -------
		 Adapted from the survdiff.fit() R function in the survival R package
	*/
	
	/* Arguments:
		 ----------
		 - time				--> real vector providing the times to event
		 - censor			-->	real vector providing the censoring indicator
		 - x					-->	real vector providing the treatment indicator
		 - strat			--> [optional] real vector providing the stratum each
											observation belongs to
		 - rho				-->	[optional] real scalar providing the rho power to use in
											performing the G-test
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> AssociativeArray scalar containing the observed and
											expected frequencies as well as the corresponding
											variances
	*/
	
	class AssociativeArray scalar res
	real scalar n, ngroup, nstrat
	real vector ord, strat2, observed, expected, var_e, risk, kaplan
	
	res.reinit("string", 1)
	
	n = length(x)
	if ((length(time) != n) | (length(x) != n)) {
		printf("{err}data length mismatch\n")
		_error(3000)
	}

	ngroup = length(uniqrows(x))
	if (ngroup < 2) {
		printf("{err}there is only 1 group\n")
		_error(3000)
	}

	if (strat == J(1, 0, .)) {
		strat = J(n, 1, 1)
	}
	else {
		strat = as_factor(strat)
	}
	nstrat = length(uniqrows(strat))
	if (length(strat) != n) {
		printf("{err}data length mismatch\n")
		_error(3000)
	}
	
	if (rho == .) rho = 0

	ord = order((strat, time, -censor), (1, 2, 3))
	strat2 = (diff(strat[ord]) :!= 0 \ 1)

	observed = J(ngroup*nstrat, 1, 0)
	expected = J(ngroup*nstrat, 1, 0)
	var_e = J(ngroup*ngroup, 1, 0)
	risk = J(ngroup, 1, 0)
	kaplan = J(n, 1, 0)
	survdiff2(n, ngroup, nstrat, rho, time[ord], censor[ord], x[ord], strat2,
		observed, expected, var_e, risk, kaplan)
	
	if (nstrat == 1) {
		res.put("observed", observed)
		res.put("expected", expected)
		res.put("var", colshape(var_e, ngroup)')
	}
	else {
		res.put("observed", colshape(observed, ngroup)')
		res.put("expected", colshape(expected, ngroup)')
		res.put("var", colshape(var_e, ngroup)')
	}
	
	return(res)
}

end
