*!stepp_utils version 0.1.0
*!Written 11Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

real vector recode(real vector V, real vector oldvals, real vector newvals)
{
	/* Description:
		 ------------
		 Function that recodes the elements of a vector using a set of new values
	*/
	
	/* Arguments:
		 ----------
		 - V					--> real vector of values to recode
		 - oldvals		-->	real vector with the values used in V
		 - newvals		-->	real vector with the new values to use for recoding the
											elements of V
	*/
	
	/* Returned value:
		 ---------------
		 - V_new 			--> real vector containing the elements of V recoded using
											elements in newvals
	*/

	real scalar n_oldvals, n_newvals, n, i
	real vector V_vals, V_new, V_new_i
	
	n_oldvals = length(oldvals)
	n_newvals = length(newvals)
	
	if (n_oldvals != n_newvals) {
		printf("{err}the number of old and new values must match\n")
		_error(3000)
	}
	
	n = length(V)
	V_vals = uniqrows(V)
	
	V_new = J(n, 1, .)
	for (i = 1; i <= n_newvals; i++) {
		if (any(V_vals :== oldvals[i])) {
			V_new_i = selectindex(V :== oldvals[i])
			V_new[V_new_i] = J(length(V_new_i), 1, newvals[i])
		}
	}
	
	return(V_new)
}

real vector median(real matrix X)
{
	/* Description:
		 ------------
		 Function that computes the median for each column of a matrix (missing
		 values are skipped)
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix of data
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> real vector containing the columnwise medians of X
	*/

	real scalar n, p, n_nomiss, j, half
	real vector res, Xcol
	
	n = rows(X)
	p = cols(X)
	
	if ((n <= 0) | (p == 0)) {
		printf("{err}the data matrix must contain at least one row and one column\n")
		_error(3000)
	}
	
	res = J(p, 1, .)
	
	for (j = 1; j <= p; j++) {
		Xcol = X[selectindex(X[., j] :!= .), j]
		n_nomiss = length(Xcol)
		if (n_nomiss > 0) {
			half = trunc((n + 1)/2)
			if ((n/2 - trunc(n/2)) > 0) {
        res[j] = sort(Xcol, 1)[half]
			}
			else {
				res[j] = mean(sort(Xcol, 1)[|half \ (half + 1)|])
			}
		}
	}
	
	return(res)
}

struct gen_tow_struct scalar gen_tailwin(real vector covariate,
	real scalar nsub, string scalar dir)
{
	/* Description:
		 ------------
		 Function for generating a vector of covariate values for tail-oriented
		 window with each subgroup roughly about the specified percentage of the
		 total cohort
	*/
	
	/* Arguments:
		 ----------
		 - covariate					--> real vector for the (continuous) covariate of
															interest
		 - nsub								--> real scalar representing the number of
															subpopulation you would like (approximately) in
															addition to the entire cohort
		 - dir								--> string scalar for subpopulations with covariate
															values less than or equal/greater than or equal to
															the generated values (either "GE" or "LE")
	*/
	
	/* Returned value:
		 ---------------
		 - res 								--> struct gen_tow_struct scalar vector containing the
															- v --> vector of covariate values to be used in a
																steppwin
															- np --> vector of subpopulation size associated
																with each tail-oriented window defined by v
	*/

	struct gen_tow_struct scalar res
	real scalar perc, remain, lcov, grpsize, i, cov_val, cov_r, cov_l
	real vector cov_t, v, np, pro_val, mat, temp
	
  if (nsub < 1) {
		printf("{err}number of subgroups must be larger than 1\n")
	}
	perc = 1/(nsub + 1)

	cov_t = sort(covariate, 1)
	remain = length(cov_t)
	lcov = remain
	grpsize = round(remain*perc)
	v = J(remain, 1, 0)
	np = J(remain, 1, 0)

	i = 0
	if (dir == "LE") {
		while (remain >= 2*grpsize) {
			cov_val = cov_t[remain - grpsize]
			// handle the case where the covariate value is the max
			if (cov_val == max(cov_t)) {
				pro_val = selectindex(cov_t :< max(cov_t))
				cov_val = cov_t[pro_val[length(pro_val)]]
			}
			mat = selectindex(cov_val :== cov_t)
			cov_r = mat[length(mat)]
			v[lcov - i] = cov_val
			cov_t = cov_t[1..cov_r]
			remain = length(cov_t)
			np[lcov - i] = remain
			i++
		}
		v = v[(lcov - i + 1)..lcov]
		np = np[(lcov - i + 1)..lcov]
	}
	else {
		while (remain >= 2*grpsize) {
			cov_val = cov_t[grpsize + 1]
			// handle the case where the covariate value is the min.
			if (cov_val == min(cov_t)) {
				pro_val = selectindex(cov_t :> min(cov_t))
				cov_val = cov_t[pro_val[1]]
			}

			temp = selectindex(cov_val :== cov_t)
			cov_l = temp[1]
			v[i + 1] = cov_val
			cov_t = cov_t[(cov_l)..length(cov_t)]
			remain = length(cov_t)
			np[i + 1] = length(cov_t)
			i++
		}
		v = v[1..i]
		np = np[1..i]
	}
	
	res.v = v
	res.np = np
	
  return(res)
}

string scalar format(real scalar val, real scalar width, real scalar last,
	| real scalar buffer)
{
	/* Description:
		 ------------
		 Function for formatting a number to fit a given width
	*/
	
	/* Arguments:
		 ----------
		 - val						--> real scalar representing the value to format
		 - width					--> real scalar representing the width to fill
		 - last						--> real scalar indicating if the value is to be displayed
													in the last column of a table
		 - buffer					--> (optional) real scalar representing the space to leave
													
	*/
	
	/* Returned value:
		 ---------------
		 - todisp					--> string scalar containing the string to display
	*/

	string scalar todisp, trail
	real scalar len, maxlen
	
	if (args() == 3) {
		buffer = 2
	}
	
	todisp = strofreal(val)
	len = strlen(todisp)
	if (last) {
		maxlen = width - buffer
	}
	else {
		maxlen = width - 2*buffer
	}
	if (len > maxlen) {
		printf("{err}the value to display is longer than the space available\n")
		_error(3000)
	}
	trail = (maxlen - len)*char(32)
	todisp = trail + todisp
	
	return(todisp)
}

real vector as_factor(transmorphic vector x)
{
	/* Description:
		 ------------
		 Function that converts values of a vector
	*/
	
	/* Arguments:
		 ----------
		 - x					--> transmorphic colvector whose values need to be changed
	*/
	
	/* Returned value:
		 ---------------
		 - f		 			--> real colvector where the original values have been changed
	*/

	real scalar i
	real vector f, levels, sel
	
	levels = sort(uniqrows(x), 1)	
	f = J(length(x), 1, .)
	for (i = 1; i <= length(levels); i++) {
		sel = selectindex(x :== levels[i])
		f[sel] = J(length(sel), 1, i)
	}
	
	return(f)
}

struct ssigma_struct scalar ssigma(real matrix imatrix)
{
	/* Description:
		 ------------
		 Function that computes a covariance matrix and its inverse
	*/
	
	/* Arguments:
		 ----------
		 - imatrix			--> real matrix containing the data
	*/
	
	/* Returned value:
		 ---------------
		 - res		 			--> struct ssigma_struct scalar containing the covariance
												matrix and its inverse
	*/

	struct ssigma_struct scalar res
	real matrix sigma, sigmainv
	
	sigma = variance(imatrix)
	sigmainv = invsym(sigma)
	
	res.sigma = sigma
	res.sigmainv = sigmainv
	
	return(res)
}

real scalar ppv(real matrix imatrix, real matrix sigmainv, real vector estarray,
	real vector est, real scalar noPerms)
{
	/* Description:
		 ------------
		 Function that computes the pvalue for a permutation test
	*/
	
	/* Arguments:
		 ----------
		 - imatrix			--> real matrix containing the data
		 - sigmainv			--> real matrix containing the inverse of the covariance
												matrix
		 - estarray			--> real vector containing the oObs values
		 - est					--> real vector containing the sObs values
		 - noPerms			--> real scalar containing the number of permutations
	*/
	
	/* Returned value:
		 ---------------
		 - pvalue		 		--> real scalar containing the permutation test pvalue
	*/

	real scalar pvalue, i
	real vector perm, obs1
	real matrix temp, obs
	
	perm = J(noPerms, 1, .)
	
	for (i = 1; i <= noPerms; i++) {
		temp = imatrix[i, .]'
		perm[i] = temp' * sigmainv * temp
	}
	obs = vec(estarray :- est)
	obs1 = obs' * sigmainv * obs
	pvalue = sum(perm :> obs1)/noPerms
	
	return(pvalue)
}

real scalar ppv2(real matrix imatrix, real vector estarray, real vector est,
	real scalar noPerms)
{
	/* Description:
		 ------------
		 Function that computes the supremum statistics and generate the pvalue
	*/
	
	/* Arguments:
		 ----------
		 - imatrix			--> real matrix containing the data
		 - estarray			--> real vector containing the oObs values
		 - est					--> real vector containing the sObs values
		 - noPerms			--> real scalar containing the number of permutations
	*/
	
	/* Returned value:
		 ---------------
		 - pvalue		 		--> real scalar containing the permutation test pvalue
	*/

	real scalar pvalue, i
	real colvector sigma, tPerm, tObs
	real matrix stdDifferences, obsDifferences, stdObsDifferences
	
	stdDifferences = J(rows(imatrix), cols(imatrix), .)
	
	sigma = sqrt(diagonal(variance(imatrix)))
	for (i = 1; i <= rows(imatrix); i++) {
		stdDifferences[i, .] = imatrix[i, .] :/ sigma'
	}
	tPerm = rowmax(abs(stdDifferences))
	
	obsDifferences = vec(estarray :- est)'
	stdObsDifferences = J(rows(obsDifferences), cols(obsDifferences), .)
	for (i = 1; i <= rows(obsDifferences); i++) {
		stdObsDifferences[i, .] = obsDifferences[i, .] :/ sigma'
	}
	tObs = rowmax(abs(stdObsDifferences))
	pvalue = sum(tPerm :> tObs)/noPerms
	
	return(pvalue)
}

void print_matrix(real matrix imatrix, string vector rownames,
	string vector colnames, | string scalar title, real scalar cutoff,
	real scalar longest)
{
	/* Description:
		 ------------
		 Function that prints the content of a (symmetric) matrix
	*/
	
	/* Arguments:
		 ----------
		 - imatrix		--> real matrix containing containing the numbers to display
		 - rownames		--> string vector containing the table's row lables
		 - colnames		--> string vector containing the table's column lables
		 - title			--> [optional] string scalar containing the table's main title
		 - cutoff			--> [optional] real scalar for not showing values smaller than
											cutoff
		 - longest		--> [optional] real scalar indicating the longest string to
											print
	*/
	
	/* Returned value:
		 ---------------
		 - void		 		--> no output returned by this function
	*/
	
	real scalar skip0, nvar, j0, j1, j, l, i
	real matrix cj
	
	if (cutoff == .) cutoff = 0
	if (longest == .) longest = 8
	if (longest > 11) {
		printf("{txt}[Warning: labels too long to be fully displayed]\n")
		return
	}
	
	skip0 = 0
	
	if (title != "") {
		printf("\n")
		printf("{txt}{space " + strofreal(skip0) + "}" + title + "\n")
	}

	nvar = length(rownames)
	j0 = 1
	cj = J(nvar, nvar, .)
	
	while (j0 <= nvar) {
		j1 = min((j0 + (longest + 1), nvar))
		j = j0
		l = (longest + 1)*(j1 - j0 + 1)
		printf("{txt}{hline " + strofreal(longest + 2) + "}{c TT}{hline " + strofreal(l) + "}\n")
		printf("{txt}{space " + strofreal(longest + 2) + "}{c |}")
		while (j <= j1) {
			printf("{txt}%" + strofreal(longest + 1) + "s", abbrev(colnames[j], longest))
			j++
		}
		printf("{txt}\n{hline " + strofreal(longest + 2) + "}{c +}{hline " + strofreal(l) + "}\n")

		i = j0
		while (i <= nvar) {
			printf("{txt}%" + strofreal(longest + 1) + "s", abbrev(rownames[i], longest))
			printf(" {c |} ")
			j = j0
			while (j <= min((j1, i))) {
				if (!hasmissing(imatrix[i, j]) & (abs(imatrix[i, j]) < cutoff)) {
					cj[i, j] = .
				}
				else {
					cj[i, j] = imatrix[i, j]
				}
				j++
			}
			j = j0
			while (j <= min((j1, i))) {
				if (cj[i, j] < .) {
					printf("{res} %" + strofreal(longest - 1) + ".4f ", cj[i, j])
				}
				else {
					printf("{res}{space " + strofreal(longest + 1) + "}")
				}
				j++
			}
			printf("\n")
			i++
		}
		j0 = j0 + 10
		printf("{txt}{hline " + strofreal(longest + 2) + "}{c BT}{hline " + strofreal(l) + "}\n")
	}
}

void print_table(real matrix omatrix, string vector matrownames,
	string matrix matcolnames, string vector firstcolname,
	real scalar firstcolwidth, real scalar colwidth, string scalar title,
	string scalar subtitle, real vector hlines, real scalar novlines,
	real scalar digits, real scalar corr, real scalar cutoff, real scalar overall)
{
	/* Description:
		 ------------
		 Function that prints the content of a table/matrix
	*/
	
	/* Arguments:
		 ----------
		 - omatrix						--> matrix containing the numbers to display
		 - matrownames				--> string vector containing the table's row lables
		 - matcolnames				--> string matrix containing the table's column lables
		 - firstcolname				--> string vector - first column's name
		 - firstcolwidth			--> real scalar - first column's width
		 - colwidth						--> real scalar - other columns' width
		 - title							--> string scalar - table's main title
		 - subtitle						--> string scalar - table's subtitle
		 - hlines							--> real vector - rows at which to insert an
															horizontal line
		 - novlines						--> real scalar - if 1, no vertical lines between
															columns
		 - digits							--> real scalar - number of digits to display
		 - corr								--> real scalar - correlation table
		 - cutoff							--> real scalar - if 1, do not show values smaller
															than cutoff
		 - overall						--> real scalar - overall label
	*/
	
	/* Returned value:
		 ---------------
		 - void		 		--> no output returned by this function
	*/
	
	real scalar skip0, skip1, skip2, nr_first, firstcolwidth_p1, ncols, nrows,
		ncols_m1, usable, usable2, tmp_skip, i, j, num
	real vector hlinestodisp
	real matrix imatrix
	string scalar tvlines, bvlines, mvlines, vlines, firstline, todisp,
		secondline, lastline, rownametodisp, rowtodisp
	string vector title1
	
	skip0 = 0
	skip1 = 1
	skip2 = 2
	nr_first = length(firstcolname)
	title1 = J(nr_first, 1, "")

	for (i = 1; i <= nr_first; i++) {
		if (strlen(firstcolname[i]) > firstcolwidth - 2*skip1) {
			 firstcolname[i] = abbrev(firstcolname[i], firstcolwidth - 2*skip1)
		}
	}
	if (colwidth < 9) {
		printf("{err}'colwidth' must be larger than 8 to properly show the table\n")
		exit()
	}
	
	imatrix = omatrix

	firstcolwidth_p1 = firstcolwidth + 1
	ncols = cols(imatrix)
	nrows = rows(imatrix)
	ncols_m1 = ncols - 1
	usable = colwidth - 2*skip1
	usable2 = colwidth - 2*skip2
	if (digits >= usable) {
		printf("{err}the number of digits chosen is too large\n")
		exit()
	}
	if (hlines != .) {
		hlinestodisp = sort(hlines, 1)
	}
	if (!novlines) {
		tvlines = "{c TT}"
		bvlines = "{c BT}"
		mvlines = "{c +}"
		vlines = "{space " + strofreal(skip1) + "}{c |}"
	}
	else {
		tvlines = "{c -}"
		bvlines = "{c -}"
		mvlines = "{c -}"
		vlines = "{space " + strofreal(skip2) + "}"
	}
	firstline = "{hline " + strofreal(firstcolwidth) + "}{c TT}"
	for (i = 1; i <= nr_first; i++) {
		todisp = firstcolname[i]
		if (strlen(todisp) > firstcolwidth - 2*skip1) {
			todisp = abbrev(todisp, firstcolwidth - 2*skip1)
		}
		tmp_skip = firstcolwidth - strlen(todisp) - 2*skip1
		title1[i] = "{space " + strofreal(skip1 + tmp_skip) + "}" + todisp +
			"{col " + strofreal(firstcolwidth_p1) + "}{c |}"
	}
	secondline = "{hline " + strofreal(firstcolwidth) + "}{c +}"
	lastline = "{hline " + strofreal(firstcolwidth) + "}{c BT}"
	
	if (corr) {
		for (i = 1; i <= nrows; i++) {
			for (j = 1; j <= ncols; j++) {
				if ((imatrix[i, j] != .) & (abs(imatrix[i, j]) < cutoff)) {
					imatrix[i, j] = .
				}
			}
		}
	}
	
	for (j = 1; j <= ncols_m1; j++) {
		firstline = firstline + "{hline " + strofreal(colwidth) + "}" + tvlines
		for (i = 1; i <= nr_first; i++) {
			todisp = matcolnames[i, j]
			if (strlen(todisp) > usable) {
				todisp = abbrev(todisp, usable)
			}
			tmp_skip = usable - strlen(todisp)
			title1[i] = title1[i] + "{space " + strofreal(skip1 + tmp_skip) + "}" +
				abbrev(todisp, colwidth) + vlines
		}
		secondline = secondline + "{hline " + strofreal(colwidth) + "}" + mvlines
		lastline = lastline + "{hline " + strofreal(colwidth) + "}" + bvlines
	}
	firstline = firstline + "{hline " + strofreal(colwidth) + "}"
	for (i = 1; i <= nr_first; i++) {
		todisp = matcolnames[i, ncols]
		if (strlen(todisp) >= usable) {
			todisp = abbrev(todisp, usable)
		}
		tmp_skip = usable - strlen(todisp)
		title1[i] = title1[i] + "{space " + strofreal(skip1 + tmp_skip) + "}" +
			abbrev(todisp, colwidth) + " "
	}
	secondline = secondline + "{hline " + strofreal(colwidth) + "}"
	lastline = lastline + "{hline " + strofreal(colwidth) + "}"
	
	if (title != "") {
		printf("{txt}{space " + strofreal(skip0) + "}" + title + "\n")
	}
	if (subtitle != "") {
		printf("{txt}{space " + strofreal(skip0) + "}" + subtitle + "\n")
	}
	printf("{txt}" + firstline + "\n")
	for (i = 1; i <= nr_first; i++) {
		printf("{txt}" + title1[i] + "\n")
	}
	printf("{txt}" + secondline + "\n")
	
	for (i = 1; i <= nrows; i++) {
		rownametodisp = matrownames[i]
		if (rownametodisp != ".") {
			if (strlen(rownametodisp) > (firstcolwidth - 2*skip1)) {
				rownametodisp = abbrev(rownametodisp, firstcolwidth - 2*skip1)
			}
			tmp_skip = firstcolwidth - strlen(rownametodisp) - 2*skip1
			rownametodisp = "{space " + strofreal(skip1 + tmp_skip) + "}" +
				rownametodisp + "{space " + strofreal(skip1) + "}{c |}"
		}
		else {
			rownametodisp = "{col " + strofreal(firstcolwidth_p1) + "}{c |}"
		}
		printf("{txt}" + rownametodisp)

		rownametodisp = matrownames[i]
		rowtodisp = ""
		for (j = 1; j <= ncols_m1; j++) {
			if (imatrix[i, j] != .) {
				if (rownametodisp != ".") {
					todisp = strtrim(strofreal(imatrix[i, j], "%" + strofreal(usable) + 
						"." + strofreal(digits) + "f"))
					if (strlen(todisp) > usable) {
						todisp = strtrim(strofreal(imatrix[i, j], "%" + strofreal(usable) +
							".2e"))
					}
					tmp_skip = usable - strlen(todisp)
					rowtodisp = rowtodisp + "{space " + strofreal(skip1 + tmp_skip) +
						"}" + "{res}" + todisp + "{txt}" + vlines
				}
				else {
					todisp = strtrim(strofreal(imatrix[i, j], "%" + strofreal(usable2) +
						"." + strofreal(digits) + "f"))
					if (strlen(todisp) > usable2) {
						todisp = strtrim(strofreal(imatrix[i, j], "%" + strofreal(usable2) +
						".2e"))
					}
					todisp = "(todisp)"
					tmp_skip = usable - strlen(todisp)
					rowtodisp = rowtodisp + "{space " + strofreal(skip1 + tmp_skip) +
						"}" + "{res}" + todisp + "{txt}" + vlines
				}
			}
			else {
				rowtodisp = rowtodisp + "{txt}{space " + strofreal(skip1 + usable) +
					"}" + vlines
			}
		}
		if (imatrix[i, ncols] != .) {
			if (rownametodisp != ".") {
				todisp = strtrim(strofreal(imatrix[i, ncols], "%" + strofreal(usable) +
					"." + strofreal(digits) + "f"))
				if (strlen(todisp) > usable) {
					todisp = strtrim(strofreal(imatrix[i, ncols], "%" +
						strofreal(usable) + ".2e"))
				}
				tmp_skip = usable - strlen(todisp)
				if (overall == i) {
					rowtodisp = rowtodisp + "{space " + strofreal(skip1 + tmp_skip) +
						"}" + "{res}" + todisp + " (entire cohort)"
				}
				else {
					rowtodisp = rowtodisp + "{space " + strofreal(skip1 + tmp_skip) +
						"}" + "{res}" + todisp
				}
			}
			else {
				todisp = strtrim(strofreal(imatrix[i, ncols], "%" + strofreal(usable2) +
					"." + strofreal(digits) + "f"))
				if (strlen(todisp) > usable2) {
					todisp = strtrim(strofreal(imatrix[i, ncols], "%" +
					strofreal(usable2) + ".2e"))
				}
				todisp = "(todisp)"
				tmp_skip = usable - strlen(todisp)
				rowtodisp = rowtodisp + "{space " + strofreal(skip1 + tmp_skip) +
						"}{res}" + todisp
			}
		}
		printf(rowtodisp + "\n")
		
		for (num = 1; num <= length(hlinestodisp); num++) {
			if (hlinestodisp[num] == i) {
				if (hlinestodisp[num] < nrows) {
					printf("{txt}" + secondline + "\n")
				}
				else if (hlinestodisp[num] == nrows) {
					printf("{txt}" + lastline + "\n")
				}
				break
			}
		}
	}
}

real vector diff(real vector V, | real scalar lag)
{
	/* Description:
		 ------------
		 Function that computes differences of any lag
	*/
	
	/* Arguments:
		 ----------
		 - V					--> real vector of values to differentiate
		 - lag				-->	real scalar repreesnting the lag to use
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> real vector containing the differences of the elements in
											V
	*/

	real scalar n
	real vector res
	
	if (lag == .) lag = 1
	
	if (lag <= 0) {
		printf("{err}the lag must be greater than or equal to 1\n")
		_error(3000)
	}
	lag = trunc(lag)
	
	n = length(V)
	res = J(n - lag, 1, .)
	
	res = V[((1 + lag)..n)] - V[(1..(n - lag))]
	
	return(res)
}

void return_results(pointer(class steppes_class scalar) scalar stobj,
	string scalar notest)
{
	real scalar ntreff, nratios, nres, j
	class AssociativeArray scalar tmp

	ntreff = (*stobj->effect.TrtEff).N()
	
	for (j = 1; j <= ntreff; j++) {
		tmp = (*stobj->effect.TrtEff).get(j)
		st_matrix("e(sObs_" + strofreal(j) + ")", tmp.get("sObs"))
		st_matrix("e(sSE_" + strofreal(j) + ")", tmp.get("sSE"))
		st_numscalar("e(oObs_" + strofreal(j) + ")", tmp.get("oObs"))
		st_numscalar("e(oSE_" + strofreal(j) + ")", tmp.get("oSE"))
	}
	
	nratios = (*stobj->effect.Ratios).N()
	
	for (j = 1; j <= nratios; j++) {
		tmp = (*stobj->effect.Ratios).get(j)
		st_matrix("e(logHR_" + strofreal(j) + ")", tmp.get("logHR"))
		st_matrix("e(logHRSE_" + strofreal(j) + ")", tmp.get("logHRSE"))
		st_numscalar("e(skmw_" + strofreal(j) + ")", tmp.get("skmw"))
		st_numscalar("e(ologHR_" + strofreal(j) + ")", tmp.get("ologHR"))
		st_numscalar("e(ologHRSE_" + strofreal(j) + ")", tmp.get("ologHRSE"))
		st_numscalar("e(logHRw_" + strofreal(j) + ")", tmp.get("logHRw"))
	}
	
	if (notest == "") {
		nres = (*stobj->result.Res).N()		
		for (j = 1; j <= nres; j++) {
			tmp = (*stobj->result.Res).get(j)
			st_matrix("e(sigma_" + strofreal(j) + ")", tmp.get("sigma"))
			st_matrix("e(HRsigma_" + strofreal(j) + ")", tmp.get("HRsigma"))
			st_numscalar("e(pvalue_" + strofreal(j) + ")", tmp.get("pvalue"))
			st_numscalar("e(chi2pvalue_" + strofreal(j) + ")", tmp.get("chi2pvalue"))
			st_numscalar("e(HRpvalue_" + strofreal(j) + ")", tmp.get("HRpvalue"))
		}
	}
}

real vector pretty(real vector x, real scalar n)
{
	/* Description:
		 ------------
		 Function that computes a sequence of about n + 1 equally spaced 'round'
		 values which cover the range of the values in x. The values are chosen so
		 that they are 1, 2 or 5 times a power of 10
	*/
	
	/* Arguments:
		 ----------
		 - res 				--> real vector containing the data to use
		 - n					--> real scalar providing the desired number of intervals
	*/
	
	
	/* Notes:
		 ----------
		 Code is partially taken from the pretty() R function
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> real vector containing the sequence of values
	*/

	real scalar min_n, shrink_sml, high_u_bias, u5_bias, eps_correct, l, u
	real vector res
	
	min_n = trunc(n/3)
	shrink_sml = 0.75
	high_u_bias = 1.5
	u5_bias = .5 + 1.5*high_u_bias
	eps_correct = 0
	
	x = x[selectindex(x :< .)]
	if (!length(x)) {
		printf("{err}the vector provided contains only missing values\n")
		_error(3000)
	}
	
	l = min(x)
	u = max(x)
	
	if (n <= 0) {
		printf("{err}the number of intervals must be strictly positive\n")
		_error(3000)
	}
	if (min_n <= 0) {
		printf("{err}the minimum number of intervals must be strictly positive\n")
		_error(3000)
	}

	// next code is taken from R_pretty() C function defined in R source code
	real scalar rounding_eps, h, h5, dx, cell, unit, base, U, ns, nu, k, i_small,
		ndiv, return_bounds

	rounding_eps = 1e-10
	h = high_u_bias
	h5 = u5_bias
	return_bounds = 1

	dx = u - l
	if (dx == 0 & u == 0) {
		cell = 1
		i_small = 1
	}
	else {
		cell = max((abs(l), abs(u)))
		U = 1 + ((h5 >= 1.5*h+.5) ? 1/(1+h) : 1.5/(1+h5))
		U = U*max((1, n))*epsilon(1)
		i_small = dx < cell*U*3
	}

	if (i_small) {
		if (cell > 10) cell = 9 + cell/10
		cell = cell*shrink_sml
		if (min_n > 1) cell = cell/min_n
	} else {
		cell = dx
		if (n > 1) cell = cell/n
	}

	if (cell < 20*smallestdouble()) {
		printf("{txt}internal pretty(): very small range -> corrected!\n")
		cell = 20*smallestdouble()
	} else if (cell * 10 > maxdouble()) {
		printf("{txt}internal pretty(): very large range -> corrected!\n")
		cell = .1*maxdouble()
	}
	base = 10^floor(log10(cell))

	unit = base
	if ((U = 2*base) - cell <  h*(cell - unit)) {
		unit = U
		if ((U = 5*base) - cell < h5*(cell - unit)) {
			unit = U
			if ((U = 10*base) - cell < h*(cell - unit)) {
				unit = U
			}
		}
	}

	ns = floor(l/unit + rounding_eps)
	nu = ceil (u/unit - rounding_eps)
	if (eps_correct & (eps_correct > 1 | !i_small)) {
		if (l != 0) {
			l = l*(1 - epsilon(1))
		}
		else {
			l = mindouble()
		}
		if (u != 0) {
			u = u*(1 + epsilon(1))
		}
		else {
			u = smallestdouble()
		}
	}

	while (ns*unit > l + rounding_eps*unit) ns--

	while (nu*unit < u - rounding_eps*unit) nu++

	k = trunc(0.5 + nu - ns)
	if (k < min_n) {
		k = min_n - k
		if (ns >= 0) {
			nu = nu + k/2
			ns = ns - k/2 + mod(k, 2)
		}
		else {
			ns = ns - k/2
			nu = nu + k/2 + mod(k, 2)
		}
		ndiv = min_n
	}
	else {
		ndiv = k
	}
	if (return_bounds) {
		if (ns * unit < l) l = ns * unit
		if (nu * unit > u) u = nu * unit
	} else {
		l = ns
		u = nu
	}

	res = runningsum((l \ J(ndiv, 1, round((ceil(u) - floor(l))/ndiv, .5))))

	return(res)
}

void cleanup()
{
	/* Description:
		 ------------
		 Function that cleans up Mata from temporary objects 
	*/
	
	/* Arguments:
		 ----------
		 - void				--> no argument passed
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	string colvector names
	real scalar i
	
	names = direxternal("__*")

	for (i = 1; i <= rows(names); i++) {
		// delete all objectes apart from __steppes__
		if (names[i] != "__steppes__") rmexternal(names[i])
	}
}

end
