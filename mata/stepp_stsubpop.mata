*!stepp_stsubpop version 0.1.0
*!Written 12Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void stsubpop_class::new()
{
	class_name = "stsubpop"
	init = 0
}

string scalar stsubpop_class::getclass()
{
	return(class_name)
}

real scalar stsubpop_class::validity()
{
	real scalar valid

	valid = 1

	// this check is redundant
	if (eltype(colvar) != "real") {
		printf("the provided colvar is not numeric\n")
		valid = 0
// 		_error(1)
	}
	if (win->type == "sliding") {
		if (length(colvar) < win->r2) {
			printf("length of colvar must be at least equal to r2\n")
			valid = 0
// 			_error(1)
		}
	}
	
	return(valid)
}

void stsubpop_class::generate(pointer(class stsubpop_class scalar) scalar sp)
{
	class stwin_class scalar lwin   // local copy of the stwin_class instance that
																	// sp.win points to
	real scalar r1, r2, nsubpop, i, sel, k
	string scalar type	
	real vector covariate, zvals, npatsub, I0, I1
	
	lwin = *sp->win   // entire win object from sp
	
	r1 = lwin.r1
	r2 = lwin.r2
	type = lwin.type
	covariate = sp->colvar
	
	if (type == "tail-oriented") {
		zvals = uniqrows(sort(covariate, 1))
    if ((min(r1) < min(zvals)) | (max(r2) > max(zvals))) {
			printf("{err}cannot create tail-oriented window\n")
			_error(1)
		}
    nsubpop = length(r1) + length(r2) + 1
    npatsub = J(nsubpop, 1, 0)
    I0 = J(nsubpop, 1, 1)
    I1 = J(nsubpop, 1, length(zvals))
    if (length(r1) != 0) {
			for (i = 1; i <= length(r1); i++) {
				npatsub[i] = sum(covariate :<= r1[i])
				sel = selectindex(zvals :<= r1[i])
				I1[i] = sel[length(sel)]
			}
    }
    npatsub[length(r1) + 1] = length(covariate)
    if (length(r2) != 0) {
			for (i = 1; i <= length(r2); i++) {
				k = i + length(r1) + 1
				npatsub[k] = sum(covariate :>= r2[i])
				I0[k]	= selectindex(zvals :>= r2[i])[1]
			}
    }
	}
	else {
		real scalar j, stopflag, indinf, indsup
		real vector absfreq, sortedz, cumfreq
		
    zvals = J(length(covariate), 1, 0)
    absfreq = J(length(covariate), 1, 0)
    sortedz = sort(covariate, 1)
    zvals[1] = sortedz[1]
    j = 1
    absfreq[1] = 1
    for (i = 2; i <= length(covariate); i++) {
      if (sortedz[i] != zvals[j]) {
        j++
        zvals[j] = sortedz[i]
        absfreq[j] = 1
      }
      else {
        absfreq[j] = absfreq[j] + 1
      }
    }
    zvals = zvals[|1 \ j|]
    absfreq = absfreq[|1 \ j|]
    cumfreq = absfreq
    for (i = 2; i <= length(cumfreq); i++) {
			cumfreq[i] = cumfreq[i] + cumfreq[i - 1]
		}

    I0 = J(1000, 1, 0)   // max size of I0 and I1 is 1000; these limit the 
    I1 = J(1000, 1, 0)   // max number of stepp subpopulations to 1000 also
    I0[1] = 1
    I1[1] = sum(cumfreq :< r2) + 1
    stopflag = 0
    nsubpop = 2
    while (stopflag == 0) {
      indinf = I0[nsubpop - 1] + 1
      while ((cumfreq[I1[nsubpop - 1]] - cumfreq[indinf - 1]) > r1) {
        indinf++
      }
      I0[nsubpop] = indinf
      indsup = I1[nsubpop - 1]
      while (((cumfreq[indsup] - cumfreq[I0[nsubpop]] + absfreq[I0[nsubpop]]) < r2) & (stopflag == 0)) {
        indsup++
        stopflag = (indsup == length(zvals))
      }
      I1[nsubpop] = indsup
      nsubpop++
    }

    nsubpop--
    npatsub = J(nsubpop, 1, 0)
    npatsub[1] = cumfreq[I1[1]]
    for (i = 2; i <= nsubpop; i++) {
			npatsub[i] = cumfreq[I1[i]] - cumfreq[I0[i] - 1]
		}
    I0 = I0[|1 \ nsubpop|]
    I1 = I1[|1 \ nsubpop|]
	}
	
	real scalar npats
	real vector medians, minz, minc, maxz, maxc
	real matrix subpop

  medians = J(nsubpop, 1, 0)
  minz = J(nsubpop, 1, .)
  minc = J(nsubpop, 1, .)
  maxz = J(nsubpop, 1, .)
  maxc = J(nsubpop, 1, .)

  npats = length(covariate)
  subpop = J(npats, nsubpop, 0)
  for (i = 1; i <= nsubpop; i++) {
    subpop[., i] = (covariate :>= zvals[I0[i]]) :* (covariate :<= zvals[I1[i]])
    medians[i] = round(median(covariate[selectindex(subpop[., i] :== 1)]), .01)
    minz[i] =  round(zvals[I0[i]], .0001)
    minc[i] =  zvals[I0[i]]
    maxz[i] =  round(zvals[I1[i]], .0001)
    maxc[i] =  zvals[I1[i]]
  }

  // update the object
  sp->nsubpop = nsubpop
  sp->subpop = subpop
  sp->npatsub = npatsub
  sp->medianz = medians
  sp->minz = minz
  sp->minc = minc
  sp->maxz = maxz
  sp->maxc = maxc
  sp->init = 1
}

class stsubpop_class scalar stsubpop_class::merge(
	pointer(class stsubpop_class scalar) scalar sp, real vector mergevec)
{
	class stsubpop_class scalar subp_new
	real scalar len_mv, merge_ind, nsubpop_new, i, first, last
	real vector v1, run, medianz, minz, minc, maxz, maxc, covariate, subpop_cov
	real matrix subpop_matrix
	
	subp_new = *sp
	
	len_mv = length(mergevec)
	
	if (sp->nsubpop != len_mv) {
		printf("{err}invalid merge vector length\n")
		_error(1)
	}

	if (sp->nsubpop > 1) {
		v1 = (mergevec[1] \ mergevec[|1 \ (len_mv - 1)|])
		run = (1 \ selectindex(v1 :!= mergevec) \ (len_mv + 1))
		merge_ind = mergevec[1]
	
		subpop_matrix = J(rows(sp->subpop), 0, .)

		nsubpop_new = 0
		for (i = 1; i <= (length(run) - 1); i++) {
			first = run[i]
			last = run[i + 1] - 1
			if (merge_ind & (first != last)) {
				subpop_matrix = (subpop_matrix, rowmax(sp->subpop[., (first..last)]))
				nsubpop_new++
			}
			else {
				subpop_matrix = (subpop_matrix, sp->subpop[., (first..last)])
				nsubpop_new = nsubpop_new + last - first + 1
			}
			merge_ind = !merge_ind
		}

		medianz = J(nsubpop_new, 1, 0)
		minz = J(nsubpop_new, 1, 0)
		minc = J(nsubpop_new, 1, 0)
		maxz = J(nsubpop_new, 1, 0)
		maxc = J(nsubpop_new, 1, 0)

		covariate = sp->colvar
		for (i = 1; i <= nsubpop_new; i++) {
			subpop_cov = covariate[selectindex(subpop_matrix[., i] :== 1)]
			medianz[i] = round((median(subpop_cov)), .01)
			minc[i] = min(subpop_cov)
			minz[i] = round(minc[i], .0001)
			maxc[i] = max(subpop_cov)
			maxz[i] = round(maxc[i], .0001)
		}

		subp_new.nsubpop = nsubpop_new
		subp_new.npatsub = colsum(subpop_matrix)
		subp_new.subpop	= subpop_matrix
		subp_new.medianz = medianz
		subp_new.minz	= minz
		subp_new.minc	= minc
		subp_new.maxz	= maxz
		subp_new.maxc	= maxc
		subp_new.init	= 1
	}

	return(subp_new)
}

class stsubpop_class scalar stsubpop_class::edge(
	pointer(class stsubpop_class scalar) scalar sp, real scalar j,
	string scalar side)
{
	class stsubpop_class scalar subp_new
	real scalar left_edge, right_i, n, start, right_edge, i
	real vector covariate, Z, ZU, left_i, rside, medianz, minz, minc, maxz, maxc,
		subpop_cov
	real matrix subpop_matrix
	
	subp_new = *sp

	if (side == "L" & j == 1) {
		printf("{err}invalid argument j\n")
		_error(1)
	}

	// sort the covariate value
	covariate = sp->colvar
	Z = sort(covariate, 1)
	ZU = uniqrows(Z)

	if (side == "L") {
		left_edge = sp->minc[j - 1]
		left_i = selectindex(ZU :== left_edge)[1]
		rside = selectindex(ZU :== (sp->maxc[j - 1]))
		right_i = rside[length(rside)]
		n = min((right_i - left_i + 1, length(ZU) - right_i + 1))
		start = left_i
	}
	else if (side == "R") {
		right_edge = sp->maxc[j + 1]
		rside	= selectindex(ZU :== right_edge)
		right_i	= rside[length(rside)]
		left_i = selectindex(ZU :== (sp->minc[j + 1]))[1]
		n = min((right_i - left_i + 1, left_i))
		start	= left_i - n + 1
	}
	else {
		printf("{err}unknown side\n")
		_error(1)
	}

	// create a special stepp subpopulation just for this edge
	subpop_matrix = J(length(covariate), n, 0)
	medianz = J(n, 1, 0)
	minz = J(n, 1, 0)
	minc = J(n, 1, 0)
	maxz = J(n, 1, 0)
	maxc = J(n, 1, 0)

	for (i = 1; i <= n; i++) {
		subpop_matrix[., i] = (covariate :>= ZU[start + i - 1]) :* (covariate :<= ZU[start + i + n - 2])
		subpop_cov = covariate[selectindex(subpop_matrix[., i] :== 1)]
		medianz[i] = round(median(subpop_cov), .01)
		minc[i] = min(subpop_cov)
		minz[i] = round(minc[i], .0001)
		maxc[i] = max(subpop_cov)
		maxz[i] = round(maxc[i], .0001)
	}

	subp_new.win = sp->win
	subp_new.colvar = sp->colvar
	subp_new.nsubpop = n
	subp_new.subpop = subpop_matrix
	subp_new.npatsub = colsum(subpop_matrix)
	subp_new.medianz = medianz
	subp_new.minz = minz
	subp_new.minc = minc
	subp_new.maxz = maxz
	subp_new.maxc = maxc
	subp_new.init = 1

	return(subp_new)
}

void stsubpop_class::summary()
{
	real scalar i
	real vector nper
	real matrix temp
	string scalar title
	string vector fnames, cnames, fcname, rnames
	
	nper = J(nsubpop, 1, .)
	temp = J(nsubpop, 5, .)
	
	printf("\n")
	win->summary()
	printf("\n")

	if (init) {
		printf("{txt}Number of subpopulations generated: {res}%f\n", nsubpop)
		title = "Subpopulation covariate summary information (including all treatments):"
		nper = colsum(subpop)
		temp[., 1] = (1::nsubpop)
		temp[., 2] = medianz
		temp[., 3] = minz
		temp[., 4] = maxz
		temp[., 5] = nper'
		rnames = strofreal(temp[., 1])
		fcname = ("Subpopulation")
		cnames = ("Median", "Minimum", "Maximum", "Size")

		if (win->type == "tail-oriented") {
			print_table(temp[., 2..5], rnames, cnames, fcname, 15, 10, title, "", 
				nsubpop, 0, 2, 0, 0,
				selectindex(temp[., 1] :== (length(win->r1) + 1)))
		}
		else {
			print_table(temp[., 2..5], rnames, cnames, fcname, 15, 10, title, "", 
				nsubpop, 0, 2, 0, 0, 0)
		}
	}
	else {
		printf("{err}subpopulations have not been generated yet\n")
	}
	displayflush()
}

end
