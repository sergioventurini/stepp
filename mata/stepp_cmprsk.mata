*!stepp_cmprsk version 0.1.0
*!Written 12Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void cinc(real vector y, real vector ic, real vector icc, real vector x,
	real vector f, real vector v)
{
	/* Description:
		 ------------
		 Function that calculates estimates of cumulative incidence functions and
		 their variances
	*/
	
	/* Source:
		 -------
		 Adapted from the cinc() Fortran code developed by Robert Gray in the
		 cmprsk R package
	*/
	
	/* Arguments:
		 ----------
		 INPUTS
		 - y					--> real vector of failure times (assumed >= 0)
		 - ic					-->	real vector with elements equal to 1 if the case has
											failed (from any cause), 0 otherwise
		 - icc				-->	real vector with elements equal to 1 if the case has
											failed from the specific cause for which the estimate is
											being calculated, 0 otherwise
		 OUTPUTS
		 - x					--> real vector of times of the cumulative incidence function
		 - f					--> real vector of estimates of the cumulative incidence
											function
		 - v					--> real vector of variances of the cumulative incidence
											function
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

				real scalar n, rs, fk, lcnt, l, ll, nf, v1, v2, v3, ty, i, nd1, nd2, nd,
					fkn, t2, t3, t4, t5, t6
				
				n = rs = length(y)
				
				fk = lcnt = l = ll = 1
				nf = v1 = v2 = v3 = x[1] = f[1] = v[1] = 0
				ty = y[1]

	s10:	l++
				if (l > n) goto s60
				if (y[l] == ty) goto s10
  s60:	l--
				nd1 = nd2 = 0
				for (i = ll; i <= l; i++) {
					nd1 = nd1 + icc[i]
					nd2 = nd2 + ic[i] - icc[i]
				}
				nd = nd1 + nd2
				if (nd == 0) goto s40
				fkn = fk*(rs - nd)/rs
				if (nd1 > 0) {
					lcnt = lcnt + 2
					f[lcnt - 1] = f[lcnt - 2]
					f[lcnt] = f[lcnt - 1] + fk*nd1/rs
				}
				if (nd2 <= 0 | fkn <= 0) goto s30
				t5 = 1
				if (nd2 > 1) {
					t5 = 1 - (nd2 - 1)/(rs - 1)
				}
				t6 = (fk^2)*t5*nd2/(rs^2)
				t3 = 1/fkn
				t4 = f[lcnt]/fkn
				v1 = v1 + (t4^2)*t6
				v2 = v2 + t3*t4*t6
				v3 = v3 + (t3^2)*t6
	s30:	if (nd1 <= 0) goto s35
				t5 = 1
				if (nd1 > 1) {
					t5 = 1 - (nd1 - 1)/(rs - 1)
				}
				t6 = (fk^2)*t5*nd1/(rs^2)
				t3 = 0
				if (fkn > 0) {
					t3 = 1/fkn
				}
				t4 = 1 + t3*f[lcnt]
				v1 = v1 + (t4^2)*t6
				v2 = v2 + t3*t4*t6
				v3 = v3 + (t3^2)*t6
				t2 = f[lcnt]
				x[lcnt - 1] = y[l]
				x[lcnt] = y[l]
				v[lcnt - 1] = v[lcnt - 2]
				v[lcnt] = v1 + (t2^2)*v3 - 2*t2*v2
	s35:	fk = fkn
				nf = nf + nd1
	s40:	rs = n - l
				l++
				if (l > n) goto s50
				ll = l
				ty = y[l]
				goto s10
	s50:	lcnt++
				x[lcnt] = y[n]
				f[lcnt] = f[lcnt - 1]
				v[lcnt] = v[lcnt - 1]
				
				return
}

void crstm(real vector y, real vector m, real vector ig, real vector ist,
	real scalar rho, real vector s, real matrix vs)
{
	/* Description:
		 ------------
		 Function that calculates the score and variance matrix for comparing
		 cumulative incidence curves for a specific cause among groups
	*/
	
	/* Source:
		 -------
		 Adapted from the crstm() Fortran code developed by Robert Gray in the
		 cmprsk R package
	*/
	
	/* Arguments:
		 ----------
		 INPUTS
		 - y					--> real vector of failure times (sorted in increasing order)
		 - m					-->	real vector with elements coded as 0 if censored, 1 if
											failed from the cause of interest, 2 if failed from some
											other cause
		 - ig					-->	real vector denoting group membership, coded 1, 2,..., ng
		 - ist				-->	real vector denoting strata membership, coded 1, 2,...,
											nst (code all 1's if there is only 1 strata)
		 - rho				--> real scalar denoting the power used in the weight function
											in the test statistic
		 
		 OUTPUTS
		 - s					--> real vector denoting the scores for the first (ng - 1)
											groups
		 - vs					--> real matrix containing the estimated variance-covariance
											matrix of the scores in s
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	/* Notes:
		 ------
		 - ng	represents the number of groups
		 - nst represents the number of strata provided
		 - ys	is a vector storing the y values in each stratum
		 - ms is a vector storing the m values in each stratum
		 - igs is a vector storing the ig values in each stratum
		 - The test statistic is given by s' inv(vs) s, which is distributed
			 approximately as a chi-square with (ng - 1) degrees of freedom
		 - Length of ys, ms, and igs are the subvectors related to the strata and
			 thus must be as long as the size of the largest strata
		 - Length of s and st is (ng - 1)
		 - Length of v and vt is ng*(ng - 1)/2
		 - Length of vs is (ng - 1)^2
		 - v contains the same elements as in vs but stacked in a vector
		 - Length of wk is ng*(4 + 3*ng)
		 - Length of iwk is 4*ng
	*/

	real scalar ng, ng1, ng2, nst, l, i, j, ks
	
	real vector ig_vals
	
	ig_vals = sort(uniqrows(ig), 1)
	
	ng = length(ig_vals)
	if (any(ig_vals' :!= (1..ng))) {
		printf("{err}treatments must be coded using labels 1, 2, 3,..., ng\n")
		_error(3000)
	}
	
	ng1 = ng - 1
	ng2 = ng*ng1/2
	nst = length(uniqrows(ist))
	
	real vector si, ys, ms, igs, v, st, vt
	
	s = J(ng1, 1, 0)
	v = J(ng2, 1, 0)
	st = J(ng1, 1, 0)
	vt = J(ng2, 1, 0)

	for (ks = 1; ks <= nst; ks++) {
		si = selectindex(ist :== ks)
		ys = y[si]
		ms = m[si]
		igs = ig[si]
		
		crst(ys, ms, igs, ng, rho, st, vt)
    
		l = 0
		for (i = 1; i <= ng1; i++) {
      s[i] = s[i] + st[i]
			for (j = 1; j <= i; j++) {
				l++
				v[l] = v[l] + vt[l]
			}
		}
	}
	
  l = 0
	for (i = 1; i <= ng1; i++) {
		for (j = 1; j <= i; j++) {
			l++
			vs[i, j] = v[l]
			vs[j, i] = vs[i, j]
		}
	}
	
	return
}

void crst(real vector y, real vector m, real vector ig, real scalar ng,
	real scalar rho, real scalar s, real scalar v)
{
	/* Description:
		 ------------
		 Function that calculates the score and variance matrix for comparing
		 cumulative incidence curves for a specific cause among groups for a given
		 stratum
	*/
	
	/* Source:
		 -------
		 Adapted from the crstm() Fortran code developed by Robert Gray in the
		 cmprsk R package
	*/
	
	/* Arguments:
		 ----------
		 INPUTS
		 - y					--> real vector of failure times (sorted in increasing order)
		 - m					-->	real vector with elements coded as 0 if censored, 1 if
											failed from the cause of interest, 2 if failed from some
											other cause
		 - ig					-->	real vector denoting group membership, coded 1, 2, ..., ng
		 - ng					--> real scalar providing the number of groups
		 - rho				--> real scalar denoting the power used in the weight function
											in the test statistic
		 
		 OUTPUTS
		 - s					--> real vector denoting the scores for the first (ng - 1)
											groups
		 - v					--> real vector 
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

				real scalar ng1, n, i, j, l, fm, f, ll, lu, nd1, nd2, k, tr,
					tq, td, fb, t1, t2, t3, t4, t5

				real vector f1m, f1, skmm, skm, v3, rs
				real matrix d, c, a, v2
	
				ng1 = ng - 1
				n = length(y)
	
				f1m = J(ng, 1, 0)
				f1 = J(ng, 1, 0)
				skmm = J(ng, 1, 1)
				skm = J(ng, 1, 1)
				
				c = J(ng, ng, 0)
				a = J(ng, ng, 0)
				v3 = J(ng, 1, 0)
				v2 = J(ng1, ng, 0)
				
				rs = J(ng, 1, 0)
				// rs[j] will be the risk set size in group j at the current failure
				// time (initially the sample size in each group)
				for (i = 1; i <= n; i++) {
					j = ig[i]
					rs[j] = rs[j] + 1
				}
				
				fm = 0
				f = 0
				
				// begin looping over unique times
				ll = 1
				lu = ll
	s50:	lu++
				if (lu > n) goto s55
				if (y[lu] > y[ll]) goto s55
				goto s50
  s55:	lu--
				nd1 = 0
				nd2 = 0
				
				// d will contain the number censored in each group, failed from
				// cause 1, and failing from cause 2, at this time
				d = J(3, ng, 0)
				for (i = ll; i <= lu; i++) {
					j = ig[i]
					k = m[i]
					d[k + 1, j] = d[k + 1, j] + 1
				}
				for (i = 1; i <= ng; i++) {
					nd1 = nd1 + d[2, i]
					nd2 = nd2 + d[3, i]
				}
				if ((nd1 == 0) & (nd2 == 0)) goto s90
				tr = 0
				tq = 0
				for (i = 1; i <= ng; i++) {
					if (rs[i] <= 0) continue
					td = d[2, i] + d[3, i]
					// skmm is left continuous, and skm right continuous, km est.
					skm[i] = skmm[i]*(rs[i] - td)/rs[i]
					// f1m is left continuous, and f1 right continuous, cuminc est.
					f1[i] = f1m[i] + (skmm[i]*d[2, i])/rs[i]
					tr = tr + rs[i]/skmm[i]
					tq = tq + rs[i]*(1 - f1m[i])/skmm[i]
				}
				f = fm + nd1/tr
				fb = (1 - fm)^rho
				for (i = 1; i <= ng; i++) {
					for (j = i; j <= ng; j++) {
						a[i, j] = 0
					}
					if (rs[i] <= 0) continue
					t1 = rs[i]/skmm[i]
					a[i, i] = fb*t1*(1 - t1/tr)
					c[i, i] = c[i, i] + a[i, i]*nd1/(tr*(1 - fm))
					k = i + 1
					if (k > ng) continue
					for (j = k; j <= ng; j++) {
						if (rs[j] <= 0) continue
						a[i, j] = -fb*t1*rs[j]/(skmm[j]*tr)
						c[i, j] = c[i, j] + a[i, j]*nd1/(tr*(1 - fm))
					}
				}
				for (i = 2; i <= ng; i++) {
					k = i - 1
					for (j = 1; j <= k; j++) {
						a[i, j] = a[j, i]
						c[i, j] = c[j, i]
					}
				}
				for (i = 1; i <= ng1; i++) {
					if (rs[i] <= 0) continue
					s[i] = s[i] + fb*(d[2, i] - nd1*rs[i]*(1 - f1m[i])/(skmm[i]*tq))
				}
				if (nd1 <= 0) goto s77
				for (k = 1; k <= ng; k++) {
					if (rs[k] <= 0) goto s72
					t4 = 1
					if (skm[k] > 0) t4 = 1 - (1 - f)/skm[k]
					t5 = 1
					if (nd1 > 1) t5 = 1 - (nd1 - 1)/(tr*skmm[k] - 1)
					t3 = t5*skmm[k]*nd1/(tr*rs[k])
					v3[k] = v3[k] + t4*t4*t3
					for (i = 1; i <= ng1; i++) {
						t1 = a[i, k] - t4*c[i, k]
						v2[i, k] = v2[i, k] + t1*t4*t3
						for (j = 1; j <= i; j++) {
							l = i*(i - 1)/2 + j
							t2 = a[j, k] - t4*c[j, k]
							v[l] = v[l] + t1*t2*t3
						}
					}
	s72:	}
	s77:	if (nd2 == 0) goto s90
				for (k = 1; k <= ng; k++) {
					if (skm[k] <= 0 | d[3, k] <= 0) continue
					t4 = (1 - f)/skm[k]
					t5 = 1
					if (d[3, k] > 1) t5 = 1 - (d[3, k] - 1.d0)/(rs[k] - 1.d0)
					t3 = t5*((skmm[k]^2)*d[3, k])/(rs[k]^2)
					v3[k] = v3[k] + t4*t4*t3
					for (i = 1; i <= ng1; i++) {
						t1 = t4*c[i, k]
						v2[i, k] = v2[i, k] - t1*t4*t3
						for (j = 1; j <= i; j++) {
							l = i*(i - 1)/2 + j
							t2 = t4*c[j, k]
							v[l] = v[l] + t1*t2*t3
						}
					}
				}
	s90:	if (lu >= n) goto s30
				for (i = ll; i <= lu; i++) {
					j = ig[i]
					rs[j] = rs[j] - 1
				}
				fm = f
				for (i = 1; i <= ng; i++) {
					f1m[i] = f1[i]
					skmm[i] = skm[i]
				}
				ll = lu + 1
				lu = ll
				goto s50
	s30:	l = 0
				for (i = 1; i <= ng1; i++) {
					for (j = 1; j <= i; j++) {
						l++
						for (k = 1; k <= ng; k++) {
							v[l] = v[l] + c[i, k]*c[j, k]*v3[k]
							v[l] = v[l] + c[i, k]*v2[j, k]
							v[l] = v[l] + c[j, k]*v2[i, k]
						}
					}
				}
				
				return
}

end
