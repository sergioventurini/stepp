*!stepp_calc version 0.1.0
*!Written 08Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

class AssociativeArray scalar cuminc_HR(real vector ftime, real vector fstatus,
	| real vector group, real vector strata, real scalar rho, real scalar cencode,
	real vector subset, real scalar na_omit)
{
	class AssociativeArray scalar d, pf, pf_item, stat_res
	real scalar tmp1, nc, ng, ng1, l, ii, jj, n2
	real vector sub_ind, na_ind, time_idx, ugglab, censind, uclab, s, stat,
		causeind, cgind, tmp2, tmp3, tmp4, ome
	real matrix ugg, uc, v, b, omevar
	string vector cg
	
	if (rho == .) rho = 0
	if (cencode == .) cencode = 0
	if (na_omit == .) na_omit = 1   // this should always be set to 1 (i.e.
																	// remove the missing values)
	
	d.reinit("string", 1)
	pf.reinit("string", 1)
	pf_item.reinit("string", 1)
	stat_res.reinit("string", 1)
	
	d.put("time", ftime)
	d.put("cause", fstatus)
	if (group == J(1, 0, .)) {
		group = J(length(ftime), 1, 1)
	}
	d.put("group", as_factor(group))
	if (strata == J(1, 0, .)) {
		strata = J(length(ftime), 1, 1)
	}
	d.put("strata", as_factor(strata))

	tmp1 = length(d.get("time"))
	sub_ind = J(tmp1, 1, 1)
	if (subset != J(1, 0, .)) {
		sub_ind = sub_ind :* subset
		d.put("time", d.get("time")[selectindex(sub_ind)])
		d.put("cause", d.get("cause")[selectindex(sub_ind)])
		d.put("group", d.get("group")[selectindex(sub_ind)])
		d.put("strata", d.get("strata")[selectindex(sub_ind)])
	}
	
	tmp1 = length(d.get("time"))
	na_ind = J(tmp1, 1, 1)
	if (na_omit) {
		na_ind = na_ind :* (d.get("time") :!= .)
		na_ind = na_ind :* (d.get("cause") :!= .)
		na_ind = na_ind :* (d.get("group") :!= .)
		na_ind = na_ind :* (d.get("strata") :!= .)
	}
	else {
		printf("{txt}no method to deal with missing values implemented yet\n")
	}
	d.put("time", d.get("time")[selectindex(na_ind)])
	d.put("cause", d.get("cause")[selectindex(na_ind)])
	d.put("group", d.get("group")[selectindex(na_ind)])
	d.put("strata", d.get("strata")[selectindex(na_ind)])
	
	if (length(d.get("time")) != tmp1) {
		printf("{txt}%f cases omitted due to missing values\n",
			tmp1 - length(d.get("time")))
	}
	cg = J(0, 1, "")
	
	// order d with respect to time
	time_idx = order(d.get("time"), 1)
	d.put("time", d.get("time")[time_idx])
	d.put("cause", d.get("cause")[time_idx])
	d.put("group", d.get("group")[time_idx])
	d.put("strata", d.get("strata")[time_idx])

	ugg = uniqrows(d.get("group"), 1)
	ugg = ugg[selectindex(ugg[., 2] :> 0), 1]
	ugglab = recode(ugg, ugg, uniqrows(group[selectindex(sub_ind :& na_ind)]))
	censind = (d.get("cause") :!= cencode)
	uc = uniqrows(d.get("cause")[selectindex(censind)], 1)
	uclab = uc[selectindex(uc[., 2] :> 0), 1]
	nc = length(uclab)
	ng = length(ugg)
	if (ng > 1) {
		ng1 = ng - 1
		v = J(ng1, ng1, 0)
		s = J(ng1, 1, 0)
	}
	stat = J(nc, 1, 0)
	l = 0
	for (ii = 1; ii <= nc; ii++) {
		causeind = (d.get("cause") :== uclab[ii])
		for (jj = 1; jj <= length(ugg); jj++) {
			cg = (cg \ strofreal(ugglab[jj]) + " " + strofreal(uclab[ii]))
			l++
			cgind = (d.get("group") :== ugg[jj])
			n2 = length(uniqrows(d.get("time")[selectindex(cgind :& causeind)]))
			n2 = 2*n2 + 2
			tmp2 = J(n2, 1, 0)
			tmp3 = J(n2, 1, 0)
			tmp4 = J(n2, 1, 0)
			cinc(d.get("time")[selectindex(cgind)], censind[selectindex(cgind)],
				causeind[selectindex(cgind)], tmp2, tmp3, tmp4)
			pf_item.put("time", tmp2)
			pf_item.put("est", tmp3)
			pf_item.put("var", tmp4)
			pf.put(cg[l], pf_item)
		}
		if (ng > 1) {
			causeind = 2*censind - causeind
			crstm(d.get("time"), causeind, d.get("group"), d.get("strata"), rho, s,
				v)
			stat[ii] = -1
			if (rank(v) == cols(v)) {
				b = I(rows(v))
				stat[ii] = s * qrsolve(v, b) * s
			}

			if (ii == 1) {
			 ome = (-1)*s
			 omevar = v
			}
		}
	}
	if (ng > 1) {
		stat_res.put("uclab", uclab)
		stat_res.put("stat", stat)
		if (stat == J(nc, 1, 0)) {
			printf("\n")
			_error("something went wrong in the calculations; try rerunning the analysis")
		}
		stat_res.put("pv", 1 :- chi2(ng - 1, stat))
		stat_res.put("df", J(length(stat), 1, ng - 1))
		pf.put("Tests", stat_res)
		pf.put("ome", ome)
		pf.put("omevar", omevar)
	}

	return(pf)
}

class AssociativeArray scalar kmest1(real vector y, real vector m,
	real scalar ndf, real vector tpt)
{
	class AssociativeArray scalar kmest1
	real scalar f, n, kr, ntpt, ltp, var, l, i, k, k2, nd, j, t1
	real vector t, s, v, nrr, ndd
	
  f = 1
	n = length(y)
  kr = n
	ntpt = length(tpt) - 1
  nrr = J(ntpt, 1, 0)
  nrr[1] = n
  ndd = J(ntpt, 1, 0)
	t = J(2*ndf + 2, 1, 0)
	s = J(2*ndf + 2, 1, 0)
	v = J(2*ndf + 2, 1, 0)

  ltp = 1
  var = 0
  l = 1
  t[1] = 0
  s[1] = 1
  v[1] = 0
  i = 1

  while (i <= n) {
    k = i + 1
    k2 = 0

    while (k2 <= 0) {
      if (k > n) {
				k2 = 1
			}
      else {
				if (y[k] != y[i]) {
					k2 = 1
				}
				else {
					k++
				}
      }
    }

    k--
    nd = 0

    for (j = i; j <= k; j++) {
			nd = nd + m[j]
		}

    while ((ltp <= ntpt) & (y[i] > tpt[ltp + 1])) {
      ltp++
      nrr[ltp] = kr
    }

    ndd[ltp] = ndd[ltp] + nd
    if (nd > 0) {
      t1 = nd/kr
      f = f*(1 - t1)
      if (nd < kr) var = var + t1/(kr - nd)
      t[l + 1] = y[i]
      s[l + 1] = s[l]
      v[l + 1] = v[l]
      l = l + 2
      t[l] = y[i]
      s[l] = f
      v[l] = var*f*f
    }

    i = k + 1
    kr = n - k
    k = i
  }

  l++
  t[l] = y[n]
  s[l] = s[l - 1]
  v[l] = v[l - 1]
	
	kmest1.reinit("string")
	
	kmest1.put("t", t)
	kmest1.put("s", s)
	kmest1.put("v", v)
	kmest1.put("nrr", trunc(nrr))
	kmest1.put("ndd", trunc(ndd))
	
	return(kmest1)
}

struct kmest_struct vector kmest(real vector time, real vector status, | ///
	real vector group, real vector tpt, real scalar pv, real vector pv_strat,
	real vector pv_sub, real scalar rho, real vector subset, real scalar na_omit)
{
	struct kmest_struct vector z
	class AssociativeArray scalar d, a
	real scalar tmp, Tl, ntpt, i, nstep, ym, ndf
	real vector sub_ind, na_ind, subgl, o, Ty, Tm
	real matrix T1, T8
	string vector lev, tt

	if (rho == .) rho = 0
	if (na_omit == .) na_omit = 1   // this should always be set to 1 (i.e.
																	// remove the missing values)
	
	d.reinit("string", 1)
	d.put("time", time)
	d.put("status", status)
	if (group == J(1, 0, .)) {
		group = J(length(time), 1, 1)
	}
	d.put("group", as_factor(group))
	if (pv_strat == J(1, 0, .)) {
		pv_strat = J(length(time), 1, 1)
	}
	d.put("pv_strat", as_factor(pv_strat))
	if (pv_sub == J(1, 0, .)) {
		pv_sub = J(length(time), 1, 1)
	}
	d.put("pv_sub", as_factor(pv_sub))

	tmp = length(d.get("time"))
	sub_ind = J(tmp, 1, 1)
	if (subset != J(1, 0, .)) {
		sub_ind = sub_ind :* subset
		d.put("time", d.get("time")[selectindex(sub_ind)])
		d.put("status", d.get("status")[selectindex(sub_ind)])
		d.put("group", d.get("group")[selectindex(sub_ind)])
		d.put("pv_strat", d.get("pv_strat")[selectindex(sub_ind)])
		d.put("pv_sub", d.get("pv_sub")[selectindex(sub_ind)])
	}
	
	tmp = length(d.get("time"))
	na_ind = J(tmp, 1, 1)
	if (na_omit) {
		na_ind = na_ind :* (d.get("time") :!= .)
		na_ind = na_ind :* (d.get("status") :!= .)
		na_ind = na_ind :* (d.get("group") :!= .)
		na_ind = na_ind :* (d.get("pv_strat") :!= .)
		na_ind = na_ind :* (d.get("pv_sub") :!= .)
	}
	else {
		printf("{txt}no method to deal with missing values implemented yet\n")
	}
	d.put("time", d.get("time")[selectindex(na_ind)])
	d.put("status", d.get("status")[selectindex(na_ind)])
	d.put("group", d.get("group")[selectindex(na_ind)])
	d.put("pv_strat", d.get("pv_strat")[selectindex(na_ind)])
	d.put("pv_sub", d.get("pv_sub")[selectindex(na_ind)])
	
	if (length(d.get("time")) != tmp) {
		printf("{txt}%f cases omitted due to missing values\n",
			tmp - length(d.get("time")))
	}
	
	if (any(d.get("status") :!= 0 :& d.get("status") :!= 1)) {
		printf("{txt}invalid status values\n")
		_error(3000)
	}
  T1 = uniqrows(d.get("group"), 1)
  subgl = selectindex(T1[., 2] :> 0)
  T8 = d.get("group")[selectindex(d.get("status") :== 1)]
  if (rows(T8) > 0) {
    T8 = (T1[., 2], uniqrows(T8, 1)[., 2])
  } else {
    T8 = (T1[., 2], J(rows(T8), 1, 0))
  }
  T8 = T8[subgl, .]
  T1 = T1[subgl, 1]
	time = d.get("time")
	if (any(time :<= 0)) {
		time[selectindex(time :<= 0)] = J(sum(time :<= 0), 1, .00001)
	}
	nstep = 5
  if (tpt == J(1, 0, .)) {
		tpt = pretty(time, nstep)
	}
  else {
		ym = round(max(time), .01)
    tpt = (0 \ tpt \ ym)
  }
  ntpt = length(tpt) - 1
  lev = J(ntpt, 1, "")
  for (i = 1; i <= ntpt; i++) {
		lev[i] = strofreal(tpt[i]) + " - " + strofreal(tpt[i + 1])
	}
  o = order(time, 1)
  Tl = length(T1)
	z = kmest_struct(Tl)
  for (i = 1; i <= Tl; i++) {
    Ty = time[o]
    Ty = Ty[selectindex(d.get("group")[o] :== T1[i])]
    Tm = d.get("status")[o]
		Tm = Tm[selectindex(d.get("group")[o] :== T1[i])]
    ndf = length(uniqrows(Ty[selectindex(Tm :== 1)]))

    a = kmest1(Ty, Tm, ndf, tpt)

    tt = strofreal(a.get("ndd")) + J(length(a.get("ndd")), 1, "/") +
			strofreal(a.get("nrr"))
//     names(tt) = lev
    z[i].time = a.get("t")
		z[i].est = a.get("s")
		z[i].var = a.get("v")
		z[i].tint = tt
		z[i].ndd = T8[i, .]
	}

  if (pv & (Tl > 1)) {
    _error(3000)
  }

	return(z)
}

real vector tpest1(transmorphic vector x, real scalar n, real vector tp5,
	real scalar ntp)
{
	real scalar l, i, k, loop
	real vector ind
	
	ind = J(ntp, 1, .)

  l = ntp
  for (i = ntp; i >= 1; i--) {
    if (x[n] >= tp5[i]) {
			break
		}
    ind[l] = 0
    l--
  }

  if (l <= 0) {
		return(ind)
	}

  if (x[n] == tp5[l]) {
    ind[l] = n
    l--
  }

  // assuming unique values sorted in ascending order
  k = n - 1
  loop = 1
  while (loop) {
		if (l <= 0) {
			return(ind)
		}

		loop = 0
		for (i = k; i >= 1; i--) {
			if (x[k] <= tp5[l]) {
				ind[l] = k + 1
				l--
				loop = 1
				break
			}
			else {
				k--
			}
		}
	}

  for (i = 1; i <= l; i++) {
		ind[i] = 0
	}

	return(ind)
}

class AssociativeArray scalar tpest(struct kmest_struct vector w,
	real vector times)
{
	class AssociativeArray scalar res, out_nm
	real scalar ng, nt, i, j, w_i_len
	real vector times_new, slct, z1
	transmorphic vector w_i_vec
	real matrix ind, oute, outv
	
  ng = length(w)
  times_new = sort(uniqrows(times), 1)
  nt = length(times_new)
  ind = J(ng, nt, 0)
  oute = J(ng, nt, .)
  outv = J(ng, nt, .)
  slct = J(ng, 1, 1)
	for (i = 1; i <= ng; i++) {
		if (w[i].est == J(1, 0, .)) {
			slct[i] = 0
		}
		else {
			if (i == 1) {
				w_i_vec = w[i].time
			}
			else if (i == 2) {
				w_i_vec = w[i].est
			}
			else if (i == 3) {
				w_i_vec = w[i].var
			}
			else if (i == 4) {
				w_i_vec = w[i].tint
			}
			else if (i == 5) {
				w_i_vec = w[i].ndd
			}
			w_i_len = length(w_i_vec)
			z1 = trunc(tpest1(w_i_vec, w_i_len, times_new, nt))
			ind[i, .] = z1
			if (all(ind[i, .] :> 0)) {
				oute[i, selectindex(ind[i, .] :> 0)] = w[i].est[z1]
				if (w_i_len > 2) {
					outv[i, selectindex(ind[i, .] :> 0)] = w[i].var[z1]
				}
			}
		}
	}
  
	out_nm.reinit("real", 1)
  out_nm.put(1, strofreal((1..ng)))
  out_nm.put(2, strofreal(times_new))

	res.reinit("string", 1)
  res.put("est", oute[selectindex(slct), .])
  res.put("var", outv[selectindex(slct), .])
  res.put("names", out_nm)
	
	return(res)
}

end
