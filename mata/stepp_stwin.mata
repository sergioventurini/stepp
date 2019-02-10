*!stepp_stwin version 0.1.0
*!Written 11Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void stwin_class::new()
{
		class_name = "stwin"
		r1 = r2 = J(0, 1, .)
		type = "sliding"
		basedon = "all"
}

string scalar stwin_class::getclass()
{
	return(class_name)
}

real scalar stwin_class::validity()
{
	if (type == "tail-oriented") {
	  r1 = uniqrows(sort(r1, 1))
	  r2 = uniqrows(sort(r2, 1))
	}

	real scalar valid

	valid = 1

	if (type == "sliding") {
		if (any(r1 :<= 0) | any(r2 :<= 0)) {
			printf("{err}both r1 and r2 must be strictly positive\n")
			valid = 0
// 			_error(1)
		}
		if (any(r1 :> r2)) {
			printf("{err}r1 must be smaller than r2\n")
			valid = 0
// 			_error(1)
		}
	}
	if ((type != "sliding") & (r2 != "tail-oriented")) {
		printf("{err}type must be set to either 'sliding' or 'tail-oriented'\n")
		valid = 0
// 		_error(1)
	}
	if ((basedon != "event") & (basedon != "all")) {
		printf("{err}basedon must be set to either 'all' or 'event'\n")
		valid = 0
// 		_error(1)
	}
	
	return(valid)
}

void stwin_class::summary()
{
	string scalar temp1, temp2
	real scalar i
	
	temp1 = strofreal(r1[1])
	temp2 = strofreal(r2[1])
	for (i = 2; i <= length(r1); i++) {
		temp1 = temp1 + " " + strofreal(r1[i])
	}
	for (i = 2; i <= length(r2); i++) {
		temp2 = temp2 + " " + strofreal(r2[i])
	}
	
	printf("{txt}Window type: %s\n", type)
	if (type == "tail-oriented") {
		printf("{txt}Subpopulation for patients less than or equal to: {res}%f\n", temp1)
		printf("{txt}Subpopulation for patients greater than or equal to: {res}%f\n", temp2)
	}
	else if (type == "sliding") {
		printf("{txt}Number of patients per subpopulation (r2): {res}%f\n", r2)
		printf("{txt}Largest number of patients in common among consecutive subpopulations (r1): {res}%f\n", r1)
	}
	displayflush()
}

end
