*!stepp_estat version 0.3.1
*!Written 11Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program stepp_estat, rclass
	version 14.2
	gettoken subcmd rest : 0 , parse(", ")
	local lsubcmd = length("`subcmd'")
	
	if ("`subcmd'" == substr("__command__", 1, max(2, `lsubcmd'))) {
		__command__ `rest'
	}
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for stepp"
		exit
	}

	return add
end

program __command__, rclass
	version 14.2
	syntax , [ ]
	
	display as text "[for future developments]"
	
end
