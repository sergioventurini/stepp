*!stepp_estat version 0.3.2
*!Written 25Jul2020
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program stepp_estat, rclass
	version 15.1
	gettoken subcmd rest : 0 , parse(", ")
	local lsubcmd = length("`subcmd'")
	
	if ("`subcmd'" == substr("__command__", 1, max(2, `lsubcmd'))) {
		__command__ `rest'
	}
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not yet implemented for stepp"
		exit
	}

	return add
end

program __command__, rclass
	version 15.1
	syntax , [ ]
	
	display as text "[for future developments]"
	
end
