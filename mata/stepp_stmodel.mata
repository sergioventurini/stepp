*!stepp_stmodel version 0.1.0
*!Written 12Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 14.2

mata:

void stmodel_class::new()
{
	class_name = "stmodel"
}

string scalar stmodel_class::getclass()
{
	return(class_name)
}

struct steppes_effect scalar stmodel_class::estimate()
{
	printf("{txt}no estimate method available for objects of class stmodel\n")
}

struct steppes_result scalar stmodel_class::test()
{
	printf("{txt}no test method available for objects of class stmodel\n")
}

void stmodel_class::print()
{
	printf("{txt}no print method available for objects of class stmodel\n")
}

end
