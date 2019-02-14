*!stepp_classes version 0.1.0
*!Written 14Feb2019
*!Written by Sergio Venturini, Marco Bonetti and Richard D. Gelber
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

struct gen_tow_struct
{
	real vector v
	real vector np
}

struct ssigma_struct
{
	real matrix sigma
	real matrix sigmainv
}

class stwin_class {
		string scalar		type				// type of stepp window; either "sliding" or
																// "tail-oriented"
		real vector			r1					// sliding window: minimum number of patients
																// 								 allowed in overlapping
																//								 windows
																// tail-oriented window: a vector of maximum
																//											 covariate values for
																//											 each subpopulation
		real vector			r2					// sliding window: size of subpopulation in each
																//								 window
																// tail-oriented window: a vector of minimum
																//											 covariate values for
																//											 each subpopulation
		string scalar		basedon			// what the window is based on (either "event"
																// or "all" [default])
		string scalar		class_name	// class of the object
		
		void						new()				// method to create a stwin_class instance
		string scalar		getclass()	// method to get the class name
		real scalar			validity()	// method to check the validity of an
																// stwin_class instance
 		void						summary()		// method to summarize an stwin_class instance
}

class stsubpop_class
{
		pointer(class stwin_class scalar) scalar ///
																win 				// stepp window object (pointer)
		real vector									colvar			// vector of covariate of interest (V)
		real scalar									nsubpop			// number of subpopulation generated
		real matrix									subpop			// matrix of subpopulations
		real vector									npatsub			// count of each subpopulation
		real vector									medianz			// median of V for each subpopulation
		real vector									minz				// minimum of V for each subpopulation
																						// (rounded to 4 digits)
		real vector									maxz				// maximum of V for each subpopulation
																						// (rounded to 4 digits)
		real vector									minc				// minimum of V for each subpopulation
																						// (actual)
		real vector									maxc				// maximum of V for each subpopulation
																						// (actual)
		real scalar									init				// initialized (0 = no, 1 = yes)
		string scalar								class_name	// class of the object
		
		void												new()				// method to create a stsubpop_class
																						// instance
		string scalar								getclass()	// method to get the class name
		real scalar									validity()
		void 												generate()
		class stsubpop_class scalar merge()
		class stsubpop_class scalar edge()
		void 												summary()
}

class stmodel_class
{
		string scalar					class_name				// class of the object
		
		void									new()							// method to create a stmodel_class
																						// instance
		virtual string scalar	getclass()				// method to get the class name
		virtual struct steppes_effect scalar ///
													estimate()				// method to estimate the model (it
																						// calls the corresponding subclass
																						// method)
		virtual struct steppes_result scalar ///
													test()						// method to compute test statistics
																						// (it calls the corresponding
																						// subclass method)
		virtual void 					print()						// method to display the results
}

class stmodelCI_class extends stmodel_class
{
	string scalar									class_name			// class of the object
	real vector										coltrt					// treatment
	real vector										coltime					// time to event
	real vector										coltype					// competing risk type
	real vector										trts						// trt encoding
	real scalar										timePoint				// evaluated time

	void													new()						// method to create a
																								// stmodelCI_class instance
	string scalar									getclass()			// method to get the class name
	struct steppes_effect scalar	estimate()			// method to estimate the model
	struct steppes_result scalar	test()					// method to compute the test
	void													print()					// method to show model results
}

struct kmest_struct
{
	real vector time
	real vector est
	real vector var
	string vector tint
	real vector ndd
}

class stmodelKM_class extends stmodel_class
{
	string scalar									class_name			// class of the object
	real vector										coltrt					// treatment
	real vector										survTime				// time to event
	real vector										censor					// censoring indicator
	real vector										trts						// trt encoding
	real scalar										timePoint				// evaluated time

	void													new()						// method to create a
																								// stmodelCI_class instance
	string scalar									getclass()			// method to get the class name
	struct steppes_effect scalar	estimate()			// method to estimate the model
	struct steppes_result scalar	test()					// method to compute the test
	void													print()					// method to show model results
}

class stmodelGLM_class extends stmodel_class
{
	string scalar									class_name			// class of the object
	real vector										coltrt					// treatment
	real vector										colY						// outcome
	real vector										trts						// trt encoding
	string vector									MM							// vector of predictor names
	string scalar									glm							// glm model
	string scalar									link						// link function
	real scalar										debug						// debug flag

	void													new()						// method to create a
																								// stmodelCI_class instance
	string scalar									getclass()			// method to get the class name
	real scalar										validity()			// method to check the validity
																								// of an stmodelGLM class
																								// instance
	struct steppes_effect scalar	estimate()			// method to estimate the model
	struct steppes_result scalar	test()					// method to compute the test
	void													subgroup()			// method to show model results
	void													print()					// method to show model results
}

struct steppes_effect
{
	string scalar																	model
	real scalar 																	ntrts
	pointer(class AssociativeArray scalar) scalar TrtEff
	pointer(class AssociativeArray scalar) scalar	Ratios
}

struct steppes_result
{
	string scalar																	model
	real scalar 																	ntrts
	pointer(class AssociativeArray scalar) scalar	Res
}

class steppes_class
{
	pointer(class stsubpop_class scalar) scalar	///
																				subpop			// stepp subpopulations
	struct steppes_effect scalar					effect			// effect estimates
	struct steppes_result scalar					result			// test results
	pointer(class stmodel_class scalar) scalar ///
																				model				// pointer to an stmodel
	real scalar														nperm				// no. of permutations
	string scalar													class_name	// class of the object

	void																	new()				// method to create a
																										// steppes_class instance
	string scalar													getclass()	// method to get the class
																										// name	
	void																	estimate()	// method to estimate the
																										// model (it calls the
																										// corresponding stmodel
																										// method)
	void																	test()			// method to perform the
																										// permutation test
// 	void																	group()
// 	void																	reestimate()
// 	void																	cutpoint()
// 	void																	find_edge()
// 	void																	edge_boot()
	void																	summary()
	void																	plot()
	void																	print()
}

end
