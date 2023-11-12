{smcl}
{* *! version 0.1.0  12Nov2023}{...}
{vieweralsosee "stepp postestimation" "help stepp postestimation"}{...}
{vieweralsosee "steppplot" "help steppplot"}{...}
{viewerjumpto "Syntax" "stepp##syntax"}{...}
{viewerjumpto "Description" "stepp##description"}{...}
{viewerjumpto "Options" "stepp##options"}{...}
{viewerjumpto "Examples" "stepp##examples"}{...}
{viewerjumpto "Authors" "stepp##authors"}{...}
{viewerjumpto "Stored results" "stepp##results"}{...}
{viewerjumpto "References" "stepp##references"}{...}
{title:Title}

{p 4 18 2}
{hi:stepp} {hline 2} Subpopulation treatment effect pattern plot (STEPP) analysis


{marker syntax}{...}
{title:Syntax}

{pstd}
Subpopulation treatment effect pattern plot (STEPP) analysis of time to event
data using Kaplan-Meier method.

{p 8 14 8}
{cmd:stepp} {it:{help varname:responsevar}} [{it:{help varname:trtvar}}]
{ifin}{cmd:,} {cmdab:ty:pe(km)} {opth tr:ts(numlist)} {opth covs:ubpop(varname)}
{opth fail:ure(varname)} {opt ti:mepoint(#)} [{opth pat:spop(numlist)}
{opth minp:atspop(numlist)} {opt event:spop(#)} {opt minev:entspop(#)}
{opt mins:ubpops(#)} {it:{help stepp##steppopts:options}}]


{pstd}
Subpopulation treatment effect pattern plot (STEPP) analysis of competing risks
data using cumulative incidence method.

{p 8 14 8}
{cmd:stepp} {it:{help varname:responsevar}} [{it:{help varname:trtvar}}]
{ifin}{cmd:,} {cmdab:ty:pe(ci)} {opth tr:ts(numlist)} {opth covs:ubpop(varname)}
{opth comp:risk(varname)} {opt ti:mepoint(#)} [{opth pat:spop(numlist)}
{opth minp:atspop(numlist)} {opt event:spop(#)} {opt minev:entspop(#)}
{opt mins:ubpops(#)} {it:{help stepp##steppopts:options}}]


{pstd}
Subpopulation treatment effect pattern plot (STEPP) analysis for continuous,
binary and count outcomes.

{p 8 14 8}
{cmd:stepp} {it:{help varname:responsevar}} [{it:{help varname:trtvar}}]
{ifin}{cmd:,} {cmdab:ty:pe(glm)} {opt pat:spop(#)} {opth covs:ubpop(varname)}
{cmdab:fa:mily(}{it:glmfamily}{cmd:)} {cmdab:l:ink(}{it:glmlink}{cmd:)}
[{opth pat:spop(numlist)} {opth minp:atspop(numlist)} {opt event:spop(#)}
{opt minev:entspop(#)} {opt mins:ubpops(#)}
{it:{help stepp##steppopts:options}}]

{synoptset 25 tabbed}{...}
{marker steppopts}{...}
{synopthdr}
{synoptline}
{synopt:{cmdab:ty:pe(km)}}use the Kaplan-Meier method{p_end}
{synopt:{cmdab:ty:pe(ci)}}use the cumulative incidence method{p_end}
{synopt:{cmdab:ty:pe(glm)}}use the generalized linear model method{p_end}
{synopt:{opth pat:spop(numlist)}}number of patients in each subpopulation (called r2
in the literature){p_end}
{synopt:{opth min:patspop(numlist)}}largest number of patients in common among
consecutive subpopulations (called r1 in the literature){p_end}
{synopt:{opt event:spop(#)}}number of events in each subpopulation (called e2
in the literature; only relevant for event-based windows){p_end}
{synopt:{opt minev:entspop(#)}}largest number of events in common among
consecutive subpopulations (called e1 in the literature; only relevant
for event-based windows){p_end}
{synopt:{opt mins:ubpops(#)}}minimum number of subpopulations (only relevant
for event-based windows){p_end}
{synopt:{opth tr:ts(numlist)}}list of treatments included in {it:trtvar}{p_end}
{synopt:{opth covs:ubpop(varname)}}variable to use for generating the subpopulations{p_end}
{synopt:{cmdab:wint:ype(sliding)}}use the sliding window method for generating the
subpopulations; the default{p_end}
{synopt:{cmdab:wint:ype(sliding_events)}}use the event-based sliding window
method for generating the subpopulations{p_end}
{synopt:{cmdab:wint:ype(tail-oriented)}}use the tail-oriented window method for
generating the subpopulations{p_end}
{synopt:{opth fail:ure(varname)}}failure event; to use with {cmdab:type(km)}{p_end}
{synopt:{opth comp:risk(varname)}}variable with distinct codes for different
causes of failure, which must be set to 0 for censored observations, to 1 for
the event of interest and to 2 for other causes of failure; to use with
{cmdab:type(ci)}{p_end}
{synopt:{opt ti:mepoint(#)}}timepoint at which to estimate survival; to use with
{cmdab:type(km)} and {cmdab:type(ci)}{p_end}
{synopt:{opth cov:ariates(varlist)}}optional list of variables to use with
{cmdab:type(glm)}{p_end}
{synopt:{cmdab:fa:mily(gaussian)}}use when the response variable is continuous;
to use with {cmdab:type(glm)}{p_end}
{synopt:{cmdab:fa:mily(binomial)}}use when the response variable is binary;
to use with {cmdab:type(glm)}{p_end}
{synopt:{cmdab:fa:mily(poisson)}}use with a count response variable; to use with
{cmdab:type(glm)}{p_end}
{synopt:{cmdab:li:nk(identity)}}use together with {cmdab:fa:mily(gaussian)};
to use with {cmdab:type(glm)}{p_end}
{synopt:{cmdab:li:nk(logit)}}use together with {cmdab:fa:mily(binomial)};
to use with {cmdab:type(glm)}{p_end}
{synopt:{cmdab:li:nk(logarithm)}}use together with {cmdab:fa:mily(poisson)};
to use with {cmdab:type(glm)}{p_end}
{synopt:{cmdab:nocons:tant}}suppress constant term; to use with
{cmdab:type(glm)}{p_end}
{synopt:{cmdab:note:st}}suppress the permutation test{p_end}
{synopt:{opt np:erm(#)}}number of replications to use in the permutation test{p_end}
{synopt:{opt s:eed(#)}}seed number{p_end}
{synopt:{opt e:ps(#)}}value used in case of times equal to zero{p_end}
{synopt:{cmdab:noshows:ubpops}}suppress the summary of the generated
subpopulations test{p_end}
{synopt:{cmdab:noshowr:esults}}suppress the summary of the results{p_end}
{synoptline}

{p 4 6 2}
{cmd:by} is allowed with {cmd:stepp}; see {help prefix}.

{p 4 6 2}
See {helpb steppplot:steppplot} for features available after estimation.{p_end}


{marker description}{...}
{title:Description}

{pstd} {bf:stepp} implements a statistical method to explore treatment by
covariate interactions in various settings for clinical trials with up to eight
treatment arms. The method is based on constructing overlapping subpopulations
of patients with respect to a covariate of interest, and in observing the
pattern of the treatment effects estimated across subpopulations. A plot of
these treatment effects is called STEPP, or Subpopulation Treatment Effect
Pattern Plot. STEPP uses a permutation based approach for inference.

{pstd} One can explore the window parameters for generating the subpopoulations
without invoking the permutation analysis by specifying the {opt notest} option.
In this case, p-values and covariance matrices will not be produced.

{pstd} STEPP is an exploratory tool with graphical features that make it easy
for clinicians to interpret the results of the analysis. Positive results should
prompt the need for confirmation from other datasets investigating similar
treatment comparisons. Note also that STEPP is not meant to estimate specific
cutpoints in the range of values of the covariate of interest, but rather to
provide some indication on ranges of values where treatment effect might have
a particular behavior.

{pstd} STEPP considers the case in which the subpopulations are constructed
according to a sliding window pattern ({opt wintype(sliding)}). The larger
parameter ({opt patspop}) determines how many patients are included in each
subpopulation, and the smaller parameter ({opt minpatspop}) determines the
largest number of patients in common among consecutive subpopulations. A
minimum of 80-100 patients should be included in each subpopulation, but that
is not strictly necessary. The difference ({opt patspop} - {opt minpatspop}) is
the approximate number of patients replaced between any two subsequent
subpopulations, and can be used to determine the number of subpopulations once
{opt patspop} is fixed. The choice of the values of the parameters
{opt patspop} and {opt minpatspop} to be used does change the appearance of the
plot and the corresponding p-value. It is reasonable to experiment with a
few combinations to ensure that the significance (or lack of significance) is
stable with respect to that choice.

{pstd} A further possibility is to generate the subpopulations still according
to a sliding window pattern but guaranteeing that they contain a given number of
events ({opt wintype(sliding_events)}). The larger parameter ({opt eventspop})
determines how many events are included in each subpopulation, and the smaller
parameter ({opt mineventspop}) determines the largest number of events in
common among consecutive subpopulations.

{pstd} Another alternative is the generation of the subpopulations using
the so-called tail-oriented approach ({opt wintype(tail-oriented)}), whose
aim is to assess the influence of the covariate extreme values on the treatment
effect estimates.

{pstd} For best results, consider implementing 2500 permutations of the
covariate (vector of subpopulations) to obtain a detailed distribution to use
for drawing inference.

{pstd} For more details about the method see {help stepp##BonettiGelber2004:Bonetti & Gelber (2004)},
{help stepp##Bonettietal2009:Bonetti et al. (2009)},
{help stepp##Lazaretal2010:Lazar et al. (2010)}, {help stepp##Lazaretal2016:Lazar et al. (2016)} 
and {help stepp##Yipetal2016:Yip et al. (2016)}.


{marker options}{...}
{title:Options}

{phang}{opt type(modeltype)}
provides the model type to use. Alternative choices are {bf:km} (i.e. Kaplan-Meier),
{bf:ci} (i.e. cumulative incidence) or {bf:glm} (i.e. generalized linear model).

{phang}{opth patspop(numlist)}
indicates the number of patients to include in each subpopulation (usually
referred to as r2 in the literature); for tail-oriented windows it must be
a sequence of numbers corresponding to the covariate minimum values in each
subpopulation.

{phang}{opth minpatspop(numlist)}
indicates the number of patients to include in each subpopulation (usually
referred to as r1 in the literature); for tail-oriented windows it must be
a sequence of numbers corresponding to the covariate maximum values in each
subpopulation.

{phang}{opt eventspop(#)}
indicates the number of events to include in each subpopulation (usually
referred to as e2 in the literature; only relevant for event-based windows).

{phang}{opt mineventspop(#)}
indicates the number of events to include in each subpopulation (usually
referred to as e1 in the literature; only relevant for event-based windows).

{phang}{opt minsubpops(#)}
indicates the minimum number of subpopulations to generate (only relevant for
event-based windows). default to 5.

{phang}{opth trts(numlist)}
provides the list of treatments included in the treatment indicator {it:trtvar}.
Treatments must be reported in the same order as they appear in {it:trtvar}.

{phang}{opth covsubpop(varname)}
covariate variable to use for generating the subpopulations.

{phang}{opt wintype(windows_type)}
sets the type of the windows to use for generating the subpopulations.
Alternative choices are {bf:sliding} (default), {bf:sliding_events} or {bf:tail-oriented}. See the
references for further details.

{phang}{opth failure(varname)}
provides the failure event indicator to use in the analysis when {opt type(km)}
is chosen.

{phang}{opth comprisk(varname)}
provides the variable containing distinct codes for different causes of failure.
The values used in this column must be 0 (censored observations), 1 (event of
interest) or 2 (other causes of failure). To use when {opt type(ci)} is chosen.

{phang}{opt timepoint(#)}
timepoint at which to estimate the survival probability. This option is mandatory
when {opt type(km)} or {opt type(ci)} is chosen.

{phang}{opt family(glm_family)}
glm model to use when {opt type(glm)} is chosen. Alternative choices are
{bf:gaussian} for continuous outcomes, {bf:binomial} for binary outcomes or
{bf:poisson} for count data outcomes. See the references for further details,
especially {help stepp##Yipetal2016:Yip et al. (2016)}.

{phang}{opt link(glm_link)}
link function to use when {opt type(glm)} is chosen. Alternative choices are
{bf:identity} for continuous outcomes, {bf:logit} for binary outcomes or
{bf:logarithm} for count data outcomes. See the references for further details,
especially {help stepp##Yipetal2016:Yip et al. (2016)}.

{phang}{opt noconstant}
suppresses the constant term when {cmdab:type(glm)} is used.

{phang}{opt notest}
suppresses the permutation test. If provided, p-values and covariance matrices
will not be produced.
 
{phang}{opt noshowsubpops}
suppresses the subpopulations summary.
 
{phang}{opt noshowresults}
suppresses the results summary.
 
{phang}{opt nperm(#)}
number of permutation replications to use in the test. default to 25.

{phang}{opt seed(#)}
seed number to use for making the permutation test results reproducible.

{phang}{opt eps(#)}
value added to the observed times with {opt type(km)} to avoid the exclusion
of times equal to zero. default to 0.00001.


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}{it:Kaplan-Meier method}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigKM, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) failure(event) type(km) patspop(150) minpatspop(50) trts(1 2) timepoint(4.0) nperm(250)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) trtlabs(1 Letrozole 2 Tamoxifen) xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ytitle("4-year disease free survival") nopop}{p_end}
    {hline}
{pstd}{it:Cumulative incidence method}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigCI, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) comprisk(event) type(ci) patspop(150) minpatspop(50) trts(1 2) timepoint(4.0) nperm(250)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) trtlabs(1 Letrozole 2 Tamoxifen) xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ytitle("4-year disease free survival") nopop}{p_end}
    {hline}
{pstd}{it:Generalized linear model method (binary outcome)}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse aspirin, clear}{p_end}
{phang2}{cmd:. drop if missing(AD) | missing(AL)}{p_end}
{phang2}{cmd:. keep if DOSE == 0 | DOSE == 81}{p_end}
{phang2}{cmd:. generate trtA = cond(DOSE == 81, 1, 0)}{p_end}
{phang2}{cmd:. generate ADorLE = cond(AD == 1 | AL == 1, 1, 0)}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp ADorLE trtA, covsubpop(AGE) type(glm) patspop(100) minpatspop(30) trts(0 1) nperm(100) family(binomial) link(logit)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) trtlabs(1 "Placebo" 2 "81 mg aspirin") xtitle("Subpopulations by median age") ytitle(Risk) nopop}{p_end}
    {hline}
{pstd}{it:Cumulative incidence method (event-based)}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigKM, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) comprisk(event) type(ci) eventspop(20) mineventspop(10) trts(1 2) timepoint(4) nperm(250) wintype("sliding_events") minsubpops(5)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) tyscale(-13 29) dyscale(-35 35) ryscale(-1 3.2)}{p_end}
    {hline}


{marker authors}{...}
{title:Authors}

{pstd} Sergio Venturini{break}
Department of Economic and Social Sciences{break}
Università Cattolica del Sacro Cuore, Cremona, Italy{break}
{browse "mailto:sergio.venturini@unicatt.it":sergio.venturini@unicatt.it}{break}

{pstd} Marco Bonetti{break}
Carlo F. Dondena Centre for Research on Social Dynamics and Public Policy{break}
Università Bocconi, Milan, Italy{break}
{browse "mailto:marco.bonetti@unibocconi.it":marco.bonetti@unibocconi.it}{break}

{pstd} Richard D. Gelber{break}
Department of Biostatistics and Computational Biology{break}
Dana-Farber Cancer Institute, Boston, MA, USA{break}
{browse "mailto:gelber@jimmy.harvard.edu":gelber@jimmy.harvard.edu}{break}
{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:stepp} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(nperm)}}number of permutations performed{p_end}
{synopt:{cmd:e(timepoint)}}timepoint at which to estimate survival; available
when {opt type(km)} or {opt type(ci)} are chosen{p_end}
{synopt:{cmd:e(e1)}}largest number of events in common among consecutive
subpopulations (only relevant for event-based windows){p_end}
{synopt:{cmd:e(e2)}}number of events in each subpopulation (only relevant for
event-based windows){p_end}
{synopt:{cmd:e(minsubpops)}}minimum number of subpopulations (only relevant for
event-based windows){p_end}
{synopt:{cmd:e(ntrts)}}number of treatments considered{p_end}
{synopt:{cmd:e(eps)}}value added to the response in case of times equal to zero{p_end}
{synopt:{cmd:e(nsubpop)}}number of subpopulations generated{p_end}
{synopt:{cmd:e(oObs_{it:j})}}effect estimate of the entire population based on
the {it:j}-th treatment{p_end}
{synopt:{cmd:e(oSE_{it:j})}}standard error of the effect estimate of the entire
population based on the {it:j}-th treatment{p_end}
{synopt:{cmd:e(skmw_{it:j})}}Wald's statistics for the effect estimate difference
between treatment {it:j} and treatment 1{p_end}
{synopt:{cmd:e(ologHR_{it:j})}}log-hazard ratio estimate of the entire population
comparing treatment {it:j} with treatment 1{p_end}
{synopt:{cmd:e(ologHRSE_{it:j})}}standard error of the log-hazard ratio estimate
of the entire population comparing treatment {it:j} with treatment 1{p_end}
{synopt:{cmd:e(logHRw_{it:j})}}Wald's statistics for the log-hazard ratio
between treatment {it:j} and treatment 1{p_end}
{synopt:{cmd:e(pvalue_{it:j})}}supremum test p-value based on effect difference{p_end}
{synopt:{cmd:e(chi2pvalue_{it:j})}}chi-square test p-value based on effect difference{p_end}
{synopt:{cmd:e(HRpvalue_{it:j})}}supremum test p-value based on hazard ratio{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stepp}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(estat_cmd)}}program used to implement {cmd:estat}{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(properties)}}whether the test hasn't been performed{p_end}
{synopt:{cmd:e(responsevar)}}response variable{p_end}
{synopt:{cmd:e(trtvar)}}treatment indicator variable{p_end}
{synopt:{cmd:e(censorvar)}}failure event variable; available when {opt type(km)}
is chosen{p_end}
{synopt:{cmd:e(compriskvar)}}event indicator variable; available when {opt type(ci)}
is chosen{p_end}
{synopt:{cmd:e(family)}}GLM model used; available when {opt type(glm)} is
chosen{p_end}
{synopt:{cmd:e(link)}}GLM link function used; available when {opt type(glm)} is
chosen{p_end}
{synopt:{cmd:e(covariates)}}list of additional covariates used; available when
{opt type(glm)} is chosen{p_end}
{synopt:{cmd:e(noconstant)}}constant suppressed in the estimation of the GLM
model; available when {opt type(glm)} is chosen{p_end}
{synopt:{cmd:e(wintype)}}type of window used to generate the subpopulations{p_end}
{synopt:{cmd:e(covsubpop)}}covariate used to generate the subpopulations{p_end}
{synopt:{cmd:e(type)}}choice of model type (either {it:km}, {it:ci} or
{it:glm}){p_end}
{synopt:{cmd:e(r1)}}largest number of patients in common among consecutive
subpopulations; for tail-oriented windows it corresponds to the covariate
minimum values in each subpopulation{p_end}
{synopt:{cmd:e(r2)}}number of patients in each subpopulation; for tail-oriented
windows it corresponds to the covariate minimum values in each subpopulation{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(trts)}}vector of treatment names{p_end}
{synopt:{cmd:e(subpop)}}matrix of subpopulations generated based on the stepp
window and the specified covariate of interest{p_end}
{synopt:{cmd:e(npatsub)}}vector of sizes of each subpopulation{p_end}
{synopt:{cmd:e(medianz)}}vector of median value of the covariate of interest for
each subpopulation{p_end}
{synopt:{cmd:e(minc)}}vector of minimum value of the covariate of interest for
each subpopulation{p_end}
{synopt:{cmd:e(maxc)}}vector of maximum value of the covariate of interest for
each subpopulation{p_end}
{synopt:{cmd:e(sObs_{it:j})}}vector of effect estimates of all subpopulations
based on the {it:j}-th treatment{p_end}
{synopt:{cmd:e(sSE_{it:j})}}vector of standard errors of effect estimates of all
subpopulations based on the {it:j}-th treatment{p_end}
{synopt:{cmd:e(logHR_{it:j})}}vector of log-hazard ratio estimate of the
subpopulations comparing treatment {it:j} with treatment 1{p_end}
{synopt:{cmd:e(logHRSE_{it:j})}}vector of standard errors of the log-hazard ratio
estimate of the entire population comparing treatment {it:j} with treatment 1{p_end}
{synopt:{cmd:e(sigma_{it:j})}}covariance matrix for subpopulations based on effect
differences{p_end}
{synopt:{cmd:e(HRsigma_{it:j})}}covariance matrix for subpopulations based on
the hazard ratio{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References} 

{marker BonettiGelber2004}{...}
{phang}
Bonetti, M. and Gelber, R .D. 2004. Patterns of treatment effects in subsets of
patients in clinical trials. Biostatistics, 5(3):465-481.

{marker Bonettietal2009}{...}
{phang}
Bonetti, M., Zahrieh, D., Cole, B. F. and Gelber, R .D. 2009. A small sample
study of the STEPP approach to assessing treatment-covariate interactions in
survival data. Statistics in Medicine, 28(8):1255-68.

{marker Lazaretal2010}{...}
{phang}
Lazar, A. A., Cole, B. F., Bonetti, M. and Gelber, R .D. 2010. Evaluation of
treatment-effect heterogeneity using biomarkers measured on a continuous scale:
subpopulation treatment effect pattern plot. Journal of Clinical Oncology,
28(29):4539-4544.

{marker Lazaretal2016}{...}
{phang}
Lazar, A. A., Bonetti, M., Cole, B. F., Yip, W.-K. and Gelber, R .D. 2016.
Identifying treatment effect heterogeneity in clinical trials using
subpopulations of events: STEPP. Clinical Trials, 13(2):169–179.

{marker Venturinietal2023}{...}
{phang}
Venturini, S., Bonetti, M., Lazar, A. A., Cole, B. F., Wang, X.-V., Gelber, R. D., Yip, W.-K. 2023. Subpopulation treatment effect pattern plot (STEPP) methods with R and Stata. Journal of Data Science, 21(1):106-126.

{marker Yipetal2016}{...}
{phang}
Yip, W.-K., Bonetti, M., Cole, B. F., Barcella, W., Wang, X. V., Lazar, A. A.
and Gelber, R .D. 2016. Subpopulation Treatment Effect Pattern Plot (STEPP)
analysis for continuous, binary, and count outcomes. Clinical Trials, 13(4):382–390.
{p_end}
