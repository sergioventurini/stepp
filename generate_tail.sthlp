{smcl}
{* *! version 0.1.0  29Jan2021}{...}
{vieweralsosee "stepp" "help stepp"}{...}
{vieweralsosee "steppplot" "help steppplot"}{...}
{viewerjumpto "Syntax" "generate_tail##syntax"}{...}
{viewerjumpto "Description" "generate_tail##description"}{...}
{viewerjumpto "Options" "generate_tail##options"}{...}
{viewerjumpto "Examples" "generate_tail##examples"}{...}
{viewerjumpto "Authors" "generate_tail##authors"}{...}
{viewerjumpto "Stored results" "generate_tail##results"}{...}
{viewerjumpto "References" "generate_tail##references"}{...}
{title:Title}

{p 4 18 2}
{hi:generate_tail} {hline 2} Utility program to generate tail-oriented windows


{marker syntax}{...}
{title:Syntax}

{pstd}
Utility program to generate tail-oriented windows for use prior to a STEPP analaysis given the approximate number of
subpopulations desired.

{p 8 14 8}
{cmd:generate_tail} {ifin}{cmd:,} {opth cov:ariate(varname)} {opth n:sub(#)} 
[{opt d:ir(string)} {cmdab:noshowr:esults}]

{synoptset 25 tabbed}{...}
{marker generate_tailopts}{...}
{synopthdr}
{synoptline}
{synopt:{opth cov:ariate(varname)}}covariate to use for generating the subpopulations{p_end}
{synopt:{opth n:sub(#)}}approximate number of subpopulations to generate{p_end}
{synopt:{opt d:ir(string)}}direction according to which to generate the subpopulations{p_end}
{synopt:{cmdab:noshowr:esults}}suppress the output message{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd} {bf:generate_tail} implements a procedure to identify some covariate values to use for generating 
subpopulations using the tail-oriented approach. The command is intended to be used prior to a STEPP analysis.

{pstd} For more details about the STEPP approach, see {help generate_tail##BonettiGelber2004:Bonetti & Gelber (2004)},
{help generate_tail##Bonettietal2009:Bonetti et al. (2009)},
{help generate_tail##Lazaretal2010:Lazar et al. (2010)}, {help generate_tail##Lazaretal2016:Lazar et al. (2016)} 
and {help generate_tail##Yipetal2016:Yip et al. (2016)}.


{marker options}{...}
{title:Options}

{phang}{opth covariate(varname)}
covariate variable to use for generating the subpopulations.

{phang}{opt nsub(#)}
indicates the approximate number of subpopulations to generate.

{phang}{opt dir(string)}
provides the direction according to which to generate the subpopulations; allowed values are "LE" (windows all
starting from the covariate minimum value and progressively increasing up to the entire sample) or "GE" (windows that
start from the whole sample and progressively decrease but all arrive at the covariate maximum value).
 
{phang}{opt noshowresults}
suppresses the output message.


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigKM, clear}{p_end}

{pstd}Calculation{p_end}
{phang2}{cmd:. generate_tail, covariate(ki67) nsub(10) dir("LE")}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) failure(event) type(km) trts(1 2) timepoint(4.0) nperm(10) wintype("tail-oriented") patspop("`r(patspop)'") minpatspop("`r(minpatspop)'")}{p_end}
    {hline}


{marker authors}{...}
{title:Authors}

{pstd} Sergio Venturini{break}
Department of Management{break}
Università degli Studi di Torino, Turin, Italy{break}
{browse "mailto:sergio.venturini@unito.it":sergio.venturini@unito.it}{break}

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
{cmd:generate_tail} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(nsubpop)}}final number of subpopulations generated{p_end}
{synopt:{cmd:r(nobs)}}number of observations{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(dir)}}direction used in the analysis{p_end}
{synopt:{cmd:r(covariate)}}covariate variable used in the analysis{p_end}
{synopt:{cmd:r(npats)}}list of number of patients in each subpopulation{p_end}
{synopt:{cmd:r(patspop)}}list of covariate minimum values for each subpopulation{p_end}
{synopt:{cmd:r(minpatspop)}}list of covariate maximum values for each subpopulation{p_end}


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

{marker Yipetal2016}{...}
{phang}
Yip, W.-K., Bonetti, M., Cole, B. F., Barcella, W., Wang, X. V., Lazar, A. A.
and Gelber, R .D. 2016. Subpopulation Treatment Effect Pattern Plot (STEPP)
analysis for continuous, binary, and count outcomes. Clinical Trials, 13(4):382–390.
{p_end}
