{smcl}
{* *! version 0.1.0  20Jan2021}{...}
{vieweralsosee "stepp" "help stepp"}{...}
{vieweralsosee "steppplot" "help steppplot"}{...}
{viewerjumpto "Syntax" "balance_patients##syntax"}{...}
{viewerjumpto "Description" "balance_patients##description"}{...}
{viewerjumpto "Options" "balance_patients##options"}{...}
{viewerjumpto "Examples" "balance_patients##examples"}{...}
{viewerjumpto "Authors" "balance_patients##authors"}{...}
{viewerjumpto "Stored results" "balance_patients##results"}{...}
{viewerjumpto "References" "balance_patients##references"}{...}
{title:Title}

{p 4 18 2}
{hi:balance_patients} {hline 2} Balanced subpopulations determination (patients)


{marker syntax}{...}
{title:Syntax}

{pstd}
Balanced subpopulations determination for use prior to a STEPP analaysis. This
command focuses on determining subpopulations based on the number of patients
rather than events.

{p 8 14 8}
{cmd:balance_patients} {ifin}{cmd:,} {opth range_r1(numlist)}
{opth range_r2(numlist)} {opth maxn:subpops(#)} {opth cov:ar(varname)}
[{opt p:lot} {cmdab:noshowr:esults} {cmdab:st:oreresults}]

{synoptset 25 tabbed}{...}
{marker balance_patientsopts}{...}
{synopthdr}
{synoptline}
{synopt:{opth range_r1(numlist)}}range (i.e. minimum and maximum) of r1 values to
consider{p_end}
{synopt:{opth range_r2(numlist)}}range (i.e. minimum and maximum) of r2 values to
consider{p_end}
{synopt:{opth maxn:subpops(#)}}maximum number of subpopulations to
consider{p_end}
{synopt:{opth cov:ar(varname)}}covariate to use for generating the subpopulations{p_end}
{synopt:{cmdab:p:lot}}produce a graphical representation of the results (not yet
implemented){p_end}
{synopt:{cmdab:noshowr:esults}}suppress the summary of the results{p_end}
{synopt:{cmdab:st:oreresults}}store a big matrix with all results{p_end}
{synoptline}

{p 4 6 2}
{cmd:by} is allowed with {cmd:balance_patients}; see {help prefix}.


{marker description}{...}
{title:Description}

{pstd} {bf:balance_patients} implements a procedure to identify the optimal
values of r1, r2 and the number of subpopulations so that the variance of
the subpopulation sizes is minimized. The command is intended to be used
prior to a STEPP analysis and it focuses on the number of patients rather
than the number of events.

{pstd} For more details about the STEPP approach, see {help balance_patients##BonettiGelber2004:Bonetti & Gelber (2004)},
{help balance_patients##Bonettietal2009:Bonetti et al. (2009)},
{help balance_patients##Lazaretal2010:Lazar et al. (2010)}, {help balance_patients##Lazaretal2016:Lazar et al. (2016)} 
and {help balance_patients##Yipetal2016:Yip et al. (2016)}.


{marker options}{...}
{title:Options}

{phang}{opth range_r1(numlist)}
provides the range (i.e. minimum and maximum) of r1 values to consider in the
search for the optimal result.

{phang}{opth range_r2(numlist)}
provides the range (i.e. minimum and maximum) of r2 values to consider in the
search for the optimal result.

{phang}{opt maxnsubpops(#)}
indicates the maximum number of subpopulations to consider in the search for the
optimal result.

{phang}{opth covar(varname)}
covariate variable to use for generating the subpopulations.

{phang}{opt plot}
produces a graphical representation of the results (not yet implemented).
 
{phang}{opt noshowresults}
suppresses the results summary.

{phang}{opt storeresults}
stores in memory a (big) matrix with all results (warning: it slows down significantly the process).


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse balance_example, clear}{p_end}

{pstd}Calculations{p_end}
{phang2}{cmd:. balance_patients, range_r1(300 500) range_r2(950 1050) maxnsubpops(50) covar(covar) plot}{p_end}
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
{cmd:balance_patients} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(r1_best)}}overall best value for r1{p_end}
{synopt:{cmd:e(r2_best)}}overall best value for r2{p_end}
{synopt:{cmd:e(var_best)}}overall best value for the subpopulation sizes variance{p_end}
{synopt:{cmd:e(nsubpops_best)}}overall best value for the number of subpopulations{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:balance_patients}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(r1_range)}}range of r1 values used in the analysis{p_end}
{synopt:{cmd:e(r2_range)}}range of r2 values used in the analysis{p_end}
{synopt:{cmd:e(maxnsubpops)}}maximum number of subpopulations considered in the
analysis{p_end}
{synopt:{cmd:e(covariate)}}covariate variable used in the analysis{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(all_res)}}matrix containing the attained subpopulation sizes variance for
each combination of r1, r2 and number of subpopulations; available only if the
{opt storeresults} option has been specified{p_end}

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

{marker Yipetal2016}{...}
{phang}
Yip, W.-K., Bonetti, M., Cole, B. F., Barcella, W., Wang, X. V., Lazar, A. A.
and Gelber, R .D. 2016. Subpopulation Treatment Effect Pattern Plot (STEPP)
analysis for continuous, binary, and count outcomes. Clinical Trials, 13(4):382–390.
{p_end}
