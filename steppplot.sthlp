{smcl}
{* *! version 0.1.0  14Feb2019}{...}
{vieweralsosee "stepp" "help stepp"}{...}
{vieweralsosee "stepp postestimation" "help stepp postestimation"}{...}
{viewerjumpto "Syntax" "steppplot##syntax"}{...}
{viewerjumpto "Description" "steppplot##description"}{...}
{viewerjumpto "Options" "steppplot##options"}{...}
{viewerjumpto "Examples" "steppplot##examples"}{...}
{viewerjumpto "Authors" "steppplot##authors"}{...}
{viewerjumpto "References" "steppplot##references"}{...}
{title:Title}

{p 4 18 2}
{hi:steppplot} {hline 2} Graph results from {helpb stepp}{p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 21 2}{cmd:steppplot} [{cmd:,} {help steppplot##steppplotopts:options}]

{synoptset 20 tabbed}{...}
{marker steppplotopts}{...}
{synopthdr}
{synoptline}
{synopt:{cmdab:s:ubpop}}plot the generated subpopulations{p_end}
{synopt:{cmdab:t:reteff}}plot the estimated treatment effects{p_end}
{synopt:{cmdab:d:iff}}plot the difference between the estimated treatment effects estimates{p_end}
{synopt:{cmdab:r:atio}}plot the ratio estimates between any treatment and treatment 1{p_end}
{synopt:{cmdab:a:ll}}produce all the available grpahs{p_end}
{synopt:{cmdab:c:onf(#)}}set the confidence level to use in the graphs; default is 95{p_end}
{synopt:{cmdab:point:wise}}plot the pointwise confidence intervals instead of the
confidence bands{p_end}
{synopt:{cmdab:nopop:size}}do not plot the subpopulation sizes along the horizontal axis{p_end}
{synopt:{opth trt:labs(numlist)}}treatment names (optional){p_end}
{synopt:{cmdab:x:title(}{it:axis_title}{cmd:)}}specify {it:x} axis title{p_end}
{synopt:{cmdab:y:title(}{it:axis_title}{cmd:)}}specify {it:y} axis title{p_end}
{synopt:{opth tysc:ale(numlist)}}y axis scale for the treatment effects plot (optional){p_end}
{synopt:{opth dysc:ale(numlist)}}y axis scale for the treatment effects difference plot (optional){p_end}
{synopt:{opth rysc:ale(numlist)}}y axis scale for the ratio estimates plot (optional){p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:steppplot} graphs the results of the immediately preceding {helpb stepp} command.


{marker options}{...}
{title:Options}

{phang}
{opt steppplot, subpop} provides a graphical representation of the covariate values
for the generated subpopulations. Subpopulation sizes are also provided.

{phang}
{opt steppplot, trteff} produces a graphical representation of the treatment effect
estimates for each treatment.

{phang}
{opt steppplot, diff} provides a graphical representation of the treatment effect
estimate differences for each treatment versus treatment 1. Confidence intervals
are also provided.

{phang}
{opt steppplot, ratio} provides a graphical representation of the treatment effect
estimate ratios for each treatment versus treatment 1. Confidence intervals
are also provided.

{phang}
{opt steppplot, all} produces all the available graphical representations for
a STEPP analysis.

 
{marker examples}{...}
{title:Examples}

    {hline}
{pstd}{it:Kaplan-Meier method}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigKM, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) failure(event) type(km) patspop(150) minpatspop(50) trts(1 2) timepoint(4.0) nperm(250)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ytitle("4-year disease free survival") nopop}{p_end}
    {hline}
{pstd}{it:Cumulative incidence method}{p_end}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse bigCI, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. stepp time trt, covsubpop(ki67) comprisk(event) type(ci) patspop(150) minpatspop(50) trts(1 2) timepoint(4.0) nperm(250)}{p_end}

{pstd}Graphical analysis{p_end}
{phang2}{cmd:. steppplot, all conf(95) trtlabs(1 Taxmoxifen 2 Letrozole) xtitle("Median Ki-67 LI in subpopulation (% immunoreactivity)") ytitle("4-year disease free survival") nopop}{p_end}
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


{marker references}{...}
{title:References} 

{marker BonettiGelber2004}{...}
{phang}
Bonetti, M., and Gelber, R .D. 2004. Patterns of treatment effects in subsets of
patients in clinical trials. Biostatistics, 5(3):465-481.

{marker Bonettietal2009}{...}
{phang}
Bonetti, M., Zahrieh, D., Cole, B. F., and Gelber, R .D. 2009. A small sample
study of the STEPP approach to assessing treatment-covariate interactions in
survival data. Statistics in Medicine, 28(8):1255-68.

{marker Lazaretal2010}{...}
{phang}
Lazar, A. A., Cole, B. F., Bonetti, M., and Gelber, R .D. 2010. Evaluation of
treatment-effect heterogeneity using biomarkers measured on a continuous scale:
subpopulation treatment effect pattern plot. Journal of Clinical Oncology,
28(29):4539-4544.

{marker Yipetal2016}{...}
{phang}
Yip, W.-K., Bonetti, M., Cole, B. F., Barcella, W., Wang, X. V., Lazar, A. A.
and Gelber, R .D. 2016. Subpopulation Treatment Effect Pattern Plot (STEPP)
analysis for continuous, binary, and count outcomes. Clinical Trials, 13(4):382–390.
{p_end}
