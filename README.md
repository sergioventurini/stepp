# stepp
###### Current release: 0.2.2-1000
###### Stata version required: at least 17.0
Stata package for subpopulation treatment effect pattern plots (STEPP). The package is a porting to Stata of the [`stepp` R package](https://cran.r-project.org/web/packages/stepp/index.html).

# Installation note

To install `stepp` directly from GitHub you need to use the `github` Stata command. You can install the latest version of the `github` command by executing the following code in your Stata session:

    net install github, from("https://haghish.github.io/github/")

Then, you can install `stepp` simply using the following code in Stata:

    github install sergioventurini/stepp

Alternatively, you can install the package manually downloading the files from this GitHub repository and placing it in your Stata `PERSONAL` directory (if you don't know what this is, run the `adopath` command).

The `examples.do` file contains many examples taken from the literature as well as some simulated data examples. All the example data sets are downloaded and put in place together with the rest of the package.

**Note:** if you installed the package prior to 9 May 2019, you need to manually remove the previous installation. You can do it by deleting all the files used by `stepp` from your Stata `PERSONAL` directory.

# News
- 10/01/2021: the package now includes the possibility to run single group analyses
- 05/10/2020: the package now includes a procedure called `balance_patients` for determining the optimal number of subpopulations
- 01/08/2020: the package now allows to generate subpopulations using the event-based approach described in Lazar et al. (2016)

# Authors
Sergio Venturini, Department of Economic and Social Sciences, Università Cattolica del Sacro Cuore, Cremona, Italy

E-mail: sergio.venturini@unicatt.it

Marco Bonetti, Department of Social and Political Sciences, Università Bocconi, Milan, Italy

E-mail: marco.bonetti@unibocconi.it

Richard Gelber, Department of Biostatistics, Harvard T. H. Chan School of Public Health, Boston, MA, USA

E-mail: gelber@jimmy.harvard.edu

# Acknowledgements
Special thanks to:
- Kristian Romano, for testing the package and checking that it works as expected
- Wai-ki Yip, who developed the [original R package](https://cran.r-project.org/web/packages/stepp/index.html), for the suggestions and insights on the STEPP procedure

# Bugs
In case you find any bug, please send us an e-mail or open an issue on GitHub.

# Citation
To cite the `stepp` package use the following:

Venturini, S., Bonetti, M., Lazar, A. A., Cole, B. F., Wang, X.-V., Gelber, R. D., Yip, W.-K. 2023. Subpopulation treatment effect pattern plot (STEPP) methods with R and Stata. Journal of Data Science, 21(1):106-126.

GitHub repository: https://github.com/sergioventurini/stepp

# References
- Bonetti, M. and Gelber, R .D. 2004. Patterns of treatment effects in subsets of patients in clinical trials. Biostatistics, 5(3):465-481.

- Bonetti, M., Zahrieh, D., Cole, B. F. and Gelber, R .D. 2009. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine, 28(8):1255-68.

- Lazar, A. A., Cole, B. F., Bonetti, M. and Gelber, R .D. 2010. Evaluation of treatment-effect heterogeneity using biomarkers measured on a continuous scale: subpopulation treatment effect pattern plot. Journal of Clinical Oncology, 28(29):4539-4544.

- Lazar, A. A., Bonetti, M., Cole, B. F., Yip, W.-K. and Gelber, R .D. 2016. Identifying treatment effect heterogeneity in clinical trials using subpopulations of events: STEPP. Clinical Trials, 13(2):169–179.

- Venturini, S., Bonetti, M., Lazar, A. A., Cole, B. F., Wang, X.-V., Gelber, R. D., Yip, W.-K. 2023. Subpopulation treatment effect pattern plot (STEPP) methods with R and Stata. Journal of Data Science, 21(1):106-126.

- Yip, W.-K., Bonetti, M., Cole, B. F., Barcella, W., Wang, X. V., Lazar, A. A. and Gelber, R .D. 2016. Subpopulation Treatment Effect Pattern Plot (STEPP) analysis for continuous, binary, and count outcomes. Clinical Trials, 13(4):382–390.

# Copyright
This software is distributed under the GPL-3 license (see LICENSE file).
