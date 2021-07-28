# About these files

The files in this folder accompany the article **Monoclonal antibodies block transmission of genetically diverse Plasmodium falciparum strains to mosquitoes** by Jore, de Jong, Meerstein-Kessel *et al.*, which will shortly be published in **NPJ Vaccines**. The R scripts collected here were used to perform the regression analysis presented in the paper, and generate Figure 1 & Figure 3A. 

Once the files have been downloaded, we recommend double-clicking on the 'Final_analyses.Rproj' file (assuming RStudio is already installed) to open the R project. The advantage of this is that the RStudio session should start in the correct working directory (where the data files are stored). 

The regression analysis was carried out using the rethinking package- installation instructions for which can be found here: https://github.com/rmcelreath/rethinking. This package acts as a wrapper to RStan. Once a particular model has been run (e.g. model 'M') the corresponding Stan code can be viewed by running the command 'stancode(M)'.

There are 3 scripts used for the regression analysis displayed in Figure 1 (summarised in Supplementary Table 5)- 'TRA_regression_4845.R', 'TRA_regression_25', & 'TRA_regression_230'. The data for these models is prepared for analysis using the script 'ViewData.R' (there's no need to modify this script, unless you're curious). 

We use a consistent notation when naming the regression models, specifically: 

m1me: one slope for all the data

m2me: Slope correction for the assay used

m5me: Slope correction for the country for the DMFA data

m6me: Slope corrections for both the assay used, and the country for the DMFA data

Note that only model m1me is present for Pfs230. For Pfs25 & Pfs48/45, we use the 'compare' function from the rethinking package to assess the four models' goodness of fit. This information is displayed in Supplementary Table 5 in the article. The 'ensemble' function (from the same package) is then used to generate the plots displayed in Figure 1, using samples from the model ensemble.

The script 'collate_Figure1.R' can be run after the regression analyses have been performed, to join all the panels together. Specifically, run the scripts 'TRA_regression_4845.R', 'TRA_regression_25', & 'TRA_regression_230' first (in any order, no need to run 'ViewData.R' beforehand). Please note that the regression analyses may take a few minutes to run.

The script 'SMFA230_Figure3A' generates the plot shown in Figure 3A.

Note that both the ensemble weighting of the regression models, and the 95% credible intervals for the IC80s are calculated from samples from the posterior distributions. These may change slightly, compared to the published versions, when new samples are drawn.
