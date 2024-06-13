# Dynamics of SARS-CoV-2 seroassay sensitivity: a systematic review and modeling study

This repository contains the code and data required to reproduce the results
of the paper "Dynamics of SARS-CoV-2 seroassay sensitivity: a systematic
review and modelling study", Euro Surveill. 2023;28(21).
(https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2023.28.21.2200809).
The paper analyzes how the sensitivity of SARS-CoV-2 serological assays
changes over time since infection, and how this change is influenced by
assay characteristics.

## Reproducing the analysis

In the **code** directory, there are scripts that preprocess the data,
fit the models, and analyze the results. The code is written in R
and the Bayesian models are written in Stan.

The scripts in the **code** directory are led by a number,
indicating the order in which they should be run to fully
reproduce the analysis. Scripts 01 to 03 preprocess the
data, scripts 04 and 05 fit the main models in the
paper (sensitivity vs time), and scripts 06 to
11 are controls and additional analyses in the paper.

It is not necessary to run all scripts to reproduce the analysis,
as the data files in data/processed_data allow to jump
straight into running scripts 04 and 05.

The main outputs of scripts 04 and 05 are also included
in the data/analysis_results directory, so that the
figures and statistics in the paper can be reproduced
without running the analysis. The main figures of the
paper can be reproduced by running the script
**code/plotting_tabulating_scripts/plot_sensitivity_profiles.R**.


##  Analysis code

Each script includes a description at the beggining.
The scripts are numbered so as to be run in order.
Paths are built so as to be ran from the directory where
a script file is located.

Summary of the scripts in the ****code**** directory:

**01_death_dynamics_table.R**: (Preprocessing) Builds the
data/raw_data/death_dynamics.csv file by putting together different
case and death data sources.

**02_estimate_seroreversion_delays.R**: (Preprocessing) Uses the data in
data/raw_data/death_dynamics.csv to estimate the delays between
diagnosis and serology testing for the data in **PCR_to_serotest_unknown.csv**

**03_organize_seroreversion_data.R**: (Preprocessing) Tidies the
data and does som extra pre-processing. The output file of this
script is **PCR_to_serotest_all.csv**, the input to the
statistical analysis.

**04_average_sensitivity_analysis.R**: (Analysis) Fits a hierarchical
Bayesian regression to the data in PCR_to_serotest_all.csv, without
taking into account assay characteristics. Outputs files to
data/analysis_results with descriptions of the fitted model

**04_bis_average_sensitivity_analysis_CV.R**: (Analysis) Fits the
same model as the previous file, but for doing cross-validation.
Outputs the cross-validation results, but no model summary.

**05_characteristics_sensitivity_analysis.R**: (Analysis) Fits a hierarchical
Bayesian regression to the data in PCR_to_serotest_all.csv. Unlike script
04, it includes an effect for the different test characteristics. 

**05_bis_characteristics_sensitivity_analysis_CV.R**: (Analysis) Fits a
hierarchical Bayesian model like the script above, but for
doing cross-validation. Only outputs the CV results, and not a
model summary.

**06_positive_slope_analysis.R**: (Analysis) Fits a model with
two slopes, an early slope and a later slope. It does so on
a small set of tests that show positive slopes in the main analysis.

**07_manufacturer_comparison.R**: (Analysis) Compares the results from
previous model fittings to manufacturer reported sensitivities.

**08_characteristics_analysis_known_times.R**: (Analysis) Does the
same as script 05, but for excluding data points where we
estimated the time from diagnosis to testing.

**09_organize_specificity_data.R**: (Preprocessing) Prepare the
specificity data to analyze how it changes across assays.

**10_analyze_specificity_data.R**: (Analysis) Fit a Bayesian
model to the specificity data, to find effects of assay
characteristics on specificity.

**11_serotracker_analysis.R**: (Analysis) Computes how many
data points in SeroTracker, that are Unity-aligned, use
assays at high-risk of seroreversion.

**functions_auxiliary.R**: Contains miscellaneous functions for small tasks.

**functions_seroreversion_fit_analysis.R**: Contains functions related to
the Bayesian analysis fit. For example, preparing the initialization
values, extracting the posterior samples in a tidy format, etc.

Directory **plotting_tabulating_scripts** has scripts
that generate the figures for the paper. Each script includes
a description of what Figures it generates. Like the analysis
scripts, paths are made so as to have these scripts ran
from the directory where they are located.

Directory **stan_models** has the .stan files that implement
the Bayesian models that are fit in the main analysis scripts
described above.


## Data files

### Raw data files:

Data files to reproduce the analysis from scratch are found
in data/raw_data. These have to be preprocessed before
fitting the model.

**PCR_to_serotest_known.csv**: Serology testing on PCR diagnosed
individuals, with known diagnosis-to-serology time.

**PCR_to_serotest_unknown.csv**: Same as above, but with
unknown diagnosis-to-serology time.

**death_dynamics.csv**: Time series of reported cases and deaths
for different locations used in the analysis.

**assay_characteristics.csv**: Technical
characteristics of each test

### Secondary files:

Pre-processed data files for jumping straight into fitting
the analysis are in data/processed_data:

**PCR_to_serotest_estimated_times.csv**: Same as PCR_to_serotest_unknown.csv,
but with estimated times between COVID diagnosis and serology testing.

**PCR_to_serotest_all.csv**: Data table putting together all the data
to be fitted in the analysis. This is the file needed to jump straight
into model fitting.

**assay_specificity.csv**: Data table with reported specificities
for each assay

### Results data files:

Results from model fitting are in data/analysis_results/, and these
are used for making plots and statistics in the rest of the code.
Results are saved here automatically at the end of analysis scripts. 

We include some results files that have summaries of the model fits,
but we exclude the files with the posterior samples, due to their size.
Some figures of the manuscript can be done with the available files.

Files with 'XXX_sensitivity_curve.csv' name structure have the
time-dependent sensitivity curves estimated from the fitted model
specified in XXX. These files contain estimated sensitivity profiles
for each assay, as well as for the average estimated for a given
kind of assay.

Files with 'XXX_parameter_summary.csv' have the summary statistics
of each parameter fitted by the Bayesian models. 

Files with 'XXX_posterior_samples.csv' contain the posterior samples
of the model that were used to build the previous two files.
These files are not included due to their size.

Files with 'XXX_predicted_sensitivities_grouped_CV.csv', have the
cross-validation results.

Files 'XXX_statistical_significance.csv' contain statistical tests
as proportions of posterior samples that show a result of interest.

Files 'sensitivity_profile_table*.csv' contain the resulting
time-varying sensitivities for different assays and
types of assays.

### Important variables:

The data file fitted in the statistical analyses is PCR_to_serotest_all.csv.
The most important variables of this dataset are:

| Variable | Description |
| -------- | ----------- |
| **testName** | Indicates the identity of the assay used |
| **testTime** | Time elapsed from diagnosis to antibody testing |
| **citationID** | Identifier of the study reporting the data |
| **nSamples** | Number of people with COVID diagnosis that were tested for antibodies |
| **nSeropositives** | Number of positive antibody tests among the nSamples tested |
| **sensitivityMean** | Mean sensitivity for this time point |
| **sensitivityL(H)** | Binomial 95% CI lower (higher) bound of sensitivity |
| **timeKnown** | Indicates whether testTime is reported in the source, or estimated here |
| **antigenTarget** | Indicates what antigens are targeted by the assay |
| **antibodyTarget** | Indicates what antibody isotypes are targeted by the assay |
| **technique** | Indicates analytic technique (i.e. LFIA or which type of Quantitative detection) |
| **design** | Indicates the test design (indirect/direct/competitive) |
| **sampleType** | Indicates what kind of population was sampled |
| **midpointDate** | Median date of sample collection |


##  Meta-analysis summary

In directory **data/systematic_review_summary/**, several
.csv files with notes on the systematic review procedure, data
on the analyzed cohorts, a list of the PRISMA systematic review
procedure, and tables with assay characteristics are found.


## Contact

For questions, contact dherrera1911\[at\]gmail.com


