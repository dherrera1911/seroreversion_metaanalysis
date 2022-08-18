# Seroreversion meta-analysis project

## Data files

### Raw data files:

Data files to reproduce the analysis from scratch are found
in data/raw_data. These have to be preprocessed before
fitting the model.

**PCR_to_serotest_known.csv**: Serology testing on PCR diagnosed
individuals, with known test-to-serology time.

**PCR_to_serotest_unknown.csv**: Same as above, but with
unknown test-to-serology time.

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

Files with 'XXX_grouped_CV.csv', or 'XXX_random_CV.csv' have the
cross-validation results.

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
| **technique** | Indicates what analytic technique is used by the assay |
| **design** | Indicates the test design (direct/sandwich vs indirect mostly) |
| **sampleType** | Indicates what kind of population was sampled |
| **midpointDate** | Median date of sample collection |


##  Analysis code

Description of the scripts in the ****code**** directory:

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

**04_average_sensitivity_analysis_bis.R**: (Analysis) Fits the
same model as the previous file, but for doing cross-validation.
Outputs the cross-validation results, but no model summary.

**05_characteristics_sensitivity_analysis.R**: (Analysis) Fits a hierarchical
Bayesian regression to the data in PCR_to_serotest_all.csv. Unlike script
04, it includes an effect for the different test characteristics. 

**05_characteristics_sensitivity_analysis_bis.R**: (Analysis) Fits a
hierarchical Bayesian model like the script above, but for
doing cross-validation. Only outputs the CV results, and not a
model summary.

**06_positive_slope_analysis.R**: (Analysis) Fits a model with
two slopes, an early slope and a later slope. It does so on
a small set of tests that show positive slopes in the main analysis.

**07_manufacturer_comparison.R**: (Analysis) Compares the results from
previous model fittings to manufacturer reported sensitivities.

**functions_auxiliary.R**: Contains miscellaneous functions for small tasks.

**functions_seroreversion_fit_analysis.R**: Contains functions related to
the Bayesian analysis fit. For example, preparing the initialization
values, extracting the posterior samples in a tidy format, etc.

**plot_assay_specific_params_summary.R**: Generates Figure 1 of the
paper.

**plot_characteristics_effects.R**: Generates the density plots of
Figures 2-6.

**plot_later_slope_analysis.R**: Generates Figure 7.

**plot_sensitivity_profiles.R**: Generates the sensitivity profiles
in Figures 2-6.

**plot_validation_performance.R**: Computes the performance of
the cross-validation results.


