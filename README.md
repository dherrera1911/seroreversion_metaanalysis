# Seroreversion meta-analysis project

## Data files

### Raw data files:

Data files to reproduce the analysis are found in data/raw_data:

*PCR_to_serotest_known.csv*: Seroreversion data on people with
previous COVID-19 diagnosis, and where the time elapsed between
diagnosis and serology testing is known.

*PCR_to_serotest_unknown.csv*: Seroreversion data on people with
previous COVID-19 diagnosis, and where the time elapsed between
diagnosis and serology testing is unknown.

*death_dynamics.csv*: Time series of reported cases and deaths
for many locations used in the analysis. These are later used to
estimate the delay between diagnosis and serology testing
for the PCR_to_serotest_unknown.csv data.

*assay_characteristics.csv*: Table containing some technical
characteristics of each test. Mainly, the antigen used,
targeted antibody isotype, and analytic technique.

### Secondary files:

Before doing the statistical analysis, some data preprocessing
is done. These are the generated files, found in data/processed_data:

*PCR_to_serotest_estimated_times.csv*: Same as PCR_to_serotest_unknown.csv,
but with estimated times between COVID diagnosis and serology testing.

*PCR_to_serotest_all.csv*: All seroreversion data put together into
a file. Has both the PCR_to_serotest_known.csv data, and
the PCR_to_serotest_estimated_times.csv data.


### Important variables:

The data file fitted in the statistical analyses is PCR_to_serotest_all.csv.
The most important variables of this dataset are:

| Variable | Description |
| -------- | ----------- |
| *testName* | Indicates the identity of the assay used |
| *testTime* | Time elapsed from diagnosis to testing |
| *citationID* | Identifier of the study reporting the data |
| *nSamples* | Number of people with COVID diagnosis that were tested for antibodies |
| *nSeropositives* | Number of positive antibody tests among the nSamples tested |
| *timeKnown* | Indicates whether testTime is reported in the source, or estimated here |
| *antigenTarget* | Indicates what antigens are targeted by the assay |
| *antibodyTarget* | Indicates what antibody isotypes are targeted by the assay |
| *technique* | Indicates what analytic technique is used by the assay |
| *sampleType* | Indicates what kind of population was sampled |
| *midpointDate* | Median date of sample collection |



##  Analysis code

Description of the scripts in the **code** directory:

*01_death_dynamics_table.R*: Builds the data/raw_data/death_dynamics.csv file
by putting together different case and death data sources.

*02_estimate_seroreversion_delays.R*: Uses the data in data/raw_data/death_dynamics.csv 
to estimate the delays between diagnosis and serology testing for the data
in *PCR_to_serotest_unknown.csv*. Outputs the file PCR_to_serotest_estimated_times.csv.

*03_organize_seroreversion_data.R*: Puts together the information from files
*PCR_to_serotest_known.csv*, *PCR_to_serotest_estimated_tiems.csv* and
*assay_characteristics.csv* into a single file that can be fitted in
the statistical analyses. The output file with this combined data is 
*PCR_to_serotest_all.csv*.

*04_average_sensitivity_analysis.R*: Fits a hierarchical Bayesian regression
to the data in PCR_to_serotest_all.csv. It fits a model to all the data,
and the samples of the parameter posteriors are stored as output in
*sensitivity_decay.csv*. Also, a home-made cross-validation
test that is done on the model, and a summary of the performance is
stored in *data/processed_data/sensitivity_decay_validation.csv*. The model fitted is
defined in *sensitivity_change.stan*.

*05_characteristics_sensitivity_analysis.R*: Fits a hierarchical Bayesian regression
to the data in PCR_to_serotest_all.csv. Unlike script 04, it includes
an effect for the different test characteristics. It stores the fit to all
data in a series of files named
*data/processed_data/assay_characteristics_XXX.csv*, where the XXX string
indicates the characteristic analyzed. 
Also, a home-made cross-validation
test that is done on the model, and a summary of the performance is
stored in *data/processed_data/assay_characteristics_XXX_validation.csv*.

*functions_auxiliary.R*: Contains miscellaneous functions for small tasks.

*functions_seroreversion_fit_analysis.R*: Contains functions related to
the Bayesian analysis fit. For example, preparing the initialization
values, extracting the posterior samples in a tidy format,
convert posterior samples of the parameters into posterior samples
of the sensitivity curve, etc.

*plot_seroreversion_fit.R*: Generate the plots of the data and the fitted
models. The file is divided by numbered headings that indicate the kind
of data to be plotted in that part of the script.


