# RSVVOI

This R package contains code to:

- conduct cost-effectiveness analyses comparing RSV prevention strategies in Australia
- simulate hypothetical clinical trials to generate evidence to inform the cost-effectiveness analyses
- estimate the value of information for the hypothetical trial designs

# Installation

This package can be installed via *devtools::install_github("michaeldymock25/RSVVOI", dependencies = TRUE)*. The cost-effectiveness analyses use incidence estimates generated from RSV transmission models via the *RSVModels* R package which can be downloaded using *devtools::install_github("michaeldymock25/RSVModels", dependencies = TRUE)*.

## Structure

The repository is structured as follows:

### Data

Contains data files used to define the transmission model parameters and health economic model parameters (required by functions defined in *R/*). Includes data on the Australian population distribution, their life expectancy, quality-adjusted life years, and mixing dynamics.

### man

Contains helper *.Rd* files for functions defined in *R/*.

### Parameters

Contains saved parameter distributions in .rds files including *trans_parms.rds*, *cea_parms.rds* and *hyper_parms.rds*. These files are generated using functions defined in *R/*. Also contains the saved burn-in distributions from the transmission models (*y0_burn.rds*) and a script that provides additional information on the definition of some of the health economic model parameters (*cea_parameter_selection.R*).

### R

Contains R scripts that define documented functions including *trans_functions.R*, *cea_functions.R*, *trial_functions.R* and *misc_functions.R*.

### Run

Contains R scripts that use the functions defined in *R/*. The *make_mixing_matrix.R* script generates the mixing matrix and aggregated Australian population data files stored in *Data/*. The *burn_trans_model.R* script runs the base transmission model to reach a stationary solution. The *main.R* script loads the required packages, sources the required scripts and allows the user the run the models.

