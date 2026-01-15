# RSVVOI

This R package contains code to:

- conduct cost-effectiveness analyses comparing RSV prevention strategies in Australia
- simulate hypothetical clinical trials to generate evidence to inform the cost-effectiveness analyses
- estimate the value of information for the hypothetical trial designs

The cost-effectiveness analyses use incidence estimates generated from RSV transmission models via the *RSVModels* R package which can be downloaded using *devtools::install_github("michaeldymock25/RSVModels")*.

## Structure

The repository is structured as follows:

### R

Contains R scripts that define documented functions.

### Run

Contains R scripts that use the functions defined in *R/*. These scripts include *prior_distributions.R* and *parameters.R* which contain prior distributions for the statistical model parameters and the distributions of the transmission model parameters, respectively. The *main.R* script loads the required packages, sources the required scripts and allows the user the run the models.

