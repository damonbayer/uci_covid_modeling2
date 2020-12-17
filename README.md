# uci_covid_modeling

This github contains code for generating situation reports on the dynamics and future trends of COVID-19 in Orange County, CA. The manuscript associated with this work is [here](https://arxiv.org/abs/2009.02654), and the current Orange County situation report is [available here](https://damonbayer.github.io/uci_covid_modeling2/). 

## Navigation
Files used to generate website pages are stored in the [analysis folder](https://github.com/damonbayer/uci_covid_modeling2/tree/main/analysis). Each of the `yyyy-mm-dd_yyyy-mm-dd.Rmd` files follow the same format defined in the [analysis template](https://github.com/damonbayer/uci_covid_modeling2/blob/main/code/analysis_template.Rmd).

Functions used in website generation, modeling, and processing are stored in the [code folder](https://github.com/damonbayer/uci_covid_modeling2/tree/main/code). 

Aggregated data used as the model input are available in the [data folder](https://github.com/damonbayer/uci_covid_modeling2/tree/main/data)


## Data
We use data aggregated from data from provided by Orange County, California Health Care Agency. Crucially, we exclude repeat tests given to patients who have tested positive for COVID-19 (which happens when patients are hospitalized). Our data may not match publicly available sources. 

## Methodology
Our analysis relies on a mechanistic six compartment model of the COVID-19 pandemic. We then use Bayesian inference to provide inference on key disease dynamics and make predictions on future observed cases and deaths. Further descriptions of the model are available in the manuscript linked above. 

We use the [stemr](https://github.com/fintzij/stemr/) package for our Bayesian inference. See the [stemr code](https://github.com/damonbayer/uci_covid_modeling2/blob/main/code/fit_new_model.R) for modeling and inference details. 

To replicate our analysis, use the code from any of `yyyy-mm-dd_yyyy-mm-dd.Rmd` files in the [analysis folder](https://github.com/damonbayer/uci_covid_modeling2/tree/main/analysis).

To conduct your own analysis of COVID-19 trends using our model you will need the following data:

1. Daily number of reported cases
2. Daily number of tests (viral not antibody)
3. Daily number of reported deaths

Due to reporting delays, we recommend limiting the end date of your data set to well before the actual end date of available data, as past counts are often updated retroactively weeks later. In our analysis we then condensed this daily data into three day periods. You can adjust this as needed. You will also need to specify a population size. It is likely you will want to adjust the priors for initial conditions.

## Citation
Fintzi J, Bayer D, Goldstein I, Lumbard K, Ricotta E, Warner S, Busch LM, Strich JR, Chertow DS, Parker DM, Boden-Albala B, Dratch A, Chhuon R, Quick N, Zahn M, Minin VN. Using multiple data streams to estimate and forecast SARS-CoV-2 transmission dynamics, with application to the virus spread in Orange County, California, https://arxiv.org/abs/2009.02654.
