# Code for [Active sensing in the categorization of visual patterns](https://elifesciences.org/content/5/e12215)

## Note
* These code may be incomplete, as they are extracted from a bigger directory.
* Only the first layer of directory structure is kept (for now).
* Only matlab/c functions are kept (no data, fig, etc. for now).
* The folder **correlation-shuffle** contains the latest correlation analysis.
* The folder **exp-setup** contains some of the earliest code for setting up the experiment.

## The Bayesian active sensing bits
* **a_Get_BALDscoreProp.m** computes the BAS score map. Among other functions, it calls **a_Model.m** and **a_GetBALDGM_mex.c**.
* **a_Model.m** computes type posteriors, gradients of parameters, and other stuff.
* **a_GetBALDGM_mex.c** computes the BAS score of individual location.
