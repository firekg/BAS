# Code for "Active sensing in the categorization of visual patterns"

## Note
* These code may be incomplete, as they are extracted from a bigger directory.
* Only the first layer of directory structure is kept (for now).
* Only matlab functions are kept (no data, fig, etc. for now).

## The Bayesian active sensing bits
* **a_Get_BALDscoreProp.m** computes the BAS score map and calls **a_Model.m** and **a_GetBALDGM_mex.c**.
* **a_Model.m** computes type posteriors, gradients of parameters, and other stuff.
* **a_GetBALDGM_mex.c** computes the BAS score of individual location.
