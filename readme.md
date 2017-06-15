# Code for [Active sensing in the categorization of visual patterns](https://elifesciences.org/content/5/e12215)

## Note
* These code may be incomplete, as they are extracted from a bigger directory.
* Only the first layer of directory structure is kept (for now).
* Only matlab/c functions are kept (no data, fig, etc. for now).
* The folder **correlation-shuffle** contains the latest correlation analysis.
* The folder **exp-setup** contains some of the earliest code for setting up the experiment.

## The Bayesian active sensing core
* **a_Get_BALDscoreProp.m** computes the BAS score map. Among other functions, it calls **a_Model.m** and **a_GetBALDGM_mex.c**.
* **a_Model.m** computes type posteriors, gradients of parameters, and other stuff.
* **a_GetBALDGM_mex.c** computes the BAS score of individual location.

## Making Figure 3A in the paper (testing)
* Run **BAS/analysis final/simu_BASs.m**.
  * Comment out lines 65-66. Uncomment lines 47, 49-50.
  * No input needed. Output is SIM.
* Run **BAS/analysis final/revmap_sim.m**.
  * Comment out lines 138-139, 142-172, 212. Uncomment lines 107-135, 210.
  * DRemovePhase(D) is introduced during the review process.
  * Input is SIM. Outputs are [REV, RrevMap, RdrevMap].
* Run **BAS/analysis final/plot_rev_maps.m**.
  * Inputs are RrevMap, RdrevMap. Output is Figure 3A.
