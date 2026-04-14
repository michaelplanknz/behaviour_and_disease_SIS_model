# Behaviour and disease model with complex behavioural contagion 

Matlab code for the article (How complex behavioural contagion can prevent infectious diseases from becoming endemic)[https://arxiv.org/abs/2604.10995].




## Verison history

* 14 April 2026 - (arXiv preprint)[https://arxiv.org/abs/2604.10995v1] - results generated using the version of this repo tagged `v1.0`.

All analyses were run in MATALB R2022b.


## How to use this repo

* Run bifurcation_plot.m to create the two-parameter bifurcation diagrams for the behaviour-dependent transmission version of the model (Figure 1) as well as the nullcline plots in Supplementary Material.
* Run main.m to create the phase plane and time series plots for the transmission-modulated and behaviour-modulated versions of the model (Figures 2-4).
* Run runSIRS.m to create the SIRS model plot (Figure 5).

The parameter values for the 10 cases are defined in getPar().

Figures will be saved in the `figures/` folder.
