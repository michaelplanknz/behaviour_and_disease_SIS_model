# Behaviour and disease model with complex behavioural contagion 

Matlab code for the article [How complex behavioural contagion can prevent infectious diseases from becoming endemic](https://arxiv.org/abs/2604.10995).


## Verison history

* 14 April 2026 - [arXiv preprint](https://arxiv.org/abs/2604.10995v1) - results generated using the version of this repo tagged `v1.0`.

All analyses were run in MATALB R2022b.


## Abstract

Infectious disease transmission in human populations has a complex two-way interaction with changes in host behaviour. It is increasingly recognised that incorporating adaptive behavioural change into epidemic models is important for improving understanding of infectious disease dynamics and developing policy-relevant modelling tools. An important aspect of behavioural dynamics is social contagion, where people tend to adopt behaviours exhibited by others around them. In a simple behavioural contagion model, the behaviour uptake rate increases linearly with the number of contacts who have adopted a given behaviour. Here, we explore an epidemic model with complex behavioural contagion, where the behaviour uptake rate is a nonlinear function of the number of behaving contacts. We identify key bifurcation parameters of the model, which include the basic reproduction number $R_0$, the strength of the behavioural effect on disease transmission, and the speed of behaviour uptake relative to behaviour abandonment. We show that, in some regions of parameter space, the model has multiple disease-free equilibria. In this situation, the occurrence of an epidemic in a population with an initially low level of behaviour practice can trigger a self-sustaining increase in behaviour, which then causes the disease to be eliminated. In some cases, while moderate values of $R_0$ lead to the disease becoming endemic, higher values of $R_0$ may lead to behaviour-driven disease elimination. We demonstrate that this mechanism of epidemic-triggered uptake of behaviour leading to disease elimination can occur in the presence and absence of temporary post-infection immunity.

## How to use this repo

* Run bifurcation_plot.m to create the two-parameter bifurcation diagrams for the behaviour-dependent transmission version of the model (Figure 1) as well as Supplementary Figures S2 and S3.
* Run main.m to create the phase plane and time series plots for the transmission-modulated and behaviour-modulated versions of the model (Figures 2-4).
* Run runSIRS.m to create the SIRS model plot (Figure 5).
* Run check_3D_model.m to create the plots of equilibrium infection prevalence against initial conditions (Supplementary Figure S4)

The parameter values for the 10 cases are defined in getPar().

Figures will be saved in the `figures/` folder.
