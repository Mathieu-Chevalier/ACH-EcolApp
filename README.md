# ACH-EcolApp

This folder contains:
- already published abundance data for 801 species together with freely available range maps. 
- an example code where the abundant-center hypothesis is tested using three distance metrics in both environmental and geographical spaces with range and niche envelopes delineated using two different algorithms (convex hull and minimum volume ellipsoids in the environmental space, convex hull and kernel density estimator in the geographical space), totalling 12 different settings. For each setting the support for the ACH was evaluated using a traditional framework based on correlation coefficients and a hierarchical (miwed effect model) framework where the effect of sampling noise on inferences is reduced. The hierarchical framework was implemented under the Bayesian framework using JAGS but also under the frequentist framework using GAMs. 
- home-made functions created for the purpose of the study.
