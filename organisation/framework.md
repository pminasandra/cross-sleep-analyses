---
title: Future steps for this paper
author: Pranav
---

# Framework for validating metrics

Below is the description of several of the next few steps, whose results will
comprise the main results of our paper. Ideally, each of the below figures will
be made both (i) separately for each species, and (ii) one plot each with all
species there with e.g. some colour-coding.


## 0. Deviation of activity estimates (VeDBA & proportion of time active)

a.  Since VeDBA is what we are basing our entire analysis on, how does burst
    sampling affect the estimation of VeDBA? For this, let _A_ be an estimate of
    activity for a time-window (either median log VeDBA or the proportion of time
    not in a low activity level). Now, first we make a plot of time vs _Â-A_, where
    Â is the estimate of the activity metric for the same time window under a burst
    sampling regime. Each burst sampling regime can go into the plot as a separate
    line.

b.  A plot of inter-burst interval vs. the overall error in the estimate,
    captured as 1/N * sum((Â-A)²). This shows how the estimate gets slightly
    worse with the interburst interval

## I., II., and III. sleep fragmentation, total sleep time, and number of wake bouts

The same 2 plots as in 0 should be made with each of these metrics for each
species.

Further, for each of these metrics, rank the individuals for each deployment in
decreasing order. Find out how consistent the ranking is with increasing
interburst interval.
