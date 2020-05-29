# GranularityEstimation

This has the code (and test data) for estimating an adequate (regular grid) granularity to count crimes,
following the methodology presented in the article (Ramos et al. )"Too fine to be good? Issues of scale, granularity and error in spatial crime analysis". The basic idea is that grid cells that are to broad will hide significant clusters (see Ecological Fallacy), while grid cells that are too fine will show to much noise and spurious highs and lows - therefore a balance must be sought. For a given data set, the algorithm estimates internal uniformity and robustness to error using grids of various granularities, then finding the grid that maximizes both criteria.

See main_orig.r and functions_orig.r for the original code used in Chapter 2 of my PhD dissertation.

See main.r and robust_mapping.r for the updated code, more fit for use in other applications.
