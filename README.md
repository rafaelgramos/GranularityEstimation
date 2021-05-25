# GranularityEstimation

This contains code related to the paper "Too fine to be good? Issues of scale, granularity and error in spatial crime analysis." (Ramos et al, 2020), as well as Chapter 2 of my PhD dissertation "Methodological contributions to the spatial analysis of crime, with an application to residential burglaries in Belo Horizonte, Brazil. Dissertation" (Ramos, 2021). Although the original code from those studies will be kept here, check out the repository Robustmap, where the same functionalities are available as a package. See below for more detail concerning this study and the code.

## About the research

More and more studies have employed micro-geographic granularities such as street blocks and segments for mapping and analyzing crime. While the use of very fine units can reveal patterns otherwise hidden by coarser units, potentially avoiding issues like the Ecological Fallacy, the robustness of crime counts at such fine levels can be an issue. This algorithm estimates a granularity that balances these two criteria, being fine enough such as not to hide most clusters of crime (i.e. high internal uniformity per unit) while being broad enough to provide robust counts (i.e. low coefficient of variation expected for crime counts).

References:
Ramos, R. G., Silva, B. F., Clarke, K. C., & Prates, M. (2020). Too Fine to be Good? Issues of Granularity, Uniformity and Error in Spatial Crime Analysis. Journal of Quantitative Criminology, 1-25.

Chapter 2 of "RAMOS, R. G. 2019. Methodological contributions to the spatial analysis of crime, with an application to residential burglaries in Belo Horizonte, Brazil. Dissertation (PhD) - University of California Santa Barbara, Department of Geography, Santa Barbara, CA."

## About the repository

This repository has the code (and test data) for estimating an adequate (regular grid) granularity to count crimes.

See main_orig.r and functions_orig.r for the original code used in Chapter 2 of my PhD dissertation.

See main.r and robust_mapping.r for the updated code, more fit for use in other applications.

Test data: reported burglary data for Belo Horizonte (2008-2014). Source: boletins de ocorrência da Policia Militar de Minas Gerais
