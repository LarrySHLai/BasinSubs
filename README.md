# BasinSubs
This repository includes R scripts that are developed by Larry Syu-Heng Lai for sedimentary basin subsidence analyses of syn-orogenic strata in the Coastal Range of eastern Taiwan. The data and analytic results are produced by Larry S.-H. Lai and his collaborators, and are used for publications listed below.


## Scripts descriptions:
**Porosity-depth function.R**: Calculating and plotting global averaged porosity-depth function for sandy, muddy, and mean marine sediments.

**Age PB dSL-N.R** & **Age PB dSL-S.R**: Fitting linear age-thickness model for composite stratigraphic columns of northern and sourthen Coastal Range of eastern Taiwan. Linear age models are further used to predict (interpolate & extrapolate) the depositional ages for foraminifera-based paleobathymetry data. Averaged values (and their uncertainties) of estimated paleobathymetry and position of relative sea level (Miller et al., 2021, Science Advances) are calcuated using data around age constraints (unit boundaries), that are further used in subsidence analysis.

**Subsidence Analysis N.R** & **Subsidence Analysis S.R**: Conducting decompaction and backstripping procedure (with error propagation) and plotting results in geohistory diagrams for northern and sourthen Coastal Range of eastern Taiwan. Input data include global averaged porosity-depth functions (_Pofun.csv_), averaged positions of relative sea level at age constraints (_SL_Unit_N.csv_ or _SL_Unit_S.csv_), averaged paleobathymetry at age constraints (_PB_Unit_N.csv_ or _PB_Unit_S.csv_), thickness of each unit bounded by stratigraphic heights of age constraints (_Ages_N_new.csv_ or _Ages_S_new.csv_), and fractions of different lithology in each unit (_Litho_N.csv_ or _Litho_N.csv_).


## Descriptions of data/result folders:
**Thickness folder**: It includes values and uncertainty ranges of the thickness of each unit with corresponding age-constraint bounds, paleobathymetry data and averaged paleobathymetry at each unit age constraint, and fractions of different lithology in each unit.

**Sealevel folder**: It includes reorganized values of the relative sea level position from Miller et al. (2021), and calculated averaged relative sea level values and uncertainties at each age constraint (unit boundary).
 
**Backstripping folder**: It includes compiled porosity-depth functions and grain density for sandy, muddy, and mean marine sediments from previous publications. Calculated global averaged porosity-depth functions used in our analysis are also included. 

**output results folder**: It includes results of our linear regressions for age models, subsidence analysis, and calculated long-term uplift and subsidence rates of the Luzon arc crust of eastern Taiwan.


## Citation:
If you consider using these scripts and data, please cite the following references:

Lai, L.S.-H., Dorsey, R.J., Horng, C.-S., Chi, W.-R., Shea, K.-S., and Yen, J.-Y. (In Review) Extremely rapid up-and-down motions of island arc crust during arc-continent collision. Submitted to *Communication Earth & Environment*.

Lai, L.S.-H., Dorsey, R.J., Horng, C.-S., Shea, K.-S., Chi, W.-R., and Yen, J.-Y. (2021) Extreme up-and-down vertical motions of accreting arc crust in an active collisional suture, Coastal Range of eastern Taiwan. *Geological Society of America Abstracts with Programs*, v. 53, no. 6, https://doi.org/10.1130/abs/2021AM-365341.

Lai, L.S.-H., Dorsey, R.J., Horng, C.-S., Chi, W.-R., Shea, K.-S., and Yen, J.-Y. (2021) Polygenetic mélange in retrowedge foredeep of an active arc-continent collision, Coastal Range of eastern Taiwan. *Sedimentary Geology*, v. 418, 105901. https://doi.org/10.1016/j.sedgeo.2021.105901
