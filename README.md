# acp_physical: Physical climate projection code for Houser et al., (2015) and Rasmussen et al. (2016)

README file last updated by DJ Rasmussen, dmr2-at-princeton-dot-edu, Mon May 22 11:31:48 PDT 2017

## Citation

This code is intended to accompany the results of

    (1) T. Houser, R.E. Kopp, S.M. Hsiang, M. Delgado, A.S. Jina, K. Larsen,
        M. Mastrandrea, S. Mohan, R. Muir-Wood, D.J. Rasmussen, J. Rising,
        and P. Wilson. (2015). American Climate Prospectus: Economic Risks
        in the United States. Columbia University Press. ISBN: 978-0231174565

    (2) D. J. Rasmussen, M. Meinshausen, and R. E. Kopp. (2016). Probability-
        weighted ensembles of U.S. county-level climate projections for climate
        risk analysis. Journal of Applied Meteorology and Climatology.

Please cite these works when using any results generated with this code.

## Overview

To run all this code requires Python v2.6-7, MATLAB, NCAR Command Language v6.1 or greater, and a Fortran90 compiler such as gfortran.

This code is intended to help end-users who wish to work with, modify or reproduce the probabilistic physical climate projections of Houser et al. (2015) and Rasmussen et al. (2016) in greater a capacity other than the tables and figures provided therein. Key functionality these codes provide include:

1. Generate county-level probabilistic projections of temperature and precipitation using the Surrogate/Model Mixed Ensemble method (SMME)
2. Generate county-level probabilistic projections of temperature and precipitation using the Monte Carlo Pattern-Residual method (MCPR)
3. Generate county-level probabilistic projections of daily maximum wet-bulb temperature using the SMME and MCPR methods
4. Post-processing code to generate summary tables of all results at various geographic aggregations with multiple weighting schemes (e.g. population weighted and land area weighted)

Additional documentation is provided as comments in the scripts. Note that these routines generate data at the daily, monthly, and annual level. These codes exclude processing for Alaska and Hawaii, which require an additional downscaled projection data set with global coverage and a separate processing stream.

This repository includes necessary imput files such as temperature pathways from global CMIP5 models (preprocessed), SMME weights, and MAGICC6 probablistic global mean temperature pathways

In the included directories are:

###download

* **download** Codes to download the raw downscaled projections from the Bureau of Reclamation.

###preproc
* **concat**  Codes to suture each model's historical and future records
* **magicc** Codes to generate SMME weights from MAGICC6 projections and also generate MCPR pattern and residual samplings
* **patternscale** Codes to generate seasonal patterns and perform pattern scaling
* **downscale** Codes to downscale monthly means to daily means using historical weather variability
* **mcpr** Codes to generate monthly probabilistic projections using the MCPR method (random sampling of CMIP5 GCM patterns and residuals)


###postproc
* **tables** Codes to generate tables of regional, state, and county-level aggregate projections for multi-year averages
* **exceedance** Codes to count the number of daily exceedances (arbitrarily set for temperature, wetbulb or precipitation) from netCDF files
* **extract** Codes to extract monthly climate variables from netCDF files
* **wetbulb** Codes to generate historical relationships between Tavg and wetbulb temperature for future projections (requires data from the North American Regional Reanalysis data) Also codes for generating future projections of wetbulb temperature based on projections of future temperature and historical dry-bulb and wet-bulb temperature relationships.
* **stn2cnty** Codes to map U.S. counties to the nearest GHCN weather station. 
* **qc_daily_tas** Codes for assuring daily Tmin < Tavg < Tmax 
* **mcpr** Exceedance calculations and extraction codes for the MCPR method
* **bin_climate** Codes to bin daily temperature and precipitation within each year for each county and ensemble model (useful economists doing climate econometrics)



###Historical Observations
These routines require GHCN station level climate normals and are incluclude in the files. The temporal downscaling requires historical daily weather which is a large dataset. These data are a combined observation-reanalysis dataset that is available on the Rutgers' RU Core at: <https://rucore.libraries.rutgers.edu/rutgers-lib/49865/>. The construction of this observation-reanalysis dataset is described in Appendix A of Rasmussen et al. (2016).


----

    Copyright (C) 2017 by ROBERT E. KOPP AND RHODIUM GROUP LLC

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
