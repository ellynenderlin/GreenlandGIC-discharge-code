# GreenlandPeriph-code
The codes contained in this repo were used to estimate iceberg discharge for Greenland's peripheral marine-terminating glaciers. 

# Overview
The flux gate method to estimate iceberg discharge from marine-terminating glaciers calculate the flux of ice across a gate inland of the terminus as the product of the gate's width, thickness, and speed. For Greenland's peripheral glaciers, widths and speeds can easily be extracted from remote sensing datasets (we use Landsat images and ITS_LIVE velocities in the codes here) but thickness observations are scarce. These codes walk through terminus delineation, which is needed to know where to place the flux gate for multi-decadal time series, and the delineation of flux gates and extraction of velocities across the gates. Code and supporting data are also provided for the estimation of thicknesses for the peripheral glaciers, but this step could be omitted and subsequent code adjusted if thickness estimates are available for other glaciers and do not require empirical estimation using surface observations. 

# Example of Regional Discharge Outputs & Data Availability
For Greenland's peripheral glaciers, time series of discharge will be plotted as shown below. Discharge time series produced using these codes can be accessed at https://doi.org/10.18122/geo_data.5.boisestate. The flux gates for which these data were extracted can be found at https://doi.org/10.18122/cryogars_data.2.boisestate.
![GreenlandGIC-discharge-timeseries](https://user-images.githubusercontent.com/51135732/161085285-480a501e-26bb-40bb-88cd-a6c41cd5072d.png)



# Dependencies
The codes are dependent on supporting datasets that must be downloaded separately by the user. These supporting data include the following. Please notify Ellyn Enderlin (ellynenderlin@boisestate.edu) if you try to run these codes and find that you need additional data so that this list can be updated appropriately.

(1) Randolph Glacier Inventory Outlines: https://www.glims.org/RGI/rgi60_dl.html 

(2) time-averaged ArcticDEM: https://www.pgc.umn.edu/data/arcticdem/ 

(3) ITS_LIVE annual velocity time series: https://nsidc.org/apps/itslive/

(4) Regional Atmospheric Climate Model annual surface mass balance time series: https://www.projects.science.uu.nl/iceclimate/models/index.php

(5) Greenland Ice Sheet discharge and thickness data: https://dataverse01.geus.dk/dataset.xhtml?persistentId=doi:10.22008/promice/data/ice_discharge/d/v02

(6) radar-derived glacier thicknesses (GreenlandGIC_UWH_table.csv provided)

(7) manually-drawn terminus boxes for each glacier of interest (can be executed in a GIS software or Matlab)

(8: optional) individual time-stamped digital elevation models overlapping each glacier if you want to estimate thickness change over time (contact Polar Geospatial Center)

The codes are generally designed so that paths and file names are specified in the first section, but there may be places where paths need to be changed if the files are not in the same relative locations as on the code developer's computer.

Use of these data or codes should be accompanied by a reference to Bollen, K. E, E. M. Enderlin, & R. Muhlheim, in press. Dynamic mass loss from Greenland’s marine-terminating peripheral glaciers. J. Glaciol.

# Summary of codes in the order in which they should be executed
A) GreenlandGIC_terminus_centerline_gate_analysis.m: (NOT ALL SECTIONS ARE REQUIRED) This code is the workhorse for discharge estimation. It will loop through terminus boxes and identify which RGI outlines they intersect, combining the geospatial data and accompanying metadata into a single file called Greenland_GIC_centerlines.mat. If Landsat data have been downloaded for manual terminus delineations and projected into the correct coordinate reference system (EPSG 3413 for Greenland), the code will walk through the Landsat directories and allow the user to interactively manually delineate terminus positions for the glaciers, adding the data to the structure. Terminus traces can also be edited using the code. After these sections, the centerline can be traced on Landsat images overlain with velocity vectors. Centerline profiles are used to extract along-flow profiles of elevation, which are used to manually identify the grounding line along the centerline. The centerline grounding line is overlain on velocity maps so that a flux gate perpendicular to flow can be manually drawn for each glacier. Speed time series are extracted across each flux gate. 

B) GreenlandGIC_thickness_extrapolation_curves.m: This code was used to explore potential empirical relationships to estimate ice thickness from speed and width observations. It uses data from the Greenland_GIC_centerlines.mat file created with (A) along with time-average glacier velocity maps and observed ice thicknesses. For the Greenland glaciers and ice caps (GICs), thickness observations were compiled from Operation IceBridge airborne radar as saved as a csv file. In the csv file, columns 1-2 are the EPSG3413 coordinates, column 3 is the terminus box ID, column 4 is the thickness (m), column 6 is the speed (m/yr), and column 8 is the glacier width (m). The "Mankoff" data are from the Greenland Ice Sheet discharge time series presented in https://doi.org/10.5194/essd-12-1367-2020, with data links provided therein. Several plots are produced and the empirical functions are saved to H_vs_U_logarithmic-functions.mat.

C) GreenlandGIC_discharge_timeseries.m: Calculates discharge using the data saved in Greenland_GIC_centerlines.mat and the empirical thickness estimation functions in H_vs_U_logarithmic-functions.mat. Greenland-wide discharge estimates and regional estimates are plotted and can be output as tables.

D) GreenlandGIC_elevation_change_analysis.m: Use this code to estimate changes in ice thickness over time and along flow. Requires surface mass balance time series as well as terminus time series to estimate mass loss between the flux gate and glacier terminus. Digital elevation model time series can also be used to examine changes in thickness over time and potential drivers. 

# Contact information
Any questions about this repo should be directed to Ellyn Enderlin (ellynenderlin@boisestate.edu).
