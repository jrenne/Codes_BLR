
-----------------------------------------------------------
Replication package for

"The time-varying inflation-growth correlation and the yield curve"

By Tilman Bletzinger, Wolfgang Lemke, and Jean-Paul Renne.

Disclaimer: The views expressed in this paper are those of the authors and do not necessarily reflect those of the European Central Bank.

Corresponding author: Jean-Paul.renne@unil.ch
-----------------------------------------------------------

This package contains the codes and data allowing to replicate the above-mentioned paper.


------ How to use the codes?:

To generate all tables and figures of the paper, simply run "main.m" using Matlab. Set indic_estim_moments and indic_estim_MLE to 0 if you do not want to re-estimate the model. (The parameters of the model presented in the paper will then be loaded.)

------ Estimation data:

Estimation data are in the US_data_matlab.xlsx file ("Data" folder).

The output gap is log GDP minus log potential GDP (extracted from the FRED database, tickers GDPC1 and GDPPOT).
Real consumption is from the Bureau of Economic Analysis.
Population is from FRED.
CPI inflation is from the Bureau of Labor Statistics
Nominal (YIELDX) and real (TIPSX) interest rates are from the updated Gürkaynak, Sack, and Wright (2007).
Survey of Professional Forecasters (SPF) data (mean of forecasts) are from the Federal Reserve Bank of Philadelphia for average inflation over the next ten years (CPI10), the average 10-year treasury bond rate (BOND10), and the average 3-month treasury bill rate (BILL10).

------ Outputs:

Figures and tables can be found in folders of the same name.
