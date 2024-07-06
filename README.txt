
-----------------------------------------------------------
Replication package for

"Time-varying risk aversion and inflation-consumption correlation in an equilibrium term structure model"

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
Real consumption of all goods and services is from the Bureau of Economic Analysis.
Population is from FRED.
CPI all items inflation is from the Bureau of Labor Statistics
Nominal (YIELDX) and real (TIPSX) interest rates are from the updated Gürkaynak, Sack, and Wright (2007).
Backcasted real (REALRX) interest rates are from Jan Groen and Menno Middeldorp (2013, Federal Reserve Bank of New York).
The perceived target rate of inflation (PTR) and the expected federal funds rate in the long run (RTR) are from the FRB-US model of the Federal Reserve Board.
Survey of Professional Forecasters (SPF) data (mean of forecasts) are from the Federal Reserve Bank of Philadelphia for average inflation over the next ten years (CPI10), and the average 3-month treasury bill rate (BILL10).

------ Outputs:

Figures and tables can be found in folders of the same name.
