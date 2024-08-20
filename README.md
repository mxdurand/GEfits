# Fitting curves from gas-exchange data


## fitLRC 
Fit light response curve from gas exchange data.  Requires the package `minpack.lm`
  
_Arguments:_  
- **data** [no default]: Input data frame for a single curve.  
- **varnames** [default = list(Photo = "Photo", PARi = "PARi", Rd = "Rd")]: : Defines header of the data frame, defaults to Li-6400 naming convention.  
- **method** [default = "NRH"]: Equation used for light response. Support is done for _non-rectangular hyperbola_ (NRH), _rectangular hyperbola Michaelis-Menten_ (RHMM), _hyperbolic tangent_ (HT), _exponential_ (EXP).    
- **useRd** [default = FALSE]: Use measured Rd? If TRUE, provide it in data and varnames.  
- **plot** [default = TRUE]: If true, plot the fitted curve against data.  
- **add** [default = FALSE]: If true, add plot to pre-existing plot?  
- **verbose** [default = TRUE]: If true, print fitting summary.  
- Additional graphical parameters used by `plot`.  

## percentLRC
Find irradiance at which `percent` of Pgmax is reached.  
  
_Arguments:_  
- **method** [default = "NRH"]: Equation used for light response. Support is done for _non-rectangular hyperbola_ (NRH), _rectangular hyperbola Michaelis-Menten_ (RHMM), _hyperbolic tangent_ (HT), _exponential_ (EXP).    
- **percent** [default = 0]: Percent of Pgmax to look for.  
- **Pgmax** [no default]: Photosynthesis at saturating irradiance  
- **Rd** [no default]: Day respiration.  
- **phiI0** [no default]: Apparent quantum yield of photosynthesis.  
- **theta** [no default]: Curvature parameter (used only is method = NRH).

## fitCO2gm and fitCO2
Fit CO2 response curve from gas exchange data. `fitCO2` considers mesophyll conductance infinite. `fitCO2gm` tries to fit gm. Requires the package `minpack.lm` and `plantecophys`. This function take the fitaci from the latter package as starting point. It is thus heavily inspired by / copied from  it. Work still in development so use at your own risks. 



