# `convertR`: A climate4R package for unit conversion and variable derivation

`convertR` is an R package performing unit conversion and variable derivation using the [**climate4R**](http://meteo.unican.es/climate4r) Common Data Model.

Unit conversion is done via the Unidata's [UDUNITS-2](https://www.unidata.ucar.edu/software/udunits/) software libraries, tailoring the functionalities of the R package `udunits2` to the **climate4R** framework.

## Unit conversion

The function `udConvertGrid` is used for unit conversion. An example follows:

```R
library(convertR)
data("ps.iberia")
transformeR::getGridUnits(ps.iberia)
## [1] "Pascals"
ps.mmHg <- udConvertGrid(ps.iberia, new.units = "mmHg")
transformeR::getGridUnits(ps.mmHg)
## [1] "mmHg"
```

## Available variable derivations

In the following table the available variable derivations are summarized, including the input variables. The nomenclature of the input variables to calculate the derivations is according to the [climate4R vocabulary](https://github.com/SantanderMetGroup/loadeR/blob/devel/inst/vocabulary.txt) (see `loadeR::c4R.vocabulary()`)

The (*) symbol means that the function is in progress. 


| function  	| Definition                                      	| psl 	| ps 	| tas 	| dpds 	| tdps 	| zs 	| huss 	| hurs 	| rsds 	| rlds 	|
|-----------	|-------------------------------------------------	|-----	|----	|-----	|------	|------	|----	|------	|------	|------	|------	|
| `psl2ps`*    	| Sea-level pressure to surface pressure          	| X   	|    	| X   	|      	|      	| X  	|      	|      	|      	|      	|
| `ps2psl`*    	| Surface pressure to sea-level pressure          	|     	| X  	| X   	|      	|      	| X  	|      	|      	|      	|      	|
| `rad2cc`    	| Radiation to cloud cover                        	|     	|    	|     	|      	|      	|    	|      	|      	| X    	| X    	|
| `huss2hurs` 	| Specific humidity to relative humidity          	|     	| X  	| X   	|      	|      	|    	| X    	|      	|      	|      	|
| `hurs2huss` 	| Relative humidity from specific humidity        	|     	| X  	| X   	|      	|      	|    	|      	| X    	|      	|      	|
| `tas2ws`    	| Saturation vapor pressure from temperature      	|     	| X  	| X   	|      	|      	|    	|      	|      	|      	|      	|
| `hurs2w`    	| Water vapor mixing ratio from relative humidity 	|     	| X  	| X   	|      	|      	|    	|      	| X    	|      	|      	|
| `dpds2hurs`* 	| Relative humidity from dew-point depression     	|     	| X  	| X   	| X    	|      	|    	|      	|      	|      	|      	|
| `tdps2hurs`* 	| Relative humidity from dew-point temperature    	|     	|    	| X   	|      	| X    	|    	|      	|      	|      	|      	|
| `huss2pvp`  	| Partial vapor pressure from specific humidity   	|     	| X  	|     	|      	|      	|    	| X    	|      	|      	|      	|


