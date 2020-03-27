Using precip to model water level should only require one coefficient for precip
(with a lag to-be-determined (probably 0 or 1)), an MA process (which looks
likely for wetland water level data)

How do I specify moving average of errors, but also MA process of a predictor.
Precip is likely to have a MA(1-4) impact on water levels. It sounds like Box
(1994), chapter 11 goes into this.


Need to account for SWE changes due to EAB, but don’t have pre/post data. Should
be able to find enough info in literature to create a shifted probability curve
of expected swe with mean/skewness increased appropriately based on opening
size. This could also be used with the map to forecast swe/climate scenarios
based on mapped stand size, which will equate to opening size.

Models need to incorporate maximum water level and snowmelt

## Precip Generators
- CLIGEN
- LARS-WG
- AWE‐GEN‐2d

### Topics to Consider:

##### Hindcasting
  How does our study compare to 30-year normal or other periods?

##### Gas Flux

##### Hydrology
- Water levels
- Response to rain events (paired watersheds)
- Response to snow melt/rain on snow events

##### Vegetation
- Phenology shifts
  * Is black ash (or are replacement species) responsive to photoperiod or climate?
  * Will they respond?
  * How will that affect other things?


##### Ecocasting
  * These forecasts could be tracked with actual data and used to identify areas where more data could be collected to improve the model via data-model loop (Dietze, 2017)

##### Things that could be helpful
meteoland R package
chillr R package
meteo R package
