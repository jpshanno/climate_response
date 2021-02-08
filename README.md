## [Wetland Disturbance and Climate Change: Are restoration targets achievable under future climate scenarios?](./Climate_Change_Impacts/README.md)

<!-- The final target for restoration may or may not be current conditions. If future
climate conditions will not support the same wetland conditions/functions then
planning to restore them is a hopeless endeavor. -->

## Hierarchical Model

*WL<sub>t<\sub> = WL<sub>t-1<\sub> + P<sub>t-1<\sub> + ET<sub>t-1<\sub> + Q<sub>t-1<\sub>  + G<sub>t-1<\sub> *

_Precip response scaled by ecosystem specific yield__
nlf(p ~ step(wlL1 - cp) * (precip / esy) + step(cp - wlL1) * (Mp * precip)),

_PET response scaled by ecosystem specific yield and 0 when annual water__
_availability (YTD P - YTD PET) is above a threshold value__
nlf(et ~ step(cpWA - wa) * (Mpet * pet / esy)),

_Streamflow and subsurface flow losses above a water level threshold__
nlf(q ~ step(wlL1 - cp) * (Bq + Mq * wlL1^2)),

_Local subsurface inputs when annual water availability is above threshold__
nlf(g ~ step(cp - cpWA) * Bg),

_Ecosystem Specific Yield__
nlf(esy ~ Besy + Mesy ^ (step(cp - wlL1) * (wlL1 - cp)))

Are we setting achievable wetland restoration goals in light of changing climate?

#### Technical Questions
- Are hydrologic changes expected to be mitigated or amplified by climate change?

- How do control and treatment watersheds respond to projected future climate
  conditions?
- How does cliamte change impact compare to EAB impact?
- What is the cumulative impact of climate change and EAB?
- Are current goals of forested wetlands achievable?
  - Are there species that can tolerate the hydroperiod in a climate analogous setting?
- Will this lead to these sites being refugia during increased water stress, or
  will continued increased water levels mean they will be maintained as
  primarily herbaceous or shrub/scrub wetlands?
- How do these results compare with wetlands in climate-analogous locations?

#### Challenges
- How do we evaluate if we are at a new steady-state condition for wetland
  hydrology in our treatment watersheds? It may more time for a hydrologic
  regime to return to an equilibrium state following the transition from
  forested wetalnd. Comparison to other non-forested wetlands could help if data
  are available.
- How do I explain the resutls outside of specific sites? Perhaps presenting
  responses across a range of site hydrologies will help, and applying my models
  to MN data may illustrate commonalities.
- Has recovery already begun?
- Wetland recovery will be driven by both vegetative recovery and climate
  change. Vegetative response will be influenced by climate change, and the
  response to climate change will be influenced by vegetation response speed and
  'direction'
- Can current wetland function be used to evaluate future wetland recovery?

#### Limitations
- Cannot/will not be evaluating other important wetland responses to climate
  change such as shifts in flora (climate-induced or continued evolutino of EAB
  impacts), and impact of climate change on upland groundwater sources.

# Model Notes


Using precip to model water level should only require one coefficient for precip
(with a lag to-be-determined (probably 0 or 1)), an MA process (which looks
likely for wetland water level data)

# Incorporating Snow Data
- Need to account for SWE changes due to EAB, but don’t have pre/post data. Should
be able to find enough info in literature to create a shifted probability curve
of expected swe with mean/skewness increased appropriately based on opening
size. This could also be used with the map to forecast swe/climate scenarios
based on mapped stand size, which will equate to opening size.


# Currently Unaccounted for Model Processes
- Frost drought from lack of snow cover but frozen soils? Will it be
likely? Probably not by the end of the century

- Rain-on-snow events and decreased snowpack reducing the amount of upland water 
level recharge

- It may be important to do a midcentury and end of century model since the
impact of lake effect will vary between those periods. It may be that the signal
of winter snowfall is lost relatively early or not. Understanding the role of
snowfall on wetland water levels will be important to decide what impact these
shifts will have. 

### Impacted Topics

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

How do other projectsions compare to Byun and Hayhoe?
 Study           | Scenario | T          | P       | Notes  
 :---            |  :---:   |    :---:   |  :---:  | :---  
 Wang (2017)     | SRES A1B | +5.6       | +0-15%  | Annual only; probabilistic
 Notaro (2015)   | RCP 8.5  | +6 (4 - 9) | +10-15% | T is winter only
 Kutzbach (2005) | A2       | +4.4       | +11%    | Annual; seasonal trends match Byun & Hayhoe
 Kutzbach (2005) | B1       | +2.7       | +6%     | Annual; seasonal trends match Byun & Hayhoe
 Kling (2003)    | High     | +3-5       | +10-20% | Annual; seasonal matches expectations
 Kling (2003)    | Low      | +5-8       | +10-20% | Annual; seasonal matches expectations

What is the role of climate change shifting winter precip away from snow? If we 
have wetter winters & springs as rain, does that change how long wetland water
levels remain elevated? There is a max amount of surface storage before spill
out occurs. But the snow signal comes later than the rain signal and provides fully
saturated uplands and groundwater at the start of the growings eason. There may
be significant ET prior to when snowmelt generally occurs so that the starting
conditions are equivalent and the draw down is similar in scale, but there is a
shift in timing. This could be amplified if the wetlands continue to green-up at
a later date than the surrounding uplands (consider Fisher's [-@fisher-2006]
suggestion that leaf-on is inversely correlated with landscape depressions)