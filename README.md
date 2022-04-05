# Post-fire effects of rangeland restoration treatments on an exotic annual grass (*Bromus tectorum*)

Authors: Elise. S. Gornish, Jessica S. Guo, Lauren M. Porensky, Barry L. Perryman, Elizabeth A. Leger

For further questions, please contact [Elise](mailto:egornish@arizona.edu) about the manuscript, [Jessica](mailto:jessicaguo@email.arizona.edu) about the analyses, and [Lauren](mailto:lauren.porensky@usda.gov) about the data. 

## Project summary
Restoration of rangeland vegetation faces twin challenges of fire and invasion, which alter ecosystem functioning and provisioning services. Recent efforts to break the feedback between fire and invasion include restoring native plants, particularly slow-growing perennial grasses to reduce fuel loads, compete with invasive species, and tolerate grazing ([Porensky et al. 2018](#1)). To assess the efficacy of such treatments, a split-split plot design (seeding, herbicide, control) was implemented following initial herbicide applications in April 2014. These treatments were crossed with episodic disturbance in the form of seasonal grazing (fall 2015 to spring 2017), followed by a naturally-occurring fire that burned all treatments (July 2017). Our objective is to determine whether targeted grazing and native plant reseeding treatments can have enduring effects on exotic annual grass invasion following fire. What are the enduring effects of: 1) Grazing treatments on the biomass of plant functional groups? 2) Grazing and seeding treatments on cheatgrass (*Bromus tectorum*, BRTE) count? 3) Grazing and seeding treatments on cover of cheatgrass, forbs, and native grasses? 4) Grazing and seeding treatments on average plant height?

## Approach

In contrast to the challenges of transforming plant abundance data to satisfy conditions for frequentist analysis ([Damgaard 2009](#2)), Bayesian approaches allow for flexible specification of appropriate distributions, including mixture distributions that give rise to zero-inflated data ([Dagne 2004](#3)). Local factors of plant growth and limited dispersal can also lead to overdispersion, which can be accounted for with hierarchical random effects ([Dagne 2004](#3)). Here, we conduct four sets of hierarchical Bayesian ANOVA model that account for biomass, count, cover, and height data:
1) The plot-level biomass of BRTE, forbs, and native grasses were modeled with a multivariate log-normal likelihood with the fixed effect of grazing. Plots were nested within the random-effect of block. 
2) Plot-level BRTE counts were modeled as a Poisson distribution with quadrat area as a covariate. Plots were nested within the random-effect of block. 
3) Plot-level BRTE, forbs, and native grass cover were modeled as a zero-inflated beta distribution ([Bayes & Valdivieso 2016](#4)). Plots were nested within the random-effect of block.
4) Average plant height was modeled with a normal likelihood. Plots were nested within the random-effect of block.

For datasets 2-4, the data were subdivided and investigated at three hierarchical levels. 

I) At the `all` level, all observations were included and assessed for the effects of seeding treatment and grazing. 

II) At the `greenstrip` level, all observations in seeded plots were assessed for the effects of grazing, spatial arrangement, seeding rate, and seed coating. 

III) At the `mono` level, all observations in monoculture seeded plots were assessed for the effects of grazing, species, seeding rate, and seed coating. 

Each level of reference-offset ANOVA model included main effects and two-way interactions. See table below; bolding indicates reference level for each factor. 

| Models      | Factor | N   | Levels |
| ---        |  ---   |     --- |      --- |
| I      | seeding       | 3   | **Control**, herbicide, greenstrip |
| I, II, III      | grazing       | 3   | **Ungrazed**, fall, spring |
| II      | spatial       | 2   | **Monoculture**, mixture |
| II, III      | seed rate       | 2   | **Low**, high |
| II, III      | seed coat       | 2   | **Uncoated**, coated |
| III    | species       | 5   | **ELTR**, POSE, POFE, VUMI, ELEL |

## Repository description

### Data 

`raw_data/` originated by [Lauren Porensky](https://www.ars.usda.gov/plains-area/fort-collins-co/center-for-agricultural-resources-research/rangeland-resources-systems-research/people/lauren-porensky/):
  - `2019_Greenstrips_Biomass.csv` contains 54 observations of the biomass of 3 functional groups (BRTE, forbs, native grass) at the non-seeded, grazed plots. 
  - `2019_Greenstrips_Cover_and_Densities.csv` contains 460 observations of BRTE count, cover of multiple functional groups, and average plant height across all treatment combinations. 
  
`cleaned_data/` contains 7 tabular dataframes (.Rdata) produced from `raw_data/` using `scripts/01_organize.R`:
  - `biomass.Rdata` 54 observations of biomass across 3 functional groups
    - 'block', 'paddock, 'grazing', and 'control_id' are covariates
    - 'BRTE_abs', 'forbs_abs', 'grass_abs', and 'standingdead_abs' are the absolute dry biomass (g) per quadrat
    - 'id' is unique identifier for each observation
    - 'totveg' and 'liveveg' are the total biomass for all and for living classes
    - 'BRTE_rel_tot' and 'BRTE_rel_live' are the relative proportion of BRTE biomass with respect to total and living biomass 
  - `count_all.Rdata` 460 observations of BRTE count among all plots
    - 'block', 'paddock, 'grazing', and 'fuelbreak' are covariates
    - 'BRTE' is the count of cheatgrass stems
    - 'quadrat' is the area ($cm^2$) used for counting cheatgrass stems
  - `count_greenstrip.Rdata` 316 observations of BRTE count among seeding plots
    - 'block', 'paddock, 'grazing', 'spatial', 'seed_rate', 'seed_coat' are covariates
    - 'BRTE' is the count of cheatgrass stems
    - 'quadrat' is the area ($cm^2$) used for counting cheatgrass stems
  - `count_mono.Rdata` 220 observations of BRTE count among monoculture seeding plots
    - 'block', 'paddock, 'grazing', 'species', 'seed_rate', 'seed_coat' are covariates
    - 'BRTE' is the count of cheatgrass stems
    - 'quadrat' is the area ($cm^2$) used for counting cheatgrass stems
  - `cover_all.Rdata` 453 observations of cover and height across 3 functional groups among all plots
    - 'block', 'paddock, 'grazing', and 'fuelbreak' are covariates
    - 'BRTE', 'intro_forbs', 'native_grass', 'native_forbs', 'bryo', 'litter', 'rock', and 'bare' are the cover (%) of each class
    - 'height' is average plant height (cm)
    - 'totveg_all' and 'totveg_live' are calculated total cover across all or living classes
  - `cover_greenstrip.Rdata` 310 observations of cover and height across 3 functional groups among seeding plots
    - 'block', 'paddock, 'grazing', 'spatial', 'seed_rate', and 'seed_coat' are covariates
    - 'BRTE', 'intro_forbs', 'native_grass', 'native_forbs', 'bryo', 'litter', 'rock', and 'bare' are the cover (%) of each class
    - 'height' is average plant height (cm)
    - 'totveg_all' and 'totveg_live' are calculated total cover across all or living classes
  - `cover_mono.Rdata` 214 observations of cover and height across 3 functional groups among monoculture seeding plots
  - 'block', 'paddock, 'grazing', 'species', 'seed_rate', and 'seed_coat' are covariates
    - 'BRTE', 'intro_forbs', 'native_grass', 'native_forbs', 'bryo', 'litter', 'rock', and 'bare' are the cover (%) of each class
    - 'height' is average plant height (cm)
    - 'totveg_all' and 'totveg_live' are calculated total cover across all or living classes
  
  
### Models

`initial_models/` contains exploratory Bayesian models for biomass and count data

`models/` contains a nested file structure:
  - `biomass/`
    - `all/`
  - `count/`
    - `all/`
    - `greenstrip/`
    - `mono/`
  - `cover/` 
    - `BRTE/`
      - `all/`
      - `greenstrip/`
      - `mono/`
    - `forbs/`
      - `all/`
      - `greenstrip/`
      - `mono/`
    - `native/`
      - `all/`
      - `greenstrip/`
      - `mono/`
  - `height/`
    - `all/`
    - `greenstrip/`
    - `mono/`

For each of the above folders, there are:
- A control script (.R) to run the model code 
- At least one model file (.jags) that encodes the Bayesian model
- A plotting script (.R) to visualize model output
- `coda/` with posterior chains (.Rdata) for modeled parameters and replicated data
- `inits/` with initial values (.Rdata) to start the JAGS model
- `plots/` with figures (.jpg) of model results

### Scripts and plots

`scripts/` contains R scripts for organizing, plotting, and analyzing simple statistics:
  - `01_organize.R` organizes the `raw_data/` files into `cleaned_data/` files
  - `02_biomass_figures.R` plots raw data and modeled posteriors for biomass
  - `03_count_figures.R` plots raw data and modeled posteriors for BRTE count
  - `04_cover_figures_BRTE.R` plots raw data and modeled posteriors for BRTE cover
  - `05_cover_stats.R` calculates a sample statistics of cover data
  - `06_cover_figures_forbs.R` plots raw data and modeled posteriors for forb cover
  - `07_cover_figures_native.R` plots raw data and modeled posteriors for native grass cover
  - `08_height_figures.R` plots raw data and modeled posteriors for average plant height

`plots/` contains the figures (.jpg)  output from `scripts/` plotting scripts

## Citations
<a id="1">[1]</a> 
Porensky, LM, Perryman, BL, Williamson, MA, Madsen, MD, Leger, EA. (2018). Combining active restoration and targeted grazing to establish native plants and reduce fuel loads in invaded ecosystems. Ecol Evol. 8: 12533-12546. https://doi.org/10.1002/ece3.4642

<a id="2">[2]</a> 
Damgaard, C. (2019). On the distribution of plant abundance data. Ecol Inform. 4(2): 76-82. https://doi.org/10.1016/j.ecoinf.2009.02.002

<a id="3">[3]</a> 
Dagne, GA. (2004). Hierarchical Bayesian analysis of correlated zero-inflated count data. Biom J. 46: 653-663. https://doi.org/10.1002/bimj.200310077

<a id="4">[4]</a> 
Bayes, CL, Valdivieso, L. (2016). A beta inflated mean regression model for fractional response variables. J Appl Stat. 43(10): 1814-1830. https://doi.org/10.1080/02664763.2015.1120711