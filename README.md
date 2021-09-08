# Post-fire effects of rangeland restoration treatments on an exotic annual grass (Bromus tectorum)

### Project summary
Restoration of rangeland vegetation faces twin challenges of fire and invasion, which alter ecosystem functioning and provisioning services. Recent efforts to break the feedback between fire and invasion include planting 'greenstrips' of native plants, particularly slow-growing perennial grasses to reduce fuel loads, compete with invasive species, and tolerate grazing (Porensky et al. 2018). To assess the efficacy of such treatments, a split-split plot design of restoration (monoculture vs. mixes, seed coating, two kinds of controls) was implemented following herbicide applications in April 2014 and 2015. These treatments were crossed with episodic disturbance in the form of seasonal grazing (fall 2015 to spring 2017) and burning (July 2017). Our objective is to determine whether targeted grazing and native greenstrips treatments can have enduring effects on exotic annual grass invasion following fire. What are the legacy effects of: 1) Grazing treatments on the biomass of plant functional groups? 2) Grazing and fuelbreak treamtents on the count and cover of cheatgrass (Bromus tectorum, BRTE)?

### Approach

In contrast to the challenges of transforming plant abundance data to satisfy conditions for frequentist analysis (Damgaard 2009), Bayesian approaches allow for flexible specification of appropriate distributions, including mixture distributions that give rise to zero-inflated data (Dagne 2004). Local factors of plant growth and limited dispersal can also lead to overdispersion, which can be accounted for with hierarchical random effects (Dagne 2004). Here, we propose three hierarchical Bayesian ANOVA model that account for biomass, cover, and count data:
1) The plot-level biomass of BRTE, forbs, and native grasses were modeled with a multivariate log-normal likelihood with the fixed effect of grazing. Plots were nested within the random-effect of block. 
2) Plot-level BRTE counts were modeled as a Poisson distribution with quadrat area as a covariate. Plots were nested within the random-effect of block. 
3) Plot-level BRTE cover were modeled as a zero-one inflated beta distribution. 
Plots were nested within the random-effect of block.

For datasets 2 & 3, the data were further subdivided and three different ANOVA models were run (see table below). Main effects and all two-way interactions were included. 

| Models      | Factor | N   | Levels |
| ---        |  ---   |     --- |      --- |
| I      | fuelbreak       | 3   | **Control**, herbicide, greenstrip |
| I, II, III      | grazing       | 3   | **Ungrazed**, fall, spring |
| II      | spatial       | 2   | **Monoculture**, mixture |
| II, III      | seedrate       | 2   | **Low**, high |
| II, III      | seed coat       | 2   | **Uncoated**, coated |
| III    | species       | 5   | **ELTR**, POSE, POFE, VUMI, ELEL |

### Repository description
`cleaned_data/` contains 7 output csv files from `scripts/01_organize.R`, including 1 biomass file and 3 files for the 3 levels of count and cover data

`initial_models/` contains exploratory Bayesian models for biomass and count data

`models/` contains 3 subfolders for each data type:
  - `biomass/` contains 1 model:
    - `all/`
      - `BRTE_biomass.R` runs the JAGS model
      - `BRTE_biomass_mvnorm.jags` codes the multivariate log-normal model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
  - `count/` contains 3 models:
    - `all/`
      - `BRTE_counts.R` runs the JAGS model
      - `BRTE_counts_Poisson.jags` codes the Poisson model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
    - `greenstrip/`
      - `BRTE_counts.R` runs the JAGS model
      - `BRTE_counts_Poisson.jags` codes the Poisson model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
    - `mono/`
      - `BRTE_counts.R` runs the JAGS model
      - `BRTE_counts_Poisson.jags` codes the Poisson model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
  - `cover/` contains 3 models:
    - `all/`
      - `BRTE_cover.R` runs the JAGS model
      - `BRTE_cover_zoib.jags` codes the zero-one inflated beta model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
    - `greenstrip/`
      - `BRTE_cover.R` runs the JAGS model
      - `BRTE_cover_zoib.jags` codes the zero-one inflated beta model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
    - `mono/`
      - `BRTE_cover.R` runs the JAGS model
      - `BRTE_cover_zoib.jags` codes the zero-one inflated beta model
      - `BRTE_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files
  

`raw_data/` contains 2 csv files from Lauren Porensky descriting the biomass and cover/count data. 

`scripts/` contains R scripts.
  - `01_organize.R` organizes the `raw_data/` files into `cleaned_data/` files
