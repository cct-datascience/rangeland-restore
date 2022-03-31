# Post-fire effects of rangeland restoration treatments on an exotic annual grass (*Bromus tectorum*)

Authors: Elise. S. Gornish, Jessica S. Guo, Lauren M. Porensky, Barry L. Perryman, Elizabeth A. Leger

### Project summary
Restoration of rangeland vegetation faces twin challenges of fire and invasion, which alter ecosystem functioning and provisioning services. Recent efforts to break the feedback between fire and invasion include restoring native plants, particularly slow-growing perennial grasses to reduce fuel loads, compete with invasive species, and tolerate grazing ([Porensky et al. 2018](https://www.ars.usda.gov/ARSUserFiles/49263/29.%20Porensky_etal_2018_Greenstrips.pdf)). To assess the efficacy of such treatments, a split-split plot design of restoration (monoculture vs. mixes, seed coating, two kinds of controls) was implemented following herbicide applications in April 2014 and 2015. These treatments were crossed with episodic disturbance in the form of seasonal grazing (fall 2015 to spring 2017) and burning (July 2017). Our objective is to determine whether targeted grazing and native plant reseeding treatments can have enduring effects on exotic annual grass invasion following fire. What are the enduring effects of: 1) Grazing treatments on the biomass of plant functional groups? 2) Grazing and fuelbreak treatments on cheatgrass (*Bromus tectorum*, BRTE) count? 3) Grazing and fuelbreak treatments on cover of cheatgrass, forbs, and native grasses? 4) Grazing and fuelbreak treatments on average plant height?

### Approach

In contrast to the challenges of transforming plant abundance data to satisfy conditions for frequentist analysis (Damgaard 2009), Bayesian approaches allow for flexible specification of appropriate distributions, including mixture distributions that give rise to zero-inflated data (Dagne 2004). Local factors of plant growth and limited dispersal can also lead to overdispersion, which can be accounted for with hierarchical random effects (Dagne 2004). Here, we conduct four sets of hierarchical Bayesian ANOVA model that account for biomass, count, cover, and height data:
1) The plot-level biomass of BRTE, forbs, and native grasses were modeled with a multivariate log-normal likelihood with the fixed effect of grazing. Plots were nested within the random-effect of block. 
2) Plot-level BRTE counts were modeled as a Poisson distribution with quadrat area as a covariate. Plots were nested within the random-effect of block. 
3) Plot-level BRTE, forbs, and native grass cover were modeled as a zero-inflated beta distribution. Plots were nested within the random-effect of block.
4) Average plant height was modeled with a normal likelihood. Plots were nested within the random-effect of block.

For datasets 2-4, the data were subdivided and investigated at three hierarchical levels. 

I) At the `all` level, all observations were included and assessed for the effects of fuelbreak treatment and grazing. 

II) At the `greenstrip` level, all observations in seeded plots were assessed for the effects of grazing, spatial arrangement, seeding rate, and seed coating. 

III) At the `mono` level, all observations in monoculture seeded plots were assessed for the effects of grazing, species, seeding rate, and seed coating. 

Each level of reference-offset ANOVA model included main effects and two-way interactions. See table below; bolding indicates reference level for each factor. 

| Models      | Factor | N   | Levels |
| ---        |  ---   |     --- |      --- |
| I      | fuelbreak       | 3   | **Control**, herbicide, greenstrip |
| I, II, III      | grazing       | 3   | **Ungrazed**, fall, spring |
| II      | spatial       | 2   | **Monoculture**, mixture |
| II, III      | seed rate       | 2   | **Low**, high |
| II, III      | seed coat       | 2   | **Uncoated**, coated |
| III    | species       | 5   | **ELTR**, POSE, POFE, VUMI, ELEL |

### Repository description
`cleaned_data/` contains 7 .Rdata dataframes produced from `raw_data/` using `scripts/01_organize.R`:
  - `biomass.Rdata` 54 observations of biomass across 3 functional groups
  - `count_all.Rdata` 460 observations of BRTE count among all plots
  - `count_greenstrip.Rdata` 316 observations of BRTE count among seeding plots
  - `count_mono.Rdata` 220 observations of BRTE count among monoculture seeding plots
  - `cover_all.Rdata` 453 observations of cover and height across 3 functional groups among all plots
  - `cover_greenstrip.Rdata` 310 observations of cover and height across 3 functional groups among seeding plots
  - `cover_mono.Rdata` 214 observations of cover and height across 3 functional groups among monoculture seeding plots

`initial_models/` contains exploratory Bayesian models for biomass and count data

`models/` contains 4 subfolders for each data type. In most cases, multiple model forms were explored; only the final JAGS model file is listed below. 
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
  - `cover/` contains 3 subfolders for each functionl group, each containing 3 models:
    - `BRTE/`
      - `all/`
        - `BRTE_cover.R` runs the JAGS model
        - `BRTE_cover_zib.jags` codes the zero-inflated beta model
        - `BRTE_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `greenstrip/`
        - `BRTE_cover.R` runs the JAGS model
        - `BRTE_cover_zib.jags` codes the zero-inflated beta model
        - `BRTE_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `mono/`
        - `BRTE_cover.R` runs the JAGS model
        - `BRTE_cover_zib.jags` codes the zero-inflated beta model
        - `BRTE_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
    - `forbs/`
      - `all/`
        - `Forbs_cover.R` runs the JAGS model
        - `Forbs_cover_zib.jags` codes the zero-inflated beta model
        - `Forbs_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `greenstrip/`
        - `Forbs_cover.R` runs the JAGS model
        - `Forbs_cover_zib.jags` codes the zero-inflated beta model
        - `Forbs_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `mono/`
        - `Forbs_cover.R` runs the JAGS model
        - `Forbs_cover_zib.jags` codes the zero-inflated beta model
        - `Forbs_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
    - `native/`
      - `all/`
        - `Native_cover.R` runs the JAGS model
        - `Native_cover_zib.jags` codes the zero-inflated beta model
        - `Native_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `greenstrip/`
        - `Native_cover.R` runs the JAGS model
        - `Native_cover_zib.jags` codes the zero-inflated beta model
        - `Native_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files
      - `mono/`
        - `Native_cover.R` runs the JAGS model
        - `Native_cover_zib.jags` codes the zero-inflated beta model
        - `Native_plots.R` calculates and plots modeled output
        - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
        - `inits/` contains initials to run the JAGS model as .Rdata files 
        - `plots/` contains plots of model results as .jpg files        
  - `height/` contains 3 models:
    - `all/`
      - `Avg_heights.R` runs the JAGS model
      - `Avg_heights_norm.jags` codes the normal model
      - `Heights_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files        
    - `greenstrip/`
      - `Avg_heights.R` runs the JAGS model
      - `Avg_heights_norm.jags` codes the normal model
      - `Heights_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files  
    - `mono/`
      - `Avg_heights.R` runs the JAGS model
      - `Avg_heights_norm.jags` codes the normal model
      - `Heights_plots.R` calculates and plots modeled output
      - `coda/` contains posterior chains for parameters and replicated data, as .Rdata files 
      - `inits/` contains initials to run the JAGS model as .Rdata files 
      - `plots/` contains plots of model results as .jpg files  
  

`raw_data/` from [Lauren Porensky](https://www.ars.usda.gov/plains-area/fort-collins-co/center-for-agricultural-resources-research/rangeland-resources-systems-research/people/lauren-porensky/):
  - `2019_Greenstrips_Biomass.csv` contains 54 observations of the biomass of 3 functional groups (BRTE, forbs, native grass) at the non-seeded, grazed plots. 
  - `2019_Greenstrips_Cover_and_Densities.csv` contains 460 observations of BRTE count, cover of multiple functional groups, and average plant height across all treatment combinations. 

`scripts/` contains R scripts for organizing, plotting, and analyzing simple statistics:
  - `01_organize.R` organizes the `raw_data/` files into `cleaned_data/` files
  - `02_biomass_figures.R` plots raw data and modeled posteriors for biomass
  - `03_count_figures.R` plots raw data and modeled posteriors for BRTE count
  - `04_cover_figures_BRTE.R` plots raw data and modeled posteriors for BRTE cover
  - `05_cover_stats.R` calculates a sample statistics of cover data
  - `06_cover_figures_forbs.R` plots raw data and modeled posteriors for forb cover
  - `07_cover_figures_native.R` plots raw data and modeled posteriors for native grass cover
  - `08_height_figures.R` plots raw data and modeled posteriors for average plant height
  
