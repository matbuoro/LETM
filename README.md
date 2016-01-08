Latent Environmental Threshold Model
====================================

Conditional strategies are the most common form of discrete phenotypic plasticity. In a conditional strategy, the phenotype expressed by an organism is determined by the difference between an environmental cue and a threshold, both of which may vary among individuals. The Environmental Threshold (ET) model has been proposed as a mean to understand the evolution of conditional strategies, but has been surprisingly seldom applied to empirical studies. A hindrance for the application of the ET model is that often, the proximate cue triggering the phenotypic expression and the individual threshold are not measurable, and can only be assessed using a related observable cue.   

We describe a new statistical model that can be applied in this common situation (see Buoro et al. 2012). The Latent ET model (LETM) allows for a measurement error in the phenotypic expression of the individual environmental cue and a purely genetically determined threshold. We show that coupling our model with quantitative genetic methods allows an evolutionary approach including an estimation of the heritability of conditional strategies.

*Update*

We assumed that the individual thresholds covary according to the individual relatedness. The individual additive genetic effects of the threshold , a_i covary and the structure of the covariance matrix depends on the relatedness between individuals. The vector of the a_i, is multivariate normal with mean 0 and variance–covariance matrix A, where  A  is the additive genetic relationship matrix. The additive genetic relationship matrix A is built up from the pedigree. 

To handle with the complexity of the animal model, we can also use an alternative procedure modeling the variation in additive genetic values (or breeding values), which allow the use of the pedigree (see code of the model and "Earwigs" folder as an example). Then, the analysis is much faster!


Overview
--------

This folder contains:

* Simulation script : just to illustrate how the LETM is working (using R);  
	|-  SCRIPT_LETM_SIM.R : R script to simulate data and run the analysis;  
	|-  LETM.md : code of the model  
* Examples: codes and data from the *Salmo* study (Buoro et al. 2012) and *Earwigs* study (Buzatto et al. 2015). 



Literature
--------

Corey C. P., J. W. Moore, **M. Buoro**, S. A. Hayes, J.C. Garza and D. E. Pearse, 2016. Shifting thresholds: rapid evolution of migratory life histories in steelhead/rainbow trout, Oncorhynchus mykiss. *Journal of Heredity*. [Link](http://jhered.oxfordjournals.org/content/107/1/51.abstract)  

B. A. Buzatto^§^, **M. Buoro^§^**, W. N. Hazel and J.L. Tomkins, 2015. Investigating the genetic architecture of conditional strategies using the Environmental Threshold Model. *Proceedings of the Royal Society B.*,  [Link](http://rspb.royalsocietypublishing.org/content/282/1821/20152075) . (§ These authors contributed equally to the work)

**Buoro M.**, Gimenez O., and E. Prévost, 2012. Assessing adaptive phenotypic plasticity by means of conditional strategies from empirical data: the latent environmental threshold model. *Evolution*, 66(4):996-1009. [Link](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2011.01484.x/abstract)  




