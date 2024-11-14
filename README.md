This reposetory contains 
	(a) sourcecode that allows the user to determine the amounts of Pfr (15E) and Pr (15Z) in a microbial phytochrome sample using the data obtained via UV/Vis absorption spectroscopy 
 	(b) files (.F3D, .OBJ, .STL) with the 3D model used in Huber et al. (2024) to illuminate an Eppendorf tube from below.

The complex, multi-step photocycle of phytochromes can be simplified into a system of two species (Pfr and Pr) and three rates: the rate of product formation for both species, and the mono-directional rate of dark reversion. Is the illumination intensity high enough, the rate of dark reversion becomes negligible in comparison to the photochemical rates; this is the desired experimental condition. A photochemical equilibrium, or photostate, can then be understood as a mixture of Pfr and Pr. 

To describe experimental data, we have first created individual functions for the parent states (a sum of Gaussians). These are then implemented into the fit function proper, which considers the Q-band absorption to be a superposition of the parent states. The sourcecode calculates the value of alpha, the amount of Pfr in the sample, and thus characterises the photostate.

The method has been published previously:

	Darkness inhibits autokinase activity of bacterial bathy phytochromes
	Huber, Christina et al.
	Journal of Biological Chemistry, Volume 300, Issue 4, 107148
	DOI:  10.1016/j.jbc.2024.107148

A step by step guide is additionally provided in:

	WINKLER CHAPTER

The following phytochromes are already included in the sourcecode.

WT = wildtype

- Bathy:
  - Pseudomonas aeruginosa WT
  - Agrobacterium tumefaciens 2 WT
  - Agrobacterium tumefaciens D783N variant
  - Xanthomonas campestris pv. campestris WT
  - Xanthomonas campestris pv. campestris without PAS9 domain
  - Ramlibacter tataouinensis phytochrome 2 WT
  - Agrobacterium vitis WT

- Canonical:
  - Agrobacterium tumefaciens 1 
  - Pseudomonas syringae 1
