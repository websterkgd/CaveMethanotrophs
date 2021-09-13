# CaveMethanotrophs
===================

An analysis of cave methanotrophic communities from environmental 16S rRNA dna sequece data from caves in the US.
This repository contains open-source code, data, and text files to examine how methanotrophic communities 
vary in caves.

### Repo Contents

* **CaveMethanotrophData:** - Data files used in this research
	* CvMcrb03.bac.final.0.03.taxonomy - A taxonomy file of the identified sequences
	* CvMcrb03.bac.final.shared - A shared file describing the abundance of an OTU at a particular study site
	* sample_meta_data+d.txt - A txt file of the collected meta data 

* **CaveMethanotrophRcode:** - Analytical and figure generating code
	* DistDecay.r - This code examines the distance decay pattern of methanotrophs
	* MethanoCommunityMeansCore.r - Code to look at the alpha and beta diversity, and the core metanotroph microbiome of caves
	* MethanoMap.r - This code generates a spatial map of cave methanotroph abundance
	* MethanoMineralCCA - Code for cannonical correspondance analysis of cave methanotroph communities against mineral abundance
	* MethanoVsCHConc.r - This code examines how methanotroph relative abundance changes with CH4 concentration
	* MethanoVsDistEnt.r - Code for examining how methanotroph relative abundance changes with distance to a cave entrance
	* MetahnoVsSoilGrainSize.r - Code for examing how methanotroph relative abundance changes with soil grain size  
	

## Contributors

[Dr. Kevin Webster](https://websterkgd.com/): Associate Research Scientist, [Planetary Science Institute](https://www.psi.edu/about/staffpage/webster).

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Associate Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).
