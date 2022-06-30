# CaveMethanotrophs
===================

An analysis of cave methanotrophic communities from environmental 16S rRNA dna sequece data from caves in the US.
This repository contains open-source code, data, and text files to examine how methanotrophic communities 
vary in caves.

### Repo Contents

* **CaveMethanotrophData:** - Data files used in this research
	* CvMcrb03_20190225.shared - A shared file describing the abundance of an OTU at a particular study site
	* CvMcrb03_20190225.taxonomy - A taxonomy file of the identified sequences
	* MethanotrophSearch.fasta - A FASTA file of methanotroph 16S rRNA sequences from the literature
	* MinAbnMatr.txt - A mineral abundance matrix by site
	* NCBIPutMntrphIDCave_20190616.txt - Matches between OTUs from the shared file and known methanotrophs 
	* sample_meta_data+d.txt - A txt file of the collected meta data 

* **CaveMethanotrophRcode:** - Analytical and figure generating code
	* DistDecay_20220617.r - This code examines the distance decay pattern of methanotrophs
	* MethanoExploration_20220616.r - Code to look at the alpha and beta diversity, and the core metanotroph microbiome of caves
	* MethanoMap_20220617.r - This code generates a spatial map of cave methanotroph abundance
	* MethanoMineralCCA_20220617 - Code for cannonical correspondance analysis of cave methanotroph communities against mineral abundance
	* MethanoVsCHConc_20220617.r - This code examines how methanotroph relative abundance changes with CH4 concentration
	* MethanoVsDistEnt_20220620.r - Code for examining how methanotroph relative abundance changes with distance to a cave entrance
	* MetahnoVsSoilGrainSize_20220618.r - Code for examing how methanotroph relative abundance changes with soil grain size
	* Slopetest_20220623.r - Code containing fuctions for distance decay analysis  
	

## Contributors

[Dr. Kevin Webster](https://websterkgd.com/): Associate Research Scientist, [Planetary Science Institute](https://www.psi.edu/about/staffpage/webster).

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Associate Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).
