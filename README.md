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
	* CoreMethanotrophBiome.r - Code to look at the core metanotroph microbiome
	* DistanceDecay.r - This code examines the distance decay pattern of methanotrophs
	* MethanotrophRatios.r - This code examines the rates of methanotrophs in the cave communities
	* MethanotrophVCH4Concentration.r - This code examines how methanotroph relative abundance changes with CH4 concentration
	* MethanotrophVDistanceEntrance.r - Code for examining how methanotroph relative abundance changes with distance to a cave entrance
	* MetahnotrophVSoilGrainSize.r - Code for examing how methanotroph relative abundance changes with soil grain size
	* SoilCH4CCA - Code for cannonical correspondance analysis of cave methanotroph communities

## Contributors

[Dr. Kevin Webster] (http://websterkgd.com/): Associate Research Scientist, [Planetary Science Institute] (https://www.psi.edu/about/staffpage/webster)

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Principle Investigator, Associate Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).
