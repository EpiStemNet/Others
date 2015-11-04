# Others

# ChromHMM analysis

contact: ecarrillo@cnio.es Enrique Carrillo de Santa Pau  Structural Biology and BioComputing Programme, Spanish National Cancer Research Center - CNIO, Melchor Fernandez Almagro 3 28029 Madrid, Spain.

###### ChomHMM software v1.03

##Files used:
- mm9.txt --> Included in the ChromHMM package in folder CHROMSIZES
- mESC_chromHmm_CTCF.txt --> Provided (see ChromHMM manual for details http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf)
- mESC_chromMod_ChromHmm.txt --> Provided (see ChromHMM manual for details http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf)
- ESC_20_mESC_all_chrs_0-95_posterior.txt.gz -> Provided (see ChromHMM manual for details http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf)
- BAM files with alignments --> Provided in http://ubio.bioinfo.cnio.es/data/mESC_CNIO/bam_files/


1) BAM files were converted to BED files

We used bedtools v2.18.2 with the function bedtobam

2) BINARIZE HISTONE, CYTOSINE MODIFICATIONS AND CTCF DATA

java -mx20G -jar ChromHMM.jar BinarizeBed  -c ~/Documents/Histones/Inputs/ ~/.lib/ChromHMM/CHROMSIZES/mm9.txt ~/Documents/Histones/ChromHmm/ ~/Documents/Histones/mESC_chromHmm_CTCF.txt ~/Documents/Histones/ChromHmm/BINARIZED_CTCF

3) BINARIZE CHROMATIN MODIFIERS DATA

java -mx20G -jar ChromHMM.jar BinarizeBed  -c ~/Documents/Chrom_modifiers/Inputs/ ~/.lib/ChromHMM/CHROMSIZES/mm9.txt ~/Documents/Chrom_modifiers/ ~/Documents/Chrom_modifiers/mESC_chromMod_ChromHmm.txt ~/Documents/Chrom_modifiers/BINARIZED

4) LEARNMODEL WITH HISTONE, CYTOSINE MODIFICATIONS AND CTCF DATA

java -mx20G -jar ChromHMM.jar LearnModel -i mESC -l CHROMSIZES/mm9.txt -printposterior -printstatebyline  ~/Documents/Histones/ChromHmm/BINARIZED_CTCF ~/Documents/Histones/ChromHmm/RESULTS/Learnmodel_20_CTCF 20 mm9

### Only the intervals with a probability higher than 0.95 were considered for further analysis.

5) ENRICHMENT PROBABILITIES FOR HISTONE, CYTOSINE MODIFICATIONS AND CTCF DATA

java -mx20G -jar ChromHMM.jar OverlapEnrichment ~/Documents/Histones/ChromHmm/RESULTS/Learnmodel_20_CTCF/ESC_20_mESC_all_chrs_0-95_posterior.txt ~/Documents/Histones/ChromHmm/BINARIZED_CTCF/peaks ~/Documents/Histones/ChromHmm/RESULTS/Learnmodel_20_CTCF/enrichment_20_CTCF_0_95_histones

6) ENRICHMENT PROBABILITIES FOR CHROMATIN MODIFIERS DATA

java -mx20G -jar ChromHMM.jar OverlapEnrichment ~/Documents/Histones/ChromHmm/RESULTS/Learnmodel_20_CT

# ChromNets analysis

contact: dadejuan@cnio.es David de Juan. Structural Biology and BioComputing Programme, Spanish National Cancer Research Center - CNIO, Melchor Fernandez Almagro 3 28029 Madrid, Spain.

##Files used:
-   get_chromnets.R --> Provided in https://github.com/ChromatinNetwork/Others/blob/master/get_chromnets.R
-   ColocationNetwork.txt --> Provided in https://github.com/ChromatinNetwork/Others/blob/master/ColocationNetwork.txt
- 	CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz --> Provided in http://ubio.bioinfo.cnio.es/data/mESC_CNIO/others/

## Required software
-   R version 3.1.2 (or later) from https://www.r-project.org/
-   data.table R package version 1.9.4 (or later) from http://CRAN.R-project.org/package=data.table
-   plyr R package from http://www.jstatsoft.org/v40/i01/
-   pvclust R package version 1.3-2 (or later) from http://CRAN.R-project.org/package=pvclust

### Extract ChromNets from a co-location network and genome-wide location data for every epigenetic feature and chromatin states

Rscript get_chromnets.R ColocationNetwork.txt CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt \<output_file\>

or

Rscript get_chromnets.R ColocationNetwork.txt CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt \<output_file\> \<pvalue\>
