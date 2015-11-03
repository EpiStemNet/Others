library(data.table)
library(plyr)
library(pvclust)
#Read commandline arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 4){
	infile_coloc <- args[1]
	infile_peaks <- args[2]
	outfile <- args[3]
	
	chromnet_pvalue <- args[4]
}
if(length(args) == 3)
{
	infile_coloc <- args[1]
	infile_peaks <- args[2]
	outfile <- args[3]
	chromnet_pvalue = 0.05
}
if((length(args) != 4 & length(args) != 3) | chromnet_pvalue>1 | chromnet_pvalue<0 | !is.numeric(chromnet_pvalue) | !file.exists(infile_coloc) | !file.exists(infile_peaks))
	{
		write("\nError: Wrong number of of input arguments\n", stderr())
		write("get_chromnets.R is a R script to extract chromnets (clusters of interaction) according to the frequencies of colocalization of pairs of genomic features in different chromatin states along the genome\n", stderr())
		write("Usage:R --slave --args 'network_file' 'peaks_file' 'output_file' 'chromnet_pvalue' < /path/to/get_chromnets.R", stderr())
		write("Alternative: /path/to/Rscript /path/to/get_chromnets.R 'network_file' 'peaks_file' 'output_file' 'chromnet_pvalue'\n", stderr())
		write("\t\t'network_file'\t\tfile containing the colocalization network as retreived from XXX (header is required)", stderr())
		write("\t\t'peaks_file'\t\tfile containing the peak calling and chromatin states segmentation of the genome  (header is required)", stderr())
		write("\t\t'output_file'\t\toutput file where resulting chromnets are written", stderr())
		write("\t\t'chromnet_pvalue'\tp_value threashold for chromnets defintion (default value = 0.05)\n\n", stderr())
		quit()
	}

#Read Colocation network
coloc<-fread(infile_coloc)

#Get positive co-localization in any state
positive_coloc<-unique(data.table("HM"=coloc[coloc$SPCN_B==1]$HM,"CM"=coloc[coloc$SPCN_B==1]$CM))
size_positive_coloc <- dim(positive_coloc)[1]

HMs<-sort(unique(coloc$HM))
number_hms=length(HMs)
CMs<-sort(unique(coloc$CM))
number_cms=length(CMs)


#Read peaks and genome segmentation
peaks<-fread(infile_peaks)
size_peaks <- dim(peaks)[1]
for(i in 1:number_hms)
{
	
	if(!HMs[i] %in% colnames(peaks))
	{
		write("Missing genomic features in peaks file\n",stderr())
		quit();
	}
}
for(i in 1:number_cms)
{
	if(!CMs[i] %in% colnames(peaks))
	{
		write("Missing genomic features in peaks file\n",stderr())
		quit();
	}
}

#Get State labels
states<-sort(unique(coloc$STATE))
int_freqs<-data.table("state"=states)

#Count number of windows per state
states_counts<-count(peaks,'state')
empty_states<-states_counts$freq*0

#Loop for the colocalizing pairs
for(i in 1:size_positive_coloc)
{
	#Get label for colocalizing pair
	int_label<-paste(positive_coloc$HM[i],positive_coloc$CM[i],sep=".")
	
	#Get number of windows for every combination of presence/absence of the two features in the pair
	comb_counts<-count(peaks,c(positive_coloc$HM[i],positive_coloc$CM[i],'state'))
	#freqs<-comb_counts[comb_counts[,positive_coloc$HM[i]]==1 & comb_counts[,positive_coloc$CM[i]]==1,]$freq
	freqs<-comb_counts[comb_counts[,positive_coloc$HM[i]]==1 | comb_counts[,positive_coloc$CM[i]]==1,]$freq
	
	#Get frequency for the presence of both features in the pair
	size_freqs<-length(freqs)
	int_freqs[,eval(int_label):=empty_states]
	temp<-empty_states
	for(j in 1:size_freqs)
	{
		temp[comb_counts[comb_counts[,positive_coloc$HM[i]]==1 & comb_counts[,positive_coloc$CM[i]]==1,]$state[j]]<-comb_counts[comb_counts[,positive_coloc$HM[i]]==1 & comb_counts[,positive_coloc$CM[i]]==1,]$freq[j]/(comb_counts[comb_counts[,positive_coloc$HM[i]]==1 & comb_counts[,positive_coloc$CM[i]]==1,]$freq[j] + comb_counts[comb_counts[,positive_coloc$HM[i]]==1 & comb_counts[,positive_coloc$CM[i]]==0,]$freq[j] + comb_counts[comb_counts[,positive_coloc$HM[i]]==0 & comb_counts[,positive_coloc$CM[i]]==1,]$freq[j])
	}
	int_freqs[,eval(int_label):=temp]
}
rownames(int_freqs)<-int_freqs$state
int_freqs$state<-NULL
write.table(as.matrix(int_freqs),file="tmp.dist")
std_int_freqs<-t(scale(t(int_freqs)))
write.table(std_int_freqs,file="tmp_std.dist")
#Run hiearchical clustering with bootstrap
int_clustering_result<-pvclust(std_int_freqs,method.dist="correlation", method.hclust="average",nboot=1000)

#Get significant clusters of interactions = chromnets
alpha_value=1-chromnet_pvalue
chromnets<-pvpick(int_clustering_result,alpha=alpha_value)


#Format output
output<-""
for(i in 1:length(chromnets$clusters)){
	for(j in 1:length(chromnets$clusters[[i]])){
		if(j==1){
			output[i]<-sprintf("Chromnet_%i\t",i)
			output[i]<-paste(output[i],chromnets$clusters[[i]][j],sep="")
		}else{output[i]<-paste(output[i],chromnets$clusters[[i]][j],sep=",")}
	}
}


#Write chromnets to output file
write(output,file=outfile)

quit()

