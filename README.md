# WholeExomeSequencing-Pipeline: Panda-Exome-seq variant calling pieline accoring to GATK Best practices. (Somatic and Germline)	
 The script is made for a germline automated workflow that will call SNPs and INDELs
In order to make this somatic one can simply stop post Step 2 and use Mutect2 or VarScan2 with their own custom scripts and filters on the processed bam files.
For somatic Pipeline we use Mutect2 which you can find in the same pipeline.
#Call with following arguments
sh master_script.sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]
you can run the script in above mentioned way through another processing script that will log the processing with time at each step
# Putting the bwa mem
GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
ftp://ftp.broadinstitute.org/bundle/hg19/
The bundles were downloaded from the GATK 2.8 while the GATK 4.0.2 & 3.8.0 was used for the analysis
Please download all the required tools and install them in your Cluster Computing and modify the paths accordingly for effective running of the script	
