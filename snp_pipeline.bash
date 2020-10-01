#!/bin/bash

# Place your gatech userame in the below export
export NAME="arozanski3"

get_input () {
    # Function for doing your getopts
    #
    # Input: Getopts array
    # Get user input, the user is not required to add the extension to the output files
    # If the flags realign, gunzip, index, and v, are set the code will run those actions
    while getopts "a:b:r:eo:f:zvih" opt
    do
	case $opt in
            a)reads1=$OPTARG;;
            b)reads2=$OPTARG;;
            r)ref=$OPTARG;;
	    e)realign=1;;
	    o)output=$OPTARG;;
	    f)millsFile=$OPTARG;;
	    z)gunzip=1;;
	    v)v=1;;
	    i)index=1;;
	    h)echo"in order to run script need to set:
 -a <Location of first reads file>
-b <Location of seconds read file>
-r <Location of ref file>
-o <Output VCF file name>
-f <Mills file location>
optional flags:
-e <1 or 0> to preform read re-alignment 
-z  Output VCF file should be gunzipped 
-v  verbose mode; 
-i  index your output BAM file
-h  Print usage information "
 exit 0;;
esac
    done
   
}

check_files () {
    # Function for checking for presence of input files, reference genome,
    # and the output VCF file
    # 
    # Input: File locations (string)
    # Output: True, if checks pass; False, if checks fail (bool)
    # Added the millsFile check in addition to the other specified
    # If the output file the user specified already exists they can decide to write overwrite it input: "o"  or exit the script input "e"
    # If neither "o" or "e" are inserted it will also exit. 
    if [ ! -f "$ref" ]
    then
       echo "cannot find ref file"
       exit 1
    
    elif [ ! -f "$reads1" ]
    then
	echo "cannot find reads1 file"
	exit 1
    elif [ ! -f "$reads2" ]
    then
	echo "cannot find reads2 file"
	exit 1
    elif [ -f "$output" ]
    then
       echo "output exists (o)verwrite (e)xit"
       read answer
       if [ "$answer" == 'o' ]
       then
	      echo "okay"
       elif [ "$answer" == 'e' ]
       then
	      exit 1
       else
	      echo "input not recognized"
	         exit 1
       fi
    elif [ ! -f "$millsFile" ] 
    then
	echo "cannot find mills file"
	exit 1
    fi
   
    
}
prepare_temp () {
    # Preparing your temporary directory
    #
    #Makes a temporary directory within my tmp directory which can hold all the files created in the process
    #This provides easy access to the files when running the code and deletion of generated intermediate files upon exit of code
    #tmp/tmp.XXXXXX
    tmp=$(mktemp -d)
    echo "creating temp"
}


mapping () {
    # Function for the mapping step of the SNP-calling pipeline
    #
    # Input: File locations (string), Verbose flag (bool)
    # Output: File locations (string)
    # The following if statement sets verbose mode if the variable v is set in getopts 
    if [ ! -z "$v" ]
    then
	echo "Indexing the reference genome and mapping the sample file to the reference genome"
    fi
    #This command prepares the reference genome for mapping 
    bwa index "$ref"
    #This command preforms the actual mapping and outputs a file in .sam format 
    bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' "$ref" "$reads1" "$reads2" > "$tmp"/lane.sam
    #Cleans up any unusual FLAG information 
    samtools fixmate -O bam "$tmp"/lane.sam  "$tmp"/lane_fixmate.bam
    #Sort them into coordinate order output: lane_sorted.bam
    samtools sort -O bam -o "$tmp"/lane_sorted.bam -T "$tmp"/lane_temp "$tmp"/lane_fixmate.bam
    
}

improvement () {
    # Function for improving the number of miscalls
    #
    # Input: File locations (string)
    # Output: File locations (string)
   
    
    #index's the reference genome and Creates the .fai file for GATK
    samtools faidx "$ref"
    #index's the sorted mapping and creates the .bai file for GATK
    samtools index "$tmp"/lane_sorted.bam
    #creates a sequence dictionary from the file. The sed command changes the extension of the reference seq so we can save the created file as .dict 
    samtools dict "$ref" -o ./data/chr17.dict
    #realigns the files
    java -Xmx2g -jar ./lib/GenomeAnalysisTK.jar -log ./output/arozanski3.log -T RealignerTargetCreator -R "$ref" -I "$tmp"/lane_sorted.bam -o "$tmp"/lane.intervals -known "$millsFile"
    java -Xmx4g -jar ./lib/GenomeAnalysisTK.jar -log ./output/arozanski3.log -T IndelRealigner -R "$ref" -I "$tmp"/lane_sorted.bam -targetIntervals "$tmp"/lane.intervals --known "$millsFile" -o "$tmp"/lane_sorted.bam
 
   
    
}

call_variants () {
    # Function to call variants
    #
    # Input: File locations (string)
    # Ouput: None
    # if variable is set in getopts the created bam files are indexed
    if [ ! -z "$index" ]
    then
        samtools index "$tmp"/lane_sorted.bam
    fi
    #mpileup generates VCF file that contains genotype likelihoods for the alignment
    #this pipes to call which preforms SNP/indel calling
    #the -O flag determine if the output file is zipped
    bcftools mpileup -Ou -f "$ref" "$tmp"/lane_sorted.bam | bcftools call -vmO z -o "$output".vcf.gz
    
    #requires a zipped file as input
    tabix -p vcf "$output".vcf.gz
    #if the gunzip variable is not set in the geopts it will unzip the compressed file
    if [ -z "$gunzip" ]
    then
    #The -d flag decompresses the file
	gunzip -d "$output".vcf.gz
    fi
    
    
}

main() {
    # Function that defines the order in which functions will be called
    # You will see this construct and convention in a lot of structured code.
    
    # Add flow control as you see appropriate
    get_input "$@"
    check_files
    prepare_temp
    mapping
    #Only run improvement if the realign variable is set in getopts
    #If not set continue to call_variants and use the bam file generated in mapping
    #echo "$realign"
    if [ ! -z "$realign" ]
    then
	#echo "HELLO"
        improvement 
    fi
    call_variants 
}

# Calling the main function
main "$@"


# DO NOT EDIT THE BELOW FUNCTION
bats_test (){
    command -v bats
}
# DO NOT EDIT THE ABOVE FUNCTION
