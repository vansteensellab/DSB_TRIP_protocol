#!/bin/env bash
SCRIPT_DIR=$( cd $( dirname ${BASH_SOURCE[0]} ) &> /dev/null && pwd )
SCRIPT=`basename ${BASH_SOURCE[0]}`

#Initialize variables to default values.
THREADS=1
OUTDIR=$SCRIPT_DIR/example
META=$SCRIPT_DIR/config/metadata.txt
CONFIG=$SCRIPT_DIR/config/LBR_example_config.yaml
SNAKE=$SCRIPT_DIR/src/dsb_trip.snake

#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

#Help function
function HELP {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM} ${BOLD}$SCRIPT file.ext${NORM}"\\n
  echo "Command line switches are optional. The following switches are recognized."
  echo "${REV}-d${NORM}  --Sets output directory"
  echo "${REV}-t${NORM}  --threads to use"
  echo "${REV}-s${NORM}  --snakefile location of the snakemake file"
  echo "${REV}-m${NORM}  --metadata file to use (must contain SRR_id column as well as file, PCR_type and ID columns required by pipeline)"
  echo "${REV}-c${NORM}  --config file for pipeline"
  echo -e "${REV}-h${NORM}  --Displays this help message. No further functions are performed."\\n
  echo -e "Example: ${BOLD}$SCRIPT -d /path/to/data ${NORM}"\\n
  exit 1
}



while getopts :d:t:s:m:c:h FLAG; do
  case $FLAG in
    d)  #set option "a"
      OUTDIR=$OPTARG
      ;;
    t)  #set option "b"
      THREADS=$OPTARG
      ;;
    s)  #set option "b"
      SNAKE=$OPTARG
      ;;
    m)  #set option "b"
      META=$OPTARG
      ;;
    c)  #set option "b"
      CONFIG=$OPTARG
      ;;
    h)  #show help
      HELP
      ;;
    \?) #unrecognized option - show help
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      HELP
      #If you just want to display a simple error message instead of the full
      #help, remove the 2 lines above and uncomment the 2 lines below.
      #echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      #exit 2
      ;;
  esac
done

shift $((OPTIND-1))  #This tells getopts to move on to the next argument.

mkdir -p $OUTDIR/raw


srr_list=$(
	awk 'NR==1{
			for (i=1; i<=NF; i++)
	        	ix[$i] = i
		} NR > 1{
			print $ix["SRR_id"]
		}' $META
)

for srr in $srr_list;
do
	if [ ! -f $OUTDIR"/raw/"$srr"_1.fastq.gz" ]; then
	    fastq-dump --gzip --split-files -O $OUTDIR/raw $srr
	fi
done

snakemake -j $THREADS -s $SNAKE \
          --configfile $CONFIG --use-conda \
          --config indir=$OUTDIR/raw outdir=$OUTDIR/results \
