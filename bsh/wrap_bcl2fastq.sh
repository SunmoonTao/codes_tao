#! /bin/bash

#$ -cwd
#$ -q long
#$ -P regevlab

#$ -l m_mem_free=4g
#$ -e error.err
#$ -o out.log


source /broad/software/scripts/useuse 
reuse -q .bcl2fastq2-2.17.1.14

indir=$1
outdir=$2
samplesheet=$3

mkdir -p $outdir

nohup bcl2fastq --barcode-mismatches 0 --runfolder-dir $indir --output-dir $outdir --sample-sheet $samplesheet --no-lane-splitting  --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10


#bsub -n 1 -q regevlab -e $outdir/pipe.err -o $outdir/pipe.txt -N -R "span[hosts=1]rusage[mem=80]" nohup bcl2fastq --runfolder-dir $indir --output-dir $outdir --sample-sheet $samplesheet --no-lane-splitting  --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10

#nohup bcl2fastq --runfolder-dir $indir --output-dir $outdir --sample-sheet $samplesheet --no-lane-splitting  --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10

#bsub -P rnaseq -n 1 -q week -e /ahg/regev_miseq01/Data/150929_NB501164_0002_AHLNGLBGXX/pipe.err -o /ahg_regev_miseq01/Data/150929_NB501164_0002_AHLNGLBGXX/log-preproc.txt -N -R "span[hosts=1]rusage[mem=80]" nohup bcl2fastq --runfolder-dir /ahg/regev_miseq01/Data/150929_NB501164_0002_AHLNGLBGXX/ --output-dir /ahg/regev_miseq01/Data/150929_NB501164_0002_AHLNGLBGXX/Data/Intensities/BaseCalls/
