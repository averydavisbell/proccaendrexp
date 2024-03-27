#!/usr/bin/env nextflow
/*
By Avery Davis Bell
Batch downloads from SRA using SRA metadata
*/

/*
#### Set up input parameters & defaults ####
*/
// Inputs: input/output related & general
params.srainfo = "SraRunTable.txt" // Path to SraRunTable.txt containing info on SRA runs to download
params.outdir = "rnaseq-bams-PRJNA669810/bams" // output dir for FASTQs
//params.workdir = "~/scratch/run-downloadcendrrnaseqbams" // running directory - should run from in here

// Housekeeping:  create output directories
outdir = file(params.outdir)
outdir.mkdirs()

/*
#### Channel management ####
*/
//  Channel with all needed sample info as tuple.

formatChannel = Channel.fromPath(params.srainfo)
  .splitCsv(header: true, sep: ',', strip: true, quote: '"')
  .map{row->
    return [row.Run, row.'Sample Name', row.strain] // SRR run, sample name, strain
  }
  .set{sraInfo}

/*
#### Processes ####
*/

process downloadsra{
  // download 
  input:
  tuple val(SRR), val(samp), val(strain) from sraInfo
  
  output:
  tuple val(samp), path("*_1.fastq"), path ("*_2.fastq") into fastqs
  
  """
  # Get data
  prefetch ${SRR}
  fasterq-dump ${SRR}

  # Rename
  mv ${SRR}_1.fastq ${samp}_1.fastq
  mv ${SRR}_2.fastq ${samp}_2.fastq
  """
}

process dogzip{
  // gzip
  publishDir outdir, mode: 'copy'

  input:
  tuple val(samp), path(fastq1), path(fastq2) from fastqs

  output:
  path("*.fastq.gz") into outs
  
  """
  gzip ${fastq1}
  gzip ${fastq2}
  """
}


