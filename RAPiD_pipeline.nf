#!/usr/bin/env nextflow

/*  This script takes Nanopore basecalled fastq file chunks and
performes QC filtering and taxonomic profiling against a curated database
with plant pathogen genomes.
Results are then merged for each fastq chunk analysis and a normalised similarity
score is given for each alignment passing the set threshold.

Author: Stephen Knobloch, Thines Lab, Senckenberg Institute for Biodiversity
and Climate Reaserch, Germany.
*/

//  Help message
def helpMessage() {
    log.info"""
    RAPiD Pipeline v. 1.0

    Usage:
      Run RAPiD pipeline in realtime mode with basecalled data being generated in <input_folder>:
      nextflow run RAPiD_pipeline.nf --input <input_folder> --output <report_folder> --batch --realtime

      Run RAPiD pipeline in batch mode with basecalled data already in <input_folder>:
      nextflow run RAPiD_pipeline.nf --input <input_folder> --output <report_folder> --batch

    Arguments:
    --input \tpath to basecalled fastq files
    --output \tpath to save report [default: report/]
    --batch \toperates RAPiD in batch mode
    --realtime \toperates RAPiD in real-time mode

    Optional arguments:
    --index \tpath to alternative minimap2 index file
    --taxa \tpath to taxonomic look-up file
    --threads \tthreads [default: 1]
    --subspecies \tperforms analysis at sub-species level
    --sensitive \tperforms analysis in sensitive mode i.e. without cut-off values

    --help \thelp message
    """
}

params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Script parameters
params.output = "report"
params.index = "database/pathogens_index_230309.mmi"
params.taxa = "database/rapid_taxonomy_230309.txt"
params.threads = 1
params.subspecies = false
params.sensitive = false
params.realtime = false
params.batch = false

index = file(params.index)
taxa = file(params.taxa)

// Define cut-off values through params.sensitive
if (params.sensitive) {
  AS = 0
  cm = 0
  de = 1
} else {
  AS = 2000
  cm = 200
  de = 0.075
}

// Define report output through params.subspecies
if (params.subspecies) {
  sub = "\$9"
} else {
  sub = "\$10"
}

// Validate inputs
if (!params.realtime && !params.batch) {
  println "Warning: Please specify run mode: --realtime, --batch or both"
  exit 0
}

def parameters_expected = ['input', 'output', 'index', 'taxa', 'threads',
  'subspecies', 'help', 'batch', 'realtime', 'sensitive'
  ] as Set

def parameter_diff = params.keySet() - parameters_expected
  if (parameter_diff.size() != 0){
     exit 1, "[Error] Parameter(s) $parameter_diff is/are not a valid parameter!\n\
     Type --help for more information.\n"
  }

def checkPathParamList = [ params.index, params.taxa ]
for (param in checkPathParamList) {if (param) { file(param, checkIfExists: true) } }


// Remove temporary folder (if exists) from previous runs.
import groovy.io.FileType

String folderPath = "temp"
new File(folderPath).eachFile (FileType.FILES) { file ->
   file.delete()
}

println "\n\nRAPiD pipeline v 1.0. \nReport saved to $params.output/rapid_report.txt.\n\n"


// Define channels for running real-time analysis (params.realtime)
if (params.realtime) {

  query_ch = Channel
    .watchPath(params.input + '/*.fastq')
    .map { file -> tuple(file.baseName, file) }

  merge_ch = Channel
    .watchPath ( 'temp/*out' )

} else {

  query_ch = Channel.empty()
  merge_ch = Channel.empty()

}


// define channels for running batch analysis (params.batch)
if (params.batch) {

  batch_ch = Channel.fromPath(params.input + '/*.fastq')

} else {

  batch_ch = Channel.empty()

}


/*
  Process runMinimap_reatime: QC (Adapter trimming, quality filtering, size filtering,
  Lambda control removal) and run mimimap2 alignment against pathogen
  index in real-time mode
*/

process runMinimap_realtime {

  publishDir 'temp', mode: 'copy', pattern: '*.{out,stat}'

  input:
  set queryID, file(queryFile) from query_ch

  output:
  set queryID, file("${queryID}_counts.out")
  set queryID, file("${queryID}_reads.stat")

  script:
  """
  # Count reads of query file
  cat ${queryFile} | echo \$((`wc -l`/4)) > read_counts

  # Remove adapters, filter reads > 1000bp and q > 8, remove Lambda control reads
  porechop -i ${queryFile} -o adaptertrimmed -v 0 --no_split
  NanoFilt -q 8 -l 1000 adaptertrimmed > sizefilt
  cat sizefilt | NanoLyse -r $PWD/database/Lambda_ref.fasta | cat  > lambfilt

  # Count reads passed QC
  cat lambfilt | echo \$((`wc -l`/4)) | cat > read_counts_qc
  paste read_counts read_counts_qc > ${queryID}_reads.stat

  # Read alignment with minimap2
  minimap2 -ax map-ont --secondary=no -t $params.threads --split-prefix GON --sam-hit-only \
  $index lambfilt -o alignment.sam

  # Sort, index alignment, remove supplementary alignments
  samtools view -S -b -F0x900 alignment.sam | \
  samtools sort > sorted.bam

  # Filter alignments with AS, cm and de cut-off values
  samtools view -S sorted.bam | awk '{split(\$14, sf1, ":"); \
  split(\$17, sf2, ":"); split(\$20, sf3, ":"); \
  if(sf1[3]>$AS && sf2[3]>$cm && sf3[3]<$de) print \$0}' > alignment.sam

  # Extract SeqID, RefID, AS (Alignment Score), de (per-base gap-compressed sequence divergence) query sequence, length of query sequence
  awk '/^[^@]/ {split(\$14, sf, ":"); split(\$20, sd, ":"); print \$3, \$1, sf[3], sd[3], \$10, length(\$10) }' \
  OFS='\t' alignment.sam | sort > alignment_report

  # Add taxa of RefID and print values with normalised score to output
  join -a1 alignment_report $taxa -t \$'\t' | \
  awk -F \$'\t' '{print \$2, \$1, $sub, \$6, \$3/\$6*50, (1-\$4)*100, \$5}' OFS='\t' \
  > ${queryID}_counts.out

  rm read_counts read_counts_qc adaptertrimmed lambfilt \
  alignment.sam sorted.bam alignment_report sizefilt
  """
}


/*
  Process runMinimap_batch: QC (Adapter trimming, quality filtering,
  size filtering, removal of Lambda control reads) and
  run mimimap2 alignment against pathogen index in batch mode
*/

process runMinimap_batch {

  publishDir 'temp', mode: 'copy', pattern: '*.{out,stat}'

  input:
  file x from batch_ch.collect()

  output:
  file("batch_counts.out")
  file("batch_reads.stat")
  file("batch_counts.out") into report_ch


  script:
  """
  # Combined all fastq files in input folder and count reads
  cat $x > merged.fastq
  cat merged.fastq | echo \$((`wc -l`/4)) > read_counts

  # Remove adapters, filter reads > 1000bp and q > 10, remove Lambda control reads
  porechop -i merged.fastq -o adaptertrimmed -v 0 --no_split
  NanoFilt -q 8 -l 1000 adaptertrimmed > sizefilt
  cat sizefilt | NanoLyse -r $PWD/database/Lambda_ref.fasta | cat  > lambfilt

  # Counts reads passed QC
  cat lambfilt | echo \$((`wc -l`/4)) | cat > read_counts_qc
  paste read_counts read_counts_qc > batch_reads.stat

  # Read alignment with minimap2
  minimap2 -ax map-ont --secondary=no -t $params.threads --split-prefix GON --sam-hit-only \
  $index lambfilt -o alignment.sam

  # Sort, index alignment, remove supplementary alignments
  samtools view -S -b -F0x900 alignment.sam | \
  samtools sort > sorted.bam

  # Filter alignments with AS > 2000, cm > 200 and de < 0.075
  samtools view -S sorted.bam | awk '{split(\$14, sf1, ":"); \
  split(\$17, sf2, ":"); split(\$20, sf3, ":"); \
  if(sf1[3]>$AS && sf2[3]>$cm && sf3[3]<$de) print \$0}' > alignment.sam

  # Extract SeqID, RefID, AS (Alignment Score), de (per-base gap-compressed sequence divergence) query sequence, length of query sequence
  awk '/^[^@]/ {split(\$14, sf, ":"); split(\$20, sd, ":"); print \$3, \$1, sf[3], sd[3], \$10, length(\$10) }' \
  OFS='\t' alignment.sam | sort > alignment_report

  # Add taxa of RefID and print values with normalised score to output
  join -a1 alignment_report $taxa -t \$'\t' | \
  awk -F \$'\t' '{print \$2, \$1, $sub, \$6, \$3/\$6*50, (1-\$4)*100, \$5}' OFS='\t' \
  > batch_counts.out

  rm read_counts adaptertrimmed lambfilt NanoLyse.log \
  read_counts_qc alignment.sam sorted.bam alignment_report merged.fastq
  """
}


/*
Process mergeCount: Merge results from output in real-time mode. Then summarize
counts, calculate normalized alignment score and print report and report summary.
*/

process mergeCount {
  publishDir params.output, mode: 'copy', pattern: '*.{txt}'

  input:
  path counts from merge_ch

  output:
  file("RAPiD_report.txt")
  file("RAPiD_summary.txt")
  stdout into result_merge

  script:
  """
  # Combine all stat and count results, summarize species count, calculate mean score value, print RAPiD_report.txt
  cat $PWD/temp/*.stat > final.stats
  cat $PWD/temp/*.out > temp.counts
  cat temp.counts > RAPiD_report.txt
  awk -F '\t' '{sc[\$3]+=\$5; sk[\$3]+=\$6; count[\$3]++} END {for (i in sc) \
  print i"\t"count[i]"\t"sc[i]/count[i]"\t"sk[i]/count[i];}' temp.counts | sort -t\$'\t' -k2 -nr > final.counts

  # Add confidence score (low, medium or high)
  cat final.counts | awk -F'\t' 'BEGIN {OFS = FS} {if (\$2 > 1  && \$3 >= 70 && \$3 < 75) conf = "medium"; else if (\$2 > 1 && \$3 >= 75) conf = "high"; else conf = "low"; print \$1, \$2, \$3, \$4, conf}' > final2.counts

  # Print results to RAPiD_summary.txt in output directory
  echo -e "RAPiD pipeline v 1.0  Time: \$(date) \n\n \$(awk '{sum2 += \$2 } END {print sum2}' final.stats) reads passed quality filtering from \$(awk '{sum += \$1 } END {print sum}' final.stats) basecalled reads.\n \$(awk -F '\t' '{sum += \$2 } END {if (sum == "") print "0"; else print sum}' final.counts) reads matched a pathogen in the database. \n\nTarget_species \tRead_counts \tNormalized_alignment_score \tMean_per_base_identity \tConfidence_score" | \
  cat - final2.counts > RAPiD_summary.txt

  # Print results to screen
  export TERM=linux
  clear
  tabs 25
  echo -e "\nRAPiD pipeline v 1.0. \nReport saved to $params.output/RAPiD_report.txt."
  echo -e "\n\n \$(awk '{sum2 += \$2 } END {print sum2}' final.stats) reads passed quality filtering from \$(awk '{sum += \$1 } END {print sum}' final.stats) basecalled reads.\n \$(awk -F '\t' '{sum += \$2 } END {if (sum == "") print "0"; else print sum}' final.counts) reads matched a pathogen in the database. \n\n-Results-\n\nTarget Species \tRead counts \tNormalized score \tMean per-base identity \tConfidence score" | \
  cat - final2.counts | head -25
  echo "(showing top 10 matches)"
  echo -e "\n ...Waiting for new files in $params.input/. Press Ctr+C to exit."
  rm temp.counts final.counts final2.counts final.stats
  """
}

result_merge.subscribe { println it }


/*
Process reportBatch: Summarize counts, calculate normalized alignment score and
print report and report summary.
*/

process reportBatch {

  publishDir params.output, mode: 'copy', pattern: '*.{txt}'

  input:
  file("batch_counts.out") from report_ch

  output:
  file("RAPiD_report.txt")
  file("RAPiD_summary.txt")
  stdout into result_batch

  script:
  """
  # Print alignment results to RAPiD_report.txt in output directory
  cat batch_counts.out > RAPiD_report.txt

  # Summarize species count, calculate mean alignment score value and mean per-base identity
  awk -F '\t' '{sc[\$3]+=\$5; sk[\$3]+=\$6; count[\$3]++} END {for (i in sc) \
  print i"\t"count[i]"\t"sc[i]/count[i]"\t"sk[i]/count[i];}' batch_counts.out | sort -t\$'\t' -k2 -nr > batch_counts2.out

  # Add confidence score (low, medium or high)
  cat batch_counts2.out | awk -F'\t' 'BEGIN {OFS = FS} {if (\$2 > 1  && \$3 >= 70 && \$3 < 75) conf = "medium"; else if (\$2 > 1 && \$3 >= 75) conf = "high"; else conf = "low"; print \$1, \$2, \$3, \$4, conf}' > batch_counts3.out

  # Print results to RAPiD_summary.txt in output directory
  echo -e "RAPiD pipeline v 1.0  Time: \$(date) \n\n \$(awk '{sum2 += \$2 } END {print sum2}' $PWD/temp/batch_reads.stat) reads passed quality filtering from \$(awk '{sum += \$1 } END {print sum}' $PWD/temp/batch_reads.stat) basecalled reads.\n \$(awk -F '\t' '{sum += \$2 } END {if (sum == "") print "0"; else print sum}' batch_counts2.out) reads matched a pathogen in the database. \n\nTarget_species \tRead_counts \tNormalized_alignment_score \tMean_per_base_identity \tConfidence_score" | \
  cat - batch_counts3.out > RAPiD_summary.txt

  # Print results to screen
  export TERM=linux
  tabs 25
  echo -e "\nRAPiD pipeline v 1.0. \nReport saved to $params.output/RAPiD_report.txt."
  echo -e "\n\n \$(awk '{sum2 += \$2 } END {print sum2}' $PWD/temp/batch_reads.stat) reads passed quality filtering from \$(awk '{sum += \$1 } END {print sum}' $PWD/temp/batch_reads.stat) basecalled reads.\n \$(awk -F '\t' '{sum += \$2 } END {if (sum == "") print "0"; else print sum}' batch_counts2.out) reads matched a pathogen in the database. \n\n-Results-\n\nTarget Species \tRead counts \tNormalized score \tMean per-base identity \tConfidence score" | \
  cat - batch_counts3.out | head -19
  echo "(showing top 10 matches)"

  rm batch_counts2.out batch_counts3.out
  """
}

result_batch.subscribe {println it}
