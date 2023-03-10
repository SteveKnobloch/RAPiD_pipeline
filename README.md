# RAPiD pipeline
A **R**apid and **A**ccurate plant **P**athogen **iD**entification pipeline for long-read sequencing.  

Available as stand-alone tool (see [Installation](#Installation)) or as web app with live basecalling at http://agrifuture.senckenberg.de

# Requirements
Linux or macOS  
Nextflow (install here: https://www.nextflow.io/docs/latest/getstarted.html)  
Docker (install here: https://docs.docker.com/engine/install/ubuntu/)  

# Installation
### Get Nextflow script and databases:
```
git clone https://github.com/SteveKnobloch/RAPiD_pipeline.git
cd RAPiD_pipline

wget https://figshare.com/ndownloader/files/39546862
tar -xvf 39546862
```

### Build docker image:  

```
sudo docker build -t rapid_pipe .
```  

or pull from steveknobloch1444/rapid_pipe with:  
```
docker pull steveknobloch1444/rapid_pipe
docker rename steveknobloch1444/rapid_pipe rapid_pipe
```  

### Run a simple command:  
```
nextflow RAPiD_pipeline.nf --help
```  
If you get a permission error for the docker container, first run: ```sudo chmod 666 /var/run/docker.sock```


# Running the pipeline

### Quick usage:
  Run RAPiD pipeline in real-time mode with basecalled data being generated in <input_folder>:  
     ```
     nextflow run RAPiD_pipeline.nf --input <input_folder> --output <report_folder> --batch --realtime
     ```

  Run RAPiD pipeline in batch mode with basecalled data already in <input_folder>:  
      ```
      nextflow run RAPiD_pipeline.nf --input <input_folder> --output <report_folder> --batch
      ```

### Main arguments:  
```--input```       path to basecalled fastq files [e.g. /data/nanopore_run/sample/fastq_pass/]  
```--output```      path to save report [default: report/]  
```--batch```       operates RAPiD in batch mode  
```--realtime```    operates RAPiD in real-time mode  

### Optional arguments:  
```--index```       path to alternative minimap2 index file  
```--taxa```        path to taxonomic look-up file  
```--threads```     threads [default: 1]  
```--subspecies```  performs analysis at sub-species level  
```--sensitive```   performs analysis in sensitive mode i.e. without cut-off values  

```--help```        help message

# Output

### RAPiD report
Provides a tab delimited report for all sequences matching a reference in the database.  
Column 1: Query sequence ID  
Column 2: Reference sequence ID  
Column 3: Scientific name of reference organism  
Column 4: Sequence length  
Column 5: Normalized alignment score  
Column 6: Mean per base identity  
Column 7: Query sequence  

### RAPiD summary
Provides information of the number of reads that passed QC, the number of reads that matched a reference pathogen in the database, along with the read count, normalized alignment score, mean per base identity and confidence score.

The confidence score is considered "high" when at least two reads matched the reference taxon and the average normalized alignment score of that taxon is equal to or above 75. A "medium" score is given for a taxon with at least two matches and an average normalized alignment score of between 70 and 75. Otherwise, the confidence score is considered "low".

# Citations
If you use this tool please consider citing following references:  

RAPiD: Knobloch S. et al. (bioRxiv) (2023).

Porechop: Wick, R. R. "Porechop. Github https://github.com/rrwick." (2017).  

NanoFilt and NanoLyse: De Coster, W. et al. "NanoPack: visualizing and processing long-read sequencing data." Bioinformatics 34.15 (2018): 2666-2669.  

minimap2: Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34 (2018):3094-3100.  

SAMtools: Li, H. et al. The sequence alignment/map format and SAMtools. Bioinformatics, 25.16 (2009): 2078-2079.  



***
THIS SOFTWARE COMES WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
