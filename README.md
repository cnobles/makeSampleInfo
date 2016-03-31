# makeSampleInfo
Generate sampleInfo.tsv for intSiteCaller from databased information given the date when the samples were run.

To run the script:
```
Rscript path/to/makeSampleInfo.R YYYY-MM-DD                                         #Typical execution
Rscript path/to/makeSampleInfo.R YYYY-MM-DD -d [specimen.database]                  #Specify which database to connect
Rscript path/to/makeSampleInfo.R YYYY-MM-DD -d [specimen.database] -g [ref_genome]  #Specify the ref_genome (default = hg18)
```
