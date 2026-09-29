# check config
nextflow config nextflow.config
echo $?

# lint
nextflow lint nextflow.config

# dry-run
nextflow run MitoHPC2.sr.nf -preview
