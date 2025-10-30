# GTEX_my_version
Building GTEX our version for Guis lab

First I wanted to get the docker image for GTEX. As docker is not supported on computecanada, I used
```
module load apptainer
 apptainer pull gtex_rnaseq_V10.sif docker://broadinstitute/gtex_rnaseq:V10
```
