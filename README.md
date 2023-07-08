# RNA-seq_pipeline_NF
The files for running my RNA-seq nextflow pipeline

The r_dev_test.def file is a singularity definition file that contains the instructions for building a singularity container image.
When on the HPC you will need only one module to load. That is the singularityCE module.
The version I have available to me is singularity-ce/3.11.0
Once this is loaded we will build a .sif file using the remote method, since we probably do not have root access on the HPC.
First we may have to login to singularity 
```
singularity remote login
```
You will need to go to the singularity website, create an account and login then generate an access key in the top right. Then copy and paste that access key into the terminal as a way to remote login.

Next we can build finally
```
singularity build --remote <name_of_container>.sif r_dev_test.def
```
This will create the .sif file with the name of your choice. And now the container can be specified in the nextflow config file and placed in the process that will use this container.

To check if the container works as intended you can just do the following
```
singularity run <name_of_container>.sif
```
This should open R in the command line, and now you can load the three packages we want to use
```
library(DESeq2)
library(Rsubread)
library(EnhancedVolcano)
```

In the nextflow pipeline, I have a process that calls the .R file that contains the R commands I want to use. Once this script runs sussessfully in nextflow, a work directory will be created that does not contain any files. If you update the .R file and add new code or change anything, you will need to go into the slurm script and find the directory where the process r_featurecounts calls the .R script and then delete that directory only. This way the nextflow -resume command will not think the process is still cashed. Ex how to delete that directory:

```
rm work/de/111b213d38af9696476b48cc3cee04/
```
This is the directory for me that contains that 'cashed' output, even though it is empty.
Now you can run the pipeline again with the updated commands in the .R file.
