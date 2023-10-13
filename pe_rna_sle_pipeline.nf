// this works to get the read pairs in a new tuple form


// run this to see how it looks again
//params.reads = "store_fq_reads/45*-Mock*polya*_{R1,R2}*.fastq.gz"

//do not need the background genomic regions for the rna pipeline

//params.bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

FASTP='fastp/intel/0.20.1'
FASTQC='fastqc/0.11.9'
STAR='star/intel/2.7.6a'
R='r/gcc/4.2.0'

//using star for the splice aware aligner in the rna seq pipeline instead
//BWA='bwa/intel/0.7.17'

SAMTOOLS='samtools/intel/1.14'
MULTIQC='multiqc/1.9'

// i want to put all the filt files into a dir so i can get them in the next process

params.outdir = 'filt_files'

params.bams = 'store_bam_files'

params.store_fq_reads = 'store_fq_reads'


// have to change this part when using the mock, iav, or bleo conditions

//params.filts = "filt_files/45*-Mock*polya*_{R1,R2}*.filt*"



params.tech_reads_r1 = "/scratch/rj931/tf_sle_project/all_sle_data/45*-Mock*polya*_{L001,L002}*R1*.fastq.gz"
params.tech_reads_r2 = "/scratch/rj931/tf_sle_project/all_sle_data/45*-Mock*polya*_{L001,L002}*R2*.fastq.gz"
//Channel
//	.watchPath('filt_files', 'create,modify')
//	.filter { it.name ==~ params.filts}
//	.until { it} 
//	.set {filt_pe}

//params.reads = "store_fq_reads/45*-Mock*polya*{R1,R2}*.fastq.gz"

//params.filts = "filt_files/45*-Mock*polya*{R1,R2}*.filt*"

nextflow.enable.dsl=2


tech_rep_reads_r1 = Channel.fromFilePairs(params.tech_reads_r1, checkIfExists: true)
tech_rep_reads_r2 = Channel.fromFilePairs(params.tech_reads_r2, checkIfExists: true)

process cat_tech_reps {
    
    publishDir "${params.store_fq_reads}", mode: 'copy', pattern: '*fastq.gz'

    input:

    tuple val(r1_name), path(tech_rep_reads_r1)

    tuple val(r2_name), path(tech_rep_reads_r2)

    output:

    path "*fastq.gz", emit: fastq_cat_reads

    //path "${out_name_r1}", emit: fastq_r1_reads

    //path "${out_name_r2}", emit: fastq_r2_reads

    script:
    out_name_r1 = "${r1_name}_R1.fastq.gz"
    out_name_r2 = "${r2_name}_R2.fastq.gz"
    """

    cat ${tech_rep_reads_r1[0]} ${tech_rep_reads_r1[1]} > ${out_name_r1}

    cat ${tech_rep_reads_r2[0]} ${tech_rep_reads_r2[1]} > ${out_name_r2}
    """
}




//PE_reads = Channel.fromFilePairs(params.reads)

process fastp {

//tag
publishDir params.outdir, mode: 'copy'

input:
tuple val(pair_id), path(fastqs)
// input in fastp would be ${reads[0]} and ${reads[1]}
// something i can do in fastp is to now concatenate the $pair_id'_R1.filt.fastq.gz' 
// same for output 2 in fastp use $pair_id'_R2.filt.fastq.gz'


output:

path "*filt.fq.gz", emit: filt_fq_files

//path "${pair_id}_R1.filt.fq.gz", emit: fastp_out_f
//path "${pair_id}_R2.filt.fq.gz", emit: fastp_out_r

"""
#!/bin/env bash

module load $FASTP
module load $FASTQC

fastp \
-i ${fastqs[0]} \
-I ${fastqs[1]} \
-o $pair_id'_R1.filt.fq.gz' \
-O $pair_id'_R2.filt.fq.gz' \
--detect_adapter_for_pe \
--correction

# trying the correction option to correct mismatched base pairs in overlapped regions of paired end reads
# no need for trimming anymore with modern aligners
#--trim_front1 7 \
#--trim_front2 7 \

fastqc $pair_id'_R1.filt.fq.gz' $pair_id'_R2.filt.fq.gz'

# since this script is running somewhere in the work directory I need to figure out how to find where all the qc files are

#find ./../.. -name *fastqc.zip > fastqc_files.txt

#module load $MULTIQC

#multiqc -force --file-list fastqc_files.txt --filename 'multiqc_report.html'

"""

}




/*params.outdir2 = "bg_filt_files"
params.bg_filts = "bg_filt_files/461-IgG*cut*_{R1,R2}*.fq.gz"

//Channel
//	.watchPath('bg_filt_files', 'create,modify')
//	.filter { it.name ==~ params.bg_filts }
//	.take(12)
//	.set {bg_reads}

bg_reads = Channel.fromFilePairs(params.bg_regions, checkIfExists: true)

process bg_fastp{

publishDir params.outdir2, mode: 'copy'

input:
tuple val(pair_id), path(bg_reads)

output:
path "${pair_id}_R1.filt.fq.gz", emit: bg_filt_out_f
path "${pair_id}_R2.filt.fq.gz", emit: bg_filt_out_r

"""
#!/bin/env bash

module load $FASTP

fastp \
-i ${bg_reads[0]} \
-I ${bg_reads[1]} \
-o $pair_id'_R1.filt.fq.gz' \
-O $pair_id'_R2.filt.fq.gz' \
--detect_adapter_for_pe \
--trim_front1 7 \
--trim_front2 7 \
"""

}*/





//filt_pe = Channel.fromFilePairs(params.filts) 
params.outdir4 = "store_normal_bam_files"

process star {
publishDir params.outdir4, mode: 'copy'
cpus 10
executor 'slurm'
memory '45 GB'

input:
val ref
//tuple val(pair_id), path(filt_pe)
val gtf

output:
//path "${pair_id}_sort.bam", emit: bam_file
path "genome_generate_finished.txt", emit: genome_gen_finished

"""
#!/bin/env bash

module load $STAR
module load $SAMTOOLS

mkdir /scratch/rj931/tf_sle_project/ref_indices

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir "../../../ref_indices" \
--genomeFastaFiles $ref \
--sjdbGTFfile $gtf \
--sjdbOverhang 49 \
--sjdbGTFfeatureExon exon 

touch genome_generate_finished.txt

"""
}

process star_align {

publishDir params.bams , mode: 'copy', pattern: "${pair_id}*.bam"

memory '45 GB'
cpus 10
executor 'slurm'

input: 
tuple val(pair_id), path(filt_pe)

path genome_finished

output:
path "${pair_id}*.bam", emit: star_bam_files

script:

"""
#!/bin/env bash

module load $STAR

STAR --runMode alignReads \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat \
--genomeDir "../../../ref_indices" \
--outFileNamePrefix "${pair_id}" \
-- readFilesIn ${filt_pe[0]} ${filt_pe[1]} 

#cd ../../../store_bam_files

#ls *bam > bam_list.txt

"""
}


process organize {

input: 
path bams

output:
path "place_holder.txt", emit: place_holder

script:

"""
touch place_holder.txt

cd ../../../store_bam_files

ls *bam > bam_list.txt
"""

}



/*process organize {
publishDir  "${params.bams}/merged_bams" , mode: 'copy' , pattern: "${output_file_name}"

input: 
tuple val(key), path(bam_pairs)

output:
path "place_holder.txt", emit: place_holder
path "${output_file_name}", emit: merged_bams

script:

output_file_name = "${key}.merged.bam"

"""

module load $SAMTOOLS

samtools merge -o "${output_file_name}"  ${bam_pairs[0]} ${bam_pairs[1]}

touch place_holder.txt

#cd ../../../store_bam_files/merged_bams

#ls *bam > bam_list.txt
"""

}

process collected_bams {
publishDir "${params.bams}/merged_bams", mode: 'copy' , pattern: "merged_bams.txt"

input:

val collected_bams_names

output:

path "merged_bams.txt" , emit: txt_merged_bams

script:
"""
echo '${collected_bams_names.join("\n")}' > merged_bams.txt



"""

}*/

/*process r_featurecounts {


input:
val gtf
path bams

output:

path "*_counts.txt", emit: feature_counts_files

script:

"""
#!/bin/env Rscript

library(featureCounts)

output_name = paste0("${bams}", "_counts.txt")

fc = featureCounts( file = "${bams}",
annot.ext = "${gtf}",
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id",
isPairedEnd = TRUE,
countReadPairs = TRUE
)



write.table(fc\$counts, file = output_name, sep = "\t", quote = FALSE) 

"""



}*/

process r_featurecounts {
cache false

input:
path place_holder

output:

script:

"""
cd ../../..
#export PATH="/usr/local/bin/:$PATH"
Rscript feature_count.R

"""

}


/*
#bwa sampe $ref ${filt_pe[0]}'reads_1.sai' ${filt_pe[1]}'reads_2.sai' ${filt_pe[0]} 
${filt_pe[1]}
#> ${pair_id}'aligned_reads.sam'
#samtools view -b -h -q 20 ${pair_id}'aligned_reads.sam' -o $pair_id'.bam'
#samtools sort $pair_id'.bam' -o $pair_id'_sort.bam' -O bam
#samtools index -b $pair_id'_sort.bam'

"""
}*/



/*bg_filt = Channel.fromFilePairs(params.bg_filts)
params.outdir3 = "bg_sort_bam_files"


process bg_bwa {

publishDir params.outdir3, mode: 'copy'

input:
val ref
tuple val(pair_id), path(bg_filt)

output:
path "${pair_id}_sort.bam", emit: bg_bam

"""
#!/bin/env bash

module load $BWA
module load $SAMTOOLS

bwa index -a bwtsw $ref

bwa aln -t 8 $ref ${bg_filt[0]} > ${bg_filt[0]}'reads_1.sai'
bwa aln -t 8 $ref ${bg_filt[1]} > ${bg_filt[1]}'reads_2.sai'

bwa sampe $ref ${bg_filt[0]}'reads_1.sai' ${bg_filt[1]}'reads_2.sai' ${bg_filt[0]} ${bg_filt[1]} \
> ${pair_id}'aligned_reads.sam'

samtools view -b -h -q 20 ${pair_id}'aligned_reads.sam' -o $pair_id'.bam'

samtools sort $pair_id'.bam' -o $pair_id'_sort.bam' -O bam
samtools index -b $pair_id'_sort.bam'

"""


}*/



// now i want to make a process for creating tag directories
// need to make the files into a list and sort by name
// then the transpose operator will make them into a tuple
// the combine operator does it for the file paths but it is not in order

//kd_bam_list = Channel.fromPath("store_normal_bam_files/454*bam").toSortedList{it.name}
//ctr_bam_list = Channel.fromPath("store_normal_bam_files/455*bam").toSortedList{it.name}

/*
bam_files = Channel.fromFilePairs("store_normal_bam_files/{454,455}-Mock-n{1,2,3}*_S{1,2,15,16,29,30}_L00{1,2}_sort.bam")
bam_files.view()

process bam_files {


"""
cd store_normal_bam_files

ls 454* > kd_bam.txt
ls 455* > ctr_bam.txt
paste kd_bam.txt ctr_bam.txt > kd_ctr_bam.txt 

"""

}*/


//bam_tuple = kd_bam_list.combine(ctr_bam_list)


//kd_bam = Channel.fromList(kd_bam_list)
//ctr_bam = Channel.fromList(ctr_bam_list)



//kd_ctr_pair = kd_bam_list.merge(ctr_bam_list)




workflow{

// look what channel.value does. it lets me use the ref multiple times
ref = Channel.value(params.ref)
gtf = Channel.value(params.gtf)



//PE_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

//PE_reads.view()
//filt_pe.view()
//bg_filt.view()
//kd_bam_list.view()
//bam_tuple.view()
main:

cat_tech_reps(tech_rep_reads_r1, tech_rep_reads_r2)



cat_tech_reps.out.fastq_cat_reads.flatten()
.map{file -> 
    def basename = file.baseName
    def key = basename.replaceAll(/_R[12]*/, '')
    return tuple(key,file)
}
.groupTuple()
.set{ fastq_tuples}

fastp(fastq_tuples)

//fastp(PE_reads)
//bg_fastp(bg_reads)

star(ref, gtf)



fastp.out.filt_fq_files.flatten()
.map{file ->
    def basename = file.baseName
    def key = basename.replaceAll(/_R[12]*/, '')
    return tuple(key,file)
}
.groupTuple()
.set { filt_files_tuple}

star_align(filt_files_tuple, star.out.genome_gen_finished)

//star_align(filt_pe, star.out.genome_gen_finished)

organize(star_align.out.star_bam_files.collect())
r_featurecounts(organize.out.place_holder)

/*star_align.out.star_bam_files
.map{ file -> 
def basename = file.baseName
def key = basename.replaceAll(/_L00[12]* /, '')
return tuple(key, file) 
}
.groupTuple()
.set { bam_files_grouped }

//bam_files_grouped.view()

//star_align.out.star_bam_files.collect()

organize( bam_files_grouped)


organize.out.merged_bams
.map{ file ->
return file.name
}
.set { only_file_names }

only_file_names.view()

collected_bams( only_file_names.collect())



//r_featurecounts(gtf, star_align.out.star_bam_files)


r_featurecounts(collected_bams.out.txt_merged_bams)

*/


// looking to see the output in r_featurecounts

//r_featurecounts.out.feature_counts_files.view()

//bg_bwa(ref, bg_filt)
//homer( bam_tuple)
//fastp.out.fastp_out_f.view()
//fastp.out.fastp_out_r.view()
//bwa.out.bam_file.view()

}
