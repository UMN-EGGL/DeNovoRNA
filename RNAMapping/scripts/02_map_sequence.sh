 docker run -it \
    --mount type=bind,source=/project/Data/,target=/data \
    --mount type=bind,source=/scratch/sshfs/data_release/hiseq/151006_D00635_0082_BC7HAHANXX/Project_McCue_Project_022,target=/fastq \
    rnamap:latest \
    --runThreadN 8 \
    --genomeDir /data/Fasta/STARIndices \
    --readFilesIn /fastq/11F_ATGTCA_L007_R1_001.fastq /fastq/11F_ATGTCA_L007_R2_001.fastq \
    --outFileNamePrefix /data/BAMs/STAR/ \
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate > 11F.bam

