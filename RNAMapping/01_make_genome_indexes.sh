docker run -it \
    --mount type=bind,source=/project/Data/,target=/data \
    rnamap:latest \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir /data/Fasta/STARIndices \
    --genomeFastaFiles /data/Fasta/EquCab2/EquCab2.fa \
    --sjdbGTFfile /data/GFFs/ref_EquCab2.0_top_level.gff3 \
    --sjdbGTFtagExonParentTranscript Parent
