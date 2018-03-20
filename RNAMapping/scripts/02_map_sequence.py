import os
import minus80 as m80
import locuspocus as lp

basepath = '/root'

if not m80.Tools.available('Cohort','EMS_Muscle_Fat'):
    ems_cohort = m80.Cohort.from_yaml(
		    'EMS_Muscle_Fat',
		    os.path.join(basepath,'data/MDB.yaml')
		)
else:
    ems_cohort = m80.Cohort('EMS_Muscle_Fat')

genome_dir = '/project/Data/Fasta/STARIndices/EquCab3' 
out_dir = "/output/"


STAR_index_cmd = '''\
    STAR \
	--runThreadN 14 \
	--runMode genomeGenerate \
	--genomeDir /project/Data/Fasta/STARIndices/EquCab3 \
	--genomeFastaFiles /project/Data/Fasta/EquCab3/EquCab3.fasta \
	--sjdbGTFfile /project/Data/GFFs/EquCab3/ref_EquCab3.0_top_level.gff3 \
	--sjdbGTFtagExonParentTranscript Parent 
'''

STAR_index_EquCab3_cmd = '''\
    STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir /project/Data/Fasta/STARIndices/EquCab3 \
    --genomeFastaFiles /project/Data/Fasta/EquCab3/EquCab3.fasta \
    --sjdbGTFfile /project/Data/GFFs/ref_EquCab2.0_top_level.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
'''

with open('commands.txt','w') as OUT:
    for sample in ems_cohort:
        if len(sample.files) % 2 != 0:
            raise ValueError('The number of FASTQ files must be the SAME!!')
        R1s = [x for x in sample.files if 'R1' in x]
        R2s = [x.replace('R1','R2') for x in R1s]

        for r1,r2 in zip(R1s,R2s):
            bam_name = os.path.basename(r1).replace('_R1','').replace('.fastq','')
            output_dir = os.path.join(out_dir,bam_name+'/')
            os.makedirs(output_dir,exist_ok=True)
            #cmd = STAR_cmd.format(r1,r2,bam_name)
            cmd = f'''\
                STAR \
                --runThreadN 3 \
                --genomeDir {genome_dir} \
                --readFilesIn {r1} {r2} \
                --outReadsUnmapped Fastx \
                --outFileNamePrefix  {output_dir} \
                --outSAMtype BAM SortedByCoordinate \
            '''
            print(cmd,file=OUT)
        

  
