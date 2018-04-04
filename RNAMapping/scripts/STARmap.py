import os
import sys
import minus80 as m80
import locuspocus as lp
import asyncio

if not m80.Tools.available('Cohort','EMS_Muscle_Fat'):
    ems = m80.Cohort.from_yaml(
	    'EMS_Muscle_Fat',
	    os.path.join('/root/data/MDB.yaml')
	)
else:
    ems = m80.Cohort('EMS_Muscle_Fat')


class wc_protocol(asyncio.SubprocessProtocol):
    '''
        IO protocol for the word cound program
    '''
    def __init__(self,exit_future):
        self.exit_future = exit_future
        self.output = bytearray()
        self.num_lines = None

    def pipe_data_received(self, fd, data):
        self.output.extend(data)
        self.num_lines = int(data.decode('ascii').rstrip().split()[0])

    def process_exited(self):
        self.exit_future.set_result(True)


class STAR_protocol(asyncio.SubprocessProtocol):
    '''
        IO protocol for the STAR mapper
    '''
    def __init__(self,fut):
        self.exit_future = fut
        self.output = bytearray()

    def pipe_data_received(self, fd, data):
        self.output.extend(data)

    def process_exited(self):
        self.exit_future.set_result(True)


class HTSEQ_protocol(asyncio.SubprocessProtocol):
    '''
        IO Protocol for HTSeq
    '''
    def __init__(self,fut):
        self.exit_future = fut
        self.output = bytearray()
        self.counts = dict()

    def pipe_data_received(self, fd, data):
        self.output.extend(data)
        # Decode the data and read in like lines
        for line in data.decode('ascii').split('\n'):
            gene,count = line.rstrip().split('\t')
            self.counts[gene] = count   

    def process_exited(self):
        self.exit_future.set_result(True)


class STARMap(object):
    '''
        Maps RNASeq based on LinkageIO Cohorts
    '''
    def __init__(self, genome_dir, out_dir):
        self.genome_dir = genome_dir
        self.out_dir = out_dir
        # Define some config stats
        self.sem = asyncio.Semaphore(4)
        self.loop = asyncio.get_event_loop()
        # Allocate this for later
        self.cohort = None


    async def create_genome_index(self):
        # Create a genome index
        async with self.sem:
            exit_future = asyncio.Future(loop=self.loop)
            print('Indexing Genome for STAR',file=sys.stdout)
            STAR_index_EquCab3_cmd = '''\
                STAR \
                --runThreadN 7 \
                --runMode genomeGenerate \
                --genomeDir /project/Data/Fasta/STARIndices/EquCab3 \
                --genomeFastaFiles /project/Data/Fasta/EquCab3/EquCab3.fasta \
                --sjdbGTFfile /project/Data/GFFs/ref_EquCab3.0_top_level.gff3 \
                --sjdbGTFtagExonParentTranscript Parent \
            '''

    async def map_sample(self,sample):
        '''
            runs wc on r1 and r2 
        '''
        # bound the number of samples being processed using a semaphore
        async with self.sem:
            # Check that FASTQ files occur in pairs 
            if len(sample.files) % 2 != 0:
                raise ValueError('The number of FASTQ files must be the SAME!!')
            # Generate the files from each paired end read 
            R1s = [x for x in sample.files if 'R1' in x and 'fastq' in x]
            R2s = [x.replace('R1','R2') for x in R1s]
            # Loop through and could the number of lines in each
            for r1,r2 in zip(R1s,R2s):
                # Get the target name
                bam_name = os.path.basename(r1).replace('_R1','').replace('.fastq','')
                # Make an output dir based on the BAM name
                target_dir = os.path.join(self.out_dir,bam_name+'/')
                os.makedirs(target_dir,exist_ok=True)
                # Check to see if we already mapped a BAM
                if not os.path.exists(os.path.join(target_dir,'Aligned.sortedByCoord.out.bam')):
                    print(f'Mapping for {sample.name}')
                    # Get the number of lines in R1
                    print(f'Counting lines for {sample.name}')
                    num_r1 = await self.count_lines(r1)
                    print(f'{r1} has {num_r1} lines')
                    # Get the number of lines in R2
                    num_r2 = await self.count_lines(r2)
                    print(f'{r2} has {num_r2} lines')
                    # Compare the number of lines in R1 and R2
                    if num_r1 != num_r2:
                        raise ValueError(f'{r1} and {r2} must have the same number of lines')
                    # Map the reads using STAR
                    await self.map_paired_end_reads(r1,r2,target_dir)
                print(f'{sample.name} MAPPED!')
                sample.add_file(os.path.join(target_dir,'Aligned.sortedByCoord.out.bam'))
                

    async def count_lines(self,filename):
        '''
            Create a task (subprocess) that counts the number of lines 
            in a file.
        '''
        wc_future = asyncio.Future(loop=self.loop)
        wc = self.loop.subprocess_exec(
            lambda: wc_protocol(wc_future),
            'wc', '-l', filename,
            stdin=None, stderr=None
        )
        transport, protocol = await wc
        await wc_future
        transport.close()
        return protocol.num_lines

    async def map_paired_end_reads(self,r1,r2,target_dir):
        cmd = f'''
            STAR 
            --runThreadN 3 
            --genomeDir {self.genome_dir} 
            --readFilesIn {r1} {r2} 
            --outReadsUnmapped Fastx 
            --outFileNamePrefix  {target_dir} 
            --outSAMtype BAM SortedByCoordinate 
        '''.split()
        STAR_future = asyncio.Future(loop=self.loop)
        STAR = self.loop.subprocess_exec(
            lambda: STAR_protocol(STAR_future),
            *cmd,
            stdin=None,stderr=None
        )
        trans,prot = await STAR
        await STAR_future
        trans.close()

    async def count_reads(self,bam_file,gff_file):
        cmd = f'''
            htseq-count 
                -i ID 
                -f bam 
                -t exon 
                {bam_file}
                {gff_file}
        '''.split()
        HTSEQ_future = asyncio.Future(loop=self.loop)
        HTSEQ = self.loop.subprocess_exec(
            lambda: HTSEQ_protocol(HTSEQ_future),
            *cmd,
            stdin=None, stderr=None
        )
        trans,prot = await HTSEQ
        await HTSEQ_future
        trans.close()
        return prot.counts

    async def count_sample(self,sample):
        '''
            Count the number of reads per gene for a sample
        '''
        # Bound the number of task we are going to do
        async with self.sem:
            bams = [x for x in sample.files if x.endswith('bam')]
            for i,bam in enumerate(bams): 
                print(f"Mapping {bam}")
                counts = await self.count_reads(
                    bam,
                    '/project/Data/Fasta/EquCab3/scripts/EquCab3_nice.gff'
                )
                import ipdb; ipdb.set_trace() 



    def run(self, cohort):
        self.cohort = cohort
        # STEP 1 - Map all the fastq files for a sample
        # =============================================
        print('STEP 1')
        tasks = []
        # Only do 4 STAR maps at a time
        self.sem = asyncio.Semaphore(4)
        samples = list(self.cohort)
        # Define a task for each sample in the cohort
        for sample in samples:
            # Create a task for each sample
            task = self.map_sample(sample)
            tasks.append(task)
        # Gather the tasks
        results = asyncio.gather(*tasks)
        # run them in the loop
        self.loop.run_until_complete(results)

        # STEP 2 - Count all the reads for a sample
        # =========================================
        print('STEP 2')
        # Clear tasks
        tasks = []
        # Do 12 HTSeq runs at a time
        self.sem = asyncio.Semaphore(12)
        for sample in samples:
            task = self.count_sample(sample)
            tasks.append(task)
        results = asyncio.gather(*tasks)
        self.loop.run_until_complete(results)

            
