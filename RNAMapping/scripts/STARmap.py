import os
import minus80 as m80
import locuspocus as lp
import asyncio

#basepath = '/root'

#if not m80.Tools.available('Cohort','EMS_Muscle_Fat'):
#    ems_cohort = m80.Cohort.from_yaml(
#		    'EMS_Muscle_Fat',
#		    os.path.join(basepath,'data/MDB.yaml')
#		)
#else:
#    ems_cohort = m80.Cohort('EMS_Muscle_Fat')
#
#genome_dir = '/project/Data/Fasta/STARIndices/EquCab3' 
#out_dir = "/output/"


class STARMap(object):
    '''
        Maps RNASeq based on LinkageIO Cohorts
    '''
    def __init__(self, basepath, genome_dir, out_dir):
        self.basepath = basepath
        self.genome_dir = genome_dir
        self.out_dir = out_dir
        # Define some config stats
        self.sem = asyncio.Semaphore(4)
        self.loop = asyncio.get_event_loop()
        # Allocate this for later
        self.cohort = None

    class wc_protocol(asyncio.SubprocessProtocol):
        '''
            IO protocol for the word cound program
        '''
        def __init__(self,exit_future):
            self.exit_future = exit_future
            self.output = bytearray()
            self.num_lines = 5

        def pipe_data_received(self, fd, data):
            self.output.extend(data)
            #print(f'{data} on {fd}')
            #self.num_lines = int(data.decode('ascii').rstrip().split()[0])

        def process_exited(self):
            self.exit_future.set_result(True)


    async def create_genome_index(self):
        # Create a genome index
        async with self.sem:
            STAR_index_EquCab3_cmd = '''\
                STAR \
                --runThreadN 8 \
                --runMode genomeGenerate \
                --genomeDir /project/Data/Fasta/STARIndices/EquCab3 \
                --genomeFastaFiles /project/Data/Fasta/EquCab3/EquCab3.fasta \
                --sjdbGTFfile /project/Data/GFFs/ref_EquCab2.0_top_level.gff3 \
                --sjdbGTFtagExonParentTranscript Parent \
            '''

    async def fastq_wc(self,f):
        '''
            runs wc on r1 and r2 
        '''
        # grab a future
        async with self.sem:
            exit_future = asyncio.Future(loop=self.loop)
            # Creat the subprocess
            print(f'counting lines for {f}')
            create = self.loop.subprocess_exec(
                lambda: self.wc_protocol(exit_future),
                #'wc', '-l', f,
                'countdown', '5',
                stdin=None,stderr=None
            )
            # Create the future 
            transport, protocol = await create
            # let it do its work
            await exit_future
            transport.close()
            nl = protocol.num_lines
            return nl



    def run(self, cohort):
        self.cohort = cohort
        # Define a task for each sample in the cohort
        tasks = []
        for sample in self.cohort:
            # Create a task for each sample
            task = self.fastq_wc(sample)
            tasks.append(task)
        # Gather the tasks
        results = asyncio.gather(*tasks)
        # run them in the loop
        self.loop.run_until_complete(results)

    def map(self):
        for r1,r2 in zip(R1s,R2s):
            bam_name = os.path.basename(r1).replace('_R1','').replace('.fastq','')
            output_dir = os.path.join(out_dir,bam_name+'/')
            os.makedirs(output_dir,exist_ok=True)
            if not os.path.exists(os.path.join(output_dir,'Aligned.sortedByCoord.out.bam')):
                cmd = f'''\
                    STAR \
                    --runThreadN 3 \
                    --genomeDir {genome_dir} \
                    --readFilesIn {r1} {r2} \
                    --outReadsUnmapped Fastx \
                    --outFileNamePrefix  {output_dir} \
                    --outSAMtype BAM SortedByCoordinate \
                '''
                print(cmd)

    async def ensure_equal_lines(self, sample):
        if len(sample.files) % 2 != 0:
            raise ValueError('The number of FASTQ files must be the SAME!!')
        R1s = [x for x in sample.files if 'R1' in x]
        R2s = [x.replace('R1','R2') for x in R1s]
        for r1,r2 in zip(R1s,R2s):
            num_r1 = await self.fastq_wc(r1)
            num_r2 = await self.fastq_wc(r2)
            if num_r1 != num_r2:
                return False
