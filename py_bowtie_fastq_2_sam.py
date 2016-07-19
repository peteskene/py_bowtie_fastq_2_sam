# %load /home/pskene/bin/py_bowtie_fastq_2_sam.py
from subprocess import check_output
import os
import glob
import shutil
from datetime import datetime
    
def py_bowtie_fastq_2_sam(input_type='fastq.gz', manual_entry=False, list_R1=[], list_R2=[],folder=None,
                          output_folder=None, data_output_name=None,  spike_output_name=None, barcode=None,
                          data_species='hg19', spike_species='dm6', align_spike=True):
    """
    Written by Pete Skene (peteskene@gmail.com). Free for academic use only.
    
    Will take group of fastq files from a single sample from paired end sequencing and perform mapping based on standard Henikoff lab parameters.
    Option for aligning fastq files to 2 species (data_species and spike_species)
    
    Expected naming convention for read 1: e.g. PS_HsDm_CTCF_1m_ATCACG_L001_R1_001.fastq.gz
                                   read 2: e.g. PS_HsDm_CTCF_1m_ATCACG_L001_R2_001.fastq.gz
                                           i.e. data_barcode_R1/2_#.fastq.gz
                                           
    The script will assign R1/R2 files based on searching for string '_R1_' or '_R2_'. If this doesn't fit your filenames, set manual_entry=True and provide
    list_R1=['A.fastq.gz', 'B.fastq.gz'], list_R2=['C.fastq.gz', 'D.fastq.gz'] (Can also take .fastq as list_R1 or list_R2). This assignment occurs after the
    input_type files in the folder have been decompressed if provided as fastq.gz. If manual_entry is used data_output_name and spike_output_name must be specified.
    
    Manual entry can also be used if the folder contains multiple fastq(.gz) files, but only a subset are to analyzed.
    
    Defaulted output sam filename: PS_HsDm_CTCF_1m.sam (i.e. data.sam )
    
    Currently supported bowtie2Index: human hg19 (type 'hg19')
                                      Drosophila dm6 (type 'dm6')
                                      Budding yeast (type 'sacCer3')
                                      mouse mm9 (type 'mm9')
    ____________________
    Parameters:
    
    input_type: type either 'fastq' or 'fastq.gz'
    
    folder: path to folder containing fastq files. e.g. '/home/pskene/test_fastq', 
            if 'None' then looks in current working directory. 
            note this changes current working directory so output files will be to 'folder'.
    
    output_folder: path to folder where to save a copy of the output files. If None specified, it will 
                   default to the current working directory (that may have been specified by 'folder')
                   this is to allow collation of outputs in one place if multiple groups of fastq files in
                   different locations. Output_folder must already exist, as script will not make a new directory
    
    data_output_name: option of specifying output filename e.g. 'data.sam'. If None it will create the default name from input filenames
    
    spike_output_name: option of specifying output filename e.g. 'data.sam.dm6'. If None it will create the default name from input filenames
    
    barcode: specify barcode used in sequencing as string e.g. 'ATCACG'. Only needs to be provided if output_name=None,
             as it is used to to create default output_name
    
    data_species: genome build for the experimental organism. This is used to find the bowtie2 index file. If the file is unavailable the script aborts.
    
    spike_species: genome build for the spike DNA organism. This is used to find the bowtie2 index file. If the file is unavailable the script aborts.
    
    align_spike: default to True. Option for aligning fastq to 2nd species (spike_species), as well as data_species
                   
    
    ____________________
    Requirements:
    
    module load bowtie2/2.2.5 (needs to be typed into terminal before invoking jupyter instance)
    
    """
    
    startTime = datetime.now()
    
    print 'WARNING: bowtie2 must be loaded into rhino terminal before invoking jupyter instance'\
    ' by typing: module load bowtie2/2.2.5'
    print
    
    #change directory as instructed
    if folder != None:
        os.chdir(folder)
        print 'Will look for files in: ' + os.getcwd()
    
    #check input type
    if input_type == 'fastq' or input_type == 'fastq.gz':
        print 'Files will be imported as .' + input_type
    else:
        return "Input_type not recognized, please type either 'fastq' or 'fastq.gz' as a string (i.e. including '')"
    
    
    #need to separate files into R1 and R2, create empty lists to be filled automatically or by manual entry
    input_files_R1 = []
    input_files_R2 = []
    
    #import input filenames as a list and print out the files
    if manual_entry==False:
        input_files = sorted(glob.glob('*.'+input_type))
    
        if len(input_files)==0:
            return 'No files imported. Check input_type. Exiting...'
    
        print '\n'    
        print 'Input files imported as ' + input_type + ':'
        print '\n'.join(input_files)
        print '\n'
    

    
        if input_type == 'fastq.gz':
            print 'Decompressing .gz files'
        
            #create fastq filename from .fastq.gz
            temp_list = [f.replace('.gz', '') for f in input_files]
        
            #run system command to decompress file using zcat
            for i in range(len(input_files)):
                temp_str  = ('zcat ' + input_files[i] + ' > ' + temp_list[i])
                check_output(temp_str, shell=True)
            
            #replace list of input_files with temp_list as future operations on decompressed files
            input_files = temp_list
        
            print 'Decompressed files:'
            print '\n'.join(input_files)
            print '\n'
    
   
        for item in input_files:
            if '_R1_' in item and '_R2_' in item :
                return "Input file: " + item + " contains both strings 'R1' and 'R2'. Exiting..."
            elif '_R1_' in item and '_R2_' not in item :
                input_files_R1.append(item)
            elif '_R1_' not in item and '_R2_' in item:
                input_files_R2.append(item)
            else:
                return "Input file: " + item + " does not contain string 'R1' or 'R2. Exiting..."
            
        print 'Files assigned as R1:'
        print '\n'.join(input_files_R1)
        print '\n'
    
        print 'Files assigned as R2:'
        print '\n'.join(input_files_R2)
        print '\n'
    
        if len(input_files_R1) != len(input_files_R2):
            return 'Unequal numbers of files assigned as R1 and R2. Check naming convention. Exiting...'
    
        if (len(input_files_R1) + len(input_files_R2)) != len(input_files):
            return 'Not all of input files assigned as R1 or R2. Exiting...'
    
    if manual_entry==True:
        print 'Manual entry set to true, will take user-specified list_R1 and list_R2. Assumes the boths lists are ordered the same so they match'
        
        input_files_R1 = list_R1
        input_files_R2 = list_R2
        
        print 'User specified list_R1: '
        print '\n'.join(input_files_R1)
        print '\n'
        print 'User specified list_R2: '
        print '\n'.join(input_files_R2)
        print '\n'
        
        if len(input_files_R1)!=len(input_files_R2):
            return 'User specified list_R1 and list_R2 are of unequal length. Exiting...'

        
        if input_type == 'fastq.gz':
            print 'Decompressing .gz files in list_R1 and list_R2'
        
            #create fastq filename from .fastq.gz
            temp_list_R1 = [f.replace('.gz', '') for f in input_files_R1]
            temp_list_R2 = [f.replace('.gz', '') for f in input_files_R2]
            
            #run system command to decompress file using zcat
            for i in range(len(input_files_R1)):
                temp_str_R1  = ('zcat ' + input_files_R1[i] + ' > ' + temp_list_R1[i])
                check_output(temp_str_R1, shell=True)
                
                temp_str_R2  = ('zcat ' + input_files_R2[i] + ' > ' + temp_list_R2[i])
                check_output(temp_str_R2, shell=True)
                
            #replace list of input_files with temp_list as future operations on decompressed files
            input_files_R1 = temp_list_R1
            input_files_R2 = temp_list_R2
        
            
            print 'Decompressed files from list_R1:'
            print '\n'.join(input_files_R1)
            print '\n'
            print 'Decompressed files from list_R2:'
            print '\n'.join(input_files_R2)
            print '\n'
            print 'Will run bowtie2 on these user-specified lists. Please ensure the orders of each list match.'
    
    #create output .sam and .sam.SPIKE filenames
    if data_output_name != None:
        sam_name = data_output_name
        
    else:
        if manual_entry==False and (barcode==None or type(barcode)!=str):
            return 'If output_name is not specified, barcode must be provided to generate ouput filename for sam file'
        if manual_entry==True:
            return 'If manual entry is set to True, output names must be specified'
            
        
    if spike_output_name != None:
        sam_name_spike = spike_output_name

    elif align_spike == True:
        if manual_entry==False and (barcode==None or type(barcode)!=str):
            return 'If output_name is not specified, barcode must be provided to generate ouput filename for sam file'
        if manual_entry==True:
            return 'If manual entry is set to True, output names must be specified'
    
    print 'Output sam filename for data_species: ' + sam_name
    
    if align_spike==True:
        print 'Output sam filename for spike_species: ' + sam_name_spike
    print '\n'
	
   
    #####################################################
    #run bowtie on the paired files for the data_species first
    #find the correct bowtie2 index
    print 'Data_species set to: ' + data_species
    if data_species == 'hg19':
        bowtie2_data_index = '/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
    
    elif data_species == 'dm6':
        bowtie2_data_index = '~/henikoff/solexa/Bowtie2/dmel_r6_06'
        
    elif data_species == 'sacCer3':
        bowtie2_data_index = '/shared/biodata/ngs/Reference/iGenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome'
        
    elif data_species == 'mm9':
            bowtie2_data_index = '/shared/biodata/ngs/Reference/iGenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome'
        
    elif data_species not in ['hg19', 'dm6', 'sacCer3', 'mm9']:
        return 'Bowtie index unavailable for this species, need to create bowtie index separately'
    
    print 'Bowtie index for data files found at:'
    print bowtie2_data_index
    print '\n'
       
    #run bowtie
    for i in range(len(input_files_R1)):
        if i == 0: #for first file in R1 list, create the sam file
            print 'count = ' +str(i)
            print '\n'
            #create string for system command
            temp_str = 'bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12' \
            + ' -x ' + bowtie2_data_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' > ' + sam_name
            
            print temp_str
            
            check_output(temp_str, shell=True)
            print '\n'
            
            
        else: #for rest of files in R1 list, append to the sam file
            print 'count = ' +str(i)
            print '\n'
            temp_str = 'bowtie2 --no-head --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12' \
            + ' -x ' + bowtie2_data_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' >> ' + sam_name
            
            
            check_output(temp_str, shell=True)
            print temp_str
            print '\n'
            
    #####################################################           
    #run bowtie on the paired files for the spike_species
    #find the correct bowtie2 index
    if align_spike==True:
        print 'Spike_species set to: ' + spike_species
        if spike_species == 'dm6':
            bowtie2_spike_index = '~/henikoff/solexa/Bowtie2/dmel_r6_06'
     
        elif spike_species == 'hg19':
            bowtie2_spike_index = '/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
    
        elif spike_species == 'sacCer3':
            bowtie2_spike_index = '/shared/biodata/ngs/Reference/iGenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome'
            
        elif spike_species == 'mm9':
            bowtie2_spike_index = '/shared/biodata/ngs/Reference/iGenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome'
        

        elif spike_species not in ['hg19', 'dm6', 'sacCer3', 'mm9']:
            return 'Bowtie index unavailable for this species, need to create bowtie index separately'
    
    
        print 'Bowtie index for spike files found at:'
        print bowtie2_spike_index
        print '\n' 
        
        #run bowtie
        for i in range(len(input_files_R1)):
            if i == 0: #for first file in R1 list, create the sam file
                print 'count = ' +str(i)
                print '\n'
                #create string for system command
                temp_str = 'bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12' \
                + ' -x ' + bowtie2_spike_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' > ' + sam_name_spike
            
                print temp_str
            
                check_output(temp_str, shell=True)
                print '\n'
            
            
            else: #for rest of files in R1 list, append to the sam file
                print 'count = ' +str(i)
                print '\n'
                temp_str = 'bowtie2 --no-head --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12' \
                + ' -x ' + bowtie2_spike_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' >> ' + sam_name_spike
            
            
                check_output(temp_str, shell=True)
                print temp_str
                print '\n'
    
    if align_spike==False:
        print 'align_spike set to false. Only aligning to data_species.'
    
    
    if output_folder != None:
        print 'Copying sam files to output folder: ' + output_folder
        shutil.copy(sam_name, output_folder + sam_name)
        if align_spike==True:
            shutil.copy(sam_name_spike, output_folder + sam_name_spike)
    
    print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
    
    return
