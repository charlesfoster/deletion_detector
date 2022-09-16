#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:05:42 2021

@author: https://github.com/charlesfoster
"""

from Bio import SeqIO
from Bio.Seq import Seq
from src import __version__
import os
import subprocess
import argparse
import sys
import multiprocessing as mp
import tqdm
import re
from random import randrange
import shutil
import time
import shlex
from operator import itemgetter

#%%
def parse_cigar(cigar,alignment_start,coordinates):
    alpha_pos = [m.start() for m in re.finditer("[A-Z]+", cigar)]
    x=0
    positions = []
    mtypes = []
    for i in range(0,len(alpha_pos)):
        y = alpha_pos[i]
        bases = cigar[x:y]
        mtype = cigar[y]
        positions.append(bases)
        mtypes.append(mtype)
        x = y+1
    deletions = []
    for x,y in zip(positions,mtypes):
        x = int(x)
        if y != 'D':
            alignment_start+=x
            continue
        else:
            if len(coordinates) > 0:
                if coordinates[0] <= alignment_start <= coordinates[1]:
                    d_s = alignment_start+1
                    d_e = d_s+x-1
                    d_len = x
                    deletion = (d_s, d_e, d_len)
                    deletions.append(deletion)
                    alignment_start+=x
                else:
                    alignment_start+=x
                    continue
            else:
                d_s = alignment_start+1
                d_e = d_s+x-1
                d_len = x
                deletion = (d_s, d_e, d_len)
                deletions.append(deletion)
                alignment_start+=x
    return(deletions)

#%%
def cigar_gaps(record, index, reference, outfile, tempdir, coordinates, deletion_of_interest, parse_gisaid):
    fname = tempdir+'/tmp'+str(index)
    record.seq = record.seq.ungap()
    fixed_name = record.description.replace(" ","_")
    record.id = fixed_name
    record.name = fixed_name
    with open(fname+'.fasta', "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
    cmd = '{0} -c -x asm20 --sam-hit-only --secondary=no --end-bonus 500 -t 1 {1} {2} -o {3}'.format('minimap2', reference, fname+'.fasta',fname+'.paf')
    c = subprocess.check_output(shlex.split(cmd), shell=False,stderr=subprocess.PIPE)
    with open(fname+'.paf','r') as f:
        lines = [line.rstrip().split("\t") for line in f]
        gap_list = []
        for paf in lines:
            wanted = itemgetter(7,-1)(paf)
            cigar = re.sub('^.*:','',wanted[1])
            alignment_start = int(wanted[0])
            gaps = parse_cigar(cigar,alignment_start,coordinates)
            gap_list.append(gaps)
        gaps = [item for sublist in gap_list for item in sublist]
        number_of_gaps = len(gaps)

        if number_of_gaps == 0:
            result = 'None'
        else:
            result = '; '.join(map(str, gaps))
    
        if deletion_of_interest is not None:
            if deletion_of_interest in gaps and number_of_gaps == 1:
                status = 'target_deletion'
            elif deletion_of_interest in gaps and number_of_gaps == 1 > 1:
                status = 'target_deletion_plus'
            elif deletion_of_interest not in gaps and number_of_gaps >= 1:
                status = 'deletion'
            else:
                status = 'no_deletion'
        else:
            if number_of_gaps == 0:
                status = 'no_deletion'
            elif number_of_gaps > 1:
                status = 'deletion_plus'
            else:
                status = 'deletion'

        # generate output result
        status_code = 0
        result_string = f'{record.id}\t{result}\t{number_of_gaps}\t{status}'
        if parse_gisaid == True:      
            try:
                # fix record name
                DOC = [x for x in fixed_name.split("|") if x.startswith("hCoV")==False and x.startswith("EPI")==False][0]
                country = [x for x in fixed_name.split("|") if x.startswith("hCoV")][0].split("/")[1]
                id = [x for x in fixed_name.split("|") if x.startswith("hCoV")][0].split("/")[2]
                result_string = result_string+f'\t{id}\t{country}\t{DOC}'
            except:
                DOC = 'unknown'
                country = 'unknown'
                id = 'unknown'
                result_string = result_string+f'\t{id}\t{country}\t{DOC}'
                status_code += 1
        final_result = result_string+'\n'

        with open(outfile, "a") as output_handle:
            output_handle.write(final_result)
        
        #os.system('rm {0} {1}'.format(fname+'.fasta',fname+'.paf'))
        return(status_code)


#%%
def all_gaps(record, index, reference, outfile, tempdir, coordinates, deletion_of_interest, parse_gisaid):
    fname = tempdir+'/tmp'+str(index)
    record.seq = record.seq.ungap()
    fixed_name = record.description.replace(" ","_")
    record.id = fixed_name
    record.name = fixed_name
    with open(fname+'.fasta', "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
    cmd = f'''
    minimap2 -a -x asm20 --sam-hit-only --secondary=no --end-bonus 500 -t 1 {reference} {fname+'.fasta'} | \
        gofasta sam toma --pad -t1 > {fname+'.aligned.fasta'}
    '''

    c = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
   
    for new_seq in SeqIO.parse(fname+'.aligned.fasta', "fasta"):
        variant = str(new_seq.seq)
        if len(coordinates) > 0:
            start = coordinates[0]-1
            end = coordinates[1]
            variant = variant[start:end]
            
        n_perc = f"{variant.count('N')/len(variant)*100:.2f}"
        if float(n_perc) == 0:
            QC = 'NO_N'
        elif float(n_perc) > 0 and float(n_perc) <10:
            QC = 'LOW_N'
        elif float(n_perc) >= 10 and float(n_perc) <25:
            QC = 'MILD_N'
        elif float(n_perc) >=25 and float(n_perc) <50:
            QC = 'MEDIUM_N'
        elif float(n_perc) >=50 and float(n_perc) <100:
            QC = 'HIGH_N'
        elif float(n_perc) == 100:
            QC = 'FAIL'
    
        gap_perc = f"{variant.count('-')/len(variant)*100:.2f}"
        
        if len(coordinates) == 0:
            stop_position = "NA"
            prot_perc = "NA"
        else:
            normal_len = (len(range(coordinates[0],coordinates[1]))+1)/3
            for_translation = ''.join([x for x in variant if x != '-'])
            if int(len(for_translation)) % 3 == 1:
                for_translation = for_translation+'NN'
                try:
                    stop_position = Seq(for_translation).translate().index('*')+1
                    prot_perc = (stop_position/normal_len) * 100
                except:
                    stop_position = 'None'
                    prot_perc = 'NA'
            elif int(len(for_translation)) % 3 == 2:
                for_translation = for_translation+'N'
                try:
                    stop_position = Seq(for_translation).translate().index('*')+1
                    prot_perc = (stop_position/normal_len) * 100
                except:
                    stop_position = 'None'
                    prot_perc = 'NA'
            else:
                try:
                    stop_position = Seq(for_translation).translate().index('*')+1
                    prot_perc = (stop_position/normal_len) * 100
                except:
                    stop_position = 'None'
                    prot_perc = 'NA'
        
        if prot_perc == 100:
            prot_trunc = False
        elif prot_perc != 100 and prot_perc != "NA":
            prot_trunc = True
        else:
            prot_trunc = "NA"
        
        if len(coordinates) != 0:
            start_bit = coordinates[0]
            end_bit = coordinates[0]-1
            gaps = [(m.start()+start_bit,m.end()+end_bit, len(m.group(0))) for m in re.finditer("-+", variant)]
        else:
            gaps = [(m.start()+1,m.end(), len(m.group(0))) for m in re.finditer("-+", variant)]
        
        number_of_gaps = len(gaps)

        if number_of_gaps == 0:
            result = 'None'
        else:
            result = '; '.join(map(str, gaps))
       
        if deletion_of_interest is not None:
            deletion_of_interest = deletion_of_interest
            if deletion_of_interest in gaps and number_of_gaps == 1:
                status = 'target_deletion'
            elif deletion_of_interest in gaps and number_of_gaps == 1 > 1:
                status = 'target_deletion_plus'
            elif deletion_of_interest not in gaps and number_of_gaps >= 1:
                status = 'deletion'
            else:
                status = 'no_deletion'
        else:
            if number_of_gaps == 0:
                status = 'no_deletion'
            elif number_of_gaps > 1:
                status = 'deletion_plus'
            else:
                status = 'deletion'

        # generate output result
        status_code = 0
        if len(coordinates) == 0:
            result_string = f'{record.id}\t{result}\t{number_of_gaps}\t{status}\t{gap_perc}\t{n_perc}\t{QC}'
        else:
            result_string = f'{record.id}\t{result}\t{number_of_gaps}\t{status}\t{gap_perc}\t{n_perc}\t{QC}\t{stop_position}\t{prot_perc}\t{prot_trunc}'
        if parse_gisaid == True:      
            try:
                # fix record name
                DOC = [x for x in new_seq.id.split("|") if x.startswith("hCoV")==False and x.startswith("EPI")==False][0]
                country = [x for x in new_seq.id.split("|") if x.startswith("hCoV")][0].split("/")[1]
                id = [x for x in new_seq.id.split("|") if x.startswith("hCoV")][0].split("/")[2]
                result_string = result_string+f'\t{id}\t{country}\t{DOC}'
            except:
                DOC = 'unknown'
                country = 'unknown'
                id = 'unknown'
                result_string = result_string+f'\t{id}\t{country}\t{DOC}'
                status_code += 1
        final_result = result_string+'\n'

        with open(outfile, "a") as output_handle:
            output_handle.write(final_result)
        
        os.system('rm {0} {1}'.format(fname+'.fasta',fname+'.aligned.fasta'))
        return(status_code)

#%%
def printc(thing, level):
    '''
    Print in colour :)
    '''
    cols = {'green':'\033[1;32m', 'blue':'\033[96m'}
    col = cols[level]
    print(f"{col}{thing}\033[0m")
    return()


#%%
def main(sysargs = sys.argv[1:]):
    epilog = '''
    Analysis quits without running if outfile already exists and --force not specified.
    Limitations: (a) gap detection is always subject to how easy/difficult a region is to align to the reference genome, (b) if there are large deletions close to the beginning/end of a sequence,
    soft clipping of the alignment might introduce false missing data (padded as 'N's prior to deletion analysis).\n
    '''
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    #Read inputs
    parser.add_argument('fasta', nargs="*", help="Fasta file containing sequences to check")
    parser.add_argument('-c','--coordinates', required=False, action="store", default=None, help='Survey subset of genome, e.g. an open reading frame. Must be specified in the form of start:end (1-based coordinates; inclusive). If specified, stats will be given for protein truncation etc. Currently assumes standard translation table.')
    parser.add_argument('-d','--deletion_of_interest', required=False, default=None, help="Specify a deletion of interest in the form of 'start:end' (1-based coordinates; inclusive). If specified, outfile can be easily filtered to find 'target_deletion'.")
    parser.add_argument('-f','--force', required=False, action="store_true", default=False, help='Force overwrite outfile')
    parser.add_argument('-p','--parse_gisaid', required=False, action="store_true", default=False, help='If analysing fasta files from GISAID, parse out the date of collection and country')
    parser.add_argument('-r','--reference', required=False, help="Reference genome in fasta format", default = os.path.join(os.getcwd(),"MN908947.3.fasta"))
    parser.add_argument('-t','--threads', required=False, action="store", default=os.cpu_count(), help='Number of threads to use in parallel processing. Defaults to all available threads.')
    parser.add_argument('-o','--outfile', required=False, action="store", default='deletion_results.tsv', help='Name of the outfile to store results')
    parser.add_argument("--quick", action="store_true", required=False, default=False, help="Run in 'quick' mode: faster than default mode, but less rich output")
    parser.add_argument("--version",action="version", version=f"v{__version__}")
    args=parser.parse_args()

    printc("\nDeletion Detector", 'green')
    printc("Version: {}\n".format(__version__),'blue')

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    if shutil.which("gofasta") is None or shutil.which("minimap2") is None:
        print('#####\n\033[91mError\033[0m: Missing dependencies\n#####\n')
        
        dependency_text = """
Please make sure that the following programs are installed and in your path:
    - gofasta
    - minimap2
            
Additionally, the following python libraries must be installed:
    - biopython
    - tqdm
        """
        
        print(dependency_text)
        sys.exit(-1)

    if not os.path.exists(''.join(args.fasta)):
        print('#####\n\033[91mError\033[0m: Fasta file does not exist or was not specified\n#####\n')
        print('Please check your input and try again\n')
        sys.exit(-1)
    
    MSA = ''.join(args.fasta)
    
    reference = args.reference

    if not os.path.exists(reference):
        print('#####\n\033[91mError\033[0m: Reference fasta file does not exist\n#####\n')
        print('Please check your input and try again\n')
        sys.exit(-1)

    #parse genome coordinates
    coordinates = []
    if args.coordinates is not None:
        coordinates = args.coordinates
        coordinates = [int(x) for x in coordinates.split(":")]
        if len(coordinates) != 2 or coordinates[0] > coordinates[1]:
            print('#####\n\033[91mError\033[0m: Genomic coordinates not valid\n#####\n')
            print('Must be specified in the form of start:end (1-based coordinates; inclusive)\n')
            sys.exit(-1)

    #parse deletion coordinates
    deletion_of_interest = None
    if args.deletion_of_interest:
        parsed = [int(x) for x in args.deletion_of_interest.split(":")]
        if len(parsed) != 2 or parsed[0] > parsed[1]:
            print('#####\n\033[91mError\033[0m: Deletion coordinates not valid\n#####\n')
            print('Must be specified in the form of start:end (1-based coordinates; inclusive)\n')
            sys.exit(-1)
        else:
            length = len(range(parsed[0],parsed[1]))+1 # account for 'off by one' errors with ranges
            deletion_of_interest = (parsed[0], parsed[1], length)

    outfile = "results.tsv"
    if args.outfile:
        outfile = args.outfile
    
    if os.path.exists(outfile) and not args.force:
        print('#####\n\033[91mError\033[0m: outfile already exists and overwriting not enabled\n#####\n')
        parser.print_help()
        sys.exit(-1)
    elif os.path.exists(outfile) and args.force:
        print('#####\n\033[1;33mWarning\033[0m: outfile already exists - overwriting\n#####\n')
        os.remove(outfile)
    
    #set up header for outfile
    if args.quick:
        header = 'full_seq_name\tgap_stretches\tnumber_of_deletions\tsummary\n'
    else:
        header = 'full_seq_name\tgap_stretches\tnumber_of_deletions\tsummary'
        if len(coordinates)>0:
            header = header + "\taa_stop_position\tprot_perc\tprot_truncated"
        if args.parse_gisaid:
            header = header + "\tshort_seq_id\tcountry\tdate_of_collection"
        header = header+"\n"

    with open(outfile, "a") as output_handle:
        output_handle.write(header)

    tempdir = f'tempdir_{randrange(10,10000):04}'
    os.makedirs(tempdir)

    start_time = time.perf_counter ()
    if args.quick:
        printc("\nSearching {} for deletions using parallel processing in \033[1;33mquick\033[0m mode.\nBe patient - might take a while for large files.\nEven if the progress bar pauses, don't panic.\n".format(''.join(args.fasta)), 'blue')
        with open(MSA) as fd, mp.Pool(mp.cpu_count()) as pool:
            result = pool.starmap(
                cigar_gaps,
                tqdm.tqdm([(record, index, reference, outfile, tempdir, coordinates, deletion_of_interest, args.parse_gisaid) for index, record in enumerate(SeqIO.parse(fd, "fasta"))]))    
        os.rmdir(tempdir)    
    else:
        printc("\nSearching {} for deletions using parallel processing in \033[1;33mrich\033[0m mode.\nBe patient - might take a while for large files.\nEven if the progress bar pauses, don't panic.\n".format(''.join(args.fasta)), 'blue')
        with open(MSA) as fd, mp.Pool(mp.cpu_count()) as pool:
            result = pool.starmap(
                all_gaps,
                tqdm.tqdm([(record, index, reference, outfile, tempdir, coordinates, deletion_of_interest, args.parse_gisaid) for index, record in enumerate(SeqIO.parse(fd, "fasta"))]))    
        os.rmdir(tempdir)
    
    printc("\nResults in: {}".format(outfile), 'blue')
    end_time = time.perf_counter ()
    time_taken = "{:.2f}".format(end_time - start_time)
    printc("Time taken: {} seconds\n".format(time_taken),'blue')

    if(any(result) >0):
        print('#####\n\033[1;33mWarning\033[0m: at least one sequence ID could not be parsed as a GISAID ID\n#####\n')
    sys.exit(0)

if __name__ == '__main__':
    main()
