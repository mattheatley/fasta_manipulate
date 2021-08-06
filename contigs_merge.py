#!/usr/bin/python
# -*- coding: utf-8 -*-
# CONDA ENV SOFTWARE: python (3.7), samtools (v1.9), gatk4 (v4.1.4.0)

import os, sys, subprocess, shutil, argparse
from core import FindSupplementaryFile, CAPTURE

pipe_script = os.path.realpath(__file__) # extract script path
script_name, *arguments = sys.argv # extract command line arguments

split_script = f'{os.path.dirname(pipe_script)}/fasta_split.sh' # specify module script
assert(os.path.exists(split_script)), f'Problem finding the module "{split_script}".'# check module exists

parser = argparse.ArgumentParser(description='Merge contigs', prog=script_name, usage='%(prog)s [options]', epilog='see the readme file for further details.')

inputs = parser.add_argument_group(description='user inputs:') # user inputs
inputs.add_argument('-fasta', metavar='</path/to/directory>', type=str, required=True, help='specify fasta path')
inputs.add_argument('-contigs', metavar='<number>', type=int, default=500, help='specify final contig number')
inputs.add_argument('-delimiter', metavar='<character>', type=str, default='N', help='specify delimiter symbol')
inputs.add_argument('-length', metavar='<number>', type=int, default=500, help='specify delimiter number')

fasta_dir, total_merged_contigs, delimiter_character, delimiter_length = vars(parser.parse_args()).values() # define user inputs

fasta_dir = f'/{fasta_dir.strip("/")}' # ensure correct path format
assert(os.path.exists(fasta_dir)), f'Problem finding the reference directory "{fasta_dir}".'# check path exists

# file info

original_fasta_file = FindSupplementaryFile(fasta_dir, '.fasta') # find reference fasta

original_fasta_basename = os.path.basename(original_fasta_file) # extract reference basename

split_name = 'split'  # specify split fasta subdir name

split_dir = f'{fasta_dir}/{split_name}' # specify split fasta subdirectory path

file_names = [ # specify output file names
    'smallest_contigs_combined'
    ,f'sorted_{original_fasta_basename}'
    ,f'merged_{original_fasta_basename}' ]

output_files = [ f'{fasta_dir}/{name}' for name in file_names ] # specify output file paths

merged_contig_name, sorted_fasta_name, rebuilt_fasta_name = file_names # extract individual file names

merged_small_contig_file, sorted_fasta_file, rebuilt_fasta_file = output_files # extract individual file paths

large_contig_subdir, small_contig_subdir = size_subdirs = [ f'{split_dir}/{label}' for label in ['large','small'] ] # specify contig size subdirectory paths

merge_cmd = ( # specify command to merge smallest contigs
f'''awk -v delimiter=`printf \'{delimiter_character}%.0s\' {{1..{delimiter_length}}}` \\
\'BEGIN {{ print \">{merged_contig_name}\" }} \\
FNR>1 {{ printf \"%s%s\", seq, $0 ; seq=delimiter }} \\
END {{ print \"\" }} \' \\
{small_contig_subdir}/*.fasta > {merged_small_contig_file}'''
            )
recombine_sorted_cmd = (f'cat {large_contig_subdir}/*.fasta {small_contig_subdir}/*.fasta > {sorted_fasta_file}') # specify command to recombine fasta with contigs sorted by size

recombine_merged_cmd = (f'cat {large_contig_subdir}/*.fasta {merged_small_contig_file} > {rebuilt_fasta_file}') # specify command to recombine fasta with largest contigs seperate & smallest contigs merged

total_contigs = int(CAPTURE(f'grep -c ">" {original_fasta_file}')) # calculate total number of contigs

print(f'\n   ORIGINAL CONTIGS: {total_contigs:,}')

print(f'\n   MERGED CONTIGS: {total_merged_contigs:,}')

if total_contigs <= total_merged_contigs: print(f'\n   REFERENCE ALREADY HAS THE REQUIRED NUMBER OF CONTIGS\n')
else:

    
    os.makedirs(split_dir, exist_ok=True) # create subdirectory as required

    print('\n   INDEXING REFERENCE...')
    subprocess.call(f'samtools faidx {original_fasta_file}', shell=True) # index reference fasta
    print('\tCOMPLETE')

    print('\n   SPLITTING REFERENCE...')
    subprocess.run(f'bash {split_script} {original_fasta_file} {split_dir}', shell=True, stdout=subprocess.DEVNULL) # split reference fasta via module
    print('\tCOMPLETE')


    contigs_to_retain = total_merged_contigs-1 # specify number of original contigs to retain i.e. not to be merged
    
    [ os.makedirs(subdir, exist_ok=True) for subdir in size_subdirs ] # create subdirectories as required

    fai_file = FindSupplementaryFile(fasta_dir, '.fai') # find reference index

    faidx_info = [ line.strip(' \n').split('\t') for line in open(fai_file, 'r').readlines() ] # extract index info

    contig_info = [ ( name, int(length) ) for name, length, *_ in faidx_info ] # extract contig names & lengths

    sorted_contig_info = sorted( contig_info, key=lambda info: info[1], reverse=True ) # sort contig info by length; largest -> smallest

    contig_names, *_ = zip(*sorted_contig_info) # discard contig lengths

    contig_split_by_size = [ # split contig info by size
        ( contig_names[ :contigs_to_retain ], large_contig_subdir ) # largest contig info
        ,( contig_names[ contigs_to_retain: ], small_contig_subdir ) # smallest contig info
                ]

    padding = len(str(total_contigs)) # calculate leading zeroes for index label

    contig_files = { contents.name.replace('.fasta',''):contents.path for contents in os.scandir(split_dir) if contents.name.endswith('.fasta') } # list contig files

    print('\n   ORGANISING CONTIGS...')
    for contigs, subdir in contig_split_by_size:

        for name in contigs:

            contig_file = contig_files[name] # extract relevant contig file

            new_name = str( contig_names.index(name)+1 ).zfill( padding ) # specify contig index label

            renamed_contig_file = f'{subdir}/{new_name}.fasta' # specify renamed contig file
            
            subprocess.call(f'mv {contig_file} {renamed_contig_file}', shell=True) # move contig to relevant subdirectory & rename with index label

    subprocess.call( merge_cmd, shell=True ) # merge smallest contigs together
    print('\tCOMPLETE')

    print('\n   REBUILDING REFERENCE...')
    subprocess.call( recombine_sorted_cmd, shell=True ) # recombine sorted reference fasta

    subprocess.call( recombine_merged_cmd, shell=True ) # recombine merged reference fasta
    print('\tCOMPLETE')

    print(f'\n   ALL STAGES COMPLETE\n')
