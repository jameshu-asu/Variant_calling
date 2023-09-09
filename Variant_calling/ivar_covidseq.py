import shutil
import glob
import subprocess
import os
import pandas as pd

'''by James C. Hu
This script will:
1) peform ivar variant calling on all bam files within a directory.
2) Pull and combine variant data based on locations given by input file.
'''

# pandas terminal output options
pd.options.display.max_columns = 15
pd.options.display.width = 1000


def run_iVar_varaints(in_directory: str, out_directory: str, ref_gene_fa: str, reg_gene_gff3: str) -> None:
    '''
    This function will peform ivar variant analysis on all bam files within a given directory.
    '''
    for file in os.listdir(in_directory):
        if file.endswith('.bam'):
            subprocess.call(
                f'samtools mpileup -d 0 -A -aa -q 0 -Q 0 -R {in_directory}/{file} | ivar variants -p {out_directory}/{file[:7]}_iVAR -q 20 -t 0 -r {ref_gene_fa} -g {reg_gene_gff3}', shell=True)
    return None


def iVar_variant_search(in_file: str, target_nts: list) -> pd.DataFrame:
    '''
    This function will take nt positional arguments as an input and consolidate the ivar outputs for those locations.
    '''
    seq_id = in_file.split('/')[2].split('_')[0]
    df_in = pd.read_csv(in_file, sep='\t')
    df_out = df_in[df_in['POS'].isin(target_nts)]
    df_out['REGION'] = seq_id
    if len(df_out) > 0:
        return df_out


os.mkdir('../bam_files')
os.chdir('../')
bam_files = glob.glob('**/I*/*sortedTrimmed.bam')
for file in bam_files:
    shutil.copy(file, 'bam_files')
os.chdir('ivar_output')

os.mkdir('../ivar_output')
run_iVar_varaints('../bam_files', '../ivar_output', 'MN908947.3.fasta', 'MN908947.3.gff3')

target_nts = pd.read_csv('nt_positions_infile.csv')['NT_position'].copy().to_list()
ivar_output_files = [file for file in os.listdir('../ivar_output') if file.endswith('iVAR.tsv')]
ivar_output_files.sort()

df_out = pd.DataFrame()

log_file = '../ivar_output/ivar_log.txt'

for file in ivar_output_files:
    print('\n============================')
    print(f'Searching for mutations in: {file}')
    print('============================\n')
    temp_df = iVar_variant_search(f'../ivar_output/{file}', target_nts)
    df_out = pd.concat([df_out, temp_df])

df_out = df_out.set_index('REGION')
df_out = df_out.sort_values('POS').sort_values('REGION')

df_out.to_csv('../ivar_output/combined_ivar_output.csv')
