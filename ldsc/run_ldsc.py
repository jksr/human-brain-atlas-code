# -*- coding: future_fstrings -*-

import sys
import os
from pathlib2 import Path
import multiprocessing
import glob
import pandas as pd
#import subprocess


'''run_ldsc bedlist outdir'''

ctdf = pd.read_csv(sys.argv[1], names=['celltype','path'], sep='\t')

outdir = Path(sys.argv[2])
outdir.mkdir(exist_ok=True)
outdir = outdir.resolve()
outbeddir = outdir/'hg19bed'
outbeddir.mkdir(exist_ok=True)


outrundir = outdir/'run'
outrundir.mkdir(exist_ok=True)
outrltdir = outdir/'rlt'
outrltdir.mkdir(exist_ok=True)

out_cmds_annot_file = None
out_cmds_ldsc_file = None
out_cmds_seg_file = None
out_cmds_sep_h2_file = None
if len(sys.argv)>3:
    out_cmds_annot_file = str(outdir/'annot.cmds')
    out_cmds_ldsc_file = str(outdir/'ldsc.cmds')
    out_cmds_seg_file = str(outdir/'seg.cmds')
    out_cmds_sep_h2_file = str(outdir/'seph2.cmds')


hg19chain = '/gale/netapp/home/wtian/refs/ldsc/hg38ToHg19.over.chain.gz'
ldscbin = '/gale/netapp/home/wtian/local/ldsc/'
baseline = '/gale/netapp/home/wtian/refs/ldsc/1000G_EUR_Phase3_baseline/baseline.'
weights = '/gale/netapp/home/wtian/refs/ldsc/weights_hm3_no_hla/weights.'
freqs = '/gale/netapp/home/wtian/refs/ldsc/1000G_Phase3_frq/1000G.EUR.QC.'
chroms = list(range(1,23))


ldctsfile = f'{outdir}/ldcts.ldcts'
sumstats = {x.split('/')[-1].replace('.sumstats.gz',''):x for x in glob.glob('/gale/netapp/home/wtian/refs/ldsc/GWAStraits/*.sumstats.gz')}


def runcmd(cmd):
    print(cmd)
    os.system(cmd)


## prepare signal lifted bed
hg19beds = []
for _,(ct,hg38bed) in ctdf.iterrows():
    hg19bed = outbeddir/f'{ct}.bed'
    cmd_lift = f'liftOver {hg38bed} {hg19chain} {hg19bed} /dev/null'
    runcmd(cmd_lift)
    hg19beds.append(str(hg19bed))

## merge all signal beds for background bed
cmd_prepbg = f"cat {' '.join(hg19beds)} | bedtools sort -i - | bedtools merge -i - > {outbeddir}/_BACKGROUND.bed"
runcmd(cmd_prepbg)


## compute ld score for signals and background
all_cmds_ldsc = []
all_cmds_annot = []
cts = [ct for _,(ct,_) in ctdf.iterrows()] + ['_BACKGROUND']
for ct in cts:
    hg19bed = outbeddir/f'{ct}.bed'
    for chrom in chroms:
        bimfile = f'/gale/netapp/home/wtian/refs/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim'
        snpfile = f'/gale/netapp/home/wtian/refs/ldsc/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.{chrom}.snp'
        annotfile = f'{outrundir}/{ct}.{chrom}.annot.gz'
        outfile = f'{outrundir}/{ct}.{chrom}'

        cmd_annot = f'python {ldscbin}/make_annot.py --bed-file {hg19bed} --bimfile {bimfile} --annot-file {annotfile}'
        cmd_ldsc = f'python {ldscbin}/ldsc.py --l2 --bfile {bimfile.replace(".bim","")} --ld-wind-cm 1 --annot {annotfile} --thin-annot --out {outfile} --print-snps {snpfile}'
        #runcmd(cmd_annot)
        #runcmd(cmd_ldsc)
        all_cmds_annot.append(cmd_annot)
        all_cmds_ldsc.append(cmd_ldsc)



with open(ldctsfile,'w') as f:
    for _,(ct,_) in ctdf.iterrows():
        f.write(f'{ct}\t{outrundir}/{ct}.,{outrundir}/_BACKGROUND.\n')

all_cmds_seg = []
for tr,sspath in sumstats.items():
    cmd_ldscseg = f'python {ldscbin}/ldsc.py --h2-cts {sspath} --ref-ld-chr {baseline} '\
            f'--out {outrltdir}/{tr} --ref-ld-chr-cts {ldctsfile} --w-ld-chr {weights}'
    #runcmd(cmd_ldscseg)
    all_cmds_seg.append(cmd_ldscseg)

all_cmds_sep_h2 = []
for _,(ct,_) in ctdf.iterrows():
    for tr,sspath in sumstats.items():
        cmd_ldsch2 = f'python {ldscbin}/ldsc.py --h2 {sspath}  --ref-ld-chr {outrundir}/{ct}.,{baseline}  --out {outrltdir}/sep__{ct}.{tr}  --w-ld-chr {weights} --overlap-annot --frqfile-chr {freqs}'
        #runcmd(cmd_ldsch2)
        all_cmds_sep_h2.append(cmd_ldsch2)


if out_cmds_ldsc_file is not None and out_cmds_ldsc_file is not None and out_cmds_seg_file is not None and out_cmds_sep_h2_file is not None:
    with open(out_cmds_annot_file,'w') as f:
        f.writelines('\n'.join(all_cmds_annot))
    with open(out_cmds_ldsc_file,'w') as f:
        f.writelines('\n'.join(all_cmds_ldsc))
    with open(out_cmds_seg_file,'w') as f:
        f.writelines('\n'.join(all_cmds_seg))
    with open(out_cmds_sep_h2_file,'w') as f:
        f.writelines('\n'.join(all_cmds_sep_h2))
else:
    for cmd in all_cmds_annot + all_cmds_ldsc + all_cmds_seg + all_cmds_sep_h2:
        runcmd(cmd)
