#!/usr/bin/env python3

from pathlib import Path
from cyvcf2 import VCF
import multiprocessing

# globals

calls = 'data/calls.vcf.gz'
hvr_exon = 'data/hvr_exon.txt'
ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

outdir = 'output'
logdir = Path(outdir, 'logs')

biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
clustalo = 'shub://TomHarrop/singularity-containers:clustalo_1.2.4'
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'       # fixme
transindel = 'shub://TomHarrop/variant-utils:transindel_7098bd6'

# main

# parse individuals from the VCF
all_indivs = sorted(set(VCF(calls).samples))
all_indivs.append('ref')

rule target:
    input:
        Path(outdir, 'align_consensus/consensus.faa'),
        Path(outdir, 'align_consensus/consensus.dist')


rule align_consensus:
    input:
        Path(outdir, 'translate_consensus/consensus.fa')
    output:
        aln = Path(outdir, 'align_consensus/consensus.faa'),
        dist = Path(outdir, 'align_consensus/consensus.dist')
    log:
        Path(logdir, 'align_consensus.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        clustalo
    shell:
        'clustalo '
        '-i {input} '
        '--threads {threads} '
        '--dealign '
        '--full '
        '--out {output.aln} '
        '--distmat-out {output.dist} '
        '&> {log}'


rule translate_consensus:
    input:
        Path(outdir, 'combine_cds/consensus.fa')
    output:
        Path(outdir, 'translate_consensus/consensus.fa')
    singularity:
        biopython
    script:
        'src/translate_consensus.py'

rule combine_cds:
    input:
        expand(Path(outdir, 'condense_cds/{indiv}.fa').as_posix(),
               indiv=all_indivs)
    output:
        Path(outdir, 'combine_cds/consensus.fa')
    singularity:
        samtools
    shell:
        'cat {input} > {output}'

rule condense_cds:
    input:
        Path(outdir, 'extract_derived_cds/{indiv}.fa')
    output:
        Path(outdir, 'condense_cds/{indiv}.fa')
    params:
        header = '>{indiv}'
    singularity:
        samtools
    shell:
        'echo "{params.header}" > {output} ; '
        'grep -v "^>" {input} >> {output}'


rule extract_derived_cds:
    input:
        fa = ref,
        regions = hvr_exon,
        vcf = calls
    output:
        Path(outdir, 'extract_derived_cds/{indiv}.fa')
    log:
        Path(logdir, 'extract_dmax_del_lengtherived_cds.{indiv}.log')
    # wildcard_constraints:
    #     indiv = '(?!ref)'
    singularity:
        samtools
    shell:
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '2> {log} '
        '| '
        'bcftools consensus '
        '-i\'QUAL>30\' '       # HARDCODED FILTERING!
        '-s {wildcards.indiv} '
        '-H 1 '
        '{input.vcf} '
        '> {output} '
        '2>> {log}'


rule ref_cds:
    input:
        fa = ref,
        regions = hvr_exon,
    output:
        Path(outdir, 'extract_derived_cds/ref.fa')
    singularity:
        samtools
    shell:
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '> {output}'


# experimental, extract HVR from BB18 and run transindel on hvr chr
map_data = 'data/merged.bam'
csd_region = 'data/csd_region_1kb.txt'

rule transindel_target:
    input:
        bam = Path(outdir, 'extract_csd_region/csd.bam')

rule transindel_call:
    input:
        Path(outdir, 'transindel/{indiv}.bam')
    output:
        Path(outdir, 'transindel/{indiv}.indel.vcf')
    params:
        prefix = lambda wildcards:
            Path(outdir, f'transindel/{wildcards.indiv}')
    log:
        Path(logdir, 'transindel_call.{indiv}.log')
    singularity:
        transindel
    shell:
        'transIndel_call.py  '
        '-i {input} '
        '-o {params.prefix} '
        '&> {log}'

rule transindel_build:
    input:
        bam = Path(outdir, 'extract_hvr_chr/{indiv}.bam'),
        bai = Path(outdir, 'extract_hvr_chr/{indiv}.bam.bai')
    output:
        Path(outdir, 'transindel/{indiv}.bam')
    params:
        max_del = int(1e3)
    log:
        Path(logdir, 'transindel_build.{indiv}.log')
    singularity:
        transindel
    shell:
        'transIndel_build_DNA.py '
        '-i {input.bam} '
        '-o {output} '
        '--max_del_length {params.max_del} '
        '&> {log}'

rule extract_csd_region:
    input:
        bam = Path(map_data),
        regions = csd_region
    output:
        bam = Path(outdir, 'extract_csd_region/csd.bam'),
        bai = Path(outdir, 'extract_csd_region/csd.bam.bai')
    singularity:
        samtools
    shell:
        'samtools view '
        '-O BAM '
        '{input.bam} '
        '$(cat {input.regions}) '
        '> {output.bam} '
        '; '
        'samtools index {output.bam}'



