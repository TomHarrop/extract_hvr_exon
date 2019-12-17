#!/usr/bin/env python3

from pathlib import Path
import pysam
import multiprocessing

# globals

hvr_exon = 'data/hvr_exon.txt'
ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
map_data = 'data/merged.bam'
csd_region = 'data/csd_region_1kb.txt'


outdir = 'output'
logdir = Path(outdir, 'logs')


biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
clustalo = 'shub://TomHarrop/singularity-containers:clustalo_1.2.4'
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'       # fixme
transindel = 'shub://TomHarrop/variant-utils:transindel_7098bd6'
freebayes = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
muscle = 'shub://TomHarrop/align-utils:muscle_3.8.1551'
r = 'shub://TomHarrop/r-containers:r_3.6.1'

# main

# parse individuals from the BAM
h = pysam.AlignmentFile(map_data, 'rb').header
all_indivs = sorted(set(x['ID'] for x in h.to_dict()['RG']))

# parse individuals from the VCF
# all_indivs = sorted(set(VCF(calls).samples))

# generate the reference allele
all_indivs.append('ref')

rule target:
    input:
        Path(outdir, 'align_consensus', 'consensus.faa'),
        Path(outdir, 'align_consensus', 'consensus.dist'),
        Path(outdir, 'muscle_consensus', 'tree.pdf')


rule plot_muscle_output:
    input:
        tree = Path(outdir, 'muscle_consensus', 'consensus.tree')
    output:
        plot = Path(outdir, 'muscle_consensus', 'tree.pdf')
    params:
        alleles = 5
    log:
        Path(logdir, 'plot_muscle_output.log')
    singularity:
        r
    script:
        'src/plot_muscle_output.R'


rule muscle_consensus:
    input:
        Path(outdir, 'translate_consensus', 'consensus.fa')
    output:
        aln = Path(outdir, 'muscle_consensus', 'consensus.faa'),
        tree = Path(outdir, 'muscle_consensus', 'consensus.tree')
    log:
        Path(logdir, 'muscle_consensus.log')
    singularity:
        muscle
    shell:
        'muscle '
        '-in {input} '
        '-out {output.aln} '
        '-tree2 {output.tree} '
        '&> {log}'


rule align_consensus:
    input:
        Path(outdir, 'translate_consensus', 'consensus.fa')
    output:
        aln = Path(outdir, 'align_consensus', 'consensus.faa'),
        dist = Path(outdir, 'align_consensus', 'consensus.dist')
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
        Path(outdir, 'combine_cds', 'consensus.fa')
    output:
        Path(outdir, 'translate_consensus', 'consensus.fa')
    singularity:
        biopython
    script:
        'src/translate_consensus.py'

rule combine_cds:
    input:
        expand(Path(outdir, 'condense_cds', '{indiv}.fa').as_posix(),
               indiv=all_indivs)
    output:
        Path(outdir, 'combine_cds', 'consensus.fa')
    singularity:
        samtools
    shell:
        'cat {input} > {output}'

rule condense_cds:
    input:
        Path(outdir, 'extract_derived_cds', '{indiv}.fa')
    output:
        Path(outdir, 'condense_cds', '{indiv}.fa')
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
        vcf = Path(outdir, 'extract_csd_region', 'freebayes.vcf.gz')
    output:
        Path(outdir, 'extract_derived_cds', '{indiv}.fa')
    log:
        Path(logdir, 'extract_derived_cds.{indiv}.log')
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

rule transindel_call_freebayes:
    input:
        bam = Path(outdir, 'transindel', 'csd.bam'),
        fa = ref,
        regions = csd_region
    output:
        Path(outdir, 'extract_csd_region', 'freebayes.vcf')
    params:
        ploidy = 1
    log:
        Path(logdir, 'transindel_call_freebayes.log')
    singularity:
        freebayes
    shell:
        'freebayes '
        '--region  "$(cat {input.regions})" '
        '--ploidy {params.ploidy} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'

rule transindel_build:
    input:
        bam = Path(outdir, 'extract_csd_region/csd.bam'),
        bai = Path(outdir, 'extract_csd_region/csd.bam.bai')
    output:
        Path(outdir, 'transindel/csd.bam')
    params:
        max_del = int(1e3)
    log:
        Path(logdir, 'transindel_build.log')
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

# generic index rule
rule index_vcf:
    input:
        Path(outdir, '{folder}', '{file}.vcf')
    output:
        gz = Path(outdir, '{folder}', '{file}.vcf.gz'),
        tbi = Path(outdir, '{folder}', '{file}.vcf.gz.tbi')
    log:
        Path(logdir, 'index_vcf.{folder}_{file}.log')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'
