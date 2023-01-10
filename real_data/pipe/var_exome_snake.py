rule samtools_mpileup_single_cell:
    input:
        bams = getFinalBams,
        fileNames = FINALBAMOUTDIR + '{experiment}_{type}_bamFileNames.txt',
        ref = config['resources'][ORGANISM]['reference'],
        chrRegion = REGIONSOUT + config['resources'][ORGANISM]['regions'].strip().split("/")[-1].replace('.bed','') + '_{chrom}.bed'
    output:
        mpileup = temp(MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup')
    params:
        lsfoutfile = MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.lsferr.log',
        params = config['tools']['samtools']['mpileup']['params'],
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
    conda:
        'envs/samtools.yaml'
    benchmark:
        MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.benchmark'
    shell:
        'samtools mpileup -f {input.ref} {params.params} -b {input.fileNames} -l {input.chrRegion} > {output.mpileup}'

rule samtoolsCombineMpileup:
    input:
        mpileup = expand(MPILEUPOUT + '{{experiment}}_{{type}}-{chrom}.mpileup', chrom=getContigNames())
    output:
        mpileup = temp(MPILEUPOUT + '{experiment}_{type}_complete.mpileup'),
        gz = MPILEUPOUT + '{experiment}_{type}.mpileup.gz'
    params:
        lsfoutfile = MPILEUPOUT + '{experiment}_{type}_complete.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{experiment}_{type}_complete.mpileup.lsferr.log',
        params = config['tools']['samtools']['mpileup']['params'],
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
    conda:
        'envs/samtools.yaml'
    benchmark:
        MPILEUPOUT + '{experiment}_{type}_complete.mpileup.benchmark'
    shell:
        'cat {input.mpileup} > {output.mpileup}; gzip < {output.mpileup} > {output.gz}'

ruleorder: createBamFileSummaryScite > createBamFileSummary
localrules: createBamFileSummaryScite
rule createBamFileSummaryScite:
    input:
        bams = getFinalMpileupBams,
    output:
        OUTDIR + '{location}/{experiment}_all_bamFileNames.txt'
    run:
        sampleMappingFile = open(SAMPLEMAPPING, 'r')
        sampleMapping = {}
        for line in sampleMappingFile:
            sampleMapping[line.strip().split('\t')[1]] = line.strip().split('\t')[2]

        outfile = open(str(output), "w")
        for entry in input.bams:
            sample = entry.split('/')[-1].replace('.bam','')
            outfile.write(entry + '\t' + sampleMapping[sample] + '\n')
        outfile.close()

rule sciphin:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    output:
        tsv = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}_mut2Sample.tsv',
        probs = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.probs',
        gv = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.gv',
        params = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.params.txt',
        vcf = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.vcf',
        bParams = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}/best_index/nuc.tsv',
        bTree = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}/best_index/tree.gv' 
    params:
        lsfoutfile = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.lsfout.log',
        lsferrfile = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.lsferr.log',
        scratch = config['tools']['sciphi']['scratch'],
        mem = config['tools']['sciphi']['mem'],
        time = config['tools']['sciphi']['time'],
        out = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}',
        params = config['tools']['sciphi']['params'],
        outIndex = SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}/index'
    benchmark:
        SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.benchmark'
    threads:
        1
    log:
        SCIPHIOUT + '{run}/{experiment}_llp-{llp}_lpp-{lpp}.log'
    shell:
        ('{config[tools][sciphi][call]} ' +
        '-o {params.out} ' +
        '--ol {params.outIndex} ' + 
        '--in {input.fileNames} ' +
        '--ll 1 ' +
        '--lp 1 ' +
        '--llp {wildcards.llp} ' +
        '--lpp {wildcards.lpp} ' +
        '--lz 1 ' +
        '--seed {wildcards.run} ' +
        '{params.params} ' +
        '--im {input.mpileup}')

ruleorder: sciphin_max > sciphin
rule sciphin_max:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    output:
        tsv = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}_mut2Sample.tsv',
        probs = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.probs',
        gv = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.gv',
        params = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.params.txt',
        vcf = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.vcf',
        bParams = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}/best_index/nuc.tsv',
        bTree = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}/best_index/tree.gv' 
    params:
        lsfoutfile = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.lsfout.log',
        lsferrfile = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.lsferr.log',
        scratch = config['tools']['sciphi']['scratch'],
        mem = config['tools']['sciphi']['mem'],
        time = config['tools']['sciphi']['time'],
        out = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}',
        params = config['tools']['sciphi']['params'],
        outIndex = SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}/index'
    benchmark:
        SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.benchmark'
    threads:
        1
    log:
        SCIPHIOUT + '{run}/{experiment}_max_llp-{llp}_lpp-{lpp}.log'
    shell:
        ('{config[tools][sciphi][call]} ' +
        '-o {params.out} ' +
        '--ol {params.outIndex} ' + 
        '--in {input.fileNames} ' +
        '--ll 1 ' +
        '--lp 1 ' +
        '--llp {wildcards.llp} ' +
        '--lpp {wildcards.lpp} ' +
        '--lz 1 ' +
        '--mlm 1 ' +
        '--seed {wildcards.run} ' +
        '{params.params} ' +
        '--im {input.mpileup}')
rule sciphin_chi:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    output:
        tsv = SCIPHIOUT + '{run}/{experiment}_chi_mut2Sample.tsv',
        probs = SCIPHIOUT + '{run}/{experiment}_chi.probs',
        gv = SCIPHIOUT + '{run}/{experiment}_chi.gv',
        params = SCIPHIOUT + '{run}/{experiment}_chi.params.txt',
        vcf = SCIPHIOUT + '{run}/{experiment}_chi.vcf',
        bParams = SCIPHIOUT + '{run}/{experiment}_chi/best_index/nuc.tsv',
        bTree = SCIPHIOUT + '{run}/{experiment}_chi/best_index/tree.gv' 
    params:
        lsfoutfile = SCIPHIOUT + '{run}/{experiment}_chi.lsfout.log',
        lsferrfile = SCIPHIOUT + '{run}/{experiment}_chi.lsferr.log',
        scratch = config['tools']['sciphi']['scratch'],
        mem = config['tools']['sciphi']['mem'],
        time = config['tools']['sciphi']['time'],
        out = SCIPHIOUT + '{run}/{experiment}_chi',
        params = config['tools']['sciphi']['params'],
        outIndex = SCIPHIOUT + '{run}/{experiment}_chi/index'
    benchmark:
        SCIPHIOUT + '{run}/{experiment}_chi.benchmark'
    threads:
        1
    log:
        SCIPHIOUT + '{run}/{experiment}_chi.log'
    shell:
        ('{config[tools][sciphi][call]} ' +
        '-o {params.out} ' +
        '--ol {params.outIndex} ' + 
        '--in {input.fileNames} ' +
        '--ll 1 ' +
        '--lp 1 ' +
        '--chi 1 ' + 
        '--lz 1 ' +
        '--seed {wildcards.run} ' +
        '{params.params} ' +
        '--im {input.mpileup}')
rule sciphin_chi_max:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    output:
        tsv = SCIPHIOUT + '{run}/{experiment}_chi_max_mut2Sample.tsv',
        probs = SCIPHIOUT + '{run}/{experiment}_chi_max.probs',
        gv = SCIPHIOUT + '{run}/{experiment}_chi_max.gv',
        params = SCIPHIOUT + '{run}/{experiment}_chi_max.params.txt',
        vcf = SCIPHIOUT + '{run}/{experiment}_chi_max.vcf',
        bParams = SCIPHIOUT + '{run}/{experiment}_chi_max/best_index/nuc.tsv',
        bTree = SCIPHIOUT + '{run}/{experiment}_chi_max/best_index/tree.gv' 
    params:
        lsfoutfile = SCIPHIOUT + '{run}/{experiment}_chi_max.lsfout.log',
        lsferrfile = SCIPHIOUT + '{run}/{experiment}_chi_max.lsferr.log',
        scratch = config['tools']['sciphi']['scratch'],
        mem = config['tools']['sciphi']['mem'],
        time = config['tools']['sciphi']['time'],
        out = SCIPHIOUT + '{run}/{experiment}_chi_max',
        params = config['tools']['sciphi']['params'],
        outIndex = SCIPHIOUT + '{run}/{experiment}_chi_max/index'
    benchmark:
        SCIPHIOUT + '{run}/{experiment}_chi_max.benchmark'
    threads:
        1
    log:
        SCIPHIOUT + '{run}/{experiment}_chi_max.log'
    shell:
        ('{config[tools][sciphi][call]} ' +
        '-o {params.out} ' +
        '--ol {params.outIndex} ' + 
        '--in {input.fileNames} ' +
        '--ll 1 ' +
        '--lp 1 ' +
        '--chi 1 ' + 
        '--lz 1 ' +
        '--mlm 1 ' + 
        '--seed {wildcards.run} ' +
        '{params.params} ' +
        '--im {input.mpileup}')

if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = 'placeholder'
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = 'placeholder'
rule gatkHaplotypeCaller:
    input:
        bam = HAPLOTYPECALLERIN + '{sample}.bam',
        bai = HAPLOTYPECALLERIN + '{sample}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = HAPLOTYPECALLEROUT + '{sample}.g.vcf',
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsferr.log',
        scratch = config['tools']['GATK']['haplotypeCaller']['scratch'],
        mem = config['tools']['GATK']['haplotypeCaller']['mem'],
        time = config['tools']['GATK']['haplotypeCaller']['time'],
        params = config['tools']['GATK']['haplotypeCaller']['params'],
    threads:
        config['tools']['GATK']['haplotypeCaller']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.log'
    conda:
        'envs/gatk.yaml'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T HaplotypeCaller ' +
        '{params.params} ' + 
        '-R {input.reference} ' + 
        '-I {input.bam} ' +
        '--emitRefConfidence GVCF ' +
        '-L {input.regions} ' +
        '-mmq 40 ' +
        '-mbq 30 ' + 
        '-o {output.vcf}')

rule gatkGenotypeGVCFs:
    input:
        vcf = HAPLOTYPECALLEROUT + '{sample}.g.vcf',
        reference = config['resources'][ORGANISM]['reference']
    output:
        vcf = HAPLOTYPECALLEROUT + '{sample}.vcf'
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['GATK']['genotypeGVCFs']['scratch'],
        mem = config['tools']['GATK']['genotypeGVCFs']['mem'],
        time = config['tools']['GATK']['genotypeGVCFs']['time'],
        reference = config['resources'][ORGANISM]['reference'],
        params = config['tools']['GATK']['genotypeGVCFs']['params']
    threads:
        config['tools']['GATK']['genotypeGVCFs']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + '{sample}.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + '{sample}.vcf.log'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T GenotypeGVCFs ' +
        '{params.params}' +
        '-R {input.reference} ' +
        '--variant {input.vcf} ' +
        '-o {output.vcf}')
