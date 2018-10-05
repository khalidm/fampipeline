'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages
from utils import safe_make_dir
from collections import Counter, defaultdict
import re

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='fampipeline')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # Make directories for outputs
    safe_make_dir('results')
    safe_make_dir('results/alignments')
    safe_make_dir('results/variants')
    safe_make_dir('results/qc')

    #XXX Note: Temporary solution to treat multi- and single- lane samples differently
    # Check if multiple lanes per sample.
    print(fastq_files)
    sample_info = [re.search('.+/(FAM_[a-zA-Z0-9]+_SM_([a-zA-Z0-9]+))_(ID_[A-Za-z0-9-]+)_.+R([0-9]).fastq.gz',
                             filename) for filename in fastq_files]
    #print([samp.group(2) for samp in sample_info])
    sample_count = Counter([samp.group(1) for samp in sample_info])
    print(sample_count)
    single_lane = [samp for samp, count in sample_count.iteritems() if count == 2]
    print(single_lane)
    multi_lane = [samp for samp, count in sample_count.iteritems() if count > 2 and count % 2 == 0]
    print(multi_lane)
    assert(len(multi_lane + single_lane) == len(sample_count))

    # Define inputs for merge_bams
    merge_bams_input = ['results/alignments/{fam_sm}/{fam_sm}_{id}.sort.dedup.realn' \
                        '.recal.bam'.format(fam_sm=samp.group(1), sm=samp.group(2), id=samp.group(3))
                        for samp in sample_info if samp.group(1) in multi_lane and samp.group(4) == "1"]
    #print(merge_bams_input)

    # Define inputs for call_variants_gatk
    single_lane_processed_bams = ['results/alignments/{fam_sm}/{fam_sm}_{id}.sort.dedup.realn' \
                                  '.recal.bam'.format(fam_sm=samp.group(1), sm=samp.group(2), id=samp.group(3))
                                   for samp in sample_info if samp.group(1) in single_lane and samp.group(4) == "1"]
    multi_lane_processed_bams = ['results/alignments/{fam_sm}/{fam_sm}.merged.dedup.realn' \
                                 '.bam'.format(fam_sm=fam_sm, sm=re.search('FAM_[a-zA-Z0-9]+_SM_([a-zA-Z0-9]+)', fam_sm).group(1))
                                  for fam_sm in multi_lane]
    call_variants_gatk_input = single_lane_processed_bams + multi_lane_processed_bams

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.do_nothing,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        # IF THE READS ARE SPLIT IN LANES e.g. FAM_f2_SM_f2i5_ID_idx46-TCCCGA-L001-L002_LB_lb_PL_ILLUMINA_R2
        # filter=formatter(
            # '.+/FAM_(?P<famid>[a-zA-Z0-9]+)_SM_(?P<sample>[a-zA-Z0-9-]+)_ID_(?P<runid>[a-zA-Z0-9-]+)_(?P<lib>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_R1.fastq.gz'),
        filter=formatter('.+/FAM_(?P<fam>[a-zA-Z0-9]+)_SM_(?P<sample>[a-zA-Z0-9]+)' \
                 '_ID_(?P<id>[a-zA-Z0-9-]+)_LB_(?P<lb>[a-zA-Z0-9]+)' \
                 '_PL_(?P<pl>[a-zA-Z0-9]+)_R[12].fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        add_inputs=add_inputs('{path[0]}/FAM_{fam[0]}_SM_{sample[0]}_ID_{id[0]}' \
                              '_LB_{lb[0]}_PL_{pl[0]}_R2.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{fam[0]}', '{sample[0]}', '{id[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='results/alignments/FAM_{fam[0]}_SM_{sample[0]}/FAM_{fam[0]}_SM_{sample[0]}_ID_{id[0]}.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    # Mark duplicates in the BAM file using Picard
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_picard',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        # XXX should make metricsup an extra output?
        output=['.sort.dedup.bam', '.metricsdup'])

    # Generate chromosome intervals using GATK
    pipeline.transform(
        task_func=stages.chrom_intervals_gatk,
        name='chrom_intervals_gatk',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        output='.chr.intervals')

    # Local realignment using GATK
    (pipeline.transform(
        task_func=stages.local_realignment_gatk,
        name='local_realignment_gatk',
        input=output_from('chrom_intervals_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).chr.intervals'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.bam')
        .follows('mark_duplicates_picard'))

    # Base recalibration using GATK
    pipeline.transform(
        task_func=stages.base_recalibration_gatk,
        name='base_recalibration_gatk',
        input=output_from('local_realignment_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output=['.recal_data.csv', '.count_cov.log'])

    # Print reads using GATK
    (pipeline.transform(
        task_func=stages.print_reads_gatk,
        name='print_reads_gatk',
        input=output_from('base_recalibration_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).recal_data.csv'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.realn.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.recal.bam')
        .follows('local_realignment_gatk'))

    # Define multi_lane_bams
    (pipeline.originate(
        task_func=stages.do_nothing,
        name='define_multi_lane_bams',
        output=merge_bams_input)
        .follows('print_reads_gatk'))

    # Merge multi-lane bams
    pipeline.collate(
        task_func=stages.merge_bams,
        name='merge_bams',
        input=output_from('define_multi_lane_bams'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_ID_.+.sort.dedup.realn.recal.bam'),
        output='results/alignments/{sample[0]}/{sample[0]}.merged.bam')

    # Mark duplicates 2 on multi-lane samples using Picard
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_2',
        input=output_from('merge_bams'),
        filter=suffix('.merged.bam'),
        output=['.merged.dedup.bam', '.metricsdup'])

    # Generate chromosome intervals using GATK (2)
    pipeline.transform(
        task_func=stages.chrom_intervals_gatk,
        name='chrom_intervals_gatk_2',
        input=output_from('mark_duplicates_2'),
        filter=suffix('.merged.dedup.bam'),
        output='.chr.intervals')

    # Local realignment 2 using GATK
    (pipeline.transform(
        task_func=stages.local_realignment_gatk,
        name='local_realignment_2',
        input=output_from('chrom_intervals_gatk_2'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).chr.intervals'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.merged.dedup.bam'),
        output='{path[0]}/{sample[0]}.merged.dedup.realn.bam')
        .follows('mark_duplicates_2'))

    # Define processed bams
    (pipeline.originate(
         task_func=stages.do_nothing,
        name='define_processed_bams',
        output=call_variants_gatk_input)
        .follows('local_realignment_2')
        .follows('print_reads_gatk'))

    # Call variants using GATK
    (pipeline.transform(
        task_func=stages.call_variants_gatk,
        name='call_variants_gatk',
        input=output_from('define_processed_bams'),
        filter=formatter('.+/(?P<sample>FAM_[a-zA-Z0-9]+_SM_[a-zA-Z0-9]+)[_.].+.dedup.realn(.recal)?.bam'),
        output='results/variants/{sample[0]}.raw.snps.indels.g.vcf'))

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('call_variants_gatk'),
        output='results/variants/FAMExomes.mergegvcf.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.mergegvcf.vcf'),
        output='.genotyped.vcf')

    # SNP recalibration using GATK
    pipeline.transform(
        task_func=stages.snp_recalibrate_gatk,
        name='snp_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        output=['.snp_recal', '.snp_tranches', '.snp_plots.R'])

    # INDEL recalibration using GATK
    pipeline.transform(
        task_func=stages.indel_recalibrate_gatk,
        name='indel_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        output=['.indel_recal', '.indel_tranches', '.indel_plots.R'])

    # Apply SNP recalibration using GATK
    (pipeline.transform(
        task_func=stages.apply_snp_recalibrate_gatk,
        name='apply_snp_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        add_inputs=add_inputs(['results/variants/FAMExomes.snp_recal',
            'results/variants/FAMExomes.snp_tranches']),
        output='.recal_SNP.vcf')
        .follows('snp_recalibrate_gatk'))

    # Apply INDEL recalibration using GATK
    (pipeline.transform(
        task_func=stages.apply_indel_recalibrate_gatk,
        name='apply_indel_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        add_inputs=add_inputs(['results/variants/FAMExomes.indel_recal',
            'results/variants/FAMExomes.indel_tranches']),
        output='.recal_INDEL.vcf')
        .follows('indel_recalibrate_gatk'))

    # Combine variants using GATK
    (pipeline.transform(
        task_func=stages.combine_variants_gatk,
        name='combine_variants_gatk',
        input=output_from('apply_snp_recalibrate_gatk'),
        filter=suffix('.recal_SNP.vcf'),
        add_inputs=add_inputs(['results/variants/FAMExomes.recal_INDEL.vcf']),
        output='.combined.vcf')
        .follows('apply_indel_recalibrate_gatk'))

    # Filter variants using GATK
    pipeline.transform(
        task_func=stages.filter_variants_gatk,
        name='filter_variants_gatk',
        input=output_from('combine_variants_gatk'),
        filter=suffix('.combined.vcf'),
        output='.filtered.vcf')

    # Select variants using GATK
    pipeline.transform(
        task_func=stages.select_variants_gatk,
        name='select_variants_gatk',
        input=output_from('filter_variants_gatk'),
        filter=suffix('.filtered.vcf'),
        output='.selected.vcf')

    # Select variants using GATK (only fam samples)
    pipeline.transform(
        task_func=stages.select_fam_variants_gatk,
        name='select_fam_variants_gatk',
        input=output_from('filter_variants_gatk'),
        filter=suffix('.filtered.vcf'),
        output='.selected.fam.vcf')

    # Rare variant genotyping using FamSeq
    pipeline.transform(
        task_func=stages.rare_variants_famseq,
        name='rare_variants_famseq',
        input=output_from('select_variants_gatk'),
        filter=suffix('.selected.vcf'),
        output='.famseq.vcf')

    # Get coverage for each exon with BedTools
    pipeline.transform(
        task_func=stages.coverage_bedtools,
        name='coverage_bedtools',
        input=output_from('define_processed_bams'),
        filter=formatter('.+/(?P<sample>FAM_[a-zA-Z0-9]+_SM_[a-zA-Z0-9]+)[_.].+.dedup.realn(.recal)?.bam'),
        output='results/qc/{sample[0]}.exon_coverages.txt')

    # Get alignment statistics with BamTools
    pipeline.transform(
        task_func=stages.alignment_stats_bamtools,
        name='alignment_stats_bamtools',
        input=output_from('define_processed_bams'),
        filter=formatter('.+/(?P<sample>FAM_[a-zA-Z0-9]+_SM_[a-zA-Z0-9]+)[_.].+.dedup.realn(.recal)?.bam'),
        output='results/qc/{sample[0]}.alignment_stats.txt')

    return pipeline
