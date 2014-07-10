/*
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package com.custom;

import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.MultipleTextOutputFormat;
import org.apache.hadoop.util.*;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

/**
 * CustomMultiOutputFormat implements the VCF header for VA data.
 *
 * The code that follows is a mashup of
 * http://stackoverflow.com/questions/18541503/multiple-output-files-for-hadoop-streaming-with-python-mapper
 * and a modified version of
 * org.apache.hadoop.mapred.TextOutputFormat.
 **/

public class CustomMultiOutputFormat<K, V> extends MultipleTextOutputFormat<K, V> {

  static final String DATASET_HEADER =
      "##fileformat=VCFv4.1\n"
      + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
      + "##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description=\"Truth sensitivity tranche level for INDEL model at VQS Lod < -861.7921\">\n"
      + "##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description=\"Truth sensitivity tranche level for INDEL model at VQS Lod: -861.7921 <= x < -4.8291\">\n"
      + "##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description=\"Truth sensitivity tranche level for SNP model at VQS Lod < -39736.6503\">\n"
      + "##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description=\"Truth sensitivity tranche level for SNP model at VQS Lod: -39736.6503 <= x < 2.4578\">\n"
      + "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
      + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n"
      + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
      + "##GATKCommandLine=<ID=ApplyRecalibration,Version=2.7-4-g6f46d11,Date=\"Thu Feb 27 12:04:14 PST 2014\",Epoch=1393531454500,CommandLineOptions=\"analysis_type=ApplyRecalibration input_file=[] read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=cache/store/67224e7c805019dd76102dbcbf22fb744a291a43/ucsc.hg19.fa nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false input=[(RodBinding name=input source=f62ca808-a766-4bc8-b872-71a7f19dbac3/0.INDEL.vcf)] recal_file=(RodBinding name=recal_file source=f62ca808-a766-4bc8-b872-71a7f19dbac3/INDEL.recal.vcf) tranches_file=f62ca808-a766-4bc8-b872-71a7f19dbac3/INDEL.tranches out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub ts_filter_level=99.9 ignore_filter=null excludeFiltered=false mode=INDEL filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">\n"
      + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n"
      + "##GATKCommandLine=<ID=ApplyRecalibration,Version=2.7-4-g6f46d11,Date=\"Thu Feb 27 12:17:07 PST 2014\",Epoch=1393532227013,CommandLineOptions=\"analysis_type=ApplyRecalibration input_file=[] read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=cache/store/67224e7c805019dd76102dbcbf22fb744a291a43/ucsc.hg19.fa nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false input=[(RodBinding name=input source=f62ca808-a766-4bc8-b872-71a7f19dbac3/0.SNP.vcf)] recal_file=(RodBinding name=recal_file source=f62ca808-a766-4bc8-b872-71a7f19dbac3/SNP.recal.vcf) tranches_file=f62ca808-a766-4bc8-b872-71a7f19dbac3/SNP.tranches out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub ts_filter_level=99.9 ignore_filter=null excludeFiltered=false mode=SNP filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">\n"
      + "##GATKCommandLine=<ID=UnifiedGenotyper,Version=2.7-4-g6f46d11,Date=\"Thu Feb 27 09:11:55 PST 2014\",Epoch=1393521115436,CommandLineOptions=\"analysis_type=UnifiedGenotyper input_file=[f62ca808-a766-4bc8-b872-71a7f19dbac3/0.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/1.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/2.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/3.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/4.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/5.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/6.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/7.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/8.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/9.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/10.recalibrated.bam, f62ca808-a766-4bc8-b872-71a7f19dbac3/11.recalibrated.bam] read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=[f62ca808-a766-4bc8-b872-71a7f19dbac3/0.bed] excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=cache/store/67224e7c805019dd76102dbcbf22fb744a291a43/ucsc.hg19.fa nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=250 baq=CALCULATE_AS_NECESSARY baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false genotype_likelihoods_model=BOTH pcr_error_rate=1.0E-4 computeSLOD=false annotateNDA=false pair_hmm_implementation=LOGLESS_CACHING min_base_quality_score=17 max_deletion_fraction=0.05 allSitePLs=false min_indel_count_for_genotyping=5 min_indel_fraction_per_sample=0.25 indelGapContinuationPenalty=10 indelGapOpenPenalty=45 indelHaplotypeSize=80 indelDebug=false ignoreSNPAlleles=false allReadsSP=false ignoreLaneInfo=false reference_sample_calls=(RodBinding name= source=UNBOUND) reference_sample_name=null sample_ploidy=2 min_quality_score=1 max_quality_score=40 site_quality_prior=20 min_power_threshold_for_calling=0.95 min_reference_depth=100 exclude_filtered_reference_sites=false output_mode=EMIT_ALL_CONFIDENT_SITES heterozygosity=0.001 indel_heterozygosity=1.25E-4 genotyping_mode=DISCOVERY standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=0.1 alleles=(RodBinding name= source=UNBOUND) max_alternate_alleles=6 input_prior=[] contamination_fraction_to_filter=0.0 contamination_fraction_per_sample_file=null p_nonref_model=EXACT_INDEPENDENT exactcallslog=null dbsnp=(RodBinding name=dbsnp source=cache/store/9662b0afcefa5e2b90730ae22a58bc6c6abd4f98/dbsnp_135.hg19.vcf) comp=[] out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub onlyEmitSamples=[] debug_file=null metrics_file=null annotation=[AlleleBalance, Coverage, MappingQualityZero] excludeAnnotation=[] filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">\n"
      + "##INFO=<ID=ABHet,Number=1,Type=Float,Description=\"Allele Balance for hets (ref/(ref+alt))\">\n"
      + "##INFO=<ID=ABHom,Number=1,Type=Float,Description=\"Allele Balance for homs (A/(A+O))\">\n"
      + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n"
      + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n"
      + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
      + "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
      + "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">\n"
      + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n"
      + "##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">\n"
      + "##INFO=<ID=Dels,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">\n"
      + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n"
      + "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
      + "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">\n"
      + "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">\n"
      + "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">\n"
      + "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">\n"
      + "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n"
      + "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n"
      + "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
      + "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the negative training set of bad variants\">\n"
      + "##INFO=<ID=OND,Number=1,Type=Float,Description=\"Overall non-diploid ratio (alleles/(alleles+non-alleles))\">\n"
      + "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the positive training set of good variants\">\n"
      + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">\n"
      + "##INFO=<ID=RPA,Number=.,Type=Integer,Description=\"Number of times tandem repeat unit is repeated, for each allele (including reference)\">\n"
      + "##INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">\n"
      + "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n"
      + "##INFO=<ID=STR,Number=0,Type=Flag,Description=\"Variant is a short tandem repeat\">\n"
      + "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"Log odds ratio of being a true variant versus being false under the trained gaussian mixture model\">\n"
      + "##INFO=<ID=culprit,Number=1,Type=String,Description=\"The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out\">\n"
      + "##contig=<ID=chrM,length=16571,assembly=hg19>\n"
      + "##contig=<ID=chr1,length=249250621,assembly=hg19>\n"
      + "##contig=<ID=chr2,length=243199373,assembly=hg19>\n"
      + "##contig=<ID=chr3,length=198022430,assembly=hg19>\n"
      + "##contig=<ID=chr4,length=191154276,assembly=hg19>\n"
      + "##contig=<ID=chr5,length=180915260,assembly=hg19>\n"
      + "##contig=<ID=chr6,length=171115067,assembly=hg19>\n"
      + "##contig=<ID=chr7,length=159138663,assembly=hg19>\n"
      + "##contig=<ID=chr8,length=146364022,assembly=hg19>\n"
      + "##contig=<ID=chr9,length=141213431,assembly=hg19>\n"
      + "##contig=<ID=chr10,length=135534747,assembly=hg19>\n"
      + "##contig=<ID=chr11,length=135006516,assembly=hg19>\n"
      + "##contig=<ID=chr12,length=133851895,assembly=hg19>\n"
      + "##contig=<ID=chr13,length=115169878,assembly=hg19>\n"
      + "##contig=<ID=chr14,length=107349540,assembly=hg19>\n"
      + "##contig=<ID=chr15,length=102531392,assembly=hg19>\n"
      + "##contig=<ID=chr16,length=90354753,assembly=hg19>\n"
      + "##contig=<ID=chr17,length=81195210,assembly=hg19>\n"
      + "##contig=<ID=chr18,length=78077248,assembly=hg19>\n"
      + "##contig=<ID=chr19,length=59128983,assembly=hg19>\n"
      + "##contig=<ID=chr20,length=63025520,assembly=hg19>\n"
      + "##contig=<ID=chr21,length=48129895,assembly=hg19>\n"
      + "##contig=<ID=chr22,length=51304566,assembly=hg19>\n"
      + "##contig=<ID=chrX,length=155270560,assembly=hg19>\n"
      + "##contig=<ID=chrY,length=59373566,assembly=hg19>\n"
      + "##contig=<ID=chr1_gl000191_random,length=106433,assembly=hg19>\n"
      + "##contig=<ID=chr1_gl000192_random,length=547496,assembly=hg19>\n"
      + "##contig=<ID=chr4_ctg9_hap1,length=590426,assembly=hg19>\n"
      + "##contig=<ID=chr4_gl000193_random,length=189789,assembly=hg19>\n"
      + "##contig=<ID=chr4_gl000194_random,length=191469,assembly=hg19>\n"
      + "##contig=<ID=chr6_apd_hap1,length=4622290,assembly=hg19>\n"
      + "##contig=<ID=chr6_cox_hap2,length=4795371,assembly=hg19>\n"
      + "##contig=<ID=chr6_dbb_hap3,length=4610396,assembly=hg19>\n"
      + "##contig=<ID=chr6_mann_hap4,length=4683263,assembly=hg19>\n"
      + "##contig=<ID=chr6_mcf_hap5,length=4833398,assembly=hg19>\n"
      + "##contig=<ID=chr6_qbl_hap6,length=4611984,assembly=hg19>\n"
      + "##contig=<ID=chr6_ssto_hap7,length=4928567,assembly=hg19>\n"
      + "##contig=<ID=chr7_gl000195_random,length=182896,assembly=hg19>\n"
      + "##contig=<ID=chr8_gl000196_random,length=38914,assembly=hg19>\n"
      + "##contig=<ID=chr8_gl000197_random,length=37175,assembly=hg19>\n"
      + "##contig=<ID=chr9_gl000198_random,length=90085,assembly=hg19>\n"
      + "##contig=<ID=chr9_gl000199_random,length=169874,assembly=hg19>\n"
      + "##contig=<ID=chr9_gl000200_random,length=187035,assembly=hg19>\n"
      + "##contig=<ID=chr9_gl000201_random,length=36148,assembly=hg19>\n"
      + "##contig=<ID=chr11_gl000202_random,length=40103,assembly=hg19>\n"
      + "##contig=<ID=chr17_ctg5_hap1,length=1680828,assembly=hg19>\n"
      + "##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>\n"
      + "##contig=<ID=chr17_gl000204_random,length=81310,assembly=hg19>\n"
      + "##contig=<ID=chr17_gl000205_random,length=174588,assembly=hg19>\n"
      + "##contig=<ID=chr17_gl000206_random,length=41001,assembly=hg19>\n"
      + "##contig=<ID=chr18_gl000207_random,length=4262,assembly=hg19>\n"
      + "##contig=<ID=chr19_gl000208_random,length=92689,assembly=hg19>\n"
      + "##contig=<ID=chr19_gl000209_random,length=159169,assembly=hg19>\n"
      + "##contig=<ID=chr21_gl000210_random,length=27682,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000211,length=166566,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000212,length=186858,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000213,length=164239,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000214,length=137718,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000215,length=172545,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000216,length=172294,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000217,length=172149,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000218,length=161147,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000219,length=179198,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000220,length=161802,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000221,length=155397,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000222,length=186861,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000223,length=180455,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000224,length=179693,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000225,length=211173,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000226,length=15008,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000227,length=128374,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000228,length=129120,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000229,length=19913,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000230,length=43691,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000231,length=27386,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000232,length=40652,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000233,length=45941,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000234,length=40531,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000235,length=34474,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000236,length=41934,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000237,length=45867,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000238,length=39939,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000239,length=33824,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000240,length=41933,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000241,length=42152,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000242,length=43523,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000243,length=43341,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000244,length=39929,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000245,length=36651,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000246,length=38154,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000247,length=36422,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000248,length=39786,assembly=hg19>\n"
      + "##contig=<ID=chrUn_gl000249,length=38502,assembly=hg19>\n"
      + "##reference=file:///mnt/scratch/local/cache/store/67224e7c805019dd76102dbcbf22fb744a291a43/ucsc.hg19.fa\n";

  /**
   * Use they key as part of the path for the final output file.
   */
  @Override
      protected String generateFileNameForKeyValue(K key, V value, String leaf) {
    return new Path(key.toString(), leaf).toString();
  }

  protected static class LineRecordWriter<K, V>
      implements RecordWriter<K, V> {
    private static final String utf8 = "UTF-8";
    private static final byte[] newline;
    static {
      try {
        newline = "\n".getBytes(utf8);
      } catch (UnsupportedEncodingException uee) {
        throw new IllegalArgumentException("can't find " + utf8 + " encoding");
      }
    }

    protected DataOutputStream out;
    private final byte[] keyValueSeparator;

    public LineRecordWriter(DataOutputStream out, String keyValueSeparator) {
      this.out = out;
      try {
        this.keyValueSeparator = keyValueSeparator.getBytes(utf8);
      } catch (UnsupportedEncodingException uee) {
        throw new IllegalArgumentException("can't find " + utf8 + " encoding");
      }
    }

    public LineRecordWriter(DataOutputStream out) {
      this(out, "\t");
    }

    /**
     * Write the object to the byte stream, handling Text as a special
     * case.
     * @param o the object to print
     * @throws IOException if the write throws, we pass it on
     */
    private void writeObject(Object o) throws IOException {
      if (o instanceof Text) {
        Text to = (Text) o;
        out.write(to.getBytes(), 0, to.getLength());
      } else if (o instanceof BytesWritable) {
        BytesWritable bytes = (BytesWritable) o;
        out.write(bytes.getBytes(), 0, bytes.getLength());
      } else {
        out.write(o.toString().getBytes(utf8));
      }
    }

    public synchronized void write(K key, V value)
                                 throws IOException {

        boolean nullKey = key == null || key instanceof NullWritable;
        boolean nullValue = value == null || value instanceof NullWritable;
        if (nullKey || nullValue) {
          return;
        }
        if (0 == out.size()) {
          // This record will be written to a newly opened file.
          // Write the VCF header first.
          String header = DATASET_HEADER
              + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
              + "INFO\tFORMAT\t"
              + key.toString()
              + "\n";
          out.write(header.getBytes(), 0,
              header.getBytes().length);
        }
        writeObject(value);
        out.write(newline);
      }

    public synchronized void close(Reporter reporter) throws IOException {
        out.close();
      }
  }

  public RecordWriter<K, V> getBaseRecordWriter(FileSystem ignored,
      JobConf job,
      String name,
      Progressable progress)
      throws IOException {
    boolean isCompressed = getCompressOutput(job);
    String keyValueSeparator = job.get("mapred.textoutputformat.separator",
        "\t");
    if (!isCompressed) {
      Path file = FileOutputFormat.getTaskOutputPath(job, name);
      FileSystem fs = file.getFileSystem(job);
      FSDataOutputStream fileOut = fs.create(file, progress);
      return new LineRecordWriter<K, V>(fileOut, keyValueSeparator);
    } else {
      Class<? extends CompressionCodec> codecClass =
          getOutputCompressorClass(job, GzipCodec.class);
      // create the named codec
      CompressionCodec codec = ReflectionUtils.newInstance(codecClass, job);
      // build the filename including the extension
      Path file =
          FileOutputFormat.getTaskOutputPath(job,
              name + codec.getDefaultExtension());
      FileSystem fs = file.getFileSystem(job);
      FSDataOutputStream fileOut = fs.create(file, progress);
      return new LineRecordWriter<K, V>(new DataOutputStream
          (codec.createOutputStream(fileOut)),
          keyValueSeparator);
    }
  }
}
