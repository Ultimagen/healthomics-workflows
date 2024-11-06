# SVPipeline
Runs Structural variant pipeline
This pipeline supports germline and somatic modes
The input of that pipeline is cram files and the output is vcf file
The steps of the pipeline are as following:
-Create an assembly file out of the cram files
-Run UA alingnment on that
-Fix the UA alignment which are secondarily mapped to decoy or with low mapq
-Run gridss.IdentifyVariants and gridss.AnnotateVariants
-Run R script / GRIPSS for filtering and linkage the variants

<b>When Running in AWS HealthOmics this pipeline should run with [dynamic storage](https://docs.omics.ai/products/workbench/engines/parameters/aws-healthomics#storage_type-dynamic-or-static)</b>

## Inputs

### Required inputs
<p name="SVPipeline.base_file_name">
        <b>SVPipeline.base_file_name</b><br />
        <i>String </i> &mdash; 
         Base file name for the output files (to be used as the prefix) <br /> 
</p>
<p name="SVPipeline.references">
        <b>SVPipeline.references</b><br />
        <i>References </i> &mdash; 
         Reference files: fasta, dict and fai, recommended value set in the template <br /> 
</p>
<p name="SVPipeline.ua_references">
        <b>SVPipeline.ua_references</b><br />
        <i>UaReferences </i> &mdash; 
         UAReference files: ua_index, ref_alt, v_aware_alignment_flag and ua_extra_args, recommended value set in the template <br /> 
</p>
<p name="SVPipeline.wgs_calling_interval_list">
        <b>SVPipeline.wgs_calling_interval_list</b><br />
        <i>File </i> &mdash; 
         interval list defining the region to perform variant calling on, recommended value set in the template <br /> 
</p>
<p name="SVPipeline.min_base">
        <b>SVPipeline.min_base</b><br />
        <i>Int </i> &mdash; 
         Assembly parameter: Minimum base quality for using in DeBruijn graph construction. Default value in template <br /> 
</p>
<p name="SVPipeline.min_mapq">
        <b>SVPipeline.min_mapq</b><br />
        <i>Int </i> &mdash; 
         Assembly parameter: Minimum mapping quality. Default value in template <br /> 
</p>
<p name="SVPipeline.max_num_haps">
        <b>SVPipeline.max_num_haps</b><br />
        <i>Int? </i> &mdash; 
         Assembly parameter: Maximum number of haplotypes showing an evidence of SV to report <br /> 
</p>
<p name="SVPipeline.realign_mapq">
        <b>SVPipeline.realign_mapq</b><br />
        <i>Int </i> &mdash; 
         Realignment parameter: Below this value we skip realignment on the supplementary alignment <br /> 
</p>
<p name="SVPipeline.homopolymer_length">
        <b>SVPipeline.homopolymer_length</b><br />
        <i>Int </i> &mdash; 
         Realignment parameter: do realignment on homopolymeres longer than this value <br /> 
</p>
<p name="SVPipeline.is_somatic">
        <b>SVPipeline.is_somatic</b><br />
        <i>Boolean </i> &mdash; 
         run in somatic mode or in germline mode <br /> 
</p>
<p name="SVPipeline.reference_name">
        <b>SVPipeline.reference_name</b><br />
        <i>String </i> &mdash; 
         Can be 38 or 19 <br /> 
</p>
<p name="SVPipeline.run_ua">
        <b>SVPipeline.run_ua</b><br />
        <i>Boolean </i> &mdash; 
         Whether to run UA realignment on the output of the assembly (helps resolving some deletions) or not <br /> 
</p>
<p name="SVPipeline.num_shards">
        <b>SVPipeline.num_shards</b><br />
        <i>Int </i> &mdash; 
         Relevant for scatter tasks, which are CreateAssembly and gridss.AnnotateVariants <br /> 
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="SVPipeline.input_germline_crams">
        <b>SVPipeline.input_germline_crams</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM file for the germline or matched normal sample; optinal for supporting somatic calling tumor only, default []
</p>
<p name="SVPipeline.input_germline_crams_indexes">
        <b>SVPipeline.input_germline_crams_indexes</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM index for the germline or matched normal sample; optinal for supporting somatic calling tumor only
</p>
<p name="SVPipeline.input_tumor_crams">
        <b>SVPipeline.input_tumor_crams</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM file for the tumor (in case of matched T/N calling)
</p>
<p name="SVPipeline.input_tumor_crams_indexes">
        <b>SVPipeline.input_tumor_crams_indexes</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM index for the tumor (in case of matched T/N calling)
</p>
<p name="SVPipeline.min_indel_sc_size_to_include">
        <b>SVPipeline.min_indel_sc_size_to_include</b><br />
        <i>String? &mdash; Default: None</i><br />
        Assembly parameter: Minimum size of an indel and soft-clipping in the read to include the read in the assembly.
</p>
<p name="SVPipeline.blacklist_bed">
        <b>SVPipeline.blacklist_bed</b><br />
        <i>File? &mdash; Default: None</i><br />
        Gridss blacklist file
</p>
<p name="SVPipeline.prefilter_query">
        <b>SVPipeline.prefilter_query</b><br />
        <i>String? &mdash; Default: None</i><br />
        Expression (in bcftools view format) to filter the variants before annotation
</p>
<p name="SVPipeline.gridss_metrics_interval">
        <b>SVPipeline.gridss_metrics_interval</b><br />
        <i>String? &mdash; Default: None</i><br />
        Interval for collecting gridss metrics
</p>
<p name="SVPipeline.pon_sgl_file">
        <b>SVPipeline.pon_sgl_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Panel of normals for single end breakend (partially resolved) calls. Note that the default value is in template
</p>
<p name="SVPipeline.pon_sv_file">
        <b>SVPipeline.pon_sv_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: panel of normals for breakpoint (fully resolved) calls. Note that the default value is in template
</p>
<p name="SVPipeline.repeat_mask_file">
        <b>SVPipeline.repeat_mask_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Repeat mask file. Note that the default value is in template
</p>
<p name="SVPipeline.known_hotspot_file">
        <b>SVPipeline.known_hotspot_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Known locations that are hot spot for SVs (see https://github.com/hartwigmedical/hmftools/tree/master/linx), filtered less stringently
</p>
<p name="SVPipeline.min_normal_coverage">
        <b>SVPipeline.min_normal_coverage</b><br />
        <i>Int? &mdash; Default: None</i><br />
        gripss paramter: Minimum coverage in the normal sample to determine somatic status. Default value:8
</p>
<p name="SVPipeline.exclude_filters">
        <b>SVPipeline.exclude_filters</b><br />
        <i>String? &mdash; Default: None</i><br />
        gripss paramter: Exclude filters from the output vcf, separated by ;
</p>
<p name="SVPipeline.symbolic_vcf_format">
        <b>SVPipeline.symbolic_vcf_format</b><br />
        <i>Boolean &mdash; Default: None</i><br />
        Whether to convert the output vcf to the region format or not, default True
</p>
<p name="SVPipeline.cloud_provider_override">
        <b>SVPipeline.cloud_provider_override</b><br />
        <i>String? &mdash; Default: None</i><br />
        Cloud provider to use for the workflow. Currently supported: aws, gcp default: gcp
</p>
</details>


### Advanced inputs
<details>
<summary> Show/Hide </summary>
<p name="SVPipeline.config_file_string">
        <b>SVPipeline.config_file_string</b><br />
        <i>String &mdash; Default: None</i><br />
         Gridss config file content 
</p>
<p name="SVPipeline.single_strand_filter">
        <b>SVPipeline.single_strand_filter</b><br />
        <i>Boolean &mdash; Default: None</i><br />
         Whether to filter out non snp candidates that are on a single strand 
</p>
<p name="SVPipeline.create_assembly_memory_override">
        <b>SVPipeline.create_assembly_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for create_assembly task 
</p>
<p name="SVPipeline.annotate_variants_cpu_override">
        <b>SVPipeline.annotate_variants_cpu_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         cpu override for annotate_variants task 
</p>
<p name="SVPipeline.annotate_variants_memory_override">
        <b>SVPipeline.annotate_variants_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for annotate_variants task 
</p>
<p name="SVPipeline.convert_vcf_format_memory_override">
        <b>SVPipeline.convert_vcf_format_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for convert_vcf_format task 
</p>
<p name="SVPipeline.germline_link_variants_memory_override">
        <b>SVPipeline.germline_link_variants_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for germline_link_variants task 
</p>
<p name="SVPipeline.scatter_intervals_break">
        <b>SVPipeline.scatter_intervals_break</b><br />
        <i>Int &mdash; Default: None</i><br />
         Maximal resolution for scattering intervals 
</p>
</details>

## Outputs
<p name="SVPipeline.output_vcf">
        <b>SVPipeline.output_vcf</b><br />
        <i>File</i><br />
        Final VCF file
</p>
<p name="SVPipeline.output_vcf_index">
        <b>SVPipeline.output_vcf_index</b><br />
        <i>File</i><br />
        Final VCF index file
</p>
<p name="SVPipeline.assembly">
        <b>SVPipeline.assembly</b><br />
        <i>File</i><br />
        Assembly output before UA realingment
</p>
<p name="SVPipeline.assembly_index">
        <b>SVPipeline.assembly_index</b><br />
        <i>File</i><br />
        Assembly output index before UA realingment
</p>
<p name="SVPipeline.realigned_assembly">
        <b>SVPipeline.realigned_assembly</b><br />
        <i>File?</i><br />
        Assembly output after UA realingment
</p>
<p name="SVPipeline.realigned_assembly_index">
        <b>SVPipeline.realigned_assembly_index</b><br />
        <i>File?</i><br />
        Assembly output index after UA realingment
</p>
<p name="SVPipeline.converted_vcf">
        <b>SVPipeline.converted_vcf</b><br />
        <i>File?</i><br />
        Final VCF file in the region (non-breakend) format
</p>
<p name="SVPipeline.converted_vcf_index">
        <b>SVPipeline.converted_vcf_index</b><br />
        <i>File?</i><br />
        Final VCF index file in the region (non-breakend) format
</p>

<hr />

> Generated using WDL AID (1.0.1)
