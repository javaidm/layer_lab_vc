# Layer Lab Variant Calling Pipeline
This pipeline is developoed using Nextflow and implements the following bioinformatics workflows:
* Germline Variant Calling
  * Short/Indels:
    * Using GATK [HaplotypeCaller based](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
    * [Google Deepvaraint] (https://ai.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html)
  * Copy Number(Structural) Variants Calling
    * Using [TIDDIT] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5521161/)
    
* Somatic Variant Calling
  * Short/Indels:
    * GATK [Mutect2 based](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-)
  * Copy Number(Structural) Variants Calling 
    * [GATK based] (http://genomics.broadinstitute.org/data-sheets/PPT_Somatic_CNV_WKST_ASHG_2016.pdf)
    * SavvyCNV: genome-wide CNV calling from off-target reads. See more [here (https://www.biorxiv.org/content/10.1101/617605v1)]
