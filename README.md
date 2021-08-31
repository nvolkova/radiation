# radiation
Codes for the analysis of C. elegans samples of different DNA repair deficient backgrounds exposed to ionising radiation for the following publication which is currently under submission:

B. Meier, N.V. Volkova, et al. "C. elegans genome-wide analysis reveals DNA repair pathways that act cooperatively to preserve genome integrity upon ionizing radiation"

`Analysis_public.Rmd` contains the codes for data upload and mutation rates and signature investigation.

`Location_codes_for_IR_public.R` is a script for generating mutation location plots per chromosome per genotype and analysing clustering.

`plotting_functions_IR.R` contains the functions necessary for generating mutation signature plots.

`useful_functions.R` contains additional functions used for analysis or data upload.

Instructions for data upload:

- Download Supplementary Table 1, Sheet 1 (Sample_description) from https://www.nature.com/articles/s41467-020-15912-7#Sec24 as sample_description.csv
- Download Supplementary Table 1, Sheet 4 (Mutations_IR_samples) as mutations_ir_samples.csv
- Download Supplementary Table 2, Sheet 1 (Clustered mutations) as clustered_mutations.csv
- Download the archive with VCFs (Supplementary Data 6) from https://www.nature.com/articles/s41467-020-15912-7#Sec24 and extract into a 'Filtered_VCFS_recovered' folder
- For human data analysis: download the supplementary materials from Behjati et al 2016: https://www.nature.com/articles/ncomms12605, and download ICGC data from the ICGC portal for WGS analysed BRCA-UK and BOCA-UK to "ICGC-BRCA-SARC-BOCA" folder
