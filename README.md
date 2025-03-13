# hydrasim
Hybrid metagenome simulator

User inputs:
1. References CSV containing `accession` and `category_id` columns
2. Datasets CSV, indexed by params.dataset_index and containing a `platform` column. If `human_filtered_reads_1` and `human_filtered_reads_1` columns are provided should use the reads at those paths instead of downloading them
3. Other parameters specified at the top of the modules files

Expected outputs:
The output directory will be populated with folders labelled by `category_id`. Within each folder there will be simulated dataset files named `${dataset_accession}_${ref_accession}_${simulator}_${ref_coverage}x_${index}[_R1,_R2,].fq.gz` where the index numbers the `params.num_iterations` simulations for each ref, dataset, simulator and coverage.

Expected behaviour:
Note that this pipeline makes use of a store directory, where datasets and reference fasta are stored so they can be called on several times without needing to download again
Also note that in places maxForks has been used where having more than 1 fork led to errors either due to clashes over reading a file by rasusa, or errors when downloading if too many simultaneous calls.


## Test Run
To Run with test dataset on Bryn platform (from the repo directory):
```
nextflow run nfellaby/hydrasim -r dev --reference_csv test/test_data/test_hcid_accessions.csv --dataset_csv test/test_data/test_datasets_for_hcid_and_respiratory.csv -profile test,docker
```

## Optional Arguments

### Badread read length and standard deviation

Badread's default is --length 15000,13000 (mean=15000, stdev=13000)
these settings determine the length of the fragments spiked, not the final reads. details
on how fragment lengths are created can be found [here](https://github.com/rrwick/Badread?tab=readme-ov-file#fragment-lengths). when spiking smaller materials. shorter read lengths with lower standard deviation 

Below is an example of how the parameter --badread_length can be used to change the read length and standard deviation to 1000 and 500 nucleotides respectively 
```
nextflow run nfellaby/hydrasim -r dev --reference_csv test/test_data/test_hcid_accessions.csv --dataset_csv test/test_data/test_datasets_for_hcid_and_respiratory.csv --badread_length "1000,500" -profile test,docker
```

### Size of dataset output

The amount of the file that was spiked into, which remains after subsetting with the tool [rasusa](https://github.com/mbhall88/rasusa?tab=readme-ov-file#coverage)
is determined with the parameter --dataset_coverage. This is defaulted to 10k which determines the number of bases that reads are subsampled to, as described by the rasusa documentation under '--bases' , the option used by hydrasim

Below is an example of how the parameter --dataset_coverage can be used to change the read length and standard deviation to 100k nucleotides 
```
nextflow run nfellaby/hydrasim -r dev --reference_csv test/test_data/test_hcid_accessions.csv --dataset_csv test/test_data/test_datasets_for_hcid_and_respiratory.csv --dataset_coverage "100k" -profile test,docker
```

### Output coverage

Hydrasim can spike in at various coverages which are set by default as 0.1x, 1x, 10x and 100x 

This can be changed with the parameter --coverages and a list, at present the easiest way to change the coverages spiked in is to alter the params file itself