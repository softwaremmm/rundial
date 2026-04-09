## New

- add MIN_HQ_DP filter of 5 to clair3
- set MIN_FRS for clair3 indels to 0.9 to match SNPs
- add "caller" parameter to filter_process to name output file
- Simplify consensus calling process
- All hets genotypes are treated as if they had MIN_FRS flag, and genotype set to the first allele. Idea is just to use the notion of minor populations
- minor populations are detected in the consensus step, not the filtering step.
- Addition options about how to include/exclude minor population calls from main/support vcf
- Produces a true gvcf with all positions in consistent format.
- Reduce ram usage now that gvcf produced

## 0.7.0

- Rename final.vcf to variants.vcf, and associated channels
- Rename final.full.vcf to all_calls.vcf, and associated channels
- Rename filtered.gvcf to alternate-bcftools.vcf, and associated channels

## 0.6.0

- Mixed sites are now counted differently, to better match what Illumina counts and to be more representative.
    - Mixed sites now counts both hets (from GT=0/1 say) and having a minor population
    - Mixed indels nows counts mutations rather than sites effected. So a large het deletion is still just 1
    - Mixed clusters now uses both snps and indels

## 0.5.0

### Feat
- stop producing final.full.fasta
- Add filter and consensus params to clair3 pipeline
- add MIN_AF filter for clair3 (similar to MIN_FRS but for ref calls only)
- add dorado 5.0.0 models for clair3
- Filters from bcftools can override clair3 and set site to null (if requested in params).
- Stop using VDB filter for bcftools
- Stop using MIN_QUAL filter for bcftools (doesn't work well with minor alleles)
- Minor alleles can be removed (in filter vcf step) based on fraction of reads and strand bias
- interpret MIN_FRS with minor population to mean het
- Any row with minor allele gets output to main vcf in consensus step
- Calculate mixed site counts for indels and snps separately
- Mark mixed sites using INFO field, and count clusters

### Fix

- Fix issue with multithreaded mpileup which gave different results based on number of threads
- INVALID_INDELS now always filtered
- update to cargo 1.88
- add clair3 to rundial container
- use container prefix param
- use same container name for both registries


## 0.4.2 (2024-09-03)

### Fix

- standard indels do not impact first base
- output lines for missing sites
- het masks come from passed changes so bump priority
- filtered snps/ref beat filtered indels
- correct simplifying ref alt logic
- change fasta parser so that prints max 80 a line
- test genome creation report
- allow ploidy 1 genotypes

## 0.4.1 (2024-08-29)

### Fix

- just ignore lock file

## 0.4.0 (2024-08-29)

### Feat

- allow gzipped input and output

### Fix

- use Path as input type

## 0.3.0 (2024-08-27)

### Feat

- fully implement consensus creation
- add incomplete version of consensus

### Fix

- make test directories in test
- update ignore file
- update version number in test data
- remove pretty asserts from normal dependencies
- fix linting
- auto formatting
- refactor parameters to be submodule

## 0.2.0 (2024-08-13)

### Feat

- github ci for testing and bump
- add pre-commit and commitizen file
- implement filtering using custom vcf parser

### Fix

- use highest allele in draw
- mix of formatting, tests and minor changes
- update params file
- linting
- doc strings and refactor
- cargo files
- more linting
- clippy auto linting
- all noodles code
- initial commit
