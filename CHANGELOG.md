## New

- Fix issue with multithreaded mpileup which gave different results based on number of threads
- Add filter and consensus params to clair3 pipeline
- add MIN_AF filter for clair3 (similar to MIN_FRS but for ref calls only)
- add dorado 5.0.0 models for clair3
- Any row with minor allele gets output to main vcf
- MIN_VDB from bcftools can override clair3 and set site to null
- interpret MIN_FRS with minor population to mean het
- Calculate mixed site counts for indels and snps separately
- Mark mixed sites using INFO field, and count clusters
- use own container for clair3
- use container prefix param
- use same container name for both registries
- stop producing final.full.fasta
- INVALID_INDELS now always filtered

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
