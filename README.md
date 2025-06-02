# Rundial
A rust based program for filtering VCFs and creating consensus genomes.
Uses the output from bcftools or Clair3.

TODO:
- Clair3 container needs to be copied
- adjust params


## Dependencies
* Docker
* Nextflow
* Cargo for running rust

### Reference Data
Needs a reference fasta and also clair3 models if using.
This is included in the `data` folder. Just run `get_clair3_models.sh` to download the clair3 models.

## Running workflow

Can run with bcftools or with clair3
```bash
nextflow run . --workflow bcftools --input_dir test_data/assemblers --publish_dir results \
	--ref-fasta data/h37rv_20231215.fa.gz

nextflow run . --workflow clair3 --input_dir test_data/assemblers --publish_dir results \
	--ref-fasta data/h37rv_20231215.fa.gz --clair3_models_dir data/clair3_models \
    --basecalling_model dna_r10.4.1_e8.2_400bps_sup@v4.3.0
```

If using modification may need to rebuild container and use `-profile local_docker`
```
docker build -t test_container_rundial .
```

## Testing
Can test rust code with:
```
cargo test
```

Can test nextflow with:
```
nf-test test tests/nextflow/*.nf.test
```

You may need to rebuild the container before testing:
```
docker build -t test_container_rundial .
nf-test test tests/nextflow/*.nf.test --profile local_docker
```

## Parameters and Thresholds
### Minimap 2
No secondary alignments are output, and ont standard params are used.

### bcftools
- Q=10: Min base quality. Only bases above this threshold are counted for AD/ADF/ADR.
- h=100: homopolymer error coefficient, taken from tbpore.
- M=1000: Max read length for BAQ algorithm, mainly to stop computation being too great.
- m: multiallelic caller model is used. Also taken from tbpore.

### Rundial
Filtering params described [here](process/filter_params.yml) and consensus making params [here](process/consensus_params.yml).
Clair3 vcfs have less filtering [here](process/clair3_filter_params.yml) but similar consensus [params](process/clair3_consensus_params.yml).
Thresholds mainly choosen by trial and error.
Main filters are:
- Min read depth of 3/5, and min high quality depth of 2 (using the Q cut-off from earlier)
- No Strand bias/mismatch where there are differences between the forward and reverse strand

## Consensus with two VCFs
When running with clair3, the clair3 VCF doesn't cover the whole genome.
As such the bcftools assembly is used to fill in the blanks. It never adds any passed mutation but will add ref calls and null/filtered sites. It can be configured as to whether to output rows with minor populations or not.
This provides more explanation for the calls made.

When running in this mode, bcftools only calls snps in the mpileup command.

## Tags, Releases, and Committing
Use conventional commits. This is enforced with commitizen validate action and pre-commit hooks:
```bash
pre-commit install
```

This repo uses a standard gitflow approach, but with some changes to deal with docker containers in nextflow:
- There is a pyproject version which should be updated to match semantic version releases (from main/release branches)
- There is an `active_version` controlled by version_bumper which allows develop to have commit hash based versions.
- Every push to develop will cause an action to run `bumper bump <commit-hash> --no-tag --active`. This:
    - bumps the `active_version` in pyproject.toml
    - bumps the container tag used by nextflow processes
    - leaves the pyproject version as is

  Another workflow then builds and pushes the new container.\
  **Important: To deploy this you will need to use the hash of the bump commit, not the hash of the merge commit/container.**

- In a release branch you can create a release candidate with `bumper bump a.b.c-rcX`. This also updates the pyproject version. Pushing the changes and new tag (automatically created) will trigger a build action.
- When release branch is ready for main run `bumper bump a.b.c --no-tag`. Push these changes to main and make a release there to build the container.


---
---


# Details

### Simplifying indels
Some indels can be contracted before being applied. This is to help isolate what change the indel row is actually making.
An example is like:
```
NC_000962.3     2338194 .       A       .       244.589 PASS    DP=59;ADF=16;ADR=23;SCR=51;FS=0;MQ0F=0;AN=2;DP4=16,23,0,0;MQ=60 GT:SP:AD        0/0:0:39
NC_000962.3     2338194 .       ACCCC   ACC,A,AC        69.6435 PASS    INDEL;IDV=41;IMF=0.672131;DP=61;ADF=0,4,1,2;ADR=0,2,1,0;SCR=51;VDB=0.10087;SGB=-0.670168;RPBZ=-0.222826;SCBZ=-1.18038;FS=0;MQ0F=0;AC=2,0,0;AN=2;DP4=0,0,7,3;MQ=60   GT:PL:SP:AD     1/1:125,45,27,94,0,88,98,3,59,92:0:0,6,2,2
NC_000962.3     2338195 .       C       .       145.588 STRAND_BIAS     DP=20;ADF=2;ADR=10;SCR=51;FS=0;MQ0F=0;AN=2;DP4=2,10,0,0;MQ=60   GT:SP:AD        0/0:0:12
NC_000962.3     2338196 .       C       .       147.588 STRAND_BIAS     DP=29;ADF=1;ADR=10;SCR=51;FS=0;MQ0F=0;AN=2;DP4=1,10,0,0;MQ=60   GT:SP:AD        0/0:0:11
NC_000962.3     2338197 .       C       .       198.589 PASS    DP=41;ADF=3;ADR=11;SCR=50;FS=0;MQ0F=0;AN=2;DP4=3,11,0,0;MQ=60   GT:SP:AD        0/0:0:14
NC_000962.3     2338198 .       C       .       220.589 PASS    DP=54;ADF=4;ADR=13;SCR=50;FS=0;MQ0F=0;AN=2;DP4=4,13,0,0;MQ=60   GT:SP:AD        0/0:0:17
```
There are 2 C's deleted by the indel.

**But how to make this change?**
The preference in this code is to remove the first C's as these have fewer reads assigned by bcftools (20/29 instead of 41/54) so this better matches what bcftools is doing.

As such to simplify we remove later bases first, and then the front:

ACCCC -> ACC (pos 2338194) => ACC -> A (pos 2338194) => CC -> "" (pos 2338195)

Alternative:
ACCCC -> ACC (pos 2338194) => CC -> "" (pos 2338198)


However, the front base is special, as it should normally refer to a base which is not affected by the indel (in a normalised form). So AAAA -> AA should have the effect of A--A

## Row preference
The way in which rows of the vcf are applied depends on the priority placed on the different kind of changes. This is were all the subtle errors come about.

One quirk to be aware is that passed indels trump filtered ref calls. So if you have a potential deletion with GT 0/0 followed by a single base ref call with some filter like STRAND_BIAS. The resulting base call will be ref rather than F.

Current ordering (Note being het have effect):
1. Filtered Null call
2. Filtered indel
3. Filtered snp/ref
4. Passed Null call
5. Passed homologous ref call indel (indel with gt 0/0)
6. Passed homologous ref call
7. Passed indel
8. Passed snp


## Filters
Most filters are fairly easy to understand from the header line added.

### MIN_FRS and MIN_AL
These both get called MIN_FRS in the resulting VCF as they are checking that the called allele has sufficient support.
- MIN_AF compares the depth of the called allele to the total depth.
- MIN_FRS compares depth of the called allele to the sum of depth of all called alleles.
These are not the same! Some potential alleles have low support so never appear in the vcf but do contribute to overall depth. And in BCFTools the overall depth includes reads of insufficient quality but the alleles depths do not.

Currently MIN_AF is only used for clair3 when checking that ref calls have sufficient support.

### Low_VDB Overriding filters
VDB is Variant distance bias. A low VDB indicates that the variant (snp/indel) appears in the same position in all the alignements. Equivalently all the reads seem to start or end their alignment at the same place which is odd as we'd expect this to be random.

It works as an extra check for something funny with alignments even when mapping quality is good.

However, only bcftools provides it as a metric. So when taking the clai3 route `Low_VDB` must be set as an overriding filter in the params.yml.
This way the variant in the clair3 vcf will inherit the `Low_VDB` flag.
