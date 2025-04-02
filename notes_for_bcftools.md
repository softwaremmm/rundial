# Notes and warnings for using BCFTOOLS

## mpileup and call
The command can be run with `-I` to note calculate indels.
But even if the mpileup calculates indels, the `call` command can ignore them.
This gives identical results to having never calculated indels in the mpileup.

The rows from bcftools keeps indels and snps seperate.
So in the output there will often be a ref row followed by an indel row.
Beware that `norm` may merge this unhelpfully

bcftools can produce rows which claim to indels, but have no ALT. This is because alts not in GT get removed
Can be added in when calling.

The genotype selection can show very odd behaviour, where they go counter to AD:
NC_000962.3     1552    .       TAAA    TAA,TA,T        189.899 MIN_FRS;STRAND_BIAS;STRAND_MISMATCH     INDEL;IDV=181;IMF=0.800885;DP=226;ADF=0,46,3,2;ADR=12,19,23,7;SCR=198;VDB=0.829341;SGB=-0.693147;RPBZ=-2.10583;MQBZ=-0.498617;SCBZ=-0.712391;FS=0;MQ0F=0;AC=0,1,1;AN=2;DP4=0,12,51,49;MQ=60 GT:PL:SP:AD     2/3:250,255,177,145,99,125,107,137,0,99:34:12,65,26,9
