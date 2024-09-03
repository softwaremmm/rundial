

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