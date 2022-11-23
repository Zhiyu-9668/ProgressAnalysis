# PRS weights

Weigths to compute PRS used in this project. The `LONGscore` file is weights for general longevity. It includes SNP rsids, variants chromosome position under human genome build 36, 37 and 38.

HG37 folder includes weights for endpoints with rsids, chromosome position under human genome build 37; and HG38 folder includes weights for endpoints with chromosome position under human genome build 38 *without rsids*. Please selecct from these two based on your convinence.

To compute PRS, the easiest way is using [plink2 --score](https://www.cog-genomics.org/plink/1.9/score) command. The syntex looks like below:

```
plink2 --bfile InputGeno --score PRSweight i j k header --out OutputPrefix
```
where `InputGeno` is the prefix of your input individual level genotype file (in plink binary format. If you have other file format please check [here](https://www.cog-genomics.org/plink/1.9/input) for more details)

`PRSweight` is your input weight file;

`i` is the column index in the weight file that has the variant identifier to be match with your input genotype file;

`j` is the column index in the weight file that has the effect allele (index for A1 column in the files provided);

`k` is the column index in the weight file that has the variant weight (index for MegaPRSweight column in the files provided);

`OutputPrefix` is the prefix of your output scores. You should get a .log file and a .sscore file. The `SCORE1_AVG` column is the PRS to be used. To make it easier creating the input for the analysis, you can rename it as *`EndPtCode`*`_PRS`.
