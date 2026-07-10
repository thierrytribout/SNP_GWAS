The following paramter files are examples for GWAS on chromosome 21:

- estimate additive effect of the variants (GWAS_TYPE = add)

- GWAS model: y = overall_mean + breeding_value + g*add + e

- allelic dosages are provided in 1 unique minimac3 file named /travail/ttribout/GWAS_DBO_chr21/r46_chr21.mach.dose

- Kinship markers on a segment of 10,000,000 bp on the variants chromosome surrounding the variant are removed from equations
(SEGM_REM = 10000000 in Step2 parameter file)

The user specifies the Markers to be considered to model breeding values in Kinship marker map file /travail/ttribout/GWAS_DBO_chr21/map_SNP.
	
	
	
First lines of performance file /travail/ttribout/GWAS_DBO_chr21/perf_file:

1 -23 FR0001
1 7 FR0003
1 -182 FR0014
1 57 FR0025
1 27 FR0008
1 51 FR0012



First lines of Kinship marker genotype files /travail/ttribout/GWAS_DBO_chr21/SNP_genotypes:

FR0001          10000100000210
FR0003          20000200000220
FR0014          10000100000210
FR0025          20000200000110
FR0008          20000100000120
FR0012          10000100000210



First lines of Kinship marker map file /travail/ttribout/GWAS_DBO_chr21/map_SNP:

4th column indicates which Kinship markers have to be considered
1 776231 1 1
1 904533 2 0
1 904580 3 0
1 904610 4 0
1 904701 5 0
1 907810 6 1
1 989141 7 0
1 990551 8 0
1 990831 9 0

		
First lines of the minimac3 allelic dosage file (allelic dosage only shown for the first 10 variants):
FR0001->FR0001	DOSE	0.208	0.072	0.056	0.007	0.172	0.024	0.000	0.647	0.736	0.736
FR0003->FR0003	DOSE	0.230	0.085	0.056	0.004	0.195	0.038	0.000	0.602	0.730	0.730
FR0014->FR0014	DOSE	0.182	0.098	0.056	0.009	0.201	0.030	0.000	0.625	0.735	0.735
FR0025->FR0025	DOSE	0.200	0.066	0.056	0.007	0.170	0.017	0.000	0.620	0.735	0.735
FR0008->FR0008	DOSE	0.208	0.072	0.056	0.010	0.158	0.020	0.000	0.647	0.733	0.732




First lines of GWAS variant map file /travail/ttribout/GWAS_DBO_chr21/map_varGWAS:
21 12447 1 1
21 12450 2 1
21 12475 3 1
21 12547 4 1
21 12609 5 1

