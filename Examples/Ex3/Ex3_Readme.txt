The following paramter files are examples for GWAS on chromosome 14:

- use of genotype (.raw) and map (.map) files prepared with PLINK software

- estimate ADDITIVE AND DOMINANCE effects of the variants (GWAS_TYPE = add_dom)

- The GWAS model includes a fixed overal mean AND A FIXED REGRESSION COEFFICIENT whose covariate is in 2nd column of perfromance file

- Phenotypes are WEIGHTED : weights are in 3rd columns of performance file

-->    GWAS model: y = overall_mean + (regr_coef*covariate) + breeding_value + g*add + d*dom + e , with var(e_i) = resid_variance / weight_i


- discrete genotypes for GWAS variants (0 1 2) are provided in /seldiv_save/RUMIGEN/DATA/r12/Data_all_chr_50K.raw file prepared with PLINK

- Kinship marker genotype file and GWAS variant genotype file are the same file : /seldiv_save/RUMIGEN/DATA/r12/Data_all_chr_50K.raw

- Kinship marker map file and GWAS variant map file are the same file : /seldiv_save/RUMIGEN/DATA/r12/Data_all_chr_50K.map (prepared with PLINK)


- No Kinship marker surrounding the GWAS variant is excluded form the equations
(SEGM_REM = 0 in Step2 parameter file)

- The user asks the program to automatically select <= 20000 Kinship Markers with MAF > 0.10 to model Breeding Valuesamong the markers in Data_all_chr_50K.map 
(OPTION SelAutoSNPpar 20000 ; OPTION MAF_minimum 0.10)
	
	
	
First lines of performance file /seldiv_save/RUMIGEN/DATA/r12/lait/CAR1/PerfGWAS_L1_typ_lait_CAR1.txt:

1 0.05263 0.99 -857.2516721 -344.518007 FR000001
1 0.04978 0.99 -51.50006697 -121.333006 FR000003
1 0.07295 0.94 -2502.529843 -2082.526487 FR000012
1 0.10037 1 1576.4776648 1223.298278 FR000024
1 0 1 807.45817897 702.565007 FR000032




First lines of Kinship marker genotype files = GWAS variant genotype file /seldiv_save/RUMIGEN/DATA/r12/Data_all_chr_50K.raw:

FID IID PAT MAT SEX PHENOTYPE F0100190_A D5010001_A F0100220_A RGX1000_A RGX2000_A RGX1001_A F0133440_A ...
FR000001 FR000001 FR000123 FR000324 0 -9 2 0 1 0 0 0 0 ...
FR000003 FR000003 FR000123 FR000324 0 -9 2 0 1 0 0 0 0 ...
FR000012 FR000012 FR000123 FR000324 0 -9 2 0 1 0 0 0 0 ...
FR000024 FR000024 FR000123 FR000324 0 -9 1 0 1 0 0 0 0 ...
FR000032 FR000032 FR000123 FR000324 0 -9 1 0 2 0 0 0 0 ...



First lines of Kinship marker map file = GWAS variant map file /seldiv_save/RUMIGEN/DATA/r12/Data_all_chr_50K.map:

1	F0100190	0	776231
1	D5010001	0	904701
1	F0100220	0	907810
1	RGX1000	0	989141
1	RGX2000	0	990551
1	RGX1001	0	990831
1	F0133440	0	993060



		
