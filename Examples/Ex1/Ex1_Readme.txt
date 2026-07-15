The following paramter files are examples for GWAS on chromosome 21:

- estimate additive effect of the variants (GWAS_TYPE = add)

- use allelic dosages for GWAS variants in 8 minimac4 VCF files

- GWAS model: y = overall_mean + breeding_value + g*add + e

- allelic dosages are provided in 8 minimac4 VCF files named:
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_1.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_2.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_3.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_4.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_5.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_6.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_7.dose.vcf
	/travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_8.dose.vcf

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

		
First lines of the first VCF file /travail/ttribout/GWAS_DBO_chr21/SEQ_r46_chr21_1.dose.vcf:

##fileformat=VCFv4.1
##filedate=2026.1.7
##source=Minimac4.v1.0.2
##contig=<ID=21>
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only f
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">
##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description="Marker was genotyped but NOT imputed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1
##minimac4_Command=/big/save/BIN/MINIMAC/EXEC/minimac4 --refHaps refPanel21.m3vcf --haps genotypes_1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	FR0001	FR0003	FR0014	FR0025	FR0008
21	12447	21:12447:T:C	T	C	.	PASS	AF=0.10410;MAF=0.10410;R2=0.00000;IMPUTED	GT:DS	0|0:0.208	0|0:0.208
21	12450	21:12450:CT:C	CT	C	.	PASS	AF=0.03625;MAF=0.03625;R2=0.00000;IMPUTED	GT:DS	0|0:0.072	0|0:0.0
21	12475	21:12475:C:T	C	T	.	PASS	AF=0.02810;MAF=0.02810;R2=0.00000;IMPUTED	GT:DS	0|0:0.056	0|0:0.056


First lines of GWAS variant map file /travail/ttribout/GWAS_DBO_chr21/map_varGWAS:
21 12447 1 1
21 12450 2 1
21 12475 3 1
21 12547 4 1
21 12609 5 1

