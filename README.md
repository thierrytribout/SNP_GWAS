The SNP-GWAS software is a program developed in INRAE, enabling GWAS to be carried out on a large number of individuals with a phenotype and genotypes (up to hundreds of thousands of individuals with gentotypes or allelic dosages for variants at sequence level).
The software allows to jointly apply various options that are not found simultaneously in other available software:
- use of allelic dosages instead of discrete genotypes ;
- phenotype weighting ;
- consider a dominance effect in addition to the usual additive effect for the tested variant ;
- remove markers on the same chromosome as the tested variant, or remove markers on a segment surrounding the tested variant, the length of the segment being specified by user ;
- ….

The program SNP_GWAS is developped in Fortran90.
It requires an implementation of the BLAS and LAPACK libraries. It has been developed and tested with Intel oneMKL but should also work with other compatible implementations such as OpenBLAS.

To take into account the population structure and avoid potential confounding and false positives, most GWAS models include a polygenic effect (g). The genomic relationship matrix G built from SNP genotypes is dense and makes the computations very long or even impossible when the number of individuals is large (e.g. more than 40,000).
In SNP_GWAS program, we replace the polygenic effect g by Ms, i.e. by the sum of the SNP effects s multiplied by the centred genotypes M, for a set of a few tens of thousands of markers distributed across the genome and chosen for their informativeness.
The advantage of this option is that the size of the M’M matrix is equal to the number of SNPs, making the size of the equation system constant regardless of the number of individuals analysed. 

The software has been conceived to perform single-trait GWAS at sequence level, i.e. for a very large number of biallelic variants.
It analyses the chromosomes separately, by parallelizing analyses of groups of variants of a same chromosome. 
The companion program Cumul_Res_SG, also developped in Fortran90, gathers the results from the various variant groups in a single result file, which can then be used as input for Manhattan plots or other analyses.

Authors:
Thierry Tribout and Didier Boichard
Universite Paris Saclay, INRAE, AgroParisTech, GABI, 78350 Jouy-en-Josas, France
thierry.tribout@inrae.fr ; didier.boichard@inrae.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the LICENSE file for more details.
