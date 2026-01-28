!========================================================================
program SNP_GWAS_publi
  !=======================================================================
  
  ! !28/01/2026 : correction in log and removed unnecessary zeroing of DOSES

  ! 05/09/2025 : modif pour compter le nb de SNPpar pour lesquels un individu a genotype inconnu (5) et les lister sur une seule ligne plutot que 1 SNP par ligne (LOG trop long)

  ! 04/09/2025 : modif test sur valeur de interm_4 avant calul PEV pour sta test effet variant :  qlq variants ont interm_4 calcule < 0 car ZpZ-interm_2 = negatif tres petit --> calcul de stat_test plante.
  !              ces variants avaient StdDev(DOSE) tres proche de 0 (ex : 4.3E-8 --> on ne peut de toute maniere rien tirer de tels variants)
  !              --> on remplace test if(interm_4 ne 0) par if(interm_4 > 0) car en theorie interm 4 ne peut pas etre negatif sauf si limite precision machine
  !              + report des modifs faites par JVDP pour RUMIGEN (sans impact sur les resultats par rapport a version initiale TTR)

  ! 03/09/2025 : modif pour allouer test_EffVarGWAS_1 (...)  meme si calcVarRes=exact sinon test combine sur calcVarRes=exact or > test_EffVarGWAS_1 > seuil plante

  ! 02/09/2025 : correction bug  : on remplace test sur ecarttype(dose) diff 0 par presence d au moins 5 homozygotes xx et de 5 heterozygotes quand genotypes discrets
  !              en STEP_3 sinon on passe au variant suivant
  !
  !              Ajout d un test sur le signe de exp4_1 avant de calculer sqrt(exp4_1*  ) pour calculer le test approche de l effet estime du variant

  ! 26/06/2025 : correction bug lorsque variant GWAS positionne apres le dernier SNPparente du chromosome avec option pour retirer SNPparente d une fenetre encadrant le variantGWAS


  ! 17/01/2025 : home made function to calculate P-values from test staistics because problem with Burkardt function for very low p-values 

  ! 17/12/2024 : COJO OPTION extended to add_dom model for GWAS

  ! 27/11/2024 : NAG function G01ECF previously used to calculate p-values for Variants estimated effects replaced by free progf90 subroutines
  !              chi_square_cdf gamma_cdf normal_01_cdf and function r8_gamma_inc written by John Burkardt under MIT licence

  ! 20/11/2024 : on part du programme SNP_GWAS_rumigen_m2 du 19/11/2024 valide et dont la LOG a ete netoyee, et on retire tout ce qui est lie
  !              aux modeles rumigen_m2 et mut_rec et tout ce qui est lie aux effets ROH (donc lecture fichier .hom, ...)


  ! 07/11/2024 : on traduit les messages de la LOG d execution en anglais ; version originale du programme avec messages en francais = SNP_GWAS_rumigen_m2_french.f90

  ! 22/10/2024 = copie de superGWAS_vROH_outEffSNP.f90 renommee en SNP_GWAS_rumigen_m2 pour retirer parametres inutiles : NB_IND, NUMBER_OF_TRAITS, ...

  ! compil avec options releasedesespoir = release mais on ajout -CB

  ! sept 2024 : pour RUMIGEN : modification pour modele de GWAS RUMIGEN_M2 = eff_11vs12 + eff11Rvs12 + eff22vs12 + eff22Rvs12

  ! 03/05/2024 : on modifie le programme pour permettre la lecture du fichier typages parente et carte typage parente au format PLINK

  ! possibilite qu il y ait plus d individus dans fichier TYPAGES PARENTE que dans fichier PERF, car si format PLINK le fichier typages parente peut etre le meme que le fichier TYPAGES GWAS
  !             permettre fichier parente et fichier carte parente format typ_eval ou format plink

  ! on conditionne la présence de l effet ROH a la frequence des status ROH 0 et 1 du variant dans a population : par ex entre 5 et 95 pct pour que l effet soit estimable

  ! 28/03/2024 : semble OK quand 1 seule lecture, mais probleme probable quand plusieurs lectures sur exemple 15 ind ???

  ! 21/03/2024 : modif du programme pour etablir le statut ROH de chaque (variant x individu) a l etape 2 a partir du fichier typages_GWAS et du fichier PLINK (segmentsROH x individus) 

  ! 07/03/2024 : on part de la version /g2b/ttribout/superGWAS/SRC4/superGWAS_V8sr_SFR.f90 pour ajouter l effet statut_ROH (0/1) dans la partie add et dom du modele --> se et test pour cet effet

  ! modification du 28/07/2023 : modification pour ajouter une option permettant de realiser les GWAS sur des genotypes discrets (0 1 2) poiur les variants GWAS meme si on lit un
  !                              fichier de DOSES format MINIMAC

  ! modification du 28/07/2023 : modification de la version superGWAS_V8sr.f90 pour passer un modele comparant (homozygotes sauvages et heterozygotes) a (homozygotes mutes) au variant GWAS 
  !                              pour tester determinisme recessif complet pour anomalies mortalite S FRITZ

  ! modification du 27/07/2023 : on passe dimvec = longueur du vecteur inf_LEFT_EP en integer(kind=8) et on modifie calcul de dimvec et function TI(i,j,dimvec)
  !                              car tentative de GWAS en conservant 47061 SNP plantait car dimvec devenait negatif en integer
  !                              A priori pas d'incidence car DBO avait conserve 14205 SNPparente = une matrice de 200 000 000 elements --> integer suffisait 

  ! modification du 11/04/2023 : on verifie en ETAPE 1 que les individus dans le fichier PERF (eventuellement avec une perf codee manquante) sont identiques aux individus du fichier TYPAGES_PARENTE

  ! modification du 30/03/2023 : modif pour ajouter un variant specifie par utilisateur en effet fixe dans le modele pour voir si autres variants ont encore un effet significatif (approx COJO)
  !                              ATTENTION : pour l instant on va considerer que le variant FIXE est sur le meme chromosome que les variants GWAS et qu on supprime tous les SNPpar du chromosome
  !                                          cela simplifie les choses car on gere un seul set de SNPpar a supprimer comme quand pas COJO

  ! modification du 28/03/2023 : on conditionne la presence de l EFFET DE DOMINANCE dans le modele pour un varGWAS (si demande par utilisateur) a valeur des doses / typages

  ! modification du 20/03/2023 : on stocke les genotypes centrés aux SNPpar utiles pour tous les individus pour voir si gain de temps construction M' W M

  ! 13/03/2023 : MODIF PARTIE 2 pour lire fichier DOSES_VAR_SEQ en plusieurs fois (par paquets de N_Var_SEQ pour reduire quantite de memoire necessaire a PARTIE 2

  ! modification du 22/02/2023 : modification pour permettre de lire des fichiers de genotypes (nb alleles 2 portés = 0 1 2 ou 5=inconnu) AU LIEU DES DOSES

  ! modification du 16/02/2023 : OPTION pour que l utilisateur puisse choisir le nb de SNPparente utilises dans les calculs

  ! version V6 : on modifie pour permettre d avoir dans fichier DOSES lu dans partie 2 des individus absents du fichier PERF et du fichier TYPAGES_PARENTE
  !              (evite de creer des fichiers doses de travail)

  ! version V5_f : on ajoute une option pour calculer (1) la variance residuelle exacte pour chaque variant ou (2) une base de variance residuelle commune a tous les variants du Sous Groupe
  !                ou (3) une base de variance commune pour tous les variants PUIS la variance residuelle exacte pour les variants dont le TEST avec variance residuelle approchee est > seuil

  ! on cree des sous groupes de variants GWAS pour la partie 3 en fonction du nb de variants max par job fixe par l utilisateur

  ! V5_d : on separe partie 3 (relecture ligne par ligne des doses dans fichier binaire et GWAS) de partie 2 (creation groupes d eclusion et ecriture fichier binaire doses)

  ! dans V5_b on cree 1 fichier binaire Doses et 1 fichier parametres pour l execution de la partie 3 par groupe d exclusion de SNPparente

  ! dans V5 on modifie programme pour creer un fichier sortie de DOSES par gourpe de VarGWAS appartenant au meme segment d exclusion de SNPpar
  !     pour pouvoir en etape suivante lancer 1 job par segment d exclusion

  ! dans cette version V2_fmat on fait les calculs pour obtenir inf_LEFT_EP_act avec des fonctions matricielles de BLAS pour essayer de reduire les temps de calcul

  ! on intervertit l ordre des boucles sur eff1 et eff2 pour voir si gain de temps

  ! dans cette version on inverse les dimensions de la table de genotypes SNP_PARENTE = GENOPAR pour voir si on gagne du temps

  ! version non parallelisee ; partie inversion de LEFT_EP et partie GWAS sur variants non individualisees

  ! dans cette version on n integre pas R-1 dans les MME, mais on utilise Rho = varRes/varGenet

  ! 17/11/2022 : ATTENTION : VERIFIER LA FORMULE DE CALUL DE RANK

  ! 09/11/2022 : inversion indirecte de LEFT_EP partielle a partir des blocs de inverse(LEFT_complet) valide

  ! ATTENTION : POUR L INSTANT EN MONOCARACTERE



  ! 03 oct 2022 : debut developpement 
 
  implicit none

  character(len=100),parameter:: title=' SNP_GWAS_publi software - version with blas lapack without NAG' 

  character*15 fx1,fx2,fx22,fxcojo

  character(len=128)::jour
  integer::anlim=2046,mlim=04,jlim=01  ! date limite d executable = 01 avril 2046

  logical testficdir

  integer::i,j,k,d,ef,ef1,ef2,ind,ii
  integer::ip_snp

  !        Types of effects   
  integer,parameter::effcross=0,& !effects can be cross-classified 
       effcov=1     !or covariables 


  integer::ioperf,iobilan,iocojo

  !        Types of random effects
  integer,  parameter ::  g_fixed=1,&       ! fixed effect
       g_diag=2, &       ! diagonal
       g_AD=13,&     ! additive animal direct in Hybrid Model
       g_AM=14,&     ! additive animal maternal in Hybrid Model
       g_last=16     ! last type

  character (250)    ::   parfile, &     !name of parameter file
       datafile, &    !name of data set
       typparfile, &     !name of genotypes file used for breeding values
       fmttyppar, &      ! format du fichier typages des variants utilises pour la parente : typ_eval (sans espaces) ou plink (avec espaces)
       mapparfile, &     !name of map file for SNP used for breeding values / 4 colonnes : num_chromos / pos_sur_chrom / num_marker_dans_fich_typpar / INCLUSION-EXCLUSION PAR UTILISATEUR
       mapparfile_bis, &  !alternative name of map file for SNP used for breeding values / 4 colonnes : num_chromos / pos_sur_chrom / num_marker_dans_fich_typpar / INCLUSION-EXCLUSION PAR PROCEDURE AUTOMATIQUE
       typgwasfile, &     !name of genotypes file used for GWAS
       fmttypgwas,&       !format du fichier typages pour les variants consideres dans les GWAS (pour l instant : minimac (doses) ou typ_eval (nb alleles 2 sans espaces) ou plink (nb alleles 2 avec espaces et colonnes suppl)
       mapgwasfile, &     !name of map file for SNP used for GWAS / Si format ne plink 4 colonnes : num_chromos / pos_sur_chrom / num_marker_dans_fich_typGWAS / utile_GWAS_0_1
                                ! si formt PLINK 4 colonnes : num_chromos / nom variants GWAS (inutile) / position en cM (inutile) / position en pb
       cojofile, &       ! name of the file listing the SNP considered fixed in cojo analyses
       execpath, &       ! path for software executable to be called for step 3
       location          ! chemin des fichiers binaires inf_LEFT_EP et RIGHT_EP qui seront crees en fin de STEP1 et lus en debut de STEP2

  character (320) :: format_res

  integer::nb_champs,nb_chGWAS,INTtemp ! nb de champs dans fichiers carte variants parente et variants GWAS et indicateur (temporaire) de INCL/EXCL pour le variant GWAS en cours de lecture

  character,parameter::tab=achar(9)

  character (300) ::  chem_infos1, chem_bin0, chem_dosesCojo, rep_G_SG

  character (25) :: suppl_cojo

  real(kind=8)::mafmin=0.01d0
  real(kind=8)::x(3)
  character(20)::xc(3)
  integer::detailed_log=1

  character(len=11)::gwastype

  character(len=10)::StratElim='window'
  character(len=10)::CalcVarRes='optim'
  integer::NbMaxVarGWAS=5000 ! nombre maximum de variants GWAS par sous groupe de GWAS en partie 3 (modifiable par utilisateur par OPTION NbMaxVarGWAS)

  integer(kind=8)::segmrem

  integer::steptodo

  integer :: max_string_readline = 800 

  integer::indic_cojo=0    ! indicates if there is (1) or not (0) a cojo variant in the GWAS
  integer::neq_cojo=0      ! = number of equations associated to the Cojo Variant (1 if add, 2 if add_dom)
  integer::neff_cojo       ! = number of effects associated to the Cojo Variant (1 if add / 2 if add_dom)

  integer::VarDiffCojo=1   ! indicates, if COJO analyzis, if genomic information of current variant is different from genomic information of COJO Variant (=consider/skip Variant)

  integer :: ntrait=1                    !number of traits = TOUJOURS 1
  integer :: neff                        !number of effects excluding SNP effects for dirhyb and mathyb

  integer :: nb_threads_mkl = 12   ! number of threads for MKL, default=12

  real(kind=8) :: mis=-9999.0d0    !value of missing trait/effect

  integer,allocatable :: pos_y(:)        ! positions of observations
  integer,allocatable :: pos_weight(:)   ! positions of weight of records; zero if none
  integer::posIdchar,nbind      ! positions de l identifiant alphanumerique et nb d individus avec perf et typesdans le fichier performances
  integer::nlperf                        ! nb d elignes lues dans fichier DATA
  integer :: iprint_eff=6                ! to avoid printing known characteristics of long lists of hundreds of effects


  integer,allocatable :: pos_eff(:,:),&  !positions of effects for each trait
       nivanim(:,:),&   ! niveau d effet pour l individu et le caractere considere
       matNivanim(:,:,:),& ! table des niveaux d effets pour tous les individus
       nbcarperf(:),&   ! nb de car avec perf non manquante pour l individu en cours de traitement a la lecture de DATA
       nlev(:),&       !number of levels excluding SNP effects
       effecttype(:),& !type of effects
       nestedcov(:,:),&!position of nesting effect for each trait if the effect is nested covariable
       randomtype(:),& ! status of each effect, as above
       randomnumb(:) ! number of consecutive correlated effects


  integer :: temRand=0  ! indicator of presence of at least 1 non genetic random diagonal effect in the model
  !integer,allocatable :: adresse(:,:),&    ! position de niveau 0 de chaque effet de chaque caractere dans vecteur sol pour trouver la position de ni niveau l de l effet e du car t


  real(kind=8), allocatable :: r(:,:),&      ! residual (co)variance matrix
       rinv(:,:),&                           ! and its inverse
       g(:,:),&                              ! genetic (co)variance matrix
       ginv(:,:),&                           ! and its inverse
       rand(:,:,:),&                         ! The random (co)variance matrix for each trait and its inverse
       randinv(:,:,:),&                      ! and its inverse
       gsnpparinv(:,:),&                     ! inverse de la variance genetique des SNP = 1 / (g / 2*sumpq)
       Rho(:,:),&                            ! varRes / (varG / 2*sumpq)
       RhoDiag(:,:,:)                          ! varRes / (VarRbdDiag)


  !real(kind=8), allocatable :: varQTLm1(:,:)   ! Inverse of the genetic (co)variance matrix for QTL effects

  integer,parameter::   maxcorr=20     ! maximum number of correlated effects; used for g 

  integer,parameter::    io_d=50,&        ! unit number for data file
       io_s=30,&         ! unit number for solution file
       io_p=40,&         ! unit number for parameter file
       io_typpar=42,&    ! unit number for genotype file for markers used for relationships
       io_doses=43,&     ! unit number for dose file for markers used for GWAS
       io_mappar=44,&    ! unit number for map file for relationship markers
       io_mapparout=46,& ! unit number for alternative map file for relationship markers created by SNP_GWAS when SelAutoSNPpar OPTION selected by user
       io_mapgwas=45,&   ! unit number for map file for GWAS markers
       io_invleft=60,&   ! unit number for binary file containing inverse(LEFT_member)
       io_right=61,&     ! unit number for binary file containing RIGHT_member
       io_IdCojo=62,&    ! unit number for file containing the ID of cojo variants
       io_dosesCojo=63,& ! numero d unite pour le fichier qui va contenir les doses du variant Cojo cree en partie 2 et lu en partie 3
       io_infos0=64,&    ! numero d unite pour le fichiers informations sur tous les Groupes-SousGroupes de GWAS dans repertoire racine
       io_infos1=65,&    ! numero d unite pour les fichiers informations sur Groupes de GWAS intra sous_repertoire de sous-groupe
       io_dosesbin0=70   ! base de numero d unite pour les fichiers binaires de doses/typages variantsGWAS en sortie ou entree par groupe de SNPpar exclus


  integer::io


  logical :: testpresfic

  integer::data_len ! position de la derniere variable numerique utile du fichier DATA
  real(kind=8),allocatable::indata(:) ! one line of input data, only numeric fields
  real(kind=8),allocatable::y(:),weight_y(:),weight_cov(:,:) ! performances et poids des performances et covariables d une ligne du fichier DATA
  real(kind=8),allocatable :: matWeight_Cov(:,:),matPERF(:,:),matPOIDS(:,:),vecPERFpond(:),vecPERFpond2(:) ! tables contenant les covariables, les perf et les poids des perf de l ensemble des individus
  integer::aniperf ! identifiant numerique de l individu en cours de traitement dans le fichier performances, FORCEMENT DE 1 a NBINDIVTOT SANS TROU

  integer(kind=1),allocatable::genopar(:,:),genopartemp(:)

  real(kind=8),allocatable::SNPparUtCentr(:,:) ! genotypes centres aux SNPparente utilises dans les equations PREMULTIPLIES PAR SQRT(poids du phenotype de l individu)
  real(kind=8),allocatable::mattemp(:,:)  ! resultat de M' W M + I Rho, matrice temporaire

  character(200)::informatP1=''     ! chaine de caracteres qui recevra le format de lecture du fichier typages parente format typ_eval

  character(200)::formatPAR=''     ! chaine de caracteres qui recevra le format d ecriture du fichier alternatif typages Parente
  character(200)::informatT=''    ! chaine de caracteres qui recevra le format de lecture du fichier typages variants GWAS fomat typ_eval

  !character(200)::informatD=''     ! chaine de caracteres qui recevra le format de lecture du fichier dosages
  character(len=10000000)::linetyp
  integer::nc,nbSNP_0,nbsnp,len_linetyp,ifail,nltyppar,IFAIL2,tt
  integer :: lc   ! longueur de l identifiant alphanumerique des individus
  parameter (lc=20)    ! longueur de l identifiant (20 par defaut)
  integer(kind=4),dimension(:),allocatable::ordperf,ordinvperf
  character(len=lc),dimension(:),allocatable::animperf(:),animPar(:)
  character*lc :: workperf,workdoses,workTypPar

  character(len=lc)::NoNatTemp,NoNatTemp2

  ! declarations pour calcul des frequences allele 2 marqueurs parente
  real(kind=8),allocatable::freq(:)   ! vecteur des frequences de l allele 2 observees
  real(kind=8),allocatable::freq2(:)  ! vecteur des frequences de l allele 2 dans lequel on a mis a 0 ou 1 les frequences < ou > a MAFmin  
  real(kind=8),allocatable::frequtil(:) ! vecteur des frequences de l allele 2 pour uniquement les SNP utiles, utilise pour recentre les typages
  real(kind=8),allocatable::nballeles(:) ! vecteur contenant le nb total d alleles utilisables pour le calcul des frequences (5 = typage inconnu)
  integer(kind=1),allocatable:: mostFreqTyp(:) ! vecteur du genotype observe le plus frequent pour chaque SNP utilise pour remplacer typages manquants
  integer :: nbtypindet  ! compteur du nb de typages indetermines par individu remplaces par le genotype observe le plus frequent au SNP correspondant
  integer(kind=8) :: nbtypindet_tot  ! compteur du nb TOTAL de typages indetermines remplaces par le genotype observe le plus frequent au SNP correspondant
  integer,allocatable::liste_SNPpar_inco(:) ! vecteur temporaire des numeros de SNPparente avec typage inco d un individu
  integer::nbIndTypInco ! nb indiv avec typage inconnu a au moins 1 marqueur parente
  integer::PrintMissGeno=0 ! indique si l utilisateur veut par OPTION (=1) ou ne veut pas (=0, default) printer tous les individus qui ont des genotypes inconnus aus marqueurs parente

  integer::nbtyp,nblues
  real(kind=8)::sumpq_all   ! somme sur l ensemble des SNP de freq*(1-freq), utile pour calculer la variance genetique associee a 1 marqueur
  real(kind=8)::sumpq  ! somme sur l ensemble des SNP tq MAF >= MAFmin de freq2*(1-freq2) ET AYANT LA MEME VARIANCE, utile pour calculer la variance genetique associee a 1 marqueur   
  integer,allocatable::SNP0toSNP(:)  ! vecteur donnant la correspondance entre numero de SNP dans fichier typages et numero de SNP utile
  integer,allocatable::SNPtoSNP0(:)  ! vecteur donnant la correspondance entre numero de SNP utile et numero de SNP dans fichier typages

  ! declarations pour carte marqueurs parente
  integer,allocatable::MAPPAR(:,:) ! table carte marqueurs parente : 1ere col = num chrom / 2eme col = position en bp sur chromosome
  integer,allocatable::FL_PAR(:,:) ! table numero 1er et dernier marqueur parente de chaque chromosome
  integer::nlmappar,chrPrec,posPrec,numPrec
  integer::mappartemp(4)
  integer::iomappar
  integer,allocatable::incl_par(:)

  ! declarations pour selection des SNPparente choisis suite a OPTION 
  real(kind=8),allocatable::LongChromPar(:),MAF(:)
  integer,allocatable::nbSegmPar(:),InfoSegmPar(:,:,:)
  integer::nbchrompar ! nombre de chromosomes trouves dans le fichier CARTE SNPparente
  integer(kind=8),allocatable::PosPremLastSNPpar(:,:),PosSNPpar(:)
  integer::nbSNPpar_consut,nbSNPparAuto,numSegmPar,bta,BestSNPpar,SNPpar,nbSNPpar_retenus
  integer::SelAutoSNPpar=0    ! temoin indiquant si l utilisateur a demande (1) ou pas (0) une selection automatique de nbSNPparAuto SNPparente
  real(kind=8)::LongTotGenomePar,MAFmaxSegm,DistBest,MedPos
  integer(kind=8)::LongSegmPar



  ! declarations pour carte marqueurs GWAS
  integer,allocatable::MAPGWAS(:,:) ! table carte marqueurs GWAS : 1ere ligne = num chrom (identique pour tous les variants) / 2eme ligne = position en bp sur chromosome / 3eme ligne = incl/excl des GWAS
  integer::nlmapgwas,posgwasPrec,numgwasPrec,nbGWAS
  integer::mapgwastemp(4)
  integer::iomapgwas
  integer,allocatable::GWASuttoMAP(:) ! vecteur de nbGWAS elements = nb de variants du fichier doses pour lesquels on veut faire un GWAS = incl/excl=1 / element i = numero du variant dans MAP_GWAS
  integer,allocatable::MAPtoGWASut(:) ! vecteur de nlmapgwas elements = nb total de variants dans fichier CARTE var GWAS = utiles + eclus des GWAS / MAPtoGWASut(i) =0 si exclus = j si jeme varGWAS inclus

  ! declarations pour DOSES des variants GWAS
  real(kind=4),allocatable::DOSES(:,:),DOSEStemp(:)
  character(len=lc),dimension(:),allocatable::animdoses(:)
  integer(kind=1),allocatable::TYPgwasTemp(:)
  character(len=((2*lc)+2))::IdDosesTemp
  character(len=4)::chdose
  integer::nldoses,pos_sup
  integer(kind=4),dimension(:),allocatable::orddoses,ordinvdoses
  integer(kind=4),allocatable::DosesUtil(:),indDosUt(:)
  integer(kind=4),allocatable::TypParUtil(:),indParUt(:)
  integer(kind=4),dimension(:),allocatable::ordTypPar,ordinvTypPar
  integer::indice_temp
  integer::nbTGut
  real(kind=8),allocatable::nballelesTG(:)
  integer(kind=1),allocatable:: mostFreqTG(:) 
  real(kind=8),allocatable::freqTG(:) 
  integer::nb_varGWAS_fichtyp
  real(kind=8)::MaxMem=1d+10 ! 10Go par defaut
  integer::nb_tranches,DosUtlues,DosUtTrav,nb_lectures,lecture
  integer,allocatable::cadre_lec_doses(:,:)
  real(kind=8)::nDutLues
  integer::meth_mpm=1   ! indicateur d option =1 si genotypes CENTRES aux SNPpar non stockes (plus long mais moins de memoire) ou =2 si  genotypes CENTRES aux SNPpar stockes (plus de memoire mais plus rapide ?)
  integer::doses_to_geno=0 ! indicateur d'option = 0 si l utilisateur veur conserver les doses pour varGWAS ou si il veut convertir doses en genotypes 0 1 2 


  ! declarations pour variants Parente a eliminer des equations pour un variant GWAS donne
  integer::chrom_anal,groupelim
  integer::Felim,Lelim,LimInf,LimSup,VarGWAS,PremExcl,DernExcl
  integer,allocatable::SegmElim(:,:),Bilan_SegmElim(:,:),Bilan_SousGr(:,:)
  integer(kind=8),allocatable::lim_theo(:,:),lim_oper(:,:),lim_SNP(:,:)
  integer,allocatable::nbSegmExcl(:)
  integer::nbSGr,batch

  ! declarations pour construction matrices
  integer::neq_ep,neq_ep_2
  real(kind=8),allocatable::inf_LEFT_EP(:),vecGENO1(:),vecGENO2(:),vectemp(:)
  integer,allocatable::FIRSTNIV_EP(:)
  real(kind=8)::val


  !declarations pour inversion de LEFT_MB
  real(kind=8),allocatable::inf_LEFT_EP_act(:,:),A12A22m1(:,:),A21(:,:)
  integer(kind=8)::dimvec ! dimension du vecteur stockant LEFT_EP sous forme compacte triangulaire inferieure


  ! declarations pour la creation du bloc SNP_Parente*SNP_Parente a exclure de l inverse de LEFT_EP
  integer::nEXCL,nEXCLmat,excl1,excl2
  real(kind=8),allocatable::matEXCL(:,:)


  ! declarations pour la construction initiale de X'y M'y
  real(kind=8),allocatable::RIGHT_EP(:),RIGHT_EP_act(:,:)

  ! declarations pour la partie du modele correspondant au GWAS
  integer::neff_GWAS_lu ! nombre de covariables dans le modele pour le variant GWAS (1 si gwastype=add / 2 si gwastype=add_dom) : MODELE DEMANDE PAR UTILISATEUR DANS FICH_PAR
  integer::neff_GWAS ! nombre de covariables dans le modele pour le variant GWAS (1 si gwastype=add / 2 si gwastype=add_dom ) : NOMBRE EFFECTIF D APRES TYPAGES REPRESENTES OU DOSES LUES POUR LE varGWAS
  integer::TemDom_lu=0 ! indicateur de la presence demandee par utilisateur de l effet de dominance pour les variants GWAS dans le modele (0 = absent / 1 = present)
  integer::TemDom=0 ! indicateur de la presence effectif de l effet de dominance pour les variants GWAS dans le modele (0 = absent / 1 = present)
  integer,allocatable::vecTemDom(:)
  integer::coldom ! numero de colonne de l effet de dominance dans Z'Z
  !real(kind=8),allocatable::XpZMpZ(:,:) ! matrice [XpZ MpZ]
  real(kind=8),allocatable::ZpXZpM(:,:) ! matrice [ZpXZpM]
  real(kind=8),allocatable::ZpZ(:,:) ! matrice [ZpZ]
  real(kind=8),allocatable::detZPZ(:)   ! determinant de la matrice ZpZ pour chaque varGWAS
  real(kind=8),allocatable::ZpY(:,:) ! matrice [ZpY]
  real(kind=8),allocatable::vecDosesPond(:,:) ! vecteur (si GWAS additif) ou matrice (si GWAS additif+dominance) contenant les doses pour le GWAS = Dose*sqrt(poids)
  real(kind=8),allocatable::vecDoses(:,:) ! vecteur (si GWAS additif) ou matrice (si GWAS additif+dominance) contenant les doses pour le GWAS = Dose
  real(kind=8),allocatable::vec_stat_doses(:,:)
  integer,allocatable::tem_possible(:,:)
  real(kind=4),allocatable::vecDosesSD(:)
  real(kind=4)::dose_max,dose_min
  real(kind=8),allocatable::interm_1(:,:),interm_2(:,:),interm_3(:,:),interm_4(:,:),interm_5(:,:),interm_6(:,:),sol_VarGWAS(:,:)
  integer::DIFFELIM

  ! declarations pour partie 3 si COJO=1
  real(kind=8),allocatable::col_cojo(:,:),row_cojo(:,:)
  real(kind=4),allocatable::vecDosesSD_cojo(:,:)
  integer::numvar_cojo
  integer::nlcojo ! nombre de ligness dans le fichier des variants COJO (1 pour l instant)
  integer::temp3(3) ! vecteur temporaire
  real(kind=8),allocatable::B22(:,:)
  real(kind=8)::alpha
  real(kind=8),allocatable::A22(:,:)
  real(kind=8),allocatable::B12(:,:),B11(:,:)
  integer,allocatable::matNivanim_cojo(:,:,:)  ! table des niveaux d effets pour le ou les effets du VariantCOJO pour tous les individus
  real(kind=8),allocatable :: matWeight_Cov_cojo(:,:)  ! tables contenant les covariables pour le ou les effets COJO de l ensemble des individus
  real(kind=8),allocatable :: motif1_cojo(:,:) 

  ! declarations pour le calcul de la variance residuelle, de la variance d erreur et des tests des effets des variants GWAS
  real(kind=8),allocatable::BetaS(:,:)   ! vecteur qui va contenir les solutions des effets environnementaux et SNPparente pour le modele avec tous les SNPpar SANS VAR_GWAS
  real(kind=8),allocatable::BetaSnew(:,:) ! ! vecteur qui va contenir les solutions des effets environnementaux et SNPparente pour le modele avec tous les SNPpar AVEC VAR_GWAS
  real(kind=8)::BASE1   ! scalaire  = Y'WY - b'XW'Y 
  real(kind=8)::BASE2   ! scalaire  = Y'WY - b'XW'Y 
  real(kind=8)::exp4_1, exp4_2, invdenom,invdenom_1,invdenom_2, sumPoids, contrib_VarGWAS
  real(kind=8),allocatable::restemp(:,:) ! resultat de calcul matriciel temporaire
  integer::rank,rank_1,rank_2    ! rang de la partie effets fixes du modele pour calcul de la statistique de test
  integer::presmoy,nb_eff_fixcat,nb_eff_cov
  real(kind=8),allocatable::mat_TestEff_1(:,:),mat_TestEff_2(:,:)  ! matrice(2,nbVarGWAS) contenant test de effet additif et test de effet dominance de chaque variant GWAS avec VarRes approx et exact
  real(kind=8),allocatable::pval(:)
  real(kind=8),allocatable::mat_pval_1(:,:),mat_pval_2(:,:)

  !real(kind=8)::G01ECF
  !EXTERNAL  G01ECF
  real(kind=8)::cdf


  real(kind=8),allocatable::Test_effVarGWAS_1(:),Test_effVarGWAS_2(:)    ! test de l effet additif et le test de l effet de dominance du variant GWAS en cours avec VarRes approx et exact
  real(kind=8),allocatable::vec_VarRes(:,:)     ! vecteur contenant pour chaque variant la variance residuelle qui a ete calculee pour lui pour le calcul de son test
  real(kind=8)::seuil_test1=4.0d0   ! seuil de valeur de Test d effet de Variant GWAS au dessus duquel on calcule la variance residuelle de maniere exacte pour le variant dans CalcVarRes=optim


  ! declarations pour PARTIE 3
  integer::nbGWASgr,nbGWASexec,nbGWAStest
  integer::nbGWAS_limite=0
  integer::numGR,NumSGr
  integer,allocatable::NumVarGWASori(:),GWASuttoMAP_P3(:),SegmElim_p3(:,:)
  integer::chrom_p3,numgr_p3,numSgr_p3,Currline
  integer::StepFichOut=500  ! nombre de Variants GWAS en partie 3 apres lequel on sort les resultats de GWAS dans le fichier txt


  ! declarations pour lecture fichiers variants GWAS format plink
  integer::NUMCHR_GWAS ! numero du chromosome portant les variats GWAS pour la partie 2 et la partie 3 si fichiers typages et carte variants GWAS au format plink
  integer::firstVarGWAS,lastVarGWAS ! numeros du 1er et du dernier variant GWAS du chromosome etudie en etapes 2 et 3 a conserver dans le fichier typages .raw format plink
  integer::pospb,poscm
  integer,allocatable::vectemp5(:)
  character*lc :: sire_temp,dam_temp,sex_temp,pheno_temp
  character(len=20)::nomvar
  integer::iprec,nltot

  real(kind=8), parameter :: dx = 1.0d-5  ! Pas de l'intégration
  real(kind=8)::pval_ttr ! P-value calculated with normal_ttr function

  integer,allocatable::compttyp(:,:) ! nbre de genotypes 11 12 et 22 pour chaque variant pour determiner si on a une structure de donnees adaptee a l estimation des effets add et dom avec typages discrets


  ! on verifie que la date limite d execution n est pas depassee
  call fdate(jour)
  call vversion(jour,anlim,mlim,jlim)  


  call fdate(jour)

  print*,'****************************************'
  print*,'*                                      *'
  print*,'*           SNP_GWAS_publi             *' 
  print*,'*    Version 28 jan 2026 LOG_english   *' 
  print*,'*                                      *'
  print*,'* This version uses a home made        *'
  print*,'* function to compute pvalues          *'
  print*,'*                                      *'
  print*,'   ',jour
  print*,'*                                      *'
  print*,'****************************************'
  print*,' '
  print*,' '



  ! on lit le fichier parametres
  call read_parameters_hm

  ! change defaults if optional parameters present and check positive definitness
  call defaults_and_checks()

  call print_parameters_hm


  if(steptodo.ne.1) then ! section executed if STEP_2 or STEP3 because uses information not mendatory in STEP_1

     if(gwastype.eq.'add') neff_GWAS_lu=1
     if(gwastype.eq.'add_dom') then
        neff_GWAS_lu=2
        TemDom_lu=1
     endif
     print*,' '
     print*,'gwastype_read = ',gwastype,' --> neff_GWAS_read = ',neff_GWAS_lu,'  TemDom_read =',TemDom_lu
     print*,' '


     if((gwastype.ne.'add').and.(gwastype.ne.'add_dom')) then
        print*,'GWAS_TYPE read in parameter file is incorrect : value should be in  add  add_dom'
        STOP 10
     endif

  endif ! end of section executed if STEP_2 or STEP3

  ! on calcule l inverse de la variance residuelle
  ! POUR L INSTANT CAS MONOCARACTERE UNIQUEMENT
  rinv=0.0d0
  rinv(1,1)=1.0d0/r(1,1)

  ! on calcule l inverse de la variance genetique additive
  ! POUR L INSTANT CAS MONOCARACTERE UNIQUEMENT
  ginv(1,1)=1.0d0/g(1,1)

  ! on calcule l inverse de la variance pour les effets aleatoires environnementaux diagonaux
  randinv=0.0d0
  RhoDiag=0.0d0
  do ef=1,neff
     if(randomtype(ef).eq.g_diag) then
        randinv(ef,1,1)=1.0d0/rand(ef,1,1)
        RhoDiag(ef,1,1)=R(1,1)/rand(ef,1,1)
     endif
  enddo


  !data_len=max(maxval(pos_weight),maxval(pos_y),maxval(pos_eff),maxval(nestedcov))
  data_len=posIDchar-1
  if(max(maxval(pos_weight),maxval(pos_y),maxval(pos_eff),maxval(nestedcov)).gt.posIDchar) then
     print*,'Phenotype or weight or effect level is after (at the right of) alphanumeric Id of individual'
     print*,'The Id of individuals in Phenotype file has to be AFTER (on the right) all numeric variables read by the program'
     STOP 21
  endif
  allocate(indata(data_len),y(ntrait),weight_y(ntrait),weight_cov(neff,ntrait),nbcarperf(ntrait),nivanim(neff,ntrait))

  print*,' '
  print*,'data_len=',data_len


  call mkl_set_num_threads(nb_threads_mkl) ! Pas sur que ca soit utile : c est ce ue jai mis dans HSSGBLUP pour PARDISO mais ils n en parlent pas dans la doc pour BLAS ni LAPACK ...




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(io_d,file=datafile,status='old')

  ! on lit le fichier PERF une 1ere fois pour COMPTER LE NOMBRE DE LIGNES = NBIND utilise pour ALLOUER LES TABLES
  nbind=0
  do
     read(io_d,*,iostat=ioperf)indata,NoNatTemp
     if (ioperf .ne. 0) exit
     nbind=nbind+1
  enddo

  print*,' '
  print*,'First reading of performance file is complete: nbind = number of lines in the performance file = ',nbind

  rewind(io_d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! allocation des tables
  allocate(matNivanim(neff,ntrait,nbind),matWeight_cov(neff,nbind),matPERF(nbind,ntrait),matPOIDS(nbind,ntrait),animperf(nbind),vecPERFpond(nbind))
  allocate(vecPERFpond2(nbind))
  animperf=' '
  matNivanim=0
  matWeight_cov=0.0d0
  matPERF=0.0d0
  vecPERFpond=0.0d0
  vecPERFpond2=0.0d0
  matPOIDS=0.0d0
  nbcarperf=0
  ioperf=0


  ! 2eme lecture du fichier PERF pour remplir le vecteur des ID ANIM qu on  va trier par ordre croissant
  nlperf=0
  do
     read(io_d,*,iostat=ioperf)indata,NoNatTemp
     if (ioperf .ne. 0) exit
     nlperf=nlperf+1
     animperf(nlperf)=NoNatTemp
     !print*,'nbind=',nbind,' /  nlperf=',nlperf
  enddo

  if(detailed_log.eq.3) print*,'Second reading of performance file is complete: nlperf = number of lines in the performance file = ',nlperf

  rewind(io_d)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! on lit le fichier performances et on stocke l Id de l individu sa perf son poids et les niveaux d effets et covariables

  ! on trie les ID des individus du fichier PERF et on cree un vecteur d ordre des individus dans le fichier PERF
  allocate(ordperf(nlperf))
  ordperf=0

  if(nbind.lt.20) then
     print*,'animperf vector in the order in which the performance file is read, BEFORE sorting with HPSORT=',animperf(:)
  endif

  workperf=''
  call hpsort(animperf,nlperf,1,lc,ordperf,2,workperf,ifail)
  if(detailed_log.eq.3) print*,'sorting with HPSORT is complete'

  if(nbind.lt.20) then
     print*,'animperf vector AFTER sorting with HPSORT=',animperf(:)
     print*,'ordperf vector =',ordperf(:)
  endif

  ! We test to check there is no individual with more than 1 phenotype = individuals present several times in DATA file
  ! we work on the sorted list of individuals
  do i=2,nlperf
     if (animperf(i).eq.animperf(i-1)) then
        print*,'At least 1 individual appears several times in DATA file:',animperf(i)
        print*,'Repeated phenotypes per individual is not allowed'
        STOP 57
     endif
  enddo


  ! vecteur ordperf : element i = position dans fichier performances lu du ieme animal trie
  ! --> on cree vecteur ordinvperf : element i = classement de l animal de la ieme ligne du fichier perf lu apres tri sur IdAnim
  ! necessaire pour remplir matrices perf effets ...
  allocate(ordinvperf(nlperf))
  ordinvperf=0
  do i=1,nlperf
     ordinvperf(ordperf(i))=i
  enddo

  if(nbind.lt.20) then
     print*,' '
     print*,'ordinvperf vector=',ordinvperf(:)
     print*,' '
  endif

  ! on relit le fichier PERF une 2eme fois et on va remplir les tables de perf poids niveaux_effets et covariables selon l ordre croissant d ID ANIM
  nlperf=0
  do
     read(io_d,*,iostat=ioperf)indata,NoNatTemp
     if (ioperf .ne. 0) exit
     nlperf=nlperf+1
     ! le numero d ordre aniperf du nlperf_eme individu NoNatTemp du fichier perf dans la liste reordonnee est ordinvperf(nlperf)
     aniperf=ordinvperf(nlperf)

     call decode_record_effix ! pour performances et effets du modele d evaluation genomique

     do j=1,ntrait
        matPERF(aniperf,j)=y(j)
        if(y(j).ne.mis) then
           nbcarperf(j)=nbcarperf(j)+1
        endif
     enddo
     do j=1,ntrait
        matPOIDS(aniperf,j)=weight_y(j)
        if(y(j).eq.mis) matPOIDS(aniperf,j)=0.0d0
        vecPERFpond(aniperf)=matPERF(aniperf,j)*sqrt(matPOIDS(aniperf,j))
        vecPERFpond2(aniperf)=matPERF(aniperf,j)*matPOIDS(aniperf,j)
     enddo
     do i=1,neff
        do j=1,ntrait
           if(weight_cov(i,j).ne.0.0d0) matWeight_Cov(i,aniperf)=weight_cov(i,j) ! si effet dans modele plusieurs caract, weight_cov(i,j) a la meme valeur pour tous ces caract 
           matNivAnim(i,j,aniperf)=Nivanim(i,j)
        enddo
     enddo

  enddo

  print*,' '
  print*,'Number of NON MISSING PHENOTYPES for the trait = nbcarperf(1)=',nbcarperf(1)
  print*,' '


  ! on imprime les tables creees a partir de DATA pour verifier la lecture de DATA
  if(detailed_log.ge.2) then
     print*,' '
     print*,'For the first 10 individuals:'
     print*,'numeric Id / matNivanim / matWeight_cov / matPERF / matPOIDS / vecPERFpond / vecPERFpond2 / animperf'
     do i=1,min(10,nbind)
        print*,i,' / ',matNivanim(:,1,i),' / ',matWeight_cov(:,i),' / ',matPERF(i,1),' // ',matPOIDS(i,1),' / ',vecPERFpond(i),' / ',vecPERFpond2(i),' // ',animperf(i)
     enddo
     print*,' '
  endif

  !if(nlperf.lt.nbind) then
  !   print*,'le nb de lignes dans le fichier DATA est inferieur au nb d individus declares dans fich PARAM'
  !   STOP 21
  !endif

  close(io_d)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! modif du 06/05/2024 pour permettre d avoir des individus supplementaires dans fichier TYPAGES PARENTE par rapport a individus du fichier PERFORMANCES



  ! FICHIER TYPAGES VARIANTS PARENTE ! 

  ! modification du 06/05/2024 pour permettre d avoir des individus supplementaires dans fichier TYPAGES VARIANTS par rapport a fichier perf
  ! interet : pourvoir utiliser les memes fichier PLINK pour typages variants GWAS que pour typages variants parente


  ! on lit le fichier TYPAGES PARENTE une 1ere fois pour compter le nb de lignes pour pouvoir dimensionner les tables

  open (io_typpar,file=trim(typparfile),form='FORMATTED')     

  nltyppar=0
  io=0

  if(trim(fmttyppar).eq.'plink') then

     ! on lit une premiere fois le fichier CARTE DES VARIANTS PARENTE POUR CONNAITRE LE NB DE SNP dans le fichier TYPAGES DES VARIANTS PARENTE

     tt=0
     iomappar=0
     open(io_mappar,file=mapparfile,status='old')
     do
        read(io_mappar,*,iostat=iomappar) 
        if (iomappar .ne. 0) exit
        tt=tt+1
     enddo
     print*,' '
     print*,'Number of lines read in the PLINK-format MARKER map file (=number of markers on all chromosomes): ',tt
     print*,' '
     close(io_mappar)

     ! on lit maintenant le fichier des typages pour les VARIANTS PARENTE

     read(io_typpar,*,iostat=io) ! on lit la 1ere ligne d en tete qui ne sert a rien et qui est en alphanumerique
     if (io.lt.0) then
        print*,' '
        print*,'Marker Genotype file is empty'
        stop 30
     endif

     do
        read(io_typpar,*,iostat=io) NoNatTemp2 , NoNatTemp ! dans fichier .raw de plink l identifiant de l individu est en 2eme champ
        if (io .ne. 0) exit
        nltyppar=nltyppar+1
     enddo
     rewind(io_typpar)

     nbSNP_0=tt   ! on dit directement que le nb de SNP dans le fichier TYPAGES PARENTE est le nb de lignes du fichier CARTE VARIANTS PARENTE
  endif



  if(trim(fmttyppar).eq.'typ_eval') then
     ! on lit la 1ere ligne du fichier pour repérer la position du 1er typage et creer le format de lecture pour la suite
     read(io_typpar,'(a3000000)',iostat=io) linetyp
     if (io.lt.0) then
        print*,'Marker Genotype file is empty'
        stop 30
     endif
     len_linetyp=len_trim(linetyp)
     d=len_linetyp
     do while(linetyp(d:d).ne.' ')
        d=d-1
     enddo
     ip_snp=d+1
     nbSNP_0=len_linetyp-d

     write(informatP1,'(a,i0,a,i0,a)') '(a',ip_snp-2,',1x,',nbSNP_0,'i1)'

     print*,'informatP1= Reading format for Marker Genotype file in typ_eval format',informatP1
     print*,'Number of Markers read on the first line in the Marker genotype file =',nbSNP_0

     rewind(io_typpar)
     nltyppar=0
     io=0
     ! on lit le fichier 
     do
        read(io_typpar,*,iostat=io) NoNatTemp
        if (io .ne. 0) exit
        nltyppar=nltyppar+1
     enddo
     rewind(io_typpar)
  endif

  print*,' '
  print*,'Number of lines read in the Marker genotype file (=number of individuals with genotypes) =',nltyppar

  if (nltyppar.ne.nbind) then
     print*,' '
     print*,'INFO: number of lines read in Marker genotype file ',nltyppar,' different from number of individuals in DATA file ',nbind
     print*,'--> individuals present in Marker Genotype file but missing in phenotype file will be removed'
     print*,' '
     !STOP 30
  endif

  ! on relit le fichier TYPAGES une 2eme fois pour etablir la liste des individus avec TYPAGES PARENTE, verifier qu aucun individu avec PERF ne manque 
  ! et reperer les individus inutiles (=indiv sans performance)

  allocate(animPar(nltyppar))

  animPar=' '
  nltyppar=0
  io=0

  if(trim(fmttyppar).eq.'plink') then
     read(io_typpar,*,iostat=io) ! on lit la 1ere ligne d en tete qui ne sert a rien et qui est en alphanumerique

     do
        read(io_typpar,*,iostat=io) NoNatTemp2 , NoNatTemp ! dans fichier .raw de plink l identifiant de l individu est en 2eme champ
        if (io .ne. 0) exit
        nltyppar=nltyppar+1
        animPar(nltyppar)=NoNatTemp
     enddo
     rewind(io_typpar)
  endif

  if(trim(fmttyppar).eq.'typ_eval') then
     io=0
     ! on lit le fichier 
     do
        read(io_typpar,*,iostat=io) NoNatTemp
        if (io .ne. 0) exit
        nltyppar=nltyppar+1
        animPar(nltyppar)=NoNatTemp
     enddo
     rewind(io_typpar)
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! on trie les ID des individus du fichier TYPAGES PARENTE et on cree un vecteur d ordre des individus dans le fichier TYPAGES PARENTE
  allocate(ordTypPar(nltyppar))
  ordTypPar=0

  if(nltyppar.lt.20) then
     print*,'animPar vector containing the Id read in Marker Genotype file BEFORE sorting: ',animPar(:)
  endif

  workTypPar=''
  call hpsort(animPar,nltyppar,1,lc,ordTypPar,2,workTypPar,ifail)
  if(detailed_log.eq.3) print*,'sorting with HPSORT is complete'

  if((detailed_log.eq.3).and.(nltyppar.lt.20)) then
     print*,'animPar vector containing the Id read in Marker Genotype file AFTER sorting: ',animPar(:)
  endif


  ! vecteur ordTypPar : element i = position dans fichier TYPAGES PARENTE lu du ieme animal trie
  ! --> on cree vecteur ordinvTypPar : element i = classement de l animal de la ieme ligne du fichier TYPAGES PARENTE lu apres tri sur IdAnim
  ! necessaire pour remplir matrices TYPAGES PARENTE dans ordre coherent avec tables PERF EFFETS ...
  allocate(ordinvTypPar(nltyppar))
  ordinvTypPar=0
  do i=1,nltyppar
     ordinvTypPar(ordTypPar(i))=i
  enddo

  if(detailed_log.eq.3) then
     !if(nltyppar.lt.20) then 
     print*,'First 20 lines of sorted animPar vector = ',animPar(1:min(20,nltyppar))
     print*,'First 20 lines of sorted animperf vector = ',animperf(1:min(20,nbind))
     !endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! on compare animPar et animperf pour reperer les individus du fichier TYPAGES PARENTE inutiles dont on va supprimer les typages pour les GWAS 
  ! les 2 vecteurs sont tries donc on peut avancer pas a pas dans les 2
  allocate(TypParUtil(nltyppar),IndParUt(nbind))
  TypParUtil=0
  IndParUt=0
  j=1 !indice de la ligne dans animperf()
  do i=1,nltyppar ! indice de la ligne dans animPar()
     if(animPar(i).eq.animperf(j)) then
        TypParUtil(i)=j
        IndParUt(j)=i
        j=j+1
        if(j.gt.nbind) exit
     endif
  enddo

  print*,'After comparison of animPar and animperf, j=',(j-1)
  ! correction (j-1) because j+1 at the end of the loop

  if(j.lt.nbind) then
     print*,' ' 
     print*,'ERROR : some individuals present in Phenotype file are missing in Marker Genotype file'
     print*,'ALL individuals present in Performance file MUST BE in Marker Genotype file'
     STOP 33
  endif



  allocate(genopartemp(nbSNP_0))
  genopartemp=0

  ! on relit le fichier typages parente pour calculer les frequences alleles 2 pour calculer la MAF des SNPpar
  allocate(MAF(nbSNP_0),nballeles(nbSNP_0),freq(nbSNP_0))

  nbtyp=nltyppar

  ! on calcule la frequence de l allele 2 pour chaque SNP a partir du fichier de typages

  !sumpq_all=0d0 ! modif 19avr2024
  !sumpq=0d0 ! modif 19avr2024
  !freq=0d0 ! modif 19avr2024
  !freq2=0d0 ! modif 19avr2024
  ! on lit le fichier typages ligne par ligne et on incremente le nb d alleles 2 au fur et a mesure

  nblues=0
  nballeles=0d0
  io=0
  MAF=0.0d0
  freq=0.0d0

  if (fmttyppar.eq.'plink') then
     read(io_typpar,*,iostat=io) ! on lit la ligne d en-tetes qui ne sert a rien
  endif


  do
     if (fmttyppar.eq.'typ_eval') then
        read(io_typpar,informatP1,iostat=io)NoNatTemp,genopartemp(:)
        if(io/=0) exit
     endif

     if (fmttyppar.eq.'plink') then
        read(io_typpar,*,iostat=io) NoNatTemp2,NoNatTemp,SIRE_temp,DAM_temp,SEX_temp,PHENO_temp,genopartemp(:)
        if(io/=0) exit
     endif

     nblues=nblues+1 ! numro de la ligne du fichier TYPAGES VARIANTS PARENTE en cours de lecture

     ! on ne considere que les typages des individus utiles = individus types avec performance
     if(TypParUtil(nblues).ne.0) then   ! TypParUtil(k) = numero de lige dans la table animPerf() de l individu a la ligne k du fichier TYPAGES PARENTE
        do i=1,nbSNP_0
           if (genopartemp(i).eq.1) freq(i)=freq(i)+1d0
           if (genopartemp(i).eq.2) freq(i)=freq(i)+2d0
           if (genopartemp(i).ne.5) nballeles(i)=nballeles(i)+2d0 ! genotype=5 = code typage manquant --> IL FAUT AU PREALABLE AVOIR REMPLACE LES NA en 5
        enddo
     endif
  enddo

  do i=1,nbSNP_0
     freq(i)=freq(i)/nballeles(i)
  enddo

  allocate(mostFreqTyp(nbSNP_0))

  ! on remplit le vecteur mostFreqTyp avec le genotype observe le plus frequent pour chaque SNP
  mostFreqTyp=1 ! on initialise le genoptype observe le plus frequent a heterozygote
  do i=1,nbSNP_0
     if(freq(i).lt.(1.0d0/3.0d0)) mostFreqTyp(i)=0  ! si freq(allele2) < 0.3333 gzenotype le plus frequent = homozygote 11 
     if(freq(i).gt.(2.0d0/3.0d0)) mostFreqTyp(i)=2  ! si freq(allele2) > 0.6666 gzenotype le plus frequent = homozygote 22 
  enddo

  ! on remplit le vecteur freq2 des frequences utiles d alleles 2
  do i=1,nbSNP_0
     MAF(i)=freq(i)
     if(MAF(i).gt.0.5d0) MAF(i)=1.0d0-MAF(i)
  enddo

  print '('' Number of heterozygous Markers: '',t110,i8)',count( freq*(1d0-freq) /=0.0d0 )

  print '('' Reading of Marker Genotype file to calculate allele frequencies completed - nb of lines read :'',t110,i8)',nblues

  close(io_typpar)



  !print*,' '
  !print*,'MAF of the first 100 Markers :'
  !print*,MAF(1:100)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! A ce moment du programme, on connait :
  !         - le nb de variants parente a lire dans le fichier typages des variants parente
  !         - les individus du fichier TYPAGES PARENTE qui sont UTILES = qui sont presents dans le fichier PERFORMANCES
  !         - la MAF de chacun des VARIANTS PARENTE calculee en considerant uniquement les typages des individus utiles



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! modif du 17 fev 2023 pour le choix automatique des SNPparente si utilisateur utilise OPTION SelAutoSNPpar nbSNPparAuto
  ! principe : 1. on lit une premiere fois le fichier CARTE SNPparente pour connaitre les positions des SNPpar sur les chromosomes
  !            et lister les variants exclus volontairement pa l utilisateur
  !            2. on lit les typages SNPpar pour calculer les frequences et eliminer les SNPpar dont freq < MAFmin
  !            3. on choisit les SNPpar a retenir pour les analyses et on cree un fichier CARTE_bis en indiquant d exclure tous les autres SNPpar
  !            sur lequel on va executer la suite du programme ainsi le reste du programme ne change pas

  if( (SelAutoSNPpar.eq.1) .and. (steptodo.eq.1) ) then ! on n effectue ceci que si on est dans STEP1 car en STEP 2 et 3 on lit fich CARTE_PAR alternatif

     iomappar=0
     open(io_mappar,file=mapparfile,status='old')

     MAPPARtemp=0
     nlmappar=0
     chrPrec=0
     posPrec=0
     numPrec=0
     tt=0
     nb_champs=0

     if(fmttyppar.eq.'plink') then
        ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
        !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl

        call nb_chps_ttr(io_mappar,nb_champs)

        rewind(io_mappar)

        print*,' '
        print*,'Number of columns read in the Marker MAP file = ',nb_champs
        if(nb_champs.eq.4) print*,'no INCL/EXCL indicator in Marker MAP file'
        if(nb_champs.eq.5) print*,'INCL/EXCL indicator is present in Marker MAP file'
        print*,' '
     endif

     do
        if(fmttyppar.eq.'typ_eval') then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
           if (iomappar .ne. 0) exit
        endif

        if(fmttyppar.eq.'plink') then

           if(nb_champs.eq.4) then
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
           endif
           if(nb_champs.eq.5) then
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
                 print*,'INCL/EXCL indicator in Marker MAP file for Marker ',tt,'=',MAPPARtemp(4),' different from 0 and 1'
                 STOP 21
              endif
           endif
        endif

        if (MAPPARtemp(1).lt.chrPrec) then
           print*,'The Marker Map file is not sorted by ascending chromosome Id'
           STOP 10
        endif

        if ((MAPPARtemp(1).eq.chrPrec).and.(MAPPARtemp(2).le.posPrec)) then
           print*,'The Marker Map file is not sorted by ascending Marker position within chromosome : chr',chrPrec,posPrec,MAPPARtemp(2)
           STOP 10
        endif

        if (MAPPARtemp(3).le.numPrec) then
           print*,'The Marker Map file is not sorted by ascending Marker number :',numPrec,MAPPARtemp(3)
           STOP 10
        endif

        nlmappar=nlmappar+1
        if (MAPPARtemp(3).ne.nlmappar) then
           print*,'The number sequence of Markers from 1 to N is incorrect on line:',nlmappar
           STOP 10
        endif

        if(nlmappar.eq.1) then ! on est au 1er marqueur du 1er chromosome --> on initialise les compteurs etc...
           chrPrec=MAPPARtemp(1)
           posPrec=MAPPARtemp(2)
           numPrec=MAPPARtemp(3)
        endif

        if (MAPPARtemp(1).gt.chrPrec) then   ! on est passe a un nouveau chromosome
           chrprec=MAPPARtemp(1)
           posPrec=0
        endif

        if (MAPPARtemp(3).ne.nlmappar) then
           print*,'The number sequence of Markers from 1 to N is discontinuous :',nlmappar,MAPPARtemp(3)
           STOP 10
        endif

     enddo

     ! Actual total number of chromosomes is deduced from the last value of MAPPARtemp(1)
     nbCHROMpar=MAPPARtemp(1)

     print*,'Number of chromosomes present in the Marker map file =',nbCHROMpar
     allocate(LongChromPar(nbCHROMpar),nbSegmPar(nbCHROMpar),posPremLastSNPpar(2,nbCHROMpar))
     LongChromPar=0.0d0
     nbSegmPar=0
     posPremLastSNPpar=0

     print*,'Number of lines read in the Marker MAP file =',nlmappar


     rewind(io_mappar)

     ! on relit le fichier MAPPAR une 2eme fois pour remplir le vecteur INCL_PAR
     ! on sait que le fichier MAPPAR est trie par ordre de chromosomes croissant et par ordre de SNPpar croissant intra chromosome --> 2eme lecture simple
     allocate(INCL_PAR(nlmappar),posSNPpar(nlmappar))
     INCL_PAR=0
     posSNPpar=0
     i=0
     tt=0
     nb_champs=0


     if(fmttyppar.eq.'plink') then
        ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
        !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
        call nb_chps_ttr(io_mappar,nb_champs)
        rewind(io_mappar)
     endif

     do
        if(fmttyppar.eq.'typ_eval') then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
           if (iomappar .ne. 0) exit
        endif

        if(fmttyppar.eq.'plink') then
           if(nb_champs.eq.4) then
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
           endif
           if(nb_champs.eq.5) then
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
                 print*,'INCL/EXCL indicator for the MARKER ',tt,'=',MAPPARtemp(4),' is different from 0 and 1'
                 STOP 21
              endif
           endif
        endif

        i=i+1
        if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
           print*,'INCL/EXCL indicator different from 0 and 1 for Marker ',i
           STOP 10
        endif
        INCL_PAR(i)=MAPPARtemp(4)
        posSNPpar(i)=MAPPARtemp(2)

        if(posPremLastSNPpar(1,MAPPARtemp(1)).eq.0) posPremLastSNPpar(1,MAPPARtemp(1))=MAPPARtemp(2)
        if(MAPPARtemp(2).gt.posPremLastSNPpar(2,MAPPARtemp(1))) posPremLastSNPpar(2,MAPPARtemp(1))=MAPPARtemp(2)
     enddo

     print*,'INCL/EXCL indicator read for ',i,' Markers'

     rewind(io_mappar)

     nbSNPpar_consUt=sum(INCL_PAR(:))
     print*,'The number of markers the user has asked to be kept with INCL/ECL indicator is :', nbSNPpar_consUt
     if(nbSNPpar_consUt.lt.nbSNPparAuto) then
        print*,' '
        print*,'PROBLEM : User asked ',nbSNPparAuto,' Markers to be automatically selected but only ',nbSNPpar_consUt,'Mrkers variants asked to be kept in Marker MAP file'
        STOP 13
     endif

     print*,' '
     do i=1,nbCHROMpar
        !LongChromPar(i)=(posPremLastSNPpar(2,i)-posPremLastSNPpar(1,i))*0.000001d0
        LongChromPar(i)=posPremLastSNPpar(2,i)*0.000001d0   ! on considere que la longueur du chromosome est la position de son dernier SNPparente
        print*,'Length of chrom',i,' = ',LongChromPar(i),' / Position of 1st and last Markers on the chrom =',posPremLastSNPpar(:,i)
     enddo


     LongTotGenomePar=sum(LongChromPar(:))
     print*,' '
     print*,'Total Length of Genome (Mb) base on marker positions = ',LongTotGenomePar

     LongSegmPar=INT( ( LongTotGenomePar / (dfloat(nbSNPparAuto)*0.000001d0) ) ) ! longueur de chaque segment de chromosome pour avoir nbSNPparAuto SNPpar répartis sur le genome
     print*,'Length (bp) of elementary Marker Segments resulting from numbers of automatically selected Markers asked by user =',LongSegmPar

     do i=1,nbCHROMpar
        nbSegmPar(i)=INT( LongChromPar(i) / (dfloat(LongSegmPar)*0.000001d0) ) + 1 ! des segments de taille LongSegmPar + un segment pour finir le chromosome
     enddo
     print*,' '
     print*,'Number of segments per chromosom on which the most central marker with highest MAF will be kept:';
     do i=1,nbCHROMpar
        print*,'Chrom ',i,' Number of segments = ',nbSegmPar(i)
     enddo
     print*,' '


     allocate(InfoSegmPar(3,maxval(nbSegmPar(:)),nbCHROMpar))
     InfoSegmPar=0
     tt=0
     nb_champs=0

     if(fmttyppar.eq.'plink') then
        ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
        !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
        call nb_chps_ttr(io_mappar,nb_champs)
        rewind(io_mappar)
     endif

     ! on relit le fichier CARTE SNPpar pour recherche le numero du 1er et du dernier SNPpar non exclus par l utilisateur de chaque segment de chaque chromosome
     do
        if(fmttyppar.eq.'typ_eval') then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
           if (iomappar .ne. 0) exit
        endif

        if(fmttyppar.eq.'plink') then
           if(nb_champs.eq.4) then
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
           endif
           if(nb_champs.eq.5) then
              !read(io_mappar,*,iostat=iomappar) i,nomvar,poscM,pospb
              read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
              if (iomappar .ne. 0) exit
              tt=tt+1
              MAPPARtemp(3)=tt
              if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
                 print*,'INCL/EXCL indicator for Marker ',tt,'=',MAPPARtemp(4),' different from 0 and 1'
                 STOP 21
              endif
           endif
        endif

        i=i+1
        if(MAPPARtemp(4).eq.1) then   ! SNPpar non exclu par utilisateur

           numSegmPar = INT(dfloat(MAPPARtemp(2))/dfloat(LongSegmPar))+1 ! numero du segment sur son chromosome auquel le SNPpar appartient

           if (InfoSegmPar(1,numSegmPar,MAPPARtemp(1)).eq.0) InfoSegmPar(1,numSegmPar,MAPPARtemp(1))=MAPPARtemp(3) ! premier SNPpar du segment --> on le note
           if (InfoSegmPar(2,numSegmPar,MAPPARtemp(1)).lt.MAPPARtemp(3)) InfoSegmPar(2,numSegmPar,MAPPARtemp(1))=MAPPARtemp(3)

        endif

     enddo

     if(detailed_log.eq.3) print*,'Table InfoSegmPar(:,:,:) has been built and filled'

     rewind(io_mappar)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     do bta=1,nbCHROMpar
        do numSegmPar=1,nbSegmPar(bta)
           MafMaxSegm=0.0d0
           DistBest=LongSegmPar ! on initialise la meilleur distance a la longueur du segment pour que le 1er SNPpar teste soit le meilleur en distance
           if((InfoSegmPar(1,numSegmPar,bta).ne.0).and.(InfoSegmPar(2,numSegmPar,bta).ne.0)) then
              !MedPos = dfloat(posSNPpar(InfoSegmPar(2,numSegmPar,bta))) - dfloat(LongSegmPar)*0.5d0 ! position du milieu du segment NE FONCTIONNE PAS POUR SEGMENTS PLUS COURTS EN FIN DE CHROM
              MedPos = dfloat(posSNPpar(InfoSegmPar(2,numSegmPar,bta))) - dfloat( (posSNPpar(InfoSegmPar(2,numSegmPar,bta))-posSNPpar(InfoSegmPar(1,numSegmPar,bta)) ) )*0.5d0 ! position du milieu du segment
              BestSNPpar=InfoSegmPar(1,numSegmPar,bta)
              do SNPpar=InfoSegmPar(1,numSegmPar,bta),InfoSegmPar(2,numSegmPar,bta)
                 if(INCL_PAR(SNPpar).eq.1) then
                    INCL_PAR(SNPpar)=0 ! on re-initialise INCL-EXCL du SNPpar a 0 --> on le mettra a 1 si c est le meilleur variant du Segment
                    if((MAF(SNPpar).ge.MAFmin).and.(MAF(SNPpar).eq.MafMaxSegm)) then
                       if(ABS(dfloat(posSNPpar(SNPpar))-MedPos).lt.DistBest) then
                          DistBest=ABS(dfloat(posSNPpar(SNPpar))-MedPos)
                          InfoSegmPar(3,numSegmPar,bta)=SNPpar  ! meilleur SNPpar du segment
                       endif
                    endif
                    if((MAF(SNPpar).ge.MAFmin).and.(MAF(SNPpar).gt.MafMaxSegm)) then
                       MafMaxSegm=MAF(SNPpar)
                       DistBest=ABS(dfloat(posSNPpar(SNPpar))-MedPos)
                       InfoSegmPar(3,numSegmPar,bta)=SNPpar   ! meilleur SNPpar du segment
                    endif
                 endif
              enddo
           endif
           if(InfoSegmPar(3,numSegmPar,bta).ne.0) INCL_PAR(InfoSegmPar(3,numSegmPar,bta))=1 ! on actualise INCL-EXCL du meilleur SNPpar du segment
        enddo
     enddo


     if(detailed_log.eq.3) print*,'The best Marker on each elementary Segment has been identified'

     nbSNPpar_retenus=0
     do bta=1,nbCHROMpar
        do numSegmPar=1,nbSegmPar(bta)
           if(InfoSegmPar(3,numSegmPar,bta).ne.0) nbSNPpar_retenus=nbSNPpar_retenus+1
        enddo
     enddo

     print*,'Final number of best Markers kept on the complete segment set :',nbSNPpar_retenus


     ! on ecrit sur disque un fichier CARTE_SNPparente alternatif dans lequel seuls les SNPpar selectionnes par la procedure AUTOMATIQUE sont conserves
     ! ce fichier CARTE ALTERNATIF sera lu et utilise dans la suite du programme en etapes 1, etape 2 et Etape 3 si OPTION SelAutoSNPpar est demandee par utilisateur
     !open(io_mapparout,file=trim(mapparfile)//trim("_bis"),form="formatted") ! 20/03/2023 : on modifie le fichier CARTE_SNP_PAR_bis pour qu il soit dans repertoire RACINE
     !                                                                                       pour ne pas etre ecrase si traitement autre caractere
     open(io_mapparout,file=trim("carte_SNP_PAR_bis.txt"),form="formatted")

     write(formatPAR,'(a,i0,a)') '(i2,1x,i10,1x,i8,1x,i1)'

     i=0
     if(fmttyppar.eq.'typ_eval') then
        do
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
           if (iomappar .ne. 0) exit
           i=i+1
           !if(INCL_PAR(i).ne.2) INCL_PAR(i)=0
           !if(INCL_PAR(i).eq.2) INCL_PAR(i)=1
           write(io_mapparout,formatPAR) MAPPARtemp(1:3),INCL_PAR(i)
        enddo
     endif
     if(fmttyppar.eq.'plink') then
        do
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
           if (iomappar .ne. 0) exit
           i=i+1
           MAPPARtemp(3)=i
           ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
           write(io_mapparout,formatPAR) MAPPARtemp(1:3),INCL_PAR(i)
        enddo
     endif

     print*,' '
     !print*,i,'SNPparente ecrits dans fichier ',trim(mapparfile)//trim("_bis")
     if(detailed_log.eq.3) print*,i,'Markers written in new updated marker MAP file ',trim("carte_SNP_PAR_bis.txt")
     print*,' '

     !deallocate(INCL_PAR,nballeles,MAF,genopartemp)
     deallocate(INCL_PAR,nballeles,MAF)

     close(io_mappar)
     close(io_mapparout)
     !close(io_typpar)

  endif ! fin du test sur if( (SelAutoSNPpar.eq.1) .and. (steptodo.eq.1)) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! on lit le fichier carte pour les marqueurs parente une 1ere fois pour creer le vecteur INCL_PAR disant pour chaque SNPparente si on l exclut ou inclut a priori
  ! 4 colonnes : num_chrom / position_sur_chrom / numero_marqueur_dans_fichier_typages_parente / inclus(1) ou exclus(0) des equations parente independament de sa MAF

  if(SelAutoSNPpar.eq.0) then
     open(io_mappar,file=mapparfile,status='old')
  endif

  if(SelAutoSNPpar.eq.1) then
     !mapparfile_bis=trim(mapparfile)//trim("_bis")
     mapparfile_bis=trim("carte_SNP_PAR_bis.txt")
     open(io_mappar,file=mapparfile_bis,status='old')
  endif

  iomappar=0

  MAPPARtemp=0
  nlmappar=0
  chrPrec=0
  posPrec=0
  numPrec=0
  tt=0
  nb_champs=0

  if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
     ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
     !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
     call nb_chps_ttr(io_mappar,nb_champs)
     rewind(io_mappar)
     print*,' '
     print*,'Number of columns read in the Marker MAP file = ',nb_champs
     if(nb_champs.eq.4) print*,'No INCL/EXCL indicator in Marker MAP file'
     if(nb_champs.eq.5) print*,'INCL/EXCL indicator is present in Marker MAP file'
     print*,' '
  endif

  do
     if(fmttyppar.eq.'typ_eval') then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.1)) then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:) ! dans ce cas le fichier carte bis est au format typ_eval
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
        if(nb_champs.eq.4) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
        endif
        if(nb_champs.eq.5) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
              print*,'INCLE/EXCL indicator for Marker ',tt,'=',MAPPARtemp(4),' different from 0 and 1'
              STOP 21
           endif
        endif
     endif

     !! NB : on est oblige de faire tous les tests de coherence sur la structure du fichier CARTE VARIANTS PARENTE car si SelAutoSNPPar=0 c est la premiere fois qu on lit le fichier CARTE


     if (MAPPARtemp(1).lt.chrPrec) then
        print*,'The Marker MAP file is not sorted by ascending chromosome Id'
        STOP 10
     endif

     if ((MAPPARtemp(1).eq.chrPrec).and.(MAPPARtemp(2).le.posPrec)) then
        print*,'The Marker MAP file is not sorted by ascending Marker position within chromosome : chr',chrPrec,posPrec,MAPPARtemp(2)
        STOP 10
     endif

     if (MAPPARtemp(3).le.numPrec) then
        print*,'The Marker MAP file is not sorted by ascending Marker number:',numPrec,MAPPARtemp(3)
        STOP 10
     endif

     nlmappar=nlmappar+1
     if (MAPPARtemp(3).ne.nlmappar) then
        print*,'The Marker number sequence from 1 to N is incorrect on line:',nlmappar
        STOP 10
     endif

     if(nlmappar.eq.1) then ! on est au 1er marqueur du 1er chromosome --> on initialise les compteurs etc...
        chrPrec=MAPPARtemp(1)
        posPrec=MAPPARtemp(2)
        numPrec=MAPPARtemp(3)
     endif

     if (MAPPARtemp(1).gt.chrPrec) then   ! on est passe a un nouveau chromosome
        chrprec=MAPPARtemp(1)
        posPrec=0
     endif

     if (MAPPARtemp(3).ne.nlmappar) then
        print*,'The Marker number sequence from 1 to N  is discontinuous :',nlmappar,MAPPARtemp(3)
        STOP 10
     endif
  enddo

  print*,'Number of lines read in the Marker MAP file =',nlmappar

  nbCHROMpar=MAPPARtemp(1)
  print*,' '
  print*,'Number of chromosomes carying Markers considered to model Breeding Values =',nbCHROMpar
  print*,' '
  !if(fmttyppar.eq.'plink') nbSNP_0=nlmappar


  rewind(io_mappar)


  ! on relit le fichier MAPPAR une 2eme fois pour remplir le vecteur INCL_PAR
  ! on sait que le fichier MAPPAR est trie par ordre de chromosomes croissant et par ordre de SNPpar croissant intra chromosome --> 2eme lecture simple
  allocate(INCL_PAR(nlmappar))
  INCL_PAR=0
  i=0
  nb_champs=0

  if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
     ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
     !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
     call nb_chps_ttr(io_mappar,nb_champs)
     rewind(io_mappar)
  endif

  do

     if(fmttyppar.eq.'typ_eval') then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.1)) then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:) ! dans ce cas le fichier vcarte bis est au format typ_eval
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
        if(nb_champs.eq.4) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
        endif
        if(nb_champs.eq.5) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
              print*,'INCL/EXCL indicator for Marker ',tt,'=',MAPPARtemp(4),' different from 0 and 1'
              STOP 21
           endif
        endif

     endif

     i=i+1
     if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
        print*,'INCL/EXCL indicator different from 0 and 1 for Marker ',i
        STOP 10
     endif
     INCL_PAR(i)=MAPPARtemp(4)
  enddo
  print*,'In total INCL/EXCL indicators read for',i,'Markers'

  rewind(io_mappar)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! on lit le fichier des TYPAGES UTILISES POUR LA PARENTE
  ! les identifiants des individus sont des ALPHANUMERIQUES --> on convertit ALPHANUM en NUM de 1 a nbind

  ! ON NE CONSERVERA LES GENOTYPES AUX VARIANTS PARENTE UNIQUEMENT POUR LES INDIVIDUS UTILES = INDIVIDUS AVEC UNE PERFORMANCE


  ! on remplit le vecteur freq2 des frequences utiles d alleles 2

  allocate(freq2(nbSNP_0))
  freq2=0.0d0

  do i=1,nbSNP_0
     freq2(i)=freq(i)
     if (freq(i).lt.MAFmin) freq2(i)=0.0d0
     if (INCL_PAR(i).eq.0) freq2(i)=0.0d0
     if (freq(i).gt.(1-MAFmin)) freq2(i)=1.0d0
  enddo

  print '('' Number of heterozygous SNPs: '',t75,i8)',count( freq*(1d0-freq) /=0.0d0 )

  ! on compte le nb de SNP dont la MAF est superieure a MAFmin >= freqAl2 <= 1-MAFmin ET QUI N ONT PAS ETE VOLONTAIREMENT EXCLUS PAR UTILISATEUR
  print*,' '
  if(detailed_log.eq.3) print '('' Number of Markers such as MAF > MAFmin not excluded by user : '',t75,i8)',count( freq2*(1d0-freq2) /=0.0d0 )

  ! on compte le nb de SNP CONSERVES DANS LES EQUATIONS = CONSIDERES NON MONOMORPHES = SNP tq MAFmin >= freqAl2 <= 1-MAFmin ET ABSENTS DE l EVENTUELLE LISTE D EXCLUSION
  nbSNP=count( freq2*(1d0-freq2) /=0.0d0 )
  print*,'Nb of markers considered in the equations (MAF>=MAFmin and not in user s exclusion list) = nbSNP : ',nbSNP
  print*,' '


  ! on alloue les vecteurs SNP0toSNP et SNPtoSNP0
  allocate(SNP0toSNP(nbSNP_0),SNPtoSNP0(nbSNP))
  SNP0toSNP=0
  SNPtoSNP0=0

  ! on remplit le vecteur SNP0toSNP
  ! la ieme ligne du vecteur SNP0toSNP donne pour le ieme SNP lu dans le fichier typages son numero de SNPutile
  ! on remplit le vecteur SNPtoSNP0
  ! la jeme ligne du vecteur SNPtoSNP0 donne pour le jeme SNP utile son numero d ordre dans le fichier typages 

  j=0
  do i=1,nbSNP_0
     if ((freq2(i).ne.0).and.(freq2(i).ne.1)) then
        j=j+1
        SNP0toSNP(i)=j
        SNPtoSNP0(j)=i
     endif
  enddo

  ! on ecrit dans la log pour validation le nb de SNP utiles comptes lors du remplissage de SNP0toSNP
  ! et de SNPtoSNP0
  !print*,' '
  !print*,'Nb de SNP tq MAF >= MAFmin lors du remplissage des vecteurs correspondance SNP0toSNP et SNPtoSNP0 (j)=',j
  !print*,' '
  !print*,'10 premiers elements du vecteur SNP0toSNP :'
  !print*,SNP0toSNP(1:min(10,nbSNP_0))
  !print*,' '
  !print*,'10 premiers elements du vecteur SNPtoSNP0 :'
  !print*,SNPtoSNP0(1:min(10,nbSNP))
  !print*,' '

  ! on alloue et on remplit le vecteur temporaire frequtil qui contiendra les frequences de l allele 2 aux SNP utiles pour recentrer les typages
  allocate(frequtil(nbSNP))
  frequtil=0.0d0

  do i=1,nbsnp_0
     if(SNP0toSNP(i).ne.0) then
        frequtil(SNP0toSNP(i))=freq2(i)
     endif
  enddo

  if(detailed_log.eq.3) then
     print*,' '
     print*,'First 10 elements of vector frequtil:'
     print*,frequtil(1:min(10,nbSNP))
  endif

  ! on calcule sumpq_all et sumpq sur tous les SNP pour le cas general ou tous les SNP ont la meme variance
  sumpq_all=sum(freq*(1d0-freq))
  sumpq=sum(freq2*(1d0-freq2))
  print*,' '
  print'('' sumpq_all calculated considering all Markers (including Markers with MAF < MAFmin) ='',t110,f20.5)',sumpq_all
  print'('' sumpq calculated considering only Markers used in the equations (MAF >= MAFmin and not excluded) ='',t110,f20.5)',sumpq
  print*,' '

  allocate(gsnpparinv(1,1))
  gsnpparinv=0.0d0
  gsnpparinv(1,1)=(2.0d0*sumpq)*ginv(1,1)
  print*,' '
  print*,'gsnpparinv = inverse of prior variance of each Marker = ',gsnpparinv(1,1)
  print*,' '
  allocate(Rho(1,1))
  Rho(1,1)=R(1,1)*gsnpparinv(1,1)


  if(detailed_log.eq.3) then
     print*,' '
     print*,'OBSERVED frequency of allele 2 for the first 10 Markers :'
     do i=1,min(10,nbSNP)
        print*,'SNP ',i,'  freq all 2 = ',freq(i)
     enddo

     print*,' '
     print*,'CONSIDERED frequncy of allele 2 for the first 10 premiers Markers :'
     do i=1,min(10,nbSNP)
        print*,'SNP ',i,'  freq all 2 utile = ',freq2(i)
     enddo

     print*,' '
     print*,'Most frequent observed genotype (in nb of alleles 2) for the first 10 Markers:'
     do i=1,min(10,nbSNP)
        print*,'For Marker ',i,'  the most often observed genotype in Marker genotype file is ',mostFreqTyp(i)
     enddo
  endif
  print*,' '


  ! on lit le fichier typages parente totalement une derniere fois pour stocker en memoire les typages aux variants_parente utiles
  ! tels que MAF >= MAFmin


  open (io_typpar,file=trim(typparfile),form='FORMATTED')     


  ! on alloue le vecteur temporaire genopartemp qui contiendra les typages aux SNP PARENTE utiles et inutiles de l individu en cours de traitement
  !allocate(GENOPARtemp(nbSNP_0))
  GENOPARtemp=0


  if( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) then  ! on a besoin de GENOPAR (en partie 3 de toute maniere) et (en partie 1 si meth_mpm=1)
     allocate(genopar(nbind,nbSNP))   ! on ne considere que les typages des individus utiles = individus avec perf --> nb lignes de la table genopar = nbind
     genopar=0
  endif
  if((meth_mpm.eq.2).and.(steptodo.eq.1)) then
     allocate(SNPparUtCentr(nbSNP,nbind))   ! on ne considere que les typages des individus utiles = individus avec perf --> nb colonnes de la table SNPparUtCentr = nbind
     SNPparUtCentr=0.0d0
  endif

  nc=1
  io=0
  nbIndTypInco=0
  nbtypindet_tot=0
  print*,' '
  print*,'**********************************************************************************'
  if(steptodo.eq.1) then 
     print*,'Individuals with unknown genotype for Kinship Markers (max 100 indiv by default):'
     print*,'Add OPTION PrintMissGeno 1 in parameter file for STEP_1 to print all individuals with unknown genotypes'
  endif

  if(fmttyppar.eq.'plink') then
     read(io_typpar,*,iostat=io) ! on lit la ligne d en-tetes qui ne sert a rien
  endif

  allocate(liste_SNPpar_inco(nbSNP_0))  ! vecteur temporaire des marqueurs parente inconnus d un individu

  do 
     if(fmttyppar.eq.'typ_eval') then
        read(io_typpar,informatP1,iostat=io) NoNatTemp,GENOPARtemp(:)
        if(io.ne.0) exit
     endif

     if(fmttyppar.eq.'plink') then
        read(io_typpar,*,iostat=io) NoNatTemp2,NoNatTemp,SIRE_temp,DAM_temp,SEX_temp,PHENO_temp,GENOPARtemp(:)
        if(io.ne.0) exit
     endif

     if (mod(nc,5000)==0) print '(1x,i12," lines of Marker genotype file read")', nc

     if(TypParUtil(nc).ne.0) then ! l individu a une perf donc sera utile pour les GWAS

        nbtypindet=0 ! Nouvel individu : on reinitialise le compteur de typages indetermine aux marqueurs parente de l'individu a 0
        liste_SNPpar_inco=0 ! on reinitialise le vecteur des marqueurs parente inconnus de l individu

        do i=1,nbSNP_0
           if(SNP0toSNP(i).ne.0) then 
              if (GENOPARtemp(i).eq.0) then
                 if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) )  GENOPAR(TypParUtil(nc),SNP0toSNP(i))=0
                 if ((meth_mpm.eq.2).and.(steptodo.eq.1)) SNPparUtCentr(SNP0toSNP(i),TypParUtil(nc))= (0.0d0 - 2.0d0*freqUtil(SNP0toSNP(i)))*sqrt(matPOIDS(TypParUtil(nc),1))
              endif
              if (GENOPARtemp(i).eq.1) then
                 if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) GENOPAR(TypParUtil(nc),SNP0toSNP(i))=1
                 if ((meth_mpm.eq.2).and.(steptodo.eq.1)) SNPparUtCentr(SNP0toSNP(i),TypParUtil(nc))= (1.0d0 - 2.0d0*freqUtil(SNP0toSNP(i)))*sqrt(matPOIDS(TypParUtil(nc),1))
              endif
              if (GENOPARtemp(i).eq.2) then
                 if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) GENOPAR(TypParUtil(nc),SNP0toSNP(i))=2
                 if ((meth_mpm.eq.2).and.(steptodo.eq.1)) SNPparUtCentr(SNP0toSNP(i),TypParUtil(nc))= (2.0d0 - 2.0d0*freqUtil(SNP0toSNP(i)))*sqrt(matPOIDS(TypParUtil(nc),1))
              endif
              ! on remplace arbitrairement typage inconnu par typage observe le plus frequent
              if ((GENOPARtemp(i).ne.0).and.(GENOPARtemp(i).ne.1).and.(GENOPARtemp(i).ne.2)) then
                 nbtypindet=nbtypindet+1 ! compteur de l individu
                 nbtypindet_tot=nbtypindet_tot+1 ! compteur total sur l ensemble des individus
                 if( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) )  GENOPAR(TypParUtil(nc),SNP0toSNP(i)) = mostFreqTyp(i)
                 if((meth_mpm.eq.2).and.(steptodo.eq.1))  SNPparUtCentr(SNP0toSNP(i),TypParUtil(nc))= (dfloat(mostFreqTyp(i)) - 2.0d0*freqUtil(SNP0toSNP(i)))*sqrt(matPOIDS(TypParUtil(nc),1))

                 liste_SNPpar_inco(nbtypindet)=i
                 !print*,'Genotyped individual =',NoNatTemp,' has genotype =',GENOPARtemp(i),' at Marker ',i
              endif
           endif
        enddo
        ! on a termine la lecture des genotypes aux marquers parente de l individu --> on ecrit dans la LOG la liste des marqueurs pour lesuqels il a genotype inconnu
        if(nbtypindet.gt.0) then
           nbIndTypInco=nbIndTypInco+1
           if(steptodo.eq.1) then ! print only in STEP_1
              if((nbIndTypInco.le.100).or.(PrintMissGeno.eq.1)) then ! if user choosed to printlist markers with unknown genotypes in LOG for all animals with unknown genotypes
                 !print*,'Genotyped individual =',NoNatTemp,' has unknown genotype for ',nbtypindet,'Markers : Markers =',liste_SNPpar_inco(1:nbtypindet)
                 write(*,'(A,A,A,I0,A,*(I0:1X))') 'Genotyped indiv =',trim(NoNatTemp),' has unknown genotype for ',nbtypindet,'Markers : Markers =',liste_SNPpar_inco(1:nbtypindet)
              endif
           endif
        endif

     endif

     nc=nc+1
     if (nc.eq.(nltyppar+1)) exit
  enddo

  print*,' '
  write(*,'(A,I0,A)') 'In total ',nbIndTypInco,' individuals have at least 1 unknown genotype for kinship markers'
  write(*,'(A,I0,A)') 'In total ',nbtypindet_tot,' unknown genotypes for kinship markers have been replaced by the most frequent genotype'
  print*,'**********************************************************************************'
  print*,' '

  deallocate(liste_SNPpar_inco)

  rewind(io_typpar)


  if( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) then
     if((nbind.lt.20).and.(nbSNP_0.lt.200)) then
        print*,' '
        print*,'GENOPAR matrix built only with Markers with MAF > MAFmin:'
        do i=1,nbind
           print*,GENOPAR(i,:)
        enddo
        print*,' '
     endif
  endif

  close(io_typpar)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! on lit le fichier carte pour les marqueurs parente : 3 colonnes : num_chrom / position_sur_chrom / numero_marqueur_dans_fichier_typages_parente

  allocate(MAPPAR(nbSNP,2))  ! 1ere colonne = num_chromos / 2eme colonne = position sur chromosome
  MAPPAR=0
  allocate(FL_PAR(2,nbCHROMpar)) ! 1ere ligne = numero ligne du 1er marqueur parente du chrom dans table MAPPAR / 2eme ligne = numero ligne du dernier marqueur parente du chrom dans table MAPPAR 
  FL_PAR=0
  MAPPARtemp=0
  nlmappar=0
  chrPrec=0
  posPrec=0
  numPrec=0
  tt=0
  nb_champs=0

  if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
     ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant
     !on compte le nb de champs dans le fichier carte formt plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
     call nb_chps_ttr(io_mappar,nb_champs)
     rewind(io_mappar)
  endif


  do

     if(fmttyppar.eq.'typ_eval') then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:)
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.1)) then
        read(io_mappar,*,iostat=iomappar) MAPPARtemp(:) ! dans ce cas le fichier carte bis est au format typ_eval
        if (iomappar .ne. 0) exit
     endif

     if((fmttyppar.eq.'plink').and.(SelAutoSNPpar.eq.0)) then
        if(nb_champs.eq.4) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           MAPPARtemp(4)=1 ! si fichier carte format plink pas d indicateur inclusion/eclusion des variants parente doone par l utilisateur 
        endif
        if(nb_champs.eq.5) then
           read(io_mappar,*,iostat=iomappar) MAPPARtemp(1),nomvar,poscM,MAPPARtemp(2),MAPPARtemp(4)
           if (iomappar .ne. 0) exit
           tt=tt+1
           MAPPARtemp(3)=tt
           if((MAPPARtemp(4).ne.0).and.(MAPPARtemp(4).ne.1)) then
              print*,'INCL/EXCL indicator for Marker ',tt,'=',MAPPARtemp(4),' different from 0 and 1'
              STOP 21
           endif
        endif
     endif


     if (MAPPARtemp(1).lt.chrPrec) then
        print*,'The Marker MAP file is not sorted by ascending chromosome Id'
        STOP 10
     endif

     if ((MAPPARtemp(1).eq.chrPrec).and.(MAPPARtemp(2).le.posPrec)) then
        print*,'The Marker MAP file is not sorted by ascnding Marker position within chromosom: chr',chrPrec,posPrec,MAPPARtemp(2)
        STOP 10
     endif

     if (MAPPARtemp(3).le.numPrec) then
        print*,'The Marker MAP file is not sorted by ascending Marker number:',numPrec,MAPPARtemp(3)
        STOP 10
     endif

     nlmappar=nlmappar+1
     if(nlmappar.eq.1) then ! on est au 1er marqueur du 1er chromosome --> on initialise les compteurs etc...
        chrPrec=MAPPARtemp(1)
        posPrec=MAPPARtemp(2)
        numPrec=MAPPARtemp(3)
     endif

     if(SNP0toSNP(nlmappar).ne.0) then ! le marqueur est utiles (MAF >= MAFmin) et donc est présent dans les MME
        if (FL_PAR(1,MAPPARtemp(1)).eq.0) FL_PAR(1,MAPPARtemp(1))=SNP0toSNP(MAPPARtemp(3))
        if (SNP0toSNP(MAPPARtemp(3)).gt.FL_PAR(2,MAPPARtemp(1))) FL_PAR(2,MAPPARtemp(1))=SNP0toSNP(MAPPARtemp(3))
        MAPPAR(SNP0toSNP(MAPPARtemp(3)),1)=MAPPARtemp(1) ! chromosome
        MAPPAR(SNP0toSNP(MAPPARtemp(3)),2)=MAPPARtemp(2) ! position
     endif

     if (MAPPARtemp(1).gt.chrPrec) then   ! on est passe a un nouveau chromosome
        chrprec=MAPPARtemp(1)
        posPrec=0
     endif

     if (MAPPARtemp(3).ne.nlmappar) then
        print*,'The sequence of Marker numbers in Marker MAP file is not continuous:',nlmappar,MAPPARtemp(3)
        STOP 10
     endif
  enddo

  print*,' '
  print*,'Number of lines read in the Marker MAP file =',nlmappar
  if(nlmappar.lt.nbSNP_0) then
     print*,'The Marker MAP file contains less lines than the numebr of Markers in Marker Genotype file'
     STOP 10
  endif
  if(nlmappar.gt.nbSNP_0) then
     print*,'The Marker MAP file contains more lines than the numebr of Markers in Marker Genotype file'
     STOP 10
  endif

  !print*,' '
  !print*,'table MAPPAR : '
  !print*,'chrom / pos'
  !do i=1,nbSNP
  !   print*,MAPPAR(i,:)
  !enddo

  if(detailed_log.ge.2) then
     print*,' '
     print*,'FL_PAR table: '
     print*,'First Marker (num_prem) / Lasr Marker (num_dern) on each chromosom considered in the equations:'
     do i=1,nbCHROMpar
        print*,'chr ',i,'first Marker ',FL_PAR(1,i), 'last marker ',FL_PAR(2,i)
     enddo
  endif

  close(io_mappar)


  neq_ep = nbSNP + sum(nlev(1:neff))

  print*,' '
  print*,'Number of equations in the system of equations without the Variant effects = neq_ep = ',neq_ep


  allocate(FIRSTNIV_EP(neff+1))
  FIRSTNIV_EP=0
  FIRSTNIV_EP(1)=1  ! pour adresse 1er niveau du 1er effet
  do ef=2,neff+1
     FIRSTNIV_EP(ef)=(sum(nlev(1:ef-1)))+1
  enddo
  !print*,'FIRSTNIV_EP = ',FIRSTNIV_EP(:)

  !dimvec = (( dfloat(neq_ep)*dfloat(neq_ep+1)*0.000005d0))*100000
  dimvec = int(( dfloat(neq_ep)*dfloat(neq_ep+1)*0.5d0))

  print*,' '
  print*,'dimvec=neq_ep*(neq_ep+1)*0.5=',dimvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! DEBUT DE LA PARTIE SPECIFIQUE A STEP 1
  if(steptodo.eq.1) then

     ! on construit la partie du membre de gauche correspondant aux effets environnementaux et aux SNP_PARENTE qu on va inverser :
     ! X'X  X'M
     ! M'X  (M'M + Rho*I)


     ! =================================================================================================================================
     ! Methode 2 : on construit uniquement la triangulaire inferieure de  LEFT_EP sous forme d'un vecteur rempli par colonnes 
     !             et on inverse avec les subroutines LAPACK correspondantes dsptrf et dsptri pour matrices symetriques compacte
     ! =================================================================================================================================

     print*,' '
     call fdate(jour)
     print*,'Starting inf_LEFT_EP building process'
     print*,jour
     print*,' '


     allocate(inf_LEFT_EP(dimvec))
     inf_LEFT_EP=0.0d0

     if(detailed_log.eq.3) print*,'allocation of inf_LEFT_EP done'

     ! on remplit (X' X) pour les effets fixes, covariables et effets aleatoires diagonaux
     do ind=1,nbind
        if(matPERF(ind,1).ne.mis) then
           do ef1=1,neff
              do ef2=1,ef1
                 inf_LEFT_EP(TI((FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1),(FIRSTNIV_EP(ef2)+matNivAnim(ef2,1,ind)-1),neq_ep)) = inf_LEFT_EP(TI((FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1),(FIRSTNIV_EP(ef2)+matNivAnim(ef2,1,ind)-1),neq_ep)) + (matWeight_cov(ef1,ind)*matWeight_cov(ef2,ind)*matPoids(ind,1))
              enddo
           enddo
        endif
     enddo
     if(detailed_log.eq.3) print*,'Non Genetic Fixed effects of the model treated'
     ! cas des effets environnementaux aleatoires diagonaux
     do ef=1,neff
        if(randomtype(ef).eq.g_diag) then
           do i=1,nlev(ef)
              inf_LEFT_EP(TI((FIRSTNIV_EP(ef)+i-1),(FIRSTNIV_EP(ef)+i-1),neq_ep)) = inf_LEFT_EP(TI((FIRSTNIV_EP(ef)+i-1),(FIRSTNIV_EP(ef)+i-1),neq_ep)) + RhoDiag(ef,1,1)
           enddo
        endif
     enddo
     if(detailed_log.eq.3) print*,'Non Genetic random diagonal effects treated'

     print*,' '
     call fdate(jour)
     print*,'Building of XpX block of LEFT_EP completed'
     print*,jour
     print*,' '


     print*,' '

     if((detailed_log.eq.3).and.(sum(nlev(1:neff)).le.20)) then
        ! on ecrit le block 1:sum(nlev(1:neff)) x 1:sum(nlev(1:neff)) de LEFT_EP avant inversion
        print*,' '
        print*,'Block XpWX of LEFT_EP corresponding to non genetic effects before inversion:'
        call ecriBlock_inf(inf_left_ep,neq_ep,dimvec,sum(nlev(1:neff)))
        print*,' '
     endif

     !print*,'test pour partie marquurs pour voir si TI fonctionne bien avec des grands nombres :'
     !i=46123
     !j=47020
     !print*,'TI(46123,47020,neq_ep)=',TI(i,j,neq_ep)
     !print*,' ' 

     !print*,'test pour partie marquurs pour voir si TI fonctionne bien avec des grands nombres :'
     !i=47020
     !j=46123
     !print*,'TI(47020,46123,neq_ep)=',TI(i,j,neq_ep)
     !print*,' ' 

     !print*,'***************************************'
     !print*,'inf_LEFT_EP apres remplissage de XpX ='
     !do i=1,dimvec
     !   print*,inf_LEFT_EP(i)
     !enddo
     !print*,'***************************************'



     !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) then ! on construit M' W M sans avoir stocke les genotypes centres aux SNPpar
     if (meth_mpm.eq.1) then ! on construit M' W M sans avoir stocke les genotypes centres aux SNPpar

        ! on remplit la partie (M' M) de LEFT_EP pour les marqueurs PARENTE
        allocate(vecGENO1(nbind),vecGENO2(nbind))
        vecGENO1=0.0d0
        vecGENO2=0.0d0

        ! REMARQUE : A TESTER : intervertir ci-dessous les boucles sur ef1 et ef2 pour voir si gain de temps au remplissage de inf_LEFT_EP sur cellules contigues en memoire
        ! do ef2=1,nbSNP
        !    do ef1=1,nbSNP
        !        if(ef2.gt.ef1) exit
        !        .....
        !    enddo
        ! enddo


        do ef2=1,nbSNP
           vecGENO2=dfloat(GENOPAR(:,ef2))-2.0d0*freqUtil(ef2)

           do ef1=1,nbSNP
              ! on recentre les genotypes des nind individus pour le ef1 eme marqueur  

              if(ef2.le.ef1) then

                 vecGENO1=(dfloat(GENOPAR(:,ef1))-2.0d0*freqUtil(ef1))*matPOIDS(:,1)
                 !print*,'SNPpar',ef1,'   vecGENO1=',vecGENO1(:)

                 val = (dot_product(vecGENO1,vecGENO2))

                 if(ef1.ne.ef2) then
                    inf_LEFT_EP(TI( FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef2-1 , neq_ep) ) = inf_LEFT_EP(TI( FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef2-1 , neq_ep) ) + val
                 endif

                 if(ef1.eq.ef2) then
                    inf_LEFT_EP(TI( FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef2-1 , neq_ep) ) = inf_LEFT_EP(TI( FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef2-1 , neq_ep) ) + val + Rho(1,1)
                 endif

              endif

           enddo
        enddo
        deallocate(vecGENO1,VECGENO2)

     endif ! fin du test if (meth_mpm.eq.1)

     !print*,'table SNPparUtCentr='
     !do i=1,nbSNP
     !   print*,'SNPpar ',i,'genotypes centres = ',SNPparUtCentr(:,i)
     !enddo


     if(meth_mpm.eq.2) then ! on construit M' W M en ayant stocke les genotypes centres aux SNPpar

        allocate(mattemp(nbSNP,nbSNP))
        print*,'allocation of mattemp done'
        mattemp=0.0d0
        do i=1,nbSNP
           mattemp(i,i)= Rho(1,1)
        enddo

        call dgemm('N','T',nbSNP,nbSNP,nltyppar,1.0d0,SNPparUtCentr,nbSNP,SNPparUtCentr,nbSNP,1.0d0,mattemp,nbSNP) ! le POIDS DES PHENOTYPES est deja integre
        print*,'after dgemm call'
        do ef1=1,nbSNP
           inf_LEFT_EP( TI( FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef1-1 , neq_ep) : TI(FIRSTNIV_EP(neff+1)+ef1-1 , FIRSTNIV_EP(neff+1)+ef1-1 , neq_ep) + (nbSNP - ef1) ) = mattemp(ef1:nbSNP,ef1)
        enddo
        print*,'inf_LEFT_EP vector filled'

        !print*,'table mattemp :'
        !do i=1,nbSNP
        !   print*,mattemp(i,:)
        !enddo

        deallocate(mattemp)

     endif ! fin du test if(meth_mpm.eq.2)

     print*,' '
     call fdate(jour)
     print*,'MpM block of LEFT_EP filled'
     print*,jour
     print*,' '

     !print*,'***************************************'
     !print*,'inf_LEFT_EP apres remplissage de MpM ='
     !do i=1,dimvec
     !   print*,inf_LEFT_EP(i)
     !enddo
     !print*,'***************************************'



     ! on remplit la partie (M' W X) de LEFT_EP : les effets environnementaux (moyenne, effets fixes, cov, diag) sont AVANT les SNPparente donc on est toujours sous la diagonale --> TI(a,b):TI(a+nbSNP,b) est bien une plage continue du vecteur inf_LEFT_EP

     !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) )  then
     if (meth_mpm.eq.1) then
        allocate(vecGENO1(nbSNP))
        vecGENO1=0.0d0
     endif

     allocate(vectemp(nbSNP))
     vectemp=0.0d0

     do ind=1,nbind
        !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) then
        if (meth_mpm.eq.1) then
           vecGENO1=dfloat(GENOPAR(ind,:))-2.0d0*freqUtil(:)
        endif
        !print*,'Indiv',ind,'   vecGENO1=',vecGENO1
        do ef=1,neff
           !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) ) vectemp=matPOIDS(ind,1)*matWeight_cov(ef,ind)*vecGENO1(1:nbSNP)
           if (meth_mpm.eq.1) vectemp=matPOIDS(ind,1)*matWeight_cov(ef,ind)*vecGENO1(1:nbSNP)
           if (meth_mpm.eq.2) vectemp=sqrt(matPOIDS(ind,1))*matWeight_cov(ef,ind)*SNPparUtCentr(:,ind)

           inf_LEFT_EP( TI( FIRSTNIV_EP(neff+1), (FIRSTNIV_EP(ef)+matNivAnim(ef,1,ind)-1) ,neq_ep) : TI (FIRSTNIV_EP(neff+1)+nbSNP-1 , (FIRSTNIV_EP(ef)+matNivAnim(ef,1,ind)-1) ,neq_ep)) = inf_LEFT_EP( TI( FIRSTNIV_EP(neff+1), (FIRSTNIV_EP(ef)+matNivAnim(ef,1,ind)-1) ,neq_ep) : TI (FIRSTNIV_EP(neff+1)+nbSNP-1 , (FIRSTNIV_EP(ef)+matNivAnim(ef,1,ind)-1) ,neq_ep))   + vectemp  ! M' W X

        enddo
     enddo

     !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) )  deallocate(vecGENO1)
     if (meth_mpm.eq.1)  deallocate(vecGENO1)

     deallocate(vectemp)

     print*,'***************************************'
     if(dimvec.le.400) then
        print*,'inf_LEFT_EP at the end of the process='
        do i=1,dimvec
           print*,inf_LEFT_EP(i)
        enddo
     endif
     print*,'***************************************'

     !print*,' '
     !print*,'inf_LEFT_EP = '
     !print*,inf_LEFT_EP(:)
     !print*,' '


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! ======================================================================================================================================
     ! on construit la partie relative aux effets environnementaux et a tous les SNP_parente tq MAF > MAFmin du membre de droite = RIGHT_EP
     ! ======================================================================================================================================

     print*,' '
     call fdate(jour)
     print*,'Before building of RIGHT_EP'
     print*,jour
     print*,'*******************************'


     allocate(RIGHT_EP(neq_ep))
     RIGHT_EP=0.0d0

     ! partie X'y
     do ind=1,nbind

        if(matPERF(ind,1).ne.mis) then

           do ef1=1,neff
              RIGHT_EP((FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) = RIGHT_EP((FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) + (matWeight_cov(ef1,ind)*matPoids(ind,1)*matPERF(ind,1))
           enddo

        endif

     enddo

     ! partie M'y
     !if ( ((meth_mpm.eq.1).and.(steptodo.eq.1)) .or. (steptodo.eq.3) )  then
     if (meth_mpm.eq.1)  then
        allocate(vecGENO1(nbind))
        vecGENO1=0.0d0
        do ef1=1,nbSNP
           ! on recentre les genotypes des nind individus pour le ef1 eme marqueur 
           vecGENO1=(dfloat(GENOPAR(:,ef1))-2.0d0*freqUtil(ef1))*matPOIDS(:,1) ! un individu qui a perf manquante mais present dans fichier typages a POIDS=0 = ne contribuera pas

           val = dot_product(vecGENO1(:),matPERF(:,1))

           RIGHT_EP(FIRSTNIV_EP(neff+1)+ef1-1) = RIGHT_EP(FIRSTNIV_EP(neff+1)+ef1-1) + val
        enddo

        deallocate(vecGENO1)
     endif ! fin du test if (meth_mpm.eq.1)

     if(meth_mpm.eq.2) then
        call dgemv('N',nbSNP,nbind,1.0d0,SNPparUtCentr,nbSNP,vecPERFpond,1,0.0d0,RIGHT_EP(FIRSTNIV_EP(neff+1):FIRSTNIV_EP(neff+1)+nbSNP-1),1)
     endif  ! fin du test if(meth_mpm.eq.2)


     print*,' '
     call fdate(jour)
     print*,'Building of RIGHT_EP completed'
     print*,jour
     print*,'*******************************'

     if((detailed_log.eq.3).and.(sum(nlev(1:neff)).le.20)) then
        print*,' '
        print*,'Vector XpWY of RIGHT_EP corresponding to non genetic effects:'
        print*,RIGHT_EP(1:sum(nlev(1:neff)))
        print*,' '
     endif

     print*,'*******************************'
     if(neq_ep.le.20) then
        print*,'RIGHT_EP = '
        print*,RIGHT_EP(:)
        print*,' '
     endif
     print*,'*******************************'

     call fdate(jour)
     print*,jour
     print*,'Before writing RIGHT_EP on disc'
     print*,' '

     ! on sort sur disque la matrice RIGHT_EP en format binaire
     open(io_right,file=trim(location)//trim("/RIGHT_binfile"),form="unformatted") 
     write(io_right)RIGHT_EP(:)
     close(io_right)

     deallocate(RIGHT_EP)

     call fdate(jour)
     print*,jour
     print*,'writing of RIGHT_EP on disc completed'
     print*,' '




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! on n a plus besoin de SNPparUtCentr ou de GENOPAR

     if(meth_mpm.eq.1) then
        deallocate(GENOPAR)
     endif

     if(meth_mpm.eq.2) deallocate(SNPparUtCentr)


     print*,' '
     call fdate(jour)
     print*,'Beginning of inf_LEFT_EP inversion'
     print*,jour
     print*,'*******************************'


     call invmatlapack_inf(inf_LEFT_EP,neq_ep,dimvec)

     print*,' '
     call fdate(jour)
     print*,'inf_LEFT_EP inversion completed'
     print*,jour
     print*,' '
     !print*,inf_LEFT_EP(:)
     !print*,' '


     ! =================================================================================================================
     ! on sort la matrice inf_LEFT_EP sur disque format binaire triangulaire inferieure compactee par colonnes
     ! a la relectuere le nb d elements de inf_LEFT_EP et les dimensions de inf_LEFT_EP_act seront calcules a la 
     ! lecture du fichier parametres 

     call fdate(jour)
     print*,jour
     print*,'Before wrinting inf_LEFT_EP inverse on disc'
     print*,' '

     open(io_invleft,file=trim(location)//trim("/invLEFT_binfile"),form="unformatted") 
     write(io_invleft)inf_LEFT_EP(:)
     close(io_invleft)

     call fdate(jour)
     print*,jour
     print*,'Writing of inf_LEFT_EP inverse on disc completed'
     print*,' '

     ! on supprime inf_LEFT_EP pour tester la lecture de invLEFT_binfile et le remplissage de inf_LEFT_EP_act a partir du fichier sur disque

     deallocate(inf_LEFT_EP)


     call fdate(jour)
     print*,jour
     print*,'STEP_1 COMPLETED'


  endif ! fin de la partie executee uniquement si steptodo=1

  !==========================================================================================================
  !==========================================================================================================
  !==========================================================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
  ! DEBUT DE LA PARTIE 2 !
!!!!!!!!!!!!!!!!!!!!!!!!

  if(steptodo.eq.2) then


     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! si OPTION COJO : on lit le fichier donnant l identite du variant COJO
     if(indic_cojo.eq.1) then
        open(io_IdCojo,file=cojofile,status='old')
        do
           read(io_IdCojo,*,iostat=iocojo) numvar_cojo
           if (iocojo .ne. 0) exit
           nlcojo=nlcojo+1
        enddo
        print*,nlcojo,' lines read in file ',trim(cojofile)
        if(nlcojo.ne.1) then
           print*,'PROBLEM : 1 single COJO Variant expected in the file'
           STOP 22
        endif
        print*,'COJO variant is variant ',numvar_cojo,' in Variant Genotype file'
        print*,' '
     endif


     ! on lit le fichier CARTE pour les marqueurs GWAS pour connaitre le nb de marqueurs GWAS pour lecture fichier DOSES MINIMAC

     open(io_mapgwas,file=mapgwasfile,status='old')

     select case(fmttypgwas)

     case('plink')

        nb_chGWAS=0

        !on compte le nb de champs dans le fichier carte format plink : si n=4 il n y a pas d indicateur incl/excl ; si n=5 il y a un indicateur incl/excl
        call nb_chps_ttr(io_mapgwas,nb_chGWAS)
        rewind(io_mapgwas)
        if(nb_chGWAS.eq.4) print*,'No INCL/EXCL indicator in Variant MAP file'
        if(nb_chGWAS.eq.5) print*,'INCL/EXCl indicator present in Variant MAP GWAS'

        ! on lit une 1ere fois le fichier CARTE_GWAS plink pour connaitre le nb de variants pour GWAS sur le chromosome considere pour dimensionner vecteurs et tableaux
        ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant

        iomapgwas=0
        nlmapgwas=0
        do
           read(io_mapgwas,*,iostat=iomapgwas) i,nomvar,poscM,pospb
           if (iomapgwas .ne. 0) exit
           if(i.eq.NUMCHR_GWAS) nlmapgwas=nlmapgwas+1
        enddo

        rewind(io_mapgwas)

        print*,'Number of lines read in Variant MAP file for chromosom ',NUMCHR_GWAS,' is ',nlmapgwas

        allocate(MAPGWAS(3,nlmapGWAS))  ! 1ere ligne = num_chromos / 2eme ligne = position sur chromosome / 3eme ligne = inclusion/exclusion des GWAS
        MAPGWAS=0
        posgwasPrec=0
        nbGWAS=0
        i=0
        iprec=0
        pospb=0
        firstVarGWAS=0
        lastVarGWAS=0
        nltot=0

        ! on lit une 2eme fois le fichier CARTE GWAS pour remplir la table MAPGWAS
        nlmapgwas=0
        do

           if(nb_chGWAS.eq.4) then
              read(io_mapgwas,*,iostat=iomapgwas) i,nomvar,poscM,pospb
              nltot=nltot+1
              if (iomapgwas .ne. 0) exit
              if(i.eq.NUMCHR_GWAS) then ! le variant est sur le chromosome qu on etudie en ETAPES 2 et 3
                 nlmapgwas=nlmapgwas+1
                 if(firstVarGWAS.eq.0) firstVarGWAS=nltot ! ordre du premier variantGWAS du chromosome parmi l ensemble des variants ordonnes de l ensemble des chromosomes
                 lastVarGWAS=nltot                        ! ordre du dernier variantGWAS du chromosome parmi l ensemble des variants ordonnes de l ensemble des chromosomes
                 if(pospb.le.posgwasPrec) then
                    print*,'The GWAS variant MAP file for chromosom ',NUMCHR_GWAS,' is not sorted by ascending Variant position'
                    STOP 30
                 endif
                 MAPGWAS(1,nlmapgwas)=NUMCHR_GWAS
                 MAPGWAS(2,nlmapgwas)=pospb
                 MAPGWAS(3,nlmapgwas)=1
                 posgwasPrec=pospb
              endif
              iprec=i
              if(i.lt.iprec) then
                 print*,'The GWAS variant MAP file for chromosom ',NUMCHR_GWAS,' is not sorted by ascending chromosome Id'
                 STOP 30
              endif
              nbGWAS=nlmapgwas
           endif


           if(nb_chGWAS.eq.5) then
              intTEMP=0
              read(io_mapgwas,*,iostat=iomapgwas) i,nomvar,poscM,pospb,INTtemp

              nltot=nltot+1
              if (iomapgwas .ne. 0) exit
              if(i.eq.NUMCHR_GWAS) then ! le variant est sur le chromosome qu on etudie en ETAPES 2 et 3
                 nlmapgwas=nlmapgwas+1
                 if(firstVarGWAS.eq.0) firstVarGWAS=nltot ! ordre du premier variantGWAS du chromosome parmi l ensemble des variants ordonnes de l ensemble des chromosomes
                 lastVarGWAS=nltot                        ! ordre du dernier variantGWAS du chromosome parmi l ensemble des variants ordonnes de l ensemble des chromosomes
                 if(pospb.le.posgwasPrec) then
                    print*,'The GWAS Varant MAP file for chromosom ',NUMCHR_GWAS,' is not sorted by ascending position'
                    STOP 30
                 endif
                 MAPGWAS(1,nlmapgwas)=NUMCHR_GWAS
                 MAPGWAS(2,nlmapgwas)=pospb
                 if((INTtemp.ne.0).and.(INTtemp.ne.1)) then
                    print*,'INCL/EXCL indicator in GWAS Variant MAP file is different from 0 and 1 for GWAS variant',nlmapgwas,'/',nltot
                    STOP 30
                 endif
                 MAPGWAS(3,nlmapgwas)=INTtemp
                 if(INTtemp.eq.1) nbGWAS=nbGWAS+1
                 posgwasPrec=pospb
              endif
              iprec=i
              if(i.lt.iprec) then
                 print*,'The GWAS Variant MAP file for chromosome ',NUMCHR_GWAS,' is not sorted by ascending chromosome Id'
                 STOP 30
              endif

           endif

        enddo

        print*,'The number of GWAS Variants on the considered chromosome  = nbGWAS is ',nbGWAS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fin du cas fmttypgwas='plink'


     case default
        ! on lit une 1ere fois le fichier CARTE_GWAS pour connaitre le nb de variants pour GWAS pour dimensionner vecteurs et tableaux
        ! on ne considerera PAS POUR LES GWAS les variants qui ont leur indicateur inclusion/exclusion renseigné à 0 par l utilisateur dans fichier CARTE_GWAS

        iomapgwas=0
        nlmapgwas=0
        do
           read(io_mapgwas,*,iostat=iomapgwas)
           if (iomapgwas .ne. 0) exit
           nlmapgwas=nlmapgwas+1
        enddo

        rewind(io_mapgwas)

        print*,'Number of lines read in the Variant MAP file : ',nlmapgwas

        allocate(MAPGWAS(3,nlmapGWAS))  ! 1ere ligne = num_chromos / 2eme ligne = position sur chromosome / 3eme ligne = inclusion/exclusion des GWAS
        MAPGWAS=0
        MAPGWAStemp=0
        posgwasPrec=0
        numgwasPrec=0
        nbGWAS=0

        ! on lit une 2eme fois le fichier CARTE GWAS pour remplir la table MAPGWAS
        nlmapgwas=0
        do
           read(io_mapgwas,*,iostat=iomapgwas) MAPGWAStemp(:)
           if (iomapgwas .ne. 0) exit
           nlmapgwas=nlmapgwas+1
           if((MAPGWAStemp(2).le.posgwasPrec).or.(MAPGWAStemp(3).le.numgwasPrec).or.(MAPGWAStemp(3).ne.nlmapgwas)) then
              print*,'Variant MAP file is not sorted by ascending position / Variant number'
              STOP 30
           endif
           if((MAPGWAStemp(4).ne.0).and.(MAPGWAStemp(4).ne.1)) then
              print*,'INCLUSION/EXCLUSION indicator in Variant MAP file is different from 0 and 1 for Variant',nlmapgwas
              STOP 30
           endif
           MAPGWAS(1:2,nlmapgwas)=MAPGWAStemp(1:2)
           MAPGWAS(3,nlmapgwas)=MAPGWAStemp(4)
           if(MAPGWAStemp(4).eq.1) nbGWAS=nbGWAS+1
           posgwasPrec=MAPGWAStemp(2)
           numgwasPrec=MAPGWAStemp(3)
        enddo

     end select

     close(io_mapgwas)

     print*,' '
     print*,'The number of Variants to be considered in the GWAS is :',nbGWAS
     print*,' '
     allocate(GWASuttoMAP(nbGWAS),MAPtoGWASut(nlmapgwas))
     GWASuttoMAP=0
     MAPtoGWASut=0

     !print*,' '
     !print*,'table MAPGWAS :'
     j=0
     do i=1,nlmapgwas
        !   print*,'Variant GWAS ',i,'  chrom ',MAPGWAS(1,i),'  position ',MAPGWAS(2,i),'   excl/incl ',MAPGWAS(3,i)
        if(MAPGWAS(3,i).eq.1) then
           j=j+1
           GWASuttoMAP(j)=i
           MAPtoGWASut(i)=j
        endif
     enddo

     if(detailed_log.eq.3) then
        print*,' '
        print*,'First 10 elements in vector GWASuttoMAP =',GWASuttoMAP(1:min(10,nbGWAS))
        print*,'First 10 elements in vector MAPtoGWASut =',MAPtoGWASut(1:min(10,nlmapgwas))
        print*,' '
     endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !======================================================================================================================================
     ! modif du 08 fev 2023 pour permettre d avoir DES INDIVIDUS SUPPLEMENTAIRES DANS FICHIER DOSES QUI SONT ABSENTS DE PERF ET TYPAGES_PAR

     allocate(DOSEStemp(nlmapgwas),TypGwasTemp(nlmapgwas))
     DOSEStemp=0.0
     nldoses=0
     TypGwasTemp=0


     open(io_doses,file=typgwasfile,status='old')

     ! on lit une 1ere fois le fichier DOSES pour lire les identifiants MINIMAC et couper la chaine en 2 pour retrouver l identifiant animal

     do
        read(io_doses,*,iostat=io) 
        if (io .ne. 0) exit
        nldoses=nldoses+1
     enddo
     if(fmttypgwas.eq.'plink') nldoses=nldoses-1 ! on retire la 1ere ligne d en tetes qui est dans le fichier typages .raw de plink

     rewind(io_doses)


     allocate(animdoses(nldoses))

     nldoses=0

     if(trim(fmttypgwas).eq.'minimac') then
        do
           read(io_doses,*,iostat=io) IdDosesTemp
           if (io .ne. 0) exit
           nldoses=nldoses+1
           ! on cherche la position du signe Superieur dans la chaine de caractere et on cree l identifiant reel en lisant la chaine de 1 a pos_sup
           pos_sup=0
           pos_sup=scan(IdDosesTemp,'>')
           animdoses(nldoses)=IdDosesTemp(1:(pos_sup-2))
        enddo
        rewind(io_doses)
     endif



     if(trim(fmttypgwas).eq.'plink') then
        read(io_doses,*,iostat=io) ! on lit la 1ere ligne d en tete qui ne sert a rien et qui est en alphanumerique
        do
           read(io_doses,*,iostat=io) NoNatTemp2 , NoNatTemp ! dans fichier .raw de plink l identifiant de l individu est en 2eme champ
           if (io .ne. 0) exit
           nldoses=nldoses+1
           animdoses(nldoses)=NoNatTemp
        enddo
        rewind(io_doses)
     endif



     if(trim(fmttypgwas).eq.'typ_eval') then
        ! on lit la 1ere ligne du fichier pour repérer la position du 1er typage et creer le format de lecture pour la suite
        read(io_doses,'(a3000000)',iostat=io) linetyp

        len_linetyp=len_trim(linetyp)
        d=len_linetyp
        do while(linetyp(d:d).ne.' ')
           d=d-1
        enddo
        ip_snp=d+1
        nb_varGWAS_fichtyp=len_linetyp-d

        write(informatT,'(a,i0,a,i0,a)') '(a',ip_snp-2,',1x,',nb_varGWAS_fichtyp,'i1)'

        if (detailed_log.eq.3) print*,'informatT=',informatT

        print*,'Number of Variants read in Variant Genotype / allelic dosages file =',nb_varGWAS_fichtyp

        rewind(io_doses)
        nldoses=0
        io=0
        ! on lit le fichier 
        do
           read(io_doses,informatT,iostat=io) NoNatTemp,TypGWASTemp(:)
           if (io .ne. 0) exit
           nldoses=nldoses+1
           animdoses(nldoses)=NoNatTemp
        enddo
        rewind(io_doses)
     endif

     if (nldoses.ne.nbind) then
        print*,'WARNING : number of lines in Variant Genotype / allelic dosages is ',nldoses
        print*,' This is different from the number of individuals with Phenotypes :  ',nbind
        print*,'The unphenotyped individuals in Variant Genotype / allelic dosage file will be discarded'
        !STOP 30
     endif

     ! on trie les ID des individus du fichier DOSES et on cree un vecteur d ordre des individus dans le fichier DOSES
     allocate(orddoses(nldoses))
     orddoses=0

     if(detailed_log.eq.3) then
        if(nldoses.lt.20) then
           print*,'Vector animdoses containing the Id of Variant Genotypes / allelic dosages file BEFORE sorting by Id : ',animdoses(:)
        endif
     endif

     workdoses=''
     call hpsort(animdoses,nldoses,1,lc,orddoses,2,workdoses,ifail)
     print*,'sorting with HPSORT completed'

     if(nldoses.lt.20) then
        print*,'Vector animdoses containing the Id of Variant Genotypes / allelic dosages file AFTER sorting by Id : ',animdoses(:)
     endif


     ! vecteur orddoses : element i = position dans fichier DOSES lu du ieme animal trie
     ! --> on cree vecteur ordinvdoses : element i = classement de l animal de la ieme ligne du fichier doses lu apres tri sur IdAnim
     ! necessaire pour remplir matrices DOSES dans ordre coherent avec tables PERF EFFETS et TYPAGES_PARENTE ...
     allocate(ordinvdoses(nldoses))
     ordinvdoses=0
     do i=1,nldoses
        ordinvdoses(orddoses(i))=i
     enddo

     if(detailed_log.eq.3) then
        if(nldoses.lt.20) then 
           print*,'vector animdoses for comparaison = ',animdoses(:)
           print*,'vector animperf for comparaison = ',animperf(:)
        endif
     endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! on compare animdoses et animperf pour reperer les individus du fichier DOSES inutiles dont on va supprimer les doses pour les GWAS 
     ! les 2 vecteurs sont tries donc on peut avancer pas a pas dans les 2
     allocate(DosesUtil(nldoses),IndDosUt(nbind))
     DosesUtil=0
     IndDosUt=0
     j=1 !indice de la ligne dans animperf()
     do i=1,nldoses ! indice de la ligne dans animdoses()
        if(animdoses(i).eq.animperf(j)) then
           DosesUtil(i)=j
           IndDosUt(j)=i
           j=j+1
           if(j.gt.nbind) exit
        endif
     enddo

     if(detailed_log.eq.3) then
        ! (j-1) pour corriger nb fial car j=j+1 en fin de boucleto correct final value of j because j+1 at the end of the loop
        print*,'After the comparison between vectors animdoses and animperf j=',(j-1)
     endif

     if(j.lt.nbind) then
        print*,' ' 
        print*,'ERROR : Individuals present in PHENOTYPE FILE are missing in VARIANT GENOTYPE / ALLELIC DOSAGE FILE'
        print*,'All individuals present in PHENOTYPE FILE have to be present in Marker genotype file and in Variant genotype / allelic dosage file'
        STOP 33
     endif

     if(detailed_log.eq.3) then
     	print*,'  '
    	print*,'FL_PAR='
    	do i=1,nbCHROMpar
     	   print*,FL_PAR(:,i)
     	enddo
    	print*,'  '
     endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     ! modif V5 : segments de SNPpar recouvrant a 50pct a eliminer si option choisie par utilisateur
     allocate(nbSegmExcl(nbCHROMpar))
     nbSegmExcl=0 ! nombre de segmments de longueur segmrem chevauchants a 50pct de SNPparente a eliminer pour les GWAS pour chaque chromosome

     if(segmrem.eq.9999999999) then
        nbSegmExcl=1
     else if (segmrem.ne.0) then
        !if(StratElim.eq.'window') then
        do i=1,nbCHROMpar
           nbSegmExcl(i)=INT((2*MAPPAR(FL_PAR(2,i),2))/segmrem)
        enddo
        !endif
     else if (segmrem.eq.0) then
        nbSegmExcl=1
     endif

     print*,' '
     print*,'The number of 50pct overlapping segments of length ',segmrem,' for each chromosome is :'
     do i=1,nbCHROMpar
        print*,i,nbSegmExcl(i)
     enddo

     chrom_anal=MAPGWAS(1,1)
     print*,' '
     print*,'The GWAS Variants are on chromosome ',chrom_anal


     allocate(lim_theo(2,nbSegmExcl(chrom_anal)),lim_oper(2,nbSegmExcl(chrom_anal)),lim_SNP(2,nbSegmExcl(chrom_anal)))
     lim_theo=0
     lim_oper=0
     lim_SNP=0

     if(segmrem.eq.9999999999) then
        lim_theo(1,1)=0
        lim_theo(2,1)=segmrem
        lim_oper(1,1)=0
        lim_oper(2,1)=segmrem
        lim_SNP(1,1)=FL_PAR(1,chrom_anal)
        lim_SNP(2,1)=FL_PAR(2,chrom_anal)
     else if (segmrem.ne.0) then
        do i=1,nbSegmExcl(chrom_anal)
           lim_theo(1,i)=(i-1)*(INT(0.5d0*dfloat(segmrem)))
           lim_theo(2,i)=lim_theo(1,i)+segmrem
        enddo
        do i=1,nbSegmExcl(chrom_anal)
           lim_oper(1,i)=lim_theo(1,i) + (INT(0.25d0*dfloat(segmrem)))
           lim_oper(2,i)=lim_theo(2,i) - (INT(0.25d0*dfloat(segmrem)))
        enddo
        lim_oper(1,1)=0
        lim_oper(2,nbSegmExcl(chrom_anal))=lim_theo(2,nbSegmExcl(chrom_anal))
        do i=FL_PAR(1,chrom_anal),FL_PAR(2,chrom_anal)
           do j=1,nbSegmExcl(chrom_anal)
              if(MAPPAR(i,2).lt.lim_theo(1,j)) lim_SNP(1,j)=i
              if(MAPPAR(i,2).le.lim_theo(2,j)) lim_SNP(2,j)=i
           enddo
        enddo
        lim_SNP(1,1)=FL_PAR(1,chrom_anal)
        do j=2,nbSegmExcl(chrom_anal) ! pour le 1er SNPpar a exclure de chaque segment on ajoute 1 car on s etait arrete au SNPpar juste avant le debut du segment
           lim_SNP(1,j)=lim_SNP(1,j)+1
        enddo

        ! modif 16/12/2024 : Si sur un segment on a lim_SNP(1,j) > lim_SNP(2,j) c est que avant d ajouter 1 a lim_SNP(1,j) on avait lim_SNP(1,j) = lim_SNP(2,j)
        !                    c est a dire qu il n y a aucun SNP sur le segment considere --> on met lim_SNP(1,j)=0 et lim_SNP(2,j)=0
        do j=1,nbSegmExcl(chrom_anal)
           if(lim_SNP(1,j).gt.lim_SNP(2,j)) then ! il n y a aucun SNP a retirer pour le segment car aucun SNP sur le segment considere
              lim_SNP(1,j)=0
              lim_SNP(2,j)=0
           endif
        enddo

     else if(segmrem.eq.0) then
        lim_theo(1,1)=0
        lim_theo(2,1)=0
        lim_oper(1,1)=0
        lim_oper(2,1)=0
        lim_SNP(1,1)=0
        lim_SNP(2,1)=0
     endif ! fin du test si segmrem=9999999999 pour eliminer tout le chromosome et segmrem ne 0 pour travailler sur segments

     print*,'Theoretic limits of 50pct overlapping segments of length ',segmrem,' whose Markers are excluded for the GWAS :'
     do i=1,nbSegmExcl(chrom_anal)
        print*,lim_theo(:,i)
     enddo
     print*,' '

     print*,'Limits of segments corresponding to the different Groups of Variants having the same excluded Markers for the GWAS :'
     do i=1,nbSegmExcl(chrom_anal)
        print*,lim_oper(:,i)
     enddo
     print*,' '

     print*,'Id (among useful Markers) of the 1st and last Markers flanking the 50pct overlapping segments of length ',segmrem,' excluded for the GWAS :'
     do i=1,nbSegmExcl(chrom_anal)
        print*,lim_SNP(:,i)
     enddo
     print*,' '


     ! pour un variant GWAS donne on determine quels sont les variants PARENTE a supprimer des equations en fonction de la valeur de segmrem
     ! et on determine le groupe d exclusion dont il fait partie (3eme ligne de SegmElim)

     allocate(SegmElim(3,nbGWAS))
     SegmElim=0
     j=0

     do VarGwas=1,nlmapGWAS

        if(MAPGWAS(3,VarGwas).eq.1) then    ! si l utilisateur veut considerer le variant dans GWAS

           Felim=0 ! numero du 1er SNP_par a eliminer
           Lelim=0 ! numero du dernier SNP_par a eliminer
           LimInf=0
           LimSup=0
           j=j+1

           ! si segmrem=0 on n elimine aucune SNP_par
           if(segmrem.eq.0) then
              Felim=0
              Lelim=0
              groupelim=1
           endif

           ! si segmrem=9999999999 on elimine tous les SNP_par du chromosome du VarGwas considere
           if(segmrem.eq.9999999999) then
              Felim=FL_PAR(1,MAPGWAS(1,VarGwas))
              Lelim=FL_PAR(2,MAPGWAS(1,VarGwas))
              groupelim=1
           endif

           ! si segmrem different de 0 et de 99999 on determine le numero du 1er et du dernier SNP_par a supprimer en fonction des cartes
           if((segmrem.ne.0).and.(segmrem.ne.9999999999)) then

              if(StratElim.eq.'strict') then ! si l utilisateur veut eliminer un segment de longueur GLISSANT DE LONGUEUR FIXE autour du variant GWAS

                 ! limite inferieure du segment sur lequel eliminer les SNPpar
                 LimInf=max( min( MAPPAR(FL_PAR(1,MAPGWAS(1,VarGwas)),2) , MAPGWAS(2,VarGwas) ) ,  MAPGWAS(2,VarGwas) - INT(dfloat(segmrem)*0.5d0) )
                 LimSup=min( max( MAPPAR(FL_PAR(2,MAPGWAS(1,VarGwas)),2) , MAPGWAS(2,VarGwas)  ) ,  MAPGWAS(2,VarGwas) + INT(dfloat(segmrem)*0.5d0) )

                 ! on identifie le 1er et le dernier SNP_par a eliminer des equations
                 do i=FL_PAR(1,MAPGWAS(1,VarGwas)) , FL_PAR(2,MAPGWAS(1,VarGwas)) ! on teste sur tous les SNP_PARENTE du chromosome du variant_GWAS
                    if(MAPPAR(i,2).ge.LimInf) then
                       Felim=i
                       exit
                    endif
                 enddo
                 do i=FL_PAR(1,MAPGWAS(1,VarGwas)) , FL_PAR(2,MAPGWAS(1,VarGwas))
                    if((MAPPAR(i,2).le.LimSup).and.(MAPPAR(i,2).ge.LimInf)) Lelim=i
                    if(MAPPAR(i,2).gt.LimSup) exit
                 enddo

              endif

              if(StratElim.eq.'window') then ! si l utilisateur veut eliminer un segment FIXE PARMI DES SEGMENTS CHEVAUCHANTS A 50PCT

                 do i=1,nbSegmExcl(chrom_anal)
                    if((MAPGWAS(2,VarGwas).ge.lim_oper(1,i)).and.(MAPGWAS(2,VarGwas).le.lim_oper(2,i))) then
                       Felim=lim_SNP(1,i)
                       Lelim=lim_SNP(2,i)
                       groupelim=i
                    endif
                 enddo
                 ! modif 26jun2025
                 if(MAPGWAS(2,VarGwas).gt.lim_oper(2,nbSegmExcl(chrom_anal))) then  ! le variant GWAS est après le dernier marqueur parente du chromosome --> on le met dans le dernier groupe d elimination des SNPpar
                    Felim=lim_SNP(1,nbSegmExcl(chrom_anal))
                    Lelim=lim_SNP(2,nbSegmExcl(chrom_anal))
                    groupelim=nbSegmExcl(chrom_anal)
                 endif

              endif

           endif

           SegmElim(1,j)=Felim
           SegmElim(2,j)=Lelim
           SegmElim(3,j)=groupelim

           !print*,'Variant GWAS numero ',VarGwas,'  Felim / Lelim = ',SegmElim(:,j)
        endif

     enddo

     if(detailed_log.eq.3) print*,'remplissage de Segmelim termine pour ',j,' variants GWAS a considerer dans les GWAS'

     print*,'****************************************************************************************'
     if(j.ne.nbGWAS) then
        print*,'j different de nbGWAS, PROBLEME'
        STOP 30
     endif

     ! on compte le nombre de variants dans chaque groupe d exclusion
     allocate(Bilan_SegmElim(nbSegmExcl(chrom_anal),5))
     Bilan_SegmElim=0
     do i=1,nbSegmExcl(chrom_anal)
        Bilan_SegmElim(i,1)=i
        Bilan_SegmElim(i,2)=lim_SNP(1,i)
        Bilan_SegmElim(i,3)=lim_SNP(2,i)
     enddo
     do j=1,nbGWAS
        Bilan_SegmElim(SegmElim(3,j),4)=Bilan_SegmElim(SegmElim(3,j),4)+1
     enddo

     ! on compte pour chaque groupe combien de sous groupes il y a 
     do i=1,nbSegmExcl(chrom_anal)
        Bilan_SegmElim(i,5)=INT( dfloat(Bilan_SegmElim(i,4)) / dfloat(NbMaxVarGWAS) )
        if( (Bilan_SegmElim(i,5)*NbMaxVarGWAS) .lt. Bilan_SegmElim(i,4) ) then 
           Bilan_SegmElim(i,5)=Bilan_SegmElim(i,5)+1
        endif
     enddo

     print*,'       '
     print*,'Table Bilan_SegmElim summarizing the information on Groups Batches and excluded Marker segments = '
     print*,'Group Id / Id of 1st excludes Marker / Id of last excluded Marker / number of Variants in the Group / Number of Batches in Group'
     do i=1,nbSegmExcl(chrom_anal)
        print*,Bilan_SegmElim(i,:)
     enddo
     print*,' '

     ! on alloue et remplit la table qui contiendra les informations de chaque SousGroupe de GWAS
     allocate(Bilan_SousGr(7,sum(Bilan_Segmelim(:,5))))
     Bilan_SousGr=0
     nbSGr=0
     k=0 ! compteur du nb de variants GWAS 
     do i=1,nbSegmExcl(chrom_anal)
        numSGr=0
        do j=1,Bilan_Segmelim(i,5)
           nbSGr=nbSGr+1
           numSGr=numSGr+1
           Bilan_SousGr(1,nbSGr)=i     ! numero du Groupe
           Bilan_SousGr(2,nbSGr)=numSGr ! numero du SousGroupe (intra_groupe)
           Bilan_SousGr(3,nbSGr)=Bilan_SegmElim(i,2) ! numUT 1er SNPpar elim pour le Groupe
           Bilan_SousGr(4,nbSGr)=Bilan_SegmElim(i,3) ! numUT dernier SNPpar elim pour le Groupe

           Bilan_SousGr(5,nbSGr)=k+1
           if(j.lt.Bilan_Segmelim(i,5)) then
              Bilan_SousGr(6,nbSGr)=k+NbMaxVarGWAS
              k=Bilan_SousGr(6,nbSGr)
           endif
           if(j.eq.Bilan_Segmelim(i,5)) then
              Bilan_SousGr(6,nbSGr)= Bilan_SousGr(5,nbSGr) + (Bilan_Segmelim(i,4)-((j-1)*NbMaxVarGWAS)) -1
              k=Bilan_SousGr(6,nbSGr)
           endif
           Bilan_SousGr(7,nbSGr)=(Bilan_SousGr(6,nbSGr) - Bilan_SousGr(5,nbSGr)) + 1
        enddo
     enddo

     if(detailed_log.ge.2) then
        print*,' '
        print*,'Total number of batches = ',nbSGr
        print*,'Content of table Bilan_SousGr(,) :'
        do i=1,nbSGr
           print*,Bilan_SousGr(:,i)
        enddo
     endif

     if(detailed_log.ge.2) then
        print*,'Bilan_SousGr(6,sum(Bilan_Segmelim(:,5)))=',Bilan_SousGr(6,sum(Bilan_Segmelim(:,5)))
        print*,'nbGWAS=',nbGWAS
     endif


     if(Bilan_SousGr(6,sum(Bilan_Segmelim(:,5))).ne.nbGWAS) then
        print*,' '
        print*,'PROBLEM when filling the table Bilan_SousGr :'
        print*,'Bilan_SousGr(6,sum(Bilan_Segmelim(:,5)))=',Bilan_SousGr(6,sum(Bilan_Segmelim(:,5)))
        print*,'nbGWAS=',nbGWAS
        STOP 51
     endif



     !==========================================================================================================
     !==========================================================================================================
     !==========================================================================================================





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     ! on relit le fichier DOSES pour stocker les doses UTILES dans un fichier ordonné par numero animal trie

     ! ETAPE PRELIMINAIRE : on va eventuellement lire le fichier doses en plusieurs tranches de variants pour limiter les besoins en memoire 
     !                      pour la table DOSES



     nDutLues = MaxMem / ( dfloat(nbind)*4.0d0 )
     !nDutLues=35 ! pour test 28/03/2024
     if (nDutLues.ge.nbGWAS) nDutLues=nbGWAS
     if(detailed_log.eq.3) print*,'Number of Variants nDutLues=',nDutLues

     nb_tranches = INT( nDutLues  / dfloat(NbMaxVarGWAS)  )   ! nb_tranches = nombre de fichiers doses de NbMaxVarGWAS qu on va pouvoir remplir a chaque tour de lecture du fichier doses_varGWAS
     if(detailed_log.eq.3) print*,'nb_tranches=',nb_tranches

     DosUtLues = nb_tranches * NbMaxVarGWAS ! = nombre de doses utiles lues a chaque tour de lecture du fichier doses_varGWAS

     nb_lectures=0
     if (DosUtLues.ne.0) then
        nb_lectures = CEILING( dfloat(nbGWAS) / (dfloat(DosUtLues)) )
     endif

     if (nDutLues.ge.nbGWAS) then
        DosUtLues=nbGWAS
        nb_lectures=1
     endif

     if(detailed_log.ge.2) print*,'DosUtLues=',DosUtLues
     if(detailed_log.ge.2) print*,'nb_lectures=',nb_lectures


     allocate(cadre_lec_doses(2,nb_lectures))
     cadre_lec_doses=0

     if(nb_lectures.gt.1) then
        do lecture=1,nb_lectures - 1
           cadre_lec_doses(1,lecture) = ((lecture-1)*DosUtLues) + 1
           cadre_lec_doses(2,lecture) = lecture*DosUtLues
        enddo
     endif
     cadre_lec_doses(1,nb_lectures) = ((nb_lectures-1)*DosUtLues) + 1
     cadre_lec_doses(2,nb_lectures) = nbGWAS

     if(detailed_log.ge.2) then
        print*,' '
        print*,'Table cadre_lec_doses :'
        do i=1,nb_lectures
           print*,'Read ',i,' ',cadre_lec_doses(:,i)
        enddo
        print*,' '
     endif


     allocate(DOSES(nbind,DosUtLues))
     !DOSES=0.0 on n initialise pas la matrice DOSES a 0.0 ici car elle sera systematiquement re-initialisee a chaque debut de boucle lecture=1,nb_lectures
     ! et si la matrice est grosse cette initialisation a 0 prend enormement de temps !!!


     if(trim(fmttypgwas).eq.'plink') then
        if(detailed_log.eq.3) then
           print*,' '
           print*,'Before allocation of vector vectemp5 firstVarGWAS=',firstVarGWAS
           print*,' '
        endif

        allocate(vectemp5(firstVarGWAS))
        vectemp5=0
     endif

     if((trim(fmttypgwas).eq.'typ_eval').or.(trim(fmttypgwas).eq.'plink')) then
        allocate(freqTG(nbGWAS),nballelesTG(nbGWAS),mostfreqTG(nbGWAS))
     endif
     freqTG=0.0d0
     nballelesTG=0.0d0
     mostfreqTG=0.0d0



     ! on cree un sous repertoire par groupe d exclusion de SNPparente dans lequel on mettra le fichier binaire des doses de ces VarGWAS et le fichier parametres pour partie 3
     do i=1,nbSegmExcl(chrom_anal)
        write(fx1,'(i0)') chrom_anal
        write(fx2,'(i0)') i
        do j=1,Bilan_SegmElim(i,5)
           write(fx22,'(i0)') j

           if(indic_cojo.eq.0) then
              rep_G_SG=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
           endif
           if(indic_cojo.eq.1) then
              write(fxcojo,'(i0)') numvar_cojo
              rep_G_SG=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
           endif

!           print*,' '
!           print*,'AVANT instruction CALL SYSTEM / repertoire a creer : ',trim(rep_G_SG)

           call system("mkdir -p "//trim(rep_G_SG))  ! ATTENTION : il faudrait securiser en testant l existence des sous repertoires et des fichiers doses ...

!           print*,'APRES instruction CALL SYSTEM / repertoire a creer : ',trim(rep_G_SG)
!           print*,' '

        enddo
     enddo

     ! on ouvre le fichier qui va contenir la liste des groupes et sous groupes pour rassembler ulterieurement tous les fichiers RESULTATS en un seul
     !open(io_infos0,file=trim("Infos_all_Gr_SG_chr")//trim(fx1)//trim(".txt"),form="formatted")
     if(indic_cojo.eq.0) then
        open(io_infos0,file=trim("Infos_all_Gr_SG_chr")//trim(fx1)//trim(".txt"),form="formatted")
     endif
     if(indic_cojo.eq.1) then
        write(fxcojo,'(i0)') numvar_cojo
        call system("mkdir -p "//trim("cojo_")//trim(fxcojo))
        open(io_infos0,file=trim("cojo_")//trim(fxcojo)//trim("/Infos_all_Gr_SG_chr")//trim(fx1)//trim(".txt"),form="formatted")
     endif


     batch=1
     write(fx1,'(i0)') chrom_anal
     write(fx2,'(i0)') Bilan_SousGr(1,batch)
     write(fx22,'(i0)') Bilan_SousGr(2,batch)
     if(indic_cojo.eq.0) then
        rep_G_SG=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
     endif
     if(indic_cojo.eq.1) then
        write(fxcojo,'(i0)') numvar_cojo
        rep_G_SG=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
        chem_dosesCojo=trim("cojo_")//trim(fxcojo)//trim("/Doses_VarCojo_")//trim(fxcojo)//trim(".bin")
        open(io_dosesCojo,file=trim(chem_dosesCojo),form="unformatted")   ! si COJO on ouvre le fichier qui va recevoir les doses du variant cojo
     endif

     chem_bin0=trim(rep_G_SG)//trim("/Doses_GrVar.bin")
     chem_infos1=trim(rep_G_SG)//trim("/Infos_GrVar.txt")
     open(io_dosesbin0,file=trim(chem_bin0),form="unformatted")
     open(io_infos1,file=trim(chem_infos1),form="formatted")
     write(io_infos1,*) Bilan_SousGr(1:4,batch),Bilan_SousGr(7,batch)
     write(io_infos0,*) chrom_anal,Bilan_SousGr(1:4,batch),Bilan_SousGr(7,batch)



     !=============================================================================================================================================================
     !=============================================================================================================================================================
     !=============================================================================================================================================================
     ! DEBUT DE LA BOUCLE SUR LES nb_lectures LECTURES DU FICHIER DOSES_VARGWAS

     do lecture=1,nb_lectures

        DOSES=0.0

        nldoses=0

        if(lecture.lt.nb_lectures) DosUtTrav=DosUtLues
        if(lecture.eq.nb_lectures) DosUtTrav= nbGWAS - (DosUtLues*(nb_lectures-1))     ! si on est dans la derniere lecture le nb de VarGWAS utiles lues est le nb de varGWAS restant

        print*,' '
        print*,'Current number of reading pass of Variant Genotypes / allelic dosages doses_varGWAS = ',lecture,'  DosUtTrav = ',DosUtTrav
        PRINT*,' '

        if(trim(fmttypgwas).eq.'minimac') then
           do
              read(io_doses,*,iostat=io) IdDosesTemp,chdose,DOSEStemp
              if (io .ne. 0) exit
              nldoses=nldoses+1
              if(DosesUtil(ordinvdoses(nldoses)).ne.0) then ! l individu est dans le fichier PERF donc ses doses sont utiles
                 indice_temp=DosesUtil(ordinvdoses(nldoses))
                 !do i=1,nbGWAS
                 if(doses_to_geno.eq.0) then  ! on va faire les GWAS sur les doses
                    do i=1,DosUtTrav
                       DOSES(indice_temp,i)= DOSEStemp(GWASuttoMAP( cadre_lec_doses(1,lecture)-1 + i ))
                    enddo
                 endif
                 if(doses_to_geno.eq.1) then  ! on va faire les GWAS sur des genotypes discrets 0 1 2  reconstitues a partir des doses
                    do i=1,DosUtTrav
                       if(DOSEStemp(GWASuttoMAP( cadre_lec_doses(1,lecture)-1 + i )).le.0.5)  DOSES(indice_temp,i)= 0.0
                       if((DOSEStemp(GWASuttoMAP( cadre_lec_doses(1,lecture)-1 + i )).gt.0.5).and.(DOSEStemp(GWASuttoMAP( cadre_lec_doses(1,lecture)-1 + i )).lt.1.5))  DOSES(indice_temp,i)= 1.0
                       if(DOSEStemp(GWASuttoMAP( cadre_lec_doses(1,lecture)-1 + i )).ge.1.5)  DOSES(indice_temp,i)= 2.0
                    enddo
                 endif

              endif
           enddo

           print*,' '
           print*,'Number of lines read in Variant Genotype / allelic dosage file during this reading pass : ',nldoses
           print*,' '
        endif

        !======================================================================================================================================

        if((trim(fmttypgwas).eq.'typ_eval').or.(trim(fmttypgwas).eq.'plink')) then

           ! on lit les typages une 1ere fois pour determiner le typage le plus frequent pour remplacer les eventuels typages inconnus par le typage le plus frequent
           !allocate(freqTG(DosUtLues),nballelesTG(DosUtLues),mostfreqTG(DosUtLues))
           freqTG=0.0d0
           nballelesTG=0.0d0
           nbTGut=0
           mostfreqTG=1.0d0

           rewind(io_doses)

           if(trim(fmttypgwas).eq.'plink') then
              !allocate(vectemp5(firstVarGWAS)) ! vecteur temporaire qui va recevoir pour la ligne lue du fichier typages VarGWAS plink le phenotype et les genotypes aux variants GWAS des chromosomes precedant NUMCHR_GWAS
              vectemp5=0
              read(io_doses,*,iostat=io) ! on lit la 1ere ligne du fichier typages varGWAS de plink = en tetes
              if (io .ne. 0) exit
           endif

           do
              if(trim(fmttypgwas).eq.'typ_eval') then
                 read(io_doses,informatT,iostat=io) NoNatTemp,TypGwasTemp
                 if (io .ne. 0) exit
              endif
              if(trim(fmttypgwas).eq.'plink') then
                 read(io_doses,*,iostat=io) NoNatTemp2,NoNatTemp,SIRE_temp,DAM_temp,SEX_temp,vectemp5,TypGwasTemp ! inutile de lire apres la fin du chromosome etudie
                 if (io .ne. 0) exit
              endif


              nldoses=nldoses+1
              if(DosesUtil(ordinvdoses(nldoses)).ne.0) then ! l individu est dans le fichier PERF donc ses doses sont utiles
                 indice_temp=DosesUtil(ordinvdoses(nldoses))
                 !do i=1,nbGWAS
                 do i=1,DosUtTrav
                    j = cadre_lec_doses(1,lecture)-1 + i
                    if (TypGwasTemp(GWASuttoMAP(j)).eq.1) freqTG(j)=freqTG(j)+1d0
                    if (TypGwasTemp(GWASuttoMAP(j)).eq.2) freqTG(j)=freqTG(j)+2d0
                    if (TypGwasTemp(GWASuttoMAP(j)).ne.5) nballelesTG(j)=nballelesTG(j)+2d0
                 enddo
                 nbTGut=nbTGut+1
              endif
           enddo

           if(detailed_log.eq.3) then
              print*,'Number of lines read in Variant Genotype / allelic dosage file when calculating mostfreqTG : ',nldoses
           endif

           rewind(io_doses)

           !do i=1,nbGWAS
           do i=1,DosUtTrav
              j = cadre_lec_doses(1,lecture)-1 + i
              freqTG(j)=freqTG(j)/nballelesTG(j)
              if(freqTG(j).lt.(1.0d0/3.0d0))  mostfreqTG(j)=0 ! si freq(allele2) < 0.3333 gzenotype le plus frequent = homozygote 11 
              if(freqTG(j).gt.(2.0d0/3.0d0))  mostfreqTG(j)=2 ! si freq(allele2) > 0.6666 gzenotype le plus frequent = homozygote 22 
           enddo

           if(detailed_log.eq.3) then
              do i=1,10
                 j = cadre_lec_doses(1,lecture)-1 + i
                 print*,'freqTG variant',j,' =',freqTG(j),'  /  mostfreqTG(j) =',mostfreqTG(j)
              enddo
           endif

           if(trim(fmttypgwas).eq.'plink') then
              read(io_doses,*,iostat=io) ! on lit la 1ere ligne du fichier typages varGWAS de plink = en tetes
              if (io .ne. 0) exit
           endif

           nldoses=0

           do
              if(trim(fmttypgwas).eq.'typ_eval') then
                 read(io_doses,informatT,iostat=io) NoNatTemp,TypGwasTemp
                 if (io .ne. 0) exit
              endif
              if(trim(fmttypgwas).eq.'plink') then
                 read(io_doses,*,iostat=io) NoNatTemp2,NoNatTemp,SIRE_temp,DAM_temp,SEX_temp,vectemp5,TypGwasTemp ! inutile de lire apres la fin du chromosome etudie
                 if (io .ne. 0) exit
              endif

              nldoses=nldoses+1
              if(DosesUtil(ordinvdoses(nldoses)).ne.0) then ! l individu est dans le fichier PERF donc ses doses sont utiles
                 indice_temp=DosesUtil(ordinvdoses(nldoses))
                 !do i=1,nbGWAS
                 do i=1,DosUtTrav
                    j = cadre_lec_doses(1,lecture)-1 + i
                    if(TypGwasTemp(GWASuttoMAP(j)).eq.5) TypGwasTemp(GWASuttoMAP(j))=mostfreqTG(GWASuttoMAP(j))
                    !DOSES(indice_temp,i)= dfloat(TypGwasTemp(GWASuttoMAP(j)))  ! correction du 28/07/2023 car DOSES est REAL(kind=4) --> float au lieu de dfloat
                    DOSES(indice_temp,i)= float(TypGwasTemp(GWASuttoMAP(j)))
                 enddo
              endif
           enddo

           if(detailed_log.ge.1) then
              print*,' '
              print*,'Number of lines read in Variant Genotype / allelic dosage file : ',nldoses
              print*,' '
           endif

        endif

        if(detailed_log.eq.3) then
           print*,' '
           print*,'Allelic dosages or Genotypes for the first 100 GWAS Variants and the first 10 sorted individuals :'
           do i=1,min(10,nbind)
              print*,DOSES(i,1:min(100,DosUtLues))
           enddo
        endif

        !======================================================================================================================================

        rewind(io_doses)

        if(nldoses.ne.nltyppar) then
           print*,'For information : numbers of line read in Variant Genotype / allelic dosage file ',nldoses,' and in Marker Genotype file ',nltyppar,' are different'
           !STOP 30
        endif

        !print*,'Table DOSES :'
        !do i=1,nldoses
        !   print*,'indiv ',i,' DOSES = ',DOSES(i,:)
        !enddo

        !close(io_doses)



        !==========================================================================================================
        !==========================================================================================================
        !==========================================================================================================

        !if(indic_cojo.eq.1) then
        !   write(fxcojo,'(i0)') numvar_cojo
        !   chem_dosesCojo=trim("cojo_")//trim(fxcojo)//trim("/Doses_VarCojo_")//trim(fxcojo)//trim(".bin")
        !   open(io_dosesCojo,file=trim(chem_dosesCojo),form="unformatted")   ! si COJO on ouvre le fichier qui va recevoir les doses du variant cojo
        !endif

        do j=1,DosUtTrav
           i=cadre_lec_doses(1,lecture)-1 + j

           if((i).eq.(Bilan_SousGr(6,batch)+1)) then ! on est passe au sous groupe suivant
              close(io_dosesbin0)
              close(io_infos1)

              batch=batch+1
              if(batch.gt.nbSGr) exit
              write(fx2,'(i0)') Bilan_SousGr(1,batch)
              write(fx22,'(i0)') Bilan_SousGr(2,batch)

              if(indic_cojo.eq.0) then
                 rep_G_SG=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
              endif

              if(indic_cojo.eq.1) then
                 write(fxcojo,'(i0)') numvar_cojo
                 rep_G_SG=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)
                 !chem_dosesCojo=trim("cojo_")//trim(fxcojo)//trim("/Doses_VarCojo_")//trim(fxcojo)//trim(".bin")
                 !open(io_dosesCojo,file=trim(chem_dosesCojo),form="unformatted")   ! si COJO on ouvre le fichier qui va recevoir les doses du variant cojo
              endif

              chem_bin0=trim(rep_G_SG)//trim("/Doses_GrVar.bin")
              chem_infos1=trim(rep_G_SG)//trim("/Infos_GrVar.txt")

              open(io_dosesbin0,file=trim(chem_bin0),form="unformatted")

              open(io_infos1,file=trim(chem_infos1),form="formatted")
              write(io_infos1,*) Bilan_SousGr(1:4,batch),Bilan_SousGr(7,batch)
              write(io_infos0,*) chrom_anal,Bilan_SousGr(1:4,batch),Bilan_SousGr(7,batch)

           endif

!!!!!!!!!!!!!!!!!!!!!!!!!  SELON VALEUR DE GWAS_TYPE ECRIRE DOSES OU COVAR_M
           if((gwastype.eq.'add').or.(gwastype.eq.'add_dom'))  write(io_dosesbin0) i ,GWASuttoMAP(i) ,SegmElim(:,i) , DOSES(:,j)

           ! si le variant est le variant COJO alors on ecrit les doses pour ce variant dans un fichier dedie qui sera lu en partie 3
           if((indic_cojo.eq.1).and.(GWASuttoMAP(i).eq.numvar_cojo)) then
!!!!!!!!!!!!!!!!!!!!!!!!  SELON VALEUR DE GWAS_TYPE ECRIRE DOSES OU COVAR_M
              if((gwastype.eq.'add').or.(gwastype.eq.'add_dom')) write(io_dosesCojo) i,GWASuttoMAP(i),SegmElim(:,i), DOSES(:,j)
              print*,' '
              print*,'Genotypes / Allelic Dosages for Cojo Variant ',numvar_cojo,GWASuttoMAP(i),' written in file ',trim("cojo_")//trim(fxcojo)//trim("/Doses_VarCojo_")//trim(fxcojo)//trim(".bin")

           endif

        enddo

        rewind(io_doses) ! on se replace au debut du fichier typages GWAS pour la lecture suivante

        !if(trim(fmttypgwas).eq.'plink') deallocate(vectemp5)

     enddo ! fin de la boucle sur les nb_lectures lectures du fichier doses_varGWAS

     close(io_dosesbin0)
     close(io_infos1)
     close(io_infos0)
     close(io_doses)

     if(indic_cojo.eq.1) close(io_dosesCojo)

     ! la table DOSES(,) est desormais inutile --> on la desalloue
     deallocate(DOSES)                                                 


     !====================================================================================================================
     !====================================================================================================================
     !====================================================================================================================


     ! on ecrit dans chaque sous-repertoire du SousGroupe le fichier parametres correspondant a l execution de la PARTIE 3 de maniere autonome

     do i=1,nbSegmExcl(chrom_anal)
        do j=1,Bilan_SegmElim(i,5)
           call crea_fichpar_p3(chrom_anal,i,j)
        enddo
     enddo

     ! on cree le fichier contenant les lignes de commandes de lancements de la partie 3 pour chaque groupe de variants GWAS

     if(indic_cojo.eq.0) then
        open(io_dosesbin0,file=trim("lance_part3_allGR_chr")//trim(fx1)//trim(".sh"),form="formatted")
     endif
     if(indic_cojo.eq.1) then
        write(fxcojo,'(i0)') numvar_cojo
        open(io_dosesbin0,file=trim("lance_part3_allGR_chr")//trim(fx1)//trim("cojo_")//trim(fxcojo)//trim(".sh"),form="formatted")
     endif

     write(fx1,'(i0)') chrom_anal
     do i=1,nbSegmExcl(chrom_anal)
        write(fx2,'(i0)') i
        do j=1,Bilan_SegmElim(i,5)
           write(fx22,'(i0)') j
           if(indic_cojo.eq.0) suppl_cojo=''
           if(indic_cojo.eq.1) suppl_cojo=trim("cojo_")//trim(fxcojo)//trim("/")
           write(io_dosesbin0,"(a)") adjustl(trim(execpath)//" "//trim(suppl_cojo)//trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/fichpar_part3_GrVar.par > ")//" "//trim(suppl_cojo)//trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/SNP_GWAS_chr")//trim(fx1)//trim("_Gr")//trim(fx2)//trim(".log"))
        enddo
     enddo
     close(io_dosesbin0)
     if(indic_cojo.eq.0) call system(trim("chmod u+x lance_part3_allGR_chr")//trim(fx1)//trim(".sh"))
     if(indic_cojo.eq.1) call system(trim("chmod u+x lance_part3_allGR_chr")//trim(fx1)//trim("cojo_")//trim(fxcojo)//trim(".sh"))




     print*,'STEP_2 completed'
  endif ! fin du test if steptodo = 2 


  !==========================================================================================================
  !==========================================================================================================
  !==========================================================================================================

  if(steptodo.eq.3) then

     iobilan=0

!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! preparation des formats pour ecrire dans les fichiers resultats pour les differents modeles
     if(neff_GWAS_lu.eq.1) then
        if(CalcVarRes.eq.'optim') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif

        if(CalcVarRes.eq.'approx') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif
        if(CalcVarRes.eq.'exact') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif
     endif

     if(neff_GWAS_lu.eq.2) then
        if(CalcVarRes.eq.'optim') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif

        if(CalcVarRes.eq.'approx') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif
        if(CalcVarRes.eq.'exact') then
           write(format_res,'(a)') '(i8,1x,i8,1x,E20.12E3,1x,E20.12E3,1x,i1,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3,1x,E20.12E3)'
        endif
     endif


     print*,'STEP_3 STARTS'

     ! si COJO alors on commence par lire le fichier doses du variant COJO pour connaitre l identite du variant et connaitre le nom du repertoire contenant les fichiers a lire
     if(indic_cojo.eq.1) then
        if(gwastype.eq.'add') then
           allocate(vecDosesSD_cojo(nbind,1))
           vecDosesSD_cojo=0.0
        endif
        if(gwastype.eq.'add_dom') then
           allocate(vecDosesSD_cojo(nbind,2))
           vecDosesSD_cojo=0.0
        endif

        if(detailed_log.eq.3) print*,'Allocations COJO passees'

        print*,'cojofile=',trim(cojofile)

        open(io_dosesCojo,file=trim(cojofile),form="unformatted")

        if((GWAStype.eq.'add').or.(GWAStype.eq.'add_dom')) read(io_dosesCojo) j,numvar_cojo,temp3(1:3),vecDosesSD_cojo(1:nbind,1)

        if(GWAStype.eq.'add_dom') vecDosesSD_cojo(1:nbind,2) = (ABS((ABS(vecDosesSD_cojo(1:nbind,1)-1.0d0))-1.0d0))

        if(detailed_log.eq.3) print*,'Ligne  du fichier Doses pour variant COJO lue - numvar_cojo=',numvar_cojo

        close(io_dosesCojo)

        if(gwastype.eq.'add') then
           neff_cojo=1
           neq_cojo=1
        endif
        if(gwastype.eq.'add_dom') then
           neff_cojo=2
           neq_cojo=2
        endif

        allocate(matNivAnim_cojo(neff_cojo,1,nbind),matWeight_cov_cojo(neff_cojo,nbind))
        matNivAnim_cojo=0
        matWeight_cov_cojo=0.0d0
        if (gwastype.eq.'add') then
           matNivAnim_cojo(1,1,:)=1
           matWeight_cov_cojo(1,:)=vecDosesSD_cojo(:,1)*(sqrt(matPOIDS(:,1)))
        endif
        if (gwastype.eq.'add_dom') then
           matNivAnim_cojo(1,1,:)=1  ! effet additif
           matWeight_cov_cojo(1,:)=vecDosesSD_cojo(:,1)*(sqrt(matPOIDS(:,1)))  ! effet additif
           matNivAnim_cojo(2,1,:)=1  ! effet de dominance
           matWeight_cov_cojo(2,:)=vecDosesSD_cojo(:,2)*(sqrt(matPOIDS(:,1)))  ! effet de dominance
        endif

     endif  ! fin du test if indic_cojo=1


     ! on lit le fichier contenant les informations pour le Groupe de VarGWAS consideres qui seront utiles pour l execution de la partie 3

     write(fx1,'(i0)') chrom_p3
     write(fx2,'(i0)') numgr_p3
     write(fx22,'(i0)') numSgr_p3
     if(indic_cojo.eq.1) write(fxcojo,'(i0)') numvar_cojo

     if(indic_cojo.eq.0)  open(io_dosesbin0,file=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Infos_GrVar.txt"),form="formatted")
     if(indic_cojo.eq.1)  open(io_dosesbin0,file=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Infos_GrVar.txt"),form="formatted")

     if(indic_cojo.eq.0) print*,'Information about Group x Batch in file ',trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Infos_GrVar.txt")
     if(indic_cojo.eq.1) print*,'Information about Group x Batch in file ',trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Infos_GrVar.txt")

     do
        read(io_dosesbin0,*,iostat=iobilan) NumGr,NumSGr,Felim,Lelim,nbGWASgr
        if (iobilan .ne. 0) exit
        print*,'Group =',NumGr,' / Batch =',NumSGr,' /  1st Marker excluded =',Felim,' /  last Marker excluded =',Lelim, ' /  number of Variants in the Batch =',nbGWASgr
     enddo
     close(io_dosesbin0)

     if(numgr_p3.ne.NumGR) then
        print*,'PROBLEM : Group Id read in parameter file for Step_3 is different from Group Id read in file Infos_GrVar.txt'
        STOP 21
     endif
     if(numSgr_p3.ne.NumSGr) then
        print*,'PROBLEM : Batch Id read in parameter file for Step_3 is different from Batch Id read in file Infos_GrVar.txt'
        STOP 21
     endif


     allocate(numVarGWASori(nbGWASgr),GWASuttoMAP_p3(nbGWASgr),SegmElim_p3(3,nbGWASgr))
     numVarGWASori=0
     GWASuttoMAP_p3=0
     SegmElim_p3=0

     ! on relit le fichier invLEFT_binfile pour re-remplir inf_LEFT_EP et executer la suite du programme

     allocate (inf_LEFT_EP(dimvec))
     inf_LEFT_EP=0.0d0

     if(detailed_log.eq.3) then
        call fdate(jour)
        print*,jour
        print*,'Avant lecture de inf_LEFT_EP sur disque'
        print*,' '
     endif

     open(io_invleft,file=trim(location)//trim("/invLEFT_binfile"),form="unformatted") 
     read(io_invleft)inf_LEFT_EP
     close(io_invleft)

     call fdate(jour)
     print*,jour
     print*,'Reading of inv(inf_LEFT_EP) completed'
     print*,' '

     ! on relit le fichier RIGHT_binfile sur disque pour re-remplir RIGHT_EP et executer la suite du programme

     allocate(RIGHT_EP(neq_ep))
     RIGHT_EP=0.0d0

     if(detailed_log.eq.3) then
        print*,' '
        call fdate(jour)
        print*,jour
        print*,'Avant lecture de RIGHT_EP sur disque'
        print*,' '
     endif

     open(io_right,file=trim(location)//trim("/RIGHT_binfile"),form="unformatted") 
     read(io_right)RIGHT_EP
     close(io_right)

     print*,' '
     call fdate(jour)
     print*,jour
     print*,'Reading of RIGHT_EP completed'
     print*,' '


     !==========================================================================================================
     !==========================================================================================================
     !==========================================================================================================


     !==============================================================================================!
     !                                                                                              !
     ! ETAPE PREALABLE AU CALCUL DE LA VARIANCE D ESTIMATION DES EFFETS DU VARIANT GWAS CONSIDERE   !
     !                                                                                              !
     !==============================================================================================!

     if(detailed_log.eq.3) then
        print*,' '
        call fdate(jour)
        print*,'before computing invdenom_N'
        print*,jour
        print*,' '
     endif

     ! methode 1 : on considere que les effets environnementaux (beta) et des SNPparente ne dependent pas des effets du variant GWAS considere
     !             --> on estime ces effets dans un modele sans variant GWAS en utilisant l inverse de LEFT_EP_complet

     allocate(BetaS(neq_ep+neq_cojo,1))
     BetaS=0.0d0
     ! BetaS est solution du systeme LEFT_EP_complet*BetaS = RIGHT_EP_complet
     ! ATTENTION : la subroutine dsymm a besoin d une matrice symetrique stockee sous forme triangulaire meme si cela utilise de la memoire pour rien
     allocate(inf_LEFT_EP_act(neq_ep,neq_ep))
     inf_LEFT_EP_act=0.0d0
     do j=1,neq_ep
        inf_LEFT_EP_act(j:neq_ep,j) = inf_LEFT_EP(TI(j,j,neq_ep):TI(neq_ep,j,neq_ep))
     enddo

     sumPoids=0.0d0
     do i=1,nbind
        sumPoids=sumPoids+matPOIDS(i,1)
     enddo
     if(detailed_log.eq.3) print*,'sumPoids=',sumPoids


     ! On calcule le rang de XpX
     Rank=0
     presmoy=0 ! indicateur de presence d une moyenne generale dans le modele
     nb_eff_fixcat=0
     nb_eff_cov=0

     do ef=1,neff
        if(randomtype(ef).eq.g_fixed) then     ! effet non aleatoire = effet fixe categoriel ou covariable
           if(effecttype(ef).eq.effcross) then ! effet fixe categoriel
              if(nlev(ef).eq.1) then           ! effet fixe categoriel a 1 niveau = moyenne generale 
                 Rank = Rank + 1
                 presmoy = 1
              endif
              if(nlev(ef).gt.1) then           ! effet fixe categoriel a plusieurs niveaux
                 Rank = Rank + nlev(ef) - 1
                 nb_eff_fixcat=nb_eff_fixcat+1
              endif
           endif
           if(effecttype(ef).eq.effcov) then ! covariable
              Rank = Rank + nlev(ef)    ! si covariable intra effet fixe on a plusieurs niveaux
              nb_eff_cov=nb_eff_cov+1
           endif
        endif
        if(randomtype(ef).eq.g_diag) then
           Rank = Rank + 1
           nb_eff_fixcat=nb_eff_fixcat+1
        endif
     enddo

     rank = rank + 1 + 1 + neq_cojo  ! on ajoute l effet ADDITIF pour le variant GWAS (1) et +1 pour la valeur genetique par les SNPparente + 1 si on a l effet fixe du variant COJO
     ! on cree deux variables de rang qui seront utilisees en fonction du modele EFFECTIF pour le variant GWAS (utilisateur peut demander add_dom mais uniquement add possible pour varGWAS)
     rank_1 = rank ! rang si uniquement effet additif
     rank_2 = rank + 1 ! rang si (effet additif et effet de dominance) dans modele 

     print*,' '
     print*,'Presence of an overall mean in the model (Fixed categorial effect with 1 level) (0=no / 1=yes) :',presmoy
     print*,'nb of fixed categorial effects in the environmental part of the model (apart from any overall mean) : ',nb_eff_fixcat
     print*,'nb of fixed regression coefficents in the environmental part of the model (apart from any overall mean : ',nb_eff_cov
     print*,'Presence of an additional fixed regression coefficient for the effect of an additional fixed variant with COJO OPTION (0=no / 1=yes) : ',indic_cojo
     print*,' '

     if(detailed_log.eq.3) print*,'RANK=',Rank,'  // RANK_1=',rank_1,'  // RANK_2=',rank_2

     !invdenom = 1.0d0 / (sumPoids - dfloat(rank))
     invdenom_1 = 1.0d0 / (nbcarperf(1) - dfloat(rank_1))
     invdenom_2 = 1.0d0 / (nbcarperf(1) - dfloat(rank_2))

     if(detailed_log.eq.3) then
        print*,' '
        print*,' 1/denom_1=',invdenom_1,'  // 1/denom_2=',invdenom_2
        print*,' '

        print*,' '
        call fdate(jour)
        print*,'Calculation of invdenom_1 et invdenom_2 completed'
        print*,jour
        print*,' '
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! =============================================================================================================================
     ! DEBUT DE LA PARTIE SPECIFIQUE A UN VARIANT GWAS CONSISTANT A AJUSTER inf_LEFT_EP et RIGHT_EP EN EXCLUANT LES SNPpar PROCHES
     ! =============================================================================================================================

     !allocate(inf_LEFT_EP_act(neq_ep,neq_ep))
     allocate(RIGHT_EP_act(neq_ep+neq_cojo,1))
     !allocate(XpZMpZ(neq_ep,neff_GWAS_lu))
     allocate(ZpXZpM(neff_GWAS_lu,neq_ep+neq_cojo))
     allocate(vecDosesPond(nbind,neff_GWAS_lu))
     allocate(vecDoses(nbind,neff_GWAS_lu))
     allocate(vecDosesSD(nbind))
     allocate(vecGENO1(nbind))
     allocate(ZpY(neff_GWAS_lu,1)) ! pour la construction de Z'Y pour le variant GWAS considere
     allocate(ZpZ(neff_GWAS_lu,neff_GWAS_lu)) ! pour la construction de Z'Z pour le variant GWAS considere
     allocate(interm_1(neff_GWAS_lu,neq_ep+neq_cojo)) ! = t[X'Z M'Z] LEFT_EP_act
     allocate(interm_2(neff_GWAS_lu,neff_GWAS_lu)) ! = interm_1 [X'Z M'Z]
     allocate(interm_3(neff_GWAS_lu,1)) ! = interm_1 [X'y M'y]
     allocate(interm_4(neff_GWAS_lu,neff_GWAS_lu)) ! = Z'Z - interm_2
     allocate(interm_5(neff_GWAS_lu,1)) ! = Z'Y - interm_3
     allocate(interm_6(neq_ep+neq_cojo,1)) ! = RIGHT_EP_act - t(ZpXZpM)%*%SolVarGWAS
     allocate(sol_VarGWAS(neff_GWAS_lu,nbGWASgr)) ! solution pour effet GWAS additif et dominance si present dans le modele
     allocate(vec_stat_doses(2,nbGWASgr))    ! moyenne et ecart type des doses pour chaque variant GWAS
     vec_stat_doses=0.0d0
     allocate(tem_possible(1,nbGWASgr))    ! indicateur =1 si au moins 5 homozygotes (11 ou 22) et au moins 5 heterozygotes / = 0 sinon si genotypes discrets // =1 si vec_stat_doses(2,varGWAS).ne.0
     tem_possible=0


     allocate(vecTemDom(nbGWASgr))
     vecTemDom=0


     allocate(pval(neff_GWAS_lu))

     allocate(vec_VarRes(2,nbGWASgr))            ! variance residuelle calculee pour le variant utilisee pour calculer son Test: 1=exp4_1 2=exp4_2

     ! correction 03/09/2025 car test multiple sur exact ou valeur de test_effVarGWAS_1 qui n etait pas alloue dans cas exact
     !if((CalcVarRes.eq.'optim').or.(CalcVarRes.eq.'approx')) then
        allocate(mat_TestEff_1(neff_GWAS_lu,nbGWASgr)) ! tests pour effet GWAS additif et dominance si present dans le modele avec VarRes approx
        allocate(test_effVarGWAS_1(neff_GWAS_lu))      ! avec VarRes approx
        allocate(mat_pval_1(neff_GWAS_lu,nbGWASgr))
        mat_TestEff_1=0.0d0
        mat_pval_1=0.0d0
     !endif

     if((CalcVarRes.eq.'optim').or.(CalcVarRes.eq.'exact')) then
        allocate(mat_TestEff_2(neff_GWAS_lu,nbGWASgr)) ! tests pour effet GWAS additif et dominance si present dans le modele avec VarRes exact
        allocate(test_effVarGWAS_2(neff_GWAS_lu))      ! avec VarRes exact
        allocate(mat_pval_2(neff_GWAS_lu,nbGWASgr))
        mat_TestEff_2=0.0d0
        mat_pval_2=0.0d0
     endif

     allocate(restemp(neff_GWAS_lu,1))
     allocate(BetaSnew(neq_ep+neq_cojo,1))


     sol_VarGWAS=0.0d0
     vec_VarRes=0.0d0

     nbGWASexec=nbGWASgr
     if(nbGWAS_limite.eq.1) nbGWASexec=nbGWAStest


     ! on ouvre le fichier binaire contenant les doses pour les variants qu on va considerer dans la partie 3 du programme (GWAS)  1 ligne = 1 variant GWAS

     write(fx1,'(i0)') chrom_p3
     write(fx2,'(i0)') numgr_p3
     write(fx22,'(i0)') numSGr_p3
     if(indic_cojo.eq.1) write(fxcojo,'(i0)') numvar_cojo

     if(indic_cojo.eq.0) then
        open(io_dosesbin0,file=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Doses_GrVar.bin"),form="unformatted")
        print*,'Opening file ',trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Doses_GrVar.bin")
     endif
     if(indic_cojo.eq.1) then
        open(io_dosesbin0,file=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Doses_GrVar.bin"),form="unformatted")
        print*,'Opening file ',trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Doses_GrVar.bin")
     endif


     ! =========================================================================================================================
     ! on calcule l inverse de LEFT_EP apres suppression des marqueurs a exclure par rapport au variant considere dans le GWAS
     ! =========================================================================================================================


     print*,' '
     call fdate(jour)
     print*,'Start looping on GWAS variants'
     print*,jour
     print*,' '


     DIFFELIM=1 ! indique si le variant GWAS en cours a des SNPpar a exclure differents du variant GWAS precedent (1) ou pas (0)

     ! =========================================================================================================================
     ! =========================================================================================================================
     ! =========================================================================================================================
     ! =========================================================================================================================

     ! ouverture du fichier qui contiendra les resultats
     if(indic_cojo.eq.0) then
        open(io_s,file=trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Resultats_chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim(".txt"),form="formatted")
     endif
     if(indic_cojo.eq.1) then
        open(io_s,file=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Resultats_chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim(".txt"),form="formatted")
     endif


     print*,' '
     print*,'Content of elementary result file for the Group x Batch variants'
     print*,' '

     if(CalcVarRes.eq.'optim') then
        print*,' '
        if(neff_GWAS_lu.eq.1) then
           print*,'VarGWAS / Id_MAP_GWAS / Moy(Doses) / STD(Doses) / TemDom /  sol_ADD / Test_approx_ADD / P-val_approx_ADD / Test_exact_ADD / P-val_exact_ADD / VarResid_approx / VarResid_exact / detZPZ :'
        endif
        if(neff_GWAS_lu.eq.2) then
           print*,'VarGWAS / Id_MAP_GWAS / Moy(Doses) / STD(Doses) /  TemDom /  sol_ADD / sol_DOM / Test_approx_ADD / Test_approx_DOM / P-val_approx_ADD / P-val_approx_DOM / Test_exact_ADD / Test_exact_DOM / P-val_exact_ADD / P-val_exact_DOM / VarResid_approx / VarResid_exact / detZpZ :'
        endif
        print*,'_________________________________________'
     endif

     if(CalcVarRes.eq.'approx') then
        print*,' '
        if(neff_GWAS_lu.eq.1) then
           print*,'VarGWAS / Id_MAP_GWAS / Moy(Doses) / STD(Doses)  / TemDom /  sol_ADD / Test_approx_ADD / P-val_approx_ADD / VarResid_approx / detZPZ :'
        endif
        if(neff_GWAS_lu.eq.2) then
           print*,'VarGWAS / Id_MAP_GWAS / Moy(Doses) / STD(Doses) /  TemDom /  sol_ADD / sol_DOM / Test_approx_ADD / Test_approx_DOM / P-val_approx_ADD / P-val_approx_DOM / VarResid_approx / detZPZ :'
        endif
        print*,'_________________________________________'
     endif

     if(CalcVarRes.eq.'exact') then
        print*,' '
        if(neff_GWAS_lu.eq.1) then
           print*,'VarGWAS / Id_MAP_GWAS  / Moy(Doses) / STD(Doses) /  TemDom  / sol_ADD / Test_exact_ADD / P-val_exact_ADD / VarResid_exact / detZPZ :'
        endif
        if(neff_GWAS_lu.eq.2) then
           print*,'VarGWAS / Id_MAP_GWAS  / Moy(Doses) / STD(Doses)  / TemDom /  sol_ADD / sol_DOM / Test_exact_ADD / Test_exact_DOM / P-val_exact_ADD / P-val_exact_DOM / VarResid_exact / detZPZ :'
        endif
        print*,'_________________________________________'
     endif


     CurrLine=0


     !if(neff_GWAS_lu.eq.2) then
     allocate(detZpZ(nbGWASexec)) ! determinant de la matrice ZpZ si on a un effet additif et un effet de dominance
     detZpZ=0.0d0
     !endif



     ! DEBUT DE LA BOUCLE SUR LES VARIANTS GWAS UTILES

     allocate(compttyp(3,nbGWASexec))
     compttyp=0

     do VarGWAS=1,nbGWASexec ! on effectue ceci pour chaque variant GWAS UTILE du sous groupe ou pour le nb fixe par utilisateur en OPTION

        VarDiffCojo=1 ! on initialise a 1 le temoin iniduqnat si on doit traiter le variant GWAS ou passer au suivant (si info genomique identique a celle du variant Cojo on le saute)
        !  NB : un variant en DL complet avec le variant COJO doit etre saute
        ZpXZpM=0.0d0
        vecGENO1=0.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(detailed_log.eq.3) then
           if(varGWAS.le.10) print*,'varGWAS=',varGWAS,' GWAStype=',GWAStype
        endif

        vecDoses=0.0d0
        vecDosesSD=0.0
        vecDosesPond=0.0d0

        !vecDoses(:,1)=DOSES(:,VarGWAS)
        read(io_dosesbin0) numVarGWASori(VarGWAS) , GWASuttoMAP_p3(VarGWAS) , SegmElim_p3(:,VarGWAS) , vecDosesSD(:)
        vecDoses(:,1)=vecDosesSD

        ! SI COJO : on teste si les informations genomiques du variant sont differentes des informations genomiques pour le variant COJO
        if(indic_cojo.eq.1) then 
           if(all(vecDosesSD==vecDosesSD_cojo(:,1))) VarDiffCojo=0
        endif

        ! on calcule la moyenne des Doses pour le variant, pour avoir un equivalent de freq(Al2) pour le variant GWAS --> donnera une idee de la MAF
        call STAT_NUM(vecDoses(:,1),nbind,vec_stat_doses(1,VarGWAS),vec_stat_doses(2,VarGWAS)) ! 1ere ligne de vec_stat_doses=moy(doses) 2eme ligne de vec_stat_doses=std(doses)

        ! si on travaille sur dosages alleliques --> comptage des genotypes impossible --> on regarde si l ecart type des doses est nul = monomorphe
        if((trim(fmttypgwas).eq.'minimac').and.(doses_to_geno.eq.0)) then
           if(vec_stat_doses(2,VarGWAS).gt.0) tem_possible(1,varGWAS)=1
        endif

        ! si on travaille sur des genotypes discrets --> on compte le nb d individus des genotypes 11, 12 et 22 et on renseigne tem_possible(1,varGWAS)
        if( ((trim(fmttypgwas).eq.'minimac').and.(doses_to_geno.eq.1)).or.(trim(fmttypgwas).ne.'minimac') ) then
           !print*,'on passe dans le test des genotypes discrets'
           !print*,'vecDoses(1:50,1) = ',vecDoses(1:50,1)
           do i=1,nbind
              if(vecDoses(i,1).eq.0.0d0) compttyp(1,varGWAS)=compttyp(1,varGWAS)+1
              if(vecDoses(i,1).eq.1.0d0) compttyp(2,varGWAS)=compttyp(2,varGWAS)+1
              if(vecDoses(i,1).eq.2.0d0) compttyp(3,varGWAS)=compttyp(3,varGWAS)+1
           enddo
           if( ((compttyp(1,varGWAS).ge.5).or.(compttyp(3,varGWAS).ge.5)) .and. (compttyp(2,varGWAS).ge.5) ) tem_possible(1,varGWAS)=1
        endif

        neff_GWAS=1
        TemDom=0
        coldom=0

        if(neff_GWAS_lu.eq.1) invdenom=invdenom_1 ! l utilisateur a demande un modele avec uniquement effet additif des variants GWAS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 21 avr 2024 : on modifie la partie ci dessous
        ! ! on cherche la dose min et la dose max pour le variant pour déterminer si on peut considerer un effet de dominance (si demande par utilisateur) dans le modele pour le variant
        ! if(TemDom_lu.eq.1) then
        ! !if(neff_GWAS_lu.eq.2) then
        !    dose_max = maxval(vecDoses(:,1))
        !    dose_min = minval(vecDoses(:,1))
        !    if((dose_min.gt.0.5).or.(dose_max.lt.1.5)) then ! equivalent a n avoir que 1 ou 2 genotypes sur 3 representes donc pas d effet dominance estimable
        !       TemDom=0 ! pas d effet dominance utilisable dans le modele pour le variant
        !       vecTemDom(VarGWAS)=TemDom
        !
        !          neff_GWAS=1
        !          invdenom=invdenom_1
        !    endif
        !    if((dose_min.le.0.5).and.(dose_max.ge.1.5)) then
        !       TemDom=1 ! effet dominance utilisable dans le modele pour le variant
        !       vecTemDom(VarGWAS)=TemDom
        !          neff_GWAS=2
        !          coldom=2
        !          invdenom=invdenom_2 ! effet additif et effet de dominance
        !    endif
        ! endif
        ! !if(VarGWAS.le.100) then
        ! !   print*,'VarGWAS = ',VarGWAS,' dose_min = ',dose_min,' dose_max = ',dose_max,' --> neff_GWAS effectif = ',neff_GWAS,'   invdenom=',invdenom
        ! !endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! FIN DE L ANCIENNE PARTIE MODIFIEE

        ! NOUVELLE PARTIE REMPLACANT LA PARTIE CI DESSUS
        ! on cherche la dose min et la dose max pour le variant pour déterminer si on peut considerer un effet de dominance (si demande par utilisateur) dans le modele pour le variant
        if(TemDom_lu.eq.1) then
           !if(neff_GWAS_lu.eq.2) then
           dose_max = maxval(vecDoses(:,1))
           dose_min = minval(vecDoses(:,1))
           if((dose_min.gt.0.5).or.(dose_max.lt.1.5)) then ! equivalent a n avoir que 1 ou 2 genotypes sur 3 representes donc pas d effet dominance estimable
              TemDom=0 ! pas d effet dominance utilisable dans le modele pour le variant
              vecTemDom(VarGWAS)=TemDom
           endif
           if((dose_min.le.0.5).and.(dose_max.ge.1.5)) then
              TemDom=1 ! effet dominance utilisable dans le modele pour le variant
              vecTemDom(VarGWAS)=TemDom
              neff_GWAS=2 
           endif
        endif

        if(neff_GWAS.eq.1) invdenom=invdenom_1
        if(neff_GWAS.eq.2) invdenom=invdenom_2

        if(TemDom.eq.1) then 
           coldom=2
        endif



        !if(neff_GWAS.eq.2) then
        if(TemDom.eq.1) then
           vecDoses(:,coldom)= (ABS((ABS(vecDoses(:,1)-1.0d0))-1.0d0)) ! Formule correcte mais inutilisable sur petit exemple car Dose1=Dose2 pour tous les animaux
           !vecDoses(:,2)= (ABS(vecDoses(:,1)-1.0d0))  ! pour test sur petit exemple
           vecDosesPond(:,coldom)= vecDoses(:,coldom) * (sqrt(matPOIDS(:,1))) ! on multiplie la dose_Dom d un animal par son (poids)^0.5 pour mettre a 0 les doses des indiv sans perf
        endif

        vecDosesPond(:,1)=vecDoses(:,1)*(sqrt(matPOIDS(:,1))) ! on multiplie la dose_Add d un animal par son (poids)^0.5 maintenant qu on a calcule vecDosesPondDom --> un animal sans perf aura dose=0

        if(detailed_log.eq.3) then
           if(neq_ep.lt.50) then
              print*,' '
              print*,'VarGWAS ',varGWAS,' vecDosesPond = '
              do i=1,nbind
                 print*,vecDosesPond(i,:)
              enddo
           endif
        endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if((SegmElim_P3(2,VarGWAS).ne.0).and.(SegmElim_P3(1,VarGWAS).ne.0)) then
           nEXCL = SegmElim_P3(2,VarGWAS)-SegmElim_P3(1,VarGWAS)+1   ! nombre de SNP PARENTE a exclure pour le variant GWAS considere
        else
           nEXCL=0
        endif


        nEXCLmat=nEXCL*(nEXCL+1)/2
        !print*,'Le nb d elements dans le vecteur matEXCL est :',nEXCLmat

        if(VarGWAS.gt.1) then
           DIFFELIM=1
           if((SegmElim_P3(1,VarGWAS).eq.SegmElim_P3(1,VarGWAS-1)).and.(SegmElim_P3(2,VarGWAS).eq.SegmElim_P3(2,VarGWAS-1))) then
              DIFFELIM=0
           endif
        endif

        !print*,' '
        !print*,'Variant GWAS ',varGWAS,' /  DIFFELIM=',DIFFELIM,' / SegmElim_P3=',SegmElim_P3(:,VarGWAS)
        !print*,' '


        ! ===============================================================================================
        ! 1er cas : il n y a aucun SNPpar a exclure --> on peut directement utiliser inf_LEFT_EP complet
        ! ===============================================================================================


        if((varGWAS.eq.1).or.(DIFFELIM.eq.1)) then ! pour le 1er variant GWAS ou si le variant GWAS a des SNPpar a exclure differents de ceux du VarGWAS precedent
           ! il faut reinitialiser inf_LEFT_EP_act et RIGHT_EP_act et eventuellemnt les ajuster si nEXCL different de 0

           if(detailed_log.eq.3) then
              print*,' '
              call fdate(jour)
              print*,'Before initializing inf_LEFT_EP_act'
              print*,jour
              print*,' '
           endif

           inf_left_EP_act=0.0d0
           do j=1,neq_ep
              inf_LEFT_EP_act(j:neq_ep,j) = inf_LEFT_EP(TI(j,j,neq_ep):TI(neq_ep,j,neq_ep))
           enddo

           print*,' '
           call fdate(jour)
           print*,'inf_LEFT_EP_act initialized'
           print*,jour
           print*,' '

           RIGHT_EP_act=0.0d0
           RIGHT_EP_act(1:neq_ep,1) = RIGHT_EP ! on travaille sur une copie de RIGHT_EP pour conserver l original

           print*,' '
           call fdate(jour)
           print*,'RIGHT_EP_act initialized'
           print*,jour
           print*,' '

           if(detailed_log.eq.3) then
              if(neq_ep.le.20) then
                 print*,'RIGHT_EP as read in the binary file :'
                 print*,RIGHT_EP_act(1:neq_ep,1)
              endif
           endif



           ! =======================================================================================================================================================
           ! 1er cas : aucun SNPpr n est a exclure : on utilise directement l inverse de LEFT_EP complet et RIGHT_EP complet
           ! ======================================================================================================================================================

           if(nEXCL.eq.0) then
              print*,' '
              print*,'GWAS Variant number',VarGWAS,' no Marker to be excluded'
              print*,' '
           endif


           ! =======================================================================================================================================================
           ! 2eme cas : il y a des SNPpar a exclure --> on ajuste inv(LEFT_EP) pour les SNPpar a supprimer par la methode de l inverse d une matrice partitionnee
           ! =======================================================================================================================================================

           ! ======================================================================================================================================
           ! ETAPE 1 : on construit la matrice matEXCL symmetrique du bloc SNPparDL*SNPparDL des SNPpar a exclure de (LEFT_EP complet)-1
           !           et on l inverse           resultat de l inversion = matEXCL

           ! ON NE FAIT CECI QUE SI LES SNPpar A EXCLURE SONT DIFFERENTS DES SNPpar A EXCLURE POUR LE VARIANT GWAS PRECEDENT
           ! ---------------------------------------------------------------------------------------------------------------

           ! ======================================================================================================================================

           if(nEXCL.ne.0) then

              ! NB : if COJO variant fit in the model, COJO equations are AFTER MARKER EQUATIONS --> Equations for Markers to be removed ARE NOT IMPACTED BY COJO VARIANT

              PremExcl = FIRSTNIV_EP(neff+1) + SegmElim_P3(1,VarGWAS) -1  
              DernExcl = FIRSTNIV_EP(neff+1) + SegmElim_P3(2,VarGWAS) -1  

              print*,' '
              print*,'Equations for Markers number ',SegmElim_P3(1,VarGWAS),' to ',SegmElim_P3(2,VarGWAS),' will be removed'
              print*,'Equation numbers of 1st Marker to be excluded PremExcl = ',PremExcl,'    of Last Marker to be excluded DernExcl =',DernExcl
              print*,' '

              allocate(matEXCL(nEXCL,nEXCL)) ! vecteur qui va contenir la diag inf du bloc de l inverse de LEFT_EP SNPparEXCL*SNPparEXCL puis son inverse
              matEXCL=0.0d0

              do i=1,nEXCL
                 do j=1,i
                    matEXCL(i,j)=inf_LEFT_EP_act(PremExcl+i-1,PremExcl+j-1)
                 enddo
              enddo

              print*,' '
              print*,'Inverting block of Markers in LEFT_MB to be excluded'

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Start of inversion of bloc of Markers in LEFT_MB to be excluded with invmatlapack_inf'
                 print*,jour
                 print*,' '
              endif

              ! on inverse le bloc des SNP PARENTE a exclure
              call invmatlapack_ti(matEXCL,nEXCL)

              print*,' '
              call fdate(jour)
              print*,'Inversion of bloc of Markers in LEFT_MB to be excluded completed'
              print*,jour
              print*,' '


              ! on remplit l inverse de matEXCL car la fonction dgemm a besoin d une matrice pleine 
              do i=1,nEXCL
                 do j=1,i
                    matEXCL(j,i)=matEXCL(i,j)
                 enddo
              enddo


              ! ==========================================================================
              ! ETAPE 2_0 : on cree la matrice pleine A12 pour reduire temps de calculs !
              allocate(A21(nEXCL,neq_ep))
              A21=0.0d0
              do j=1,PremExcl-1
                 A21(1:nEXCL,j)=inf_LEFT_EP_act(PremExcl:DernExcl,j)
              enddo
              do j=DernExcl+1,neq_ep
                 i=0
                 do ii=PremExcl,DernExcl
                    i=i+1
                    A21(i,j)=inf_LEFT_EP_act(j,ii)
                 enddo
              enddo

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Matrix A21 built'
                 print*,jour
                 print*,' '
              endif


              ! ======================================================================================================================================
              ! ETAPE 2 : on met a 0 dans inv(LEFT_EP) tous les elements relatifs aux SNP_PARENTE a exclure = bloc 11 de inverse(LEFT_EP)
              ! ======================================================================================================================================

              ! ATTENTION : pour utiliser la fonction de multiplication de matrices dsymm() il faut inf_LEFT_EP_act sous forme de matrice triang inf et pas de vecteur
              ! --> on remplit la triangulaire inferieure de la matrice inf_LEFT_EP_act ; perte de memoire, mais pour l instant pas d autre solution rapide
              ! car la fonction dsymm doit etre optimisee et probablement beaucoup plus rapide qu une fonction maison
              inf_LEFT_EP_act=0.0d0
              do j=1,neq_ep
                 inf_LEFT_EP_act(j:neq_ep,j) = inf_LEFT_EP(TI(j,j,neq_ep):TI(neq_ep,j,neq_ep))
              enddo

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Step 1 of inf_LEFT_EP_act filling completed'
                 print*,jour
                 print*,' '
              endif


              do i=SegmElim_P3(1,VarGWAS),SegmElim_P3(2,VarGWAS) ! 
                 do j=1,neq_ep
                    excl1=i+FIRSTNIV_EP(neff+1)-1
                    excl2=j
                    if(j.gt.i+FIRSTNIV_EP(neff+1)-1) then   ! devient inutile car on ne remplit que la triang infer de inf_LEFT_EP_act ; on garde pour tester
                       k=excl1
                       excl1=j
                       excl2=k
                    endif
                    !print*,'Element ',excl1,',',excl2,' mis a 0'
                    inf_LEFT_EP_act(excl1,excl2) = 0.0d0
                 enddo
              enddo

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Step 2 of inf_LEFT_EP_act filling completed'
                 print*,jour
                 print*,' '
              endif


              ! ======================================================================================================================================
              ! ETAPE 3 : on calcule A12 (A22)-1 = t(A21) (A22)-1
              ! ======================================================================================================================================

              ! if(nEXCL.ne.0) then
              allocate(A12A22m1(neq_ep,nEXCL))
              A12A22m1=0.0d0

              call dgemm('T','N',neq_ep,nEXCL,nEXCL,1.0d0,A21,nEXCL,matEXCL,nEXCL,0.0d0,A12A22m1,neq_ep)

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Matrix A12A22m1 calculated with function dgemm'
                 print*,jour
                 print*,' '
              endif


              ! ======================================================================================================================================
              ! ETAPE 4 : on calcule  inf_LEFT_EP_act - ( [A12 (A22)-1] A21 )  
              ! ======================================================================================================================================

              !  if(nEXCL.ne.0) then

              call dgemm('N','N',neq_ep,neq_ep,nEXCL,-1.0d0,A12A22m1,neq_ep,A21,nEXCL,1.0d0,inf_LEFT_EP_act,neq_ep)

              if(detailed_log.eq.3) then
                 print*,' '
                 call fdate(jour)
                 print*,'Step 4 of inf_LEFT_EP_act filling completed'
                 print*,'WARNING : UPPER TRIANGLE OF inf_LEFT_EP_act INCORRECT but not used'
                 print*,jour
                 print*,' '
              endif

              ! pour validation on met a 0 la triangulaire supperieure de inf_LEFT_EP_act
              do i=2,neq_ep
                 inf_LEFT_EP_act(1:(i-1),i)=0.0d0
              enddo

              deallocate(matEXCL,A12A22m1,A21)

              if(detailed_log.eq.3) then
                 if(neq_ep.le.20) then
                    print*,' '
                    print*,'inverse of LEFT_EP after adjustment to remove the Markers to be excluded :'
                    do i=1,neq_ep
                       print*,inf_LEFT_EP_act(i,:)
                    enddo
                    print*,' '
                 endif
              endif


              !==================================================================================================
              ! on actualise la partie M'y de RIGHT_EP en supprimant la contribution des SNP_PARENTE a exclure
              !==================================================================================================


              do excl1=SegmElim_P3(1,VarGWAS),SegmElim_P3(2,VarGWAS)
                 RIGHT_EP_act((FIRSTNIV_EP(neff+1) + excl1 -1),1) = 0.0d0
              enddo

           endif    ! fin du test  if(nEXCL.ne.0)  

           print*,' '
           call fdate(jour)
           print*,'matrix inf_LEFT_EP_act and vector RIGHT_EP_act are ajusted for the Markers to be excluded'
           print*,jour
           print*,' '


           ! ======================================================================================
           ! On dispose maintenant de inf_LEFT_EP_act et de RIGHT_EP_act actualises pour les SNPpar a exclure 
           ! MAIS SI COJO ALORS IL FAUT ENCORE ACTUALISER POUR AJOUTER LES EQUATIONS POUR LE VARIANT COJO DANS INVERSE MBRE GAUCHE ET DANS MBRE DROITE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! modif 03 mars 2023 : on complete inf_LEFT_EP_act et RIGHT_EP_act avec les elements relatifs au variant fixe pour COJO

           ! ATTENTION : A PARTIR DE MAINTENANT, LE NOMBRE D EQUATIONS ENVIR+PARENTE EST neq_ep + neq_cojo
           neq_ep_2 = neq_ep + neq_cojo

           if(indic_cojo.eq.1) then

              allocate(col_cojo(neq_ep,neq_cojo),row_cojo(neq_cojo,neq_ep),A22(neq_cojo,neq_cojo),B22(neq_cojo,neq_cojo))
              col_cojo=0.0d0
              row_cojo=0.0d0 ! row_cojo constitue de dose_cojo'X puis Doses_cojo'M
              A22=0.0d0
              B22=0.0d0

              ! on remplit row_cojo pour les effets fixes, covariables et effets aleatoires diagonaux
              ! ATTENTION : si modele add_dom on a 2 covariables = 2 lignes cojo
              do ind=1,nbind
                 if(matPERF(ind,1).ne.mis) then
                    do ef1=1,neff
                       do ef2=1,neff_cojo
                          if(matNivAnim_cojo(ef2,1,ind).ne.0) then
                             row_cojo(ef2+matNivAnim_cojo(ef2,1,ind)-1,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) = row_cojo(ef2+matNivAnim_cojo(ef2,1,ind)-1,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) + (matWeight_cov(ef1,ind)*matWeight_cov_cojo(ef2,ind)*sqrt(matPOIDS(ind,1)))
                          endif
                       enddo
                    enddo
                 endif
              enddo

              ! on remplit row_cojo pour les variants parente  non exclus
              do ef1=1,nbSNP
                 if((ef1.lt.SegmElim_P3(1,VarGWAS)).or.(ef1.gt.(SegmElim_P3(2,VarGWAS)))) then  ! le SNPpar est a conserver dans les equations sinon 0
                    ! on recentre les genotypes des nind individus pour le ef1 eme marqueur 
                    vecGENO1=(dfloat(GENOPAR(:,ef1))-2.0d0*freqUtil(ef1))*(sqrt(matPOIDS(:,1))) ! on multplie par le (POIDS)^0.5 car on a deja multiplie matWeight_cov_cojo par les (POIDS)^0.5
                    do ef2=1,neff_cojo
                       do ind=1,nbind
                          if(matNivAnim_cojo(ef2,1,ind).ne.0) then
                             row_cojo(ef2+matNivAnim_cojo(ef2,1,ind)-1,(FIRSTNIV_EP(neff+1)+ef1-1)) = row_cojo(ef2+matNivAnim_cojo(ef2,1,ind)-1,(FIRSTNIV_EP(neff+1)+ef1-1)) + vecGENO1(ind)*matWeight_cov_cojo(ef2,ind)
                          endif
                       enddo
                    enddo
                 endif
              enddo

              ! on remplit col_cojo=t(row_cojo)
              do i=1,neq_ep
                 do j=1,neq_cojo
                    col_cojo(i,j)=row_cojo(j,i)
                 enddo
              enddo

              ! inf_LEFT_EP_act

              ! on va utiliser les proprietes de l inverse d une matrice partitionnee
              !B22 = (A22 - A21 (A11)-1 A12)-1 avec A22 = matrice cojo,cojo du membre de gauche , A21 = neq_cojo lignes de cojo du mbre de gauche = row_cojo, A12 = col_cojo et (A11)-1 = inf_LEFT_EP_act

              ! A22 = t(mat_incidence_cojo(1:nbind,1:neq_cojo))%*%mat_incidence_cojo(1:nbind,1:neq_cojo)
              do ind=1,nbind
                 if(matPERF(ind,1).ne.mis) then ! if individual has a phenotype
                    do ef1=1,neff_cojo
                       do ef2=1,neff_cojo
                          if(matNivAnim_cojo(ef2,1,ind).ne.0) then
                             A22(ef1+matNivAnim_cojo(ef1,1,ind)-1,ef2+matNivAnim_cojo(ef2,1,ind)-1) = A22(ef1+matNivAnim_cojo(ef1,1,ind)-1,ef2+matNivAnim_cojo(ef2,1,ind)-1) + (matWeight_cov_cojo(ef1,ind)*matWeight_cov_cojo(ef2,ind))
                          endif
                       enddo
                    enddo
                 endif
              enddo

              print*,'A22 = '
              do ef1=1,neq_cojo
                 print*,A22(ef1,:)
              enddo
              print*,' '

              !---------------------------------------------------------------
              ! B22 = bloc cojo x cojo de l inverse du membre de gauche = inverse( A22 - A21 (A11)-1 A12 )
              allocate(motif1_cojo(neq_cojo,neq_ep),B12(neq_ep,neq_cojo)) ! motif1_cojo = A21%*%((A11)-1)

              motif1_cojo = 0.0d0

              call dsymm('R','L',neq_cojo,neq_ep,1.0d0,inf_LEFT_EP_act,neq_ep,row_cojo,neq_cojo,0.0d0,motif1_cojo,neq_cojo)

              !allocate(B22(neq_cojo,neq_cojo))
              B22 = A22 ! car on va faire C = alpha A%*%B + Beta C avec C en entree = A22 et C en sortie =inv_B22 : inv_B22 = -1*motif1_cojo%*%col_cojo + A22

              call dgemm('N','N',neq_cojo,neq_cojo,neq_ep,-1.0d0,motif1_cojo,neq_cojo,col_cojo,neq_ep,1.0d0,B22,neq_cojo)

              ! on doit maintenant inverser inv_B22 pour obtenir B22 ; NB : inv_B22 est neq_cojo x neq_cojo = entre 1x1 et 4x4 selon modele = petit
              call invmatlapack(B22(1:neq_cojo,1:neq_cojo),neq_cojo)

              ! ON DISPOSE MAINTENANT DE B22
              !---------------------------------------------------------------

              ! on va calculer B12 = bloc neq_ep x neq_cojo (en haut a droite) de l inverse du memebre de gauche avec les equations Cojo
              ! B12 = -1 * ((A11)-1)%*%A12 %*% B22 = -1 * t(A21%*%((A11)-1)) %*% B22 = -1 * t(motif1_cojo) %*% B22
              B12=0.0d0

              call dgemm('T','N',neq_ep,neq_cojo,neq_cojo,-1.0d0,motif1_cojo,neq_cojo,B22,neq_cojo,0.0d0,B12,neq_ep)

              ! ON DISPOSE MAINTENANT DE B12
              !---------------------------------------------------------------

              ! on calcule maintenant B11 = le bloc neq_ep x neq_ep en haut à gauche de l'inverse du membre de gauche avec les equations Cojo

              ! B11 = (A11)-1  +  ((A11)-1)%*%A12 %*% B22 %*% A21 %*% (A11)-1 = (A11)-1 + t(motif1_cojo) %*% B22 %*% motif1_cojo  = (A11)-1 - t(motif1_cojo) %*% B21 = (A11)-1 - t(motif1_cojo) %*% t(B12) 

              allocate(B11(neq_ep,neq_ep))
              B11=0.0d0
              B11 = inf_LEFT_EP_act ! car on va utiliser dgemm pour calculer C = -1 * t(A)%*%t(B) + 1*C  donc on initialise C avec inf_LEFT_EP_act = (A11)-1

              call dgemm('T','T',neq_ep,neq_ep,neq_cojo,-1.0d0,motif1_cojo,neq_cojo,B12,neq_ep,1.0d0,B11,neq_ep)

              ! ON DISPOSE MAINTENANT DE B11
              !---------------------------------------------------------------

              ! on met a 0 la triangulaire superieure de B11 au dessus de la diagonale car on a besoin uniquement de la triangulaire inferieure
              ! je ne sais plus pourquoi j avais fait cela, c est peut etre inutile, mais bon ...
              do j=2,neq_ep
                 do i=1,j-1
                    B11(i,j)=0.0d0
                 enddo
              enddo

              ! B21=t(B12) --> inutile de le calculer

              deallocate(inf_LEFT_EP_act)

              allocate(inf_LEFT_EP_act(neq_ep_2,neq_ep_2))

              inf_LEFT_EP_act=0.0d0

              inf_LEFT_EP_act(1:neq_ep,1:neq_ep)=B11

              inf_LEFT_EP_act(neq_ep+1:neq_ep_2,1:neq_ep)=TRANSPOSE(B12)

              do i=1,neq_cojo
                 do j=1,i
                    inf_LEFT_EP_act(neq_ep+i,neq_ep+j)=B22(i,j)
                 enddo
              enddo

              print*,' '
              print*,'Actualized inverse matrix of inf_LEFT_EP_act including COJO Variant is ready'
              print*,' '
              if((detailed_log.eq.3).and.(neq_ep.le.20)) then
                 print*,'Actualized inverse matrix of inf_LEFT_EP_act including COJO Variant :'
                 do i=1,neq_ep_2
                    print*,inf_LEFT_EP_act(i,:)
                 enddo
              endif

              ! RIGHT_EP_act
              do ind=1,nbind
                 if(matPERF(ind,1).ne.mis) then
                    do ef1=1,neff_cojo
                       if(matNivAnim_cojo(ef1,1,ind).ne.0) then
                          RIGHT_EP_act(neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1,1) = RIGHT_EP_act(neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1,1) + (matWeight_cov_cojo(ef1,ind)*vecPERFpond(ind))
                       endif
                    enddo
                 endif
              enddo

              print*,' '
              print*,'Actualized RIGHT_EP_act vector  including COJO variant is ready'
              if((detailed_log.eq.3).and.(neq_ep.le.20)) then
                 print*,'Actualized RIGHT_EP_act vector  including COJO variant :'
                 print*,RIGHT_EP_act(:,1)
              endif
              print*,' '

              deallocate(col_cojo,row_cojo,A22,B22,motif1_cojo,B12,B11)
              

           endif ! fin du test if indic_cojo = 1


           ! on estime les effets environnementaux et SNPpar sans VarGWAS dans le modele uniquement si 1er VarGWAS ou si SNPpar exclus differents du VarGWAS precedent

           if(detailed_log.eq.3) then
              print*,' '
              call fdate(jour)
              print*,'Before computing solutions BetaS without any GWAS Variant using actualized LEFT_EP inverse and RIGHT_EP after exclusion of Markers on segment if needed'
              print*,jour
              print*,' '
           endif

           BetaS=0.0d0
           !call dsymm('L','L',neq_ep_2,1,1.0d0,inf_LEFT_EP_act,neq_ep_2,RIGHT_EP_act(1:neq_ep_2,1),neq_ep_2,0.0d0,BetaS,neq_ep_2)
           call dsymm('L','L',neq_ep_2,1,1.0d0,inf_LEFT_EP_act,neq_ep_2,RIGHT_EP_act(1:neq_ep_2,1:1),neq_ep_2,0.0d0,BetaS,neq_ep_2) ! modif JVDP
           ! BetaS contient les solutions pour les effets environnementaux et pour les effets SNPparente

           print*,' '
           print*,'Solutions for the non genetic effects in the reduced model (without any GWAS Variant):'
           print*,'(limited to 30 solutions)'
           print*,BetaS(1:min(sum(nlev(1:neff)),30),1)
           print*,' '


           print*,' '
           call fdate(jour)
           print*,'After computing solutions BetaS without any GWAS Variant using actualized LEFT_EP inverse and RIGHT_EP after exclusion of Markers on segment if needed'
           print*,jour
           print*,' '

           BASE1 = dot_product(vecPERFpond,vecPERFpond) - (dot_product( BetaS(:,1) , RIGHT_EP_act(1:neq_ep_2,1)))

           print*,' '
           print*,'================================================================================================================'
           print*,'VALUE OF BASE1 calculated without any GWAS Variant and after exclusion of undesired Markers (if any)  =',BASE1
           print*,'================================================================================================================'
           print*,' '


        endif  ! fin du test sur 1er variant GWAS ou sur difference de SNPpar a exclure par rapport a VarGWAS precedent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! inf_LEFT_EP_act et RIGHT_EP_act sont a jour pour le VarGWAS en cours et BASE1 a ete calcule si necessaire 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! ===========================================================================
        ! PARTIE SPECIFIQUE A CHAQUE VARIANT GWAS :  matrices ZpZ    [ZpX ZpM]  ZpY   
        ! ===========================================================================


        ! MODIF DU 17/12/2024 : on n effectue la suite que SI VarDiffCojo = 1 = tout le temps si pas OPTION COJO, et si info genomique du Variant est differente de celle du Variant COJO si OPTION COJO

        if(VarDiffCojo.eq.1) then

           ! MODIF DU 06fev2023 : on n effectue la suite que si le variant GWAS  n est pas monomorphe = si STD(Doses) different de 0
           !!!! PROBLEME : PAS ASSEZ SEVERE : pour modele add_dom avec GENOTYPES DISCRETS on a trouve des variants tq 7094 indiv 11, 0 indiv 12 et 1 indiv 22 --> cela doit generer qqch indéfini
           !!!! --> on remplace par (nb_11 ou nb_22 > 5 ET nb_12 > 5 --> on a au moins des heterozygotes et 1 des 2 genotype homozygote represente avec au moins 5 indiv. Pas beaucop mais les StatTest seront ridicules --> pas signif
           !if(vec_stat_doses(2,varGWAS).ne.0.0d0) then
           if(tem_possible(1,varGWAS).eq.1) then



              ! ======================================================================================================================================
              ! On construit la matrice (Z'X  Z'M) pour la partie du modele de GWAS pour le variant considere  
              ! ======================================================================================================================================


              ! on remplit (Z'X) pour les effets fixes, covariables et effets aleatoires diagonaux
              do ind=1,nbind
                 if(matPERF(ind,1).ne.mis) then
                    do ef1=1,neff
                       ZpXZpM(1,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) = ZpXZpM(1,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) + (matWeight_cov(ef1,ind)*vecDosesPond(ind,1)*sqrt(matPOIDS(ind,1))) ! contribution pour effet additif
                       if(TemDom.eq.1) then ! si on a un effet de Dominance estimable pour le variant considere
                          !if(neff_GWAS.eq.2) then
                          ZpXZpM(coldom,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) = ZpXZpM(coldom,(FIRSTNIV_EP(ef1)+matNivAnim(ef1,1,ind)-1)) + (matWeight_cov(ef1,ind)*vecDosesPond(ind,coldom)*sqrt(matPOIDS(ind,1))) ! contribution pour effet dominance
                       endif
                    enddo
                 endif
              enddo

              ! SI COJO ON REMPLIT L ELEMENT Z'X CORRESPONDANT A L EFFET DU VARIANT COJO
              if(indic_cojo.eq.1) then
                 do ind=1,nbind
                    if(matPERF(ind,1).ne.mis) then
                       do ef1=1,neff_cojo 
                          ZpXZpM(1,neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1) = ZpXZpM(1,neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1) + (matWeight_cov_cojo(ef1,ind)*vecDosesPond(ind,1)) ! contribution pour effet additif
                          if(TemDom.eq.1) then ! si on a un effet de Dominance estimable pour le variant considere
                             !if(neff_GWAS.eq.2) then
                             ZpXZpM(coldom,neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1) = ZpXZpM(coldom,neq_ep+ef1+matNivAnim_cojo(ef1,1,ind)-1) + (matWeight_cov_cojo(ef1,ind)*vecDosesPond(ind,coldom)) ! contribution pour effet dominance 
                          endif
                       enddo
                    endif
                 enddo
              endif


              ! on remplit (Z'M) pour les variants parente  non exclus

              do ef1=1,nbSNP
                 !if((ef1.lt.(FIRSTNIV_EP(neff+1) + SegmElim_P3(1,VarGWAS) -1)).or.(ef1.gt.(FIRSTNIV_EP(neff+1) + SegmElim_P3(2,VarGWAS) -1))) then  ! le SNPpar est a conserver dans les equations sinon 0
                 if((ef1.lt.SegmElim_P3(1,VarGWAS)).or.(ef1.gt.(SegmElim_P3(2,VarGWAS)))) then  ! le SNPpar est a conserver dans les equations sinon 0

                    ! on recentre les genotypes des nind individus pour le ef1 eme marqueur 
                    vecGENO1=(dfloat(GENOPAR(:,ef1))-2.0d0*freqUtil(ef1))*(sqrt(matPOIDS(:,1))) ! on multplie par le (POIDS)^0.5 car on a deja multiplie vecDosesPondAdd et vecDosesPondDom par les (POIDS)^0.5

                    ZpXZpM(1,(FIRSTNIV_EP(neff+1)+ef1-1))=  dot_product(vecGENO1,vecDosesPond(:,1))
                    !if(neff_GWAS.eq.2) ZpXZpM(2,(FIRSTNIV_EP(neff+1)+ef1-1))=  dot_product(vecGENO1,vecDosesPond(:,2))
                    if(TemDom.eq.1) ZpXZpM(coldom,(FIRSTNIV_EP(neff+1)+ef1-1))=  dot_product(vecGENO1,vecDosesPond(:,coldom)) ! si on a un effet de dominance estimable pour le variant GWAS considere
                 endif
              enddo

              !if(neq_ep.lt.50) then
              !   print*,' '
              !   print*,'Variant GWAS ',varGWAS,' matrice ZpXZpM = '
              !   do i=1,neff_GWAS_lu
              !      print*,ZpXZpM(i,:)
              !   enddo
              !   print*,' '
              !endif

              ! =============================================================
              ! On construit la matrice (Z'y) pour le variant GWAS considere  
              ! =============================================================

              ZpY=0.0d0

              ZpY(1,1)=dot_product(vecDosesPond(:,1),vecPERFpond)
              !if(neff_GWAS.eq.2) ZpY(2,1)=dot_product(vecDosesPond(:,2),vecPERFpond)
              if(TemDom.eq.1) ZpY(coldom,1)=dot_product(vecDosesPond(:,coldom),vecPERFpond) ! si on a un effet de dominance estimable pour le variant GWAS considere

              !if(neq_ep.lt.50) then
              !   print*,' '
              !   print*,'Variant GWAS ',varGWAS,' matrice ZpWY = '
              !   print*,ZpY(:,1)
              !   print*,' '
              !endif

              ! =============================================================
              ! On construit la matrice (Z'Z) pour le variant GWAS considere  
              ! =============================================================

              ZpZ=0.0d0

              call dgemm('T','N',neff_GWAS,neff_GWAS,nbind,1.0d0,vecDosesPond(:,1:neff_GWAS),nbind,vecDosesPond(:,1:neff_GWAS),nbind,0.0d0,ZpZ(1:neff_GWAS,1:neff_GWAS),neff_GWAS)

              !if(neq_ep.lt.50) then
              !   print*,' '
              !   print*,'Variant GWAS ',varGWAS,' matrice ZpWZ = '
              !   do i=1,neff_GWAS_lu
              !      print*,ZpZ(i,:)
              !   enddo
              !   print*,' '
              !endif

              ! on calcule le determinant de ZpZ si on a un effet additif et un effet de dominance
              !if(neff_GWAS_lu.eq.2) then
              if(TemDom_lu.eq.1) then
                 !if(neff_GWAS.eq.2) detZpZ(varGWAS)= (ZpZ(1,1)*ZpZ(2,2)) - (ZpZ(2,1)*ZpZ(1,2))
                 if(TemDom.eq.1) detZpZ(varGWAS)= (ZpZ(1,1)*ZpZ(2,2)) - (ZpZ(2,1)*ZpZ(1,2)) ! si on a un effet de dominance estimable pour le variant GWAS considere
              endif

           endif  ! fin du test  if(tem_possible(1,varGWAS).eq.1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           ! if (vec_stat_doses(2,varGWAS).ne.0.0d0) then
           if(tem_possible(1,varGWAS).eq.1) then

              ! ================================================================================================================
              ! ================================================================================================================

              ! MODIF V2_iter : on calcule les solutions add_hat et dom_hat de maniere exacte et betaSnew sachant add_hat et dom_hat

              ! ================================================================================================================
              ! ================================================================================================================


              ! **************************************************************************
              ! ETAPE 1 : on calcule les solutions add_hat et dom_hat de maniere exacte 
              ! **************************************************************************

              ! ==================================================================================
              ! On calcule interm_1 = [Z'X Z'M] inv(LEFT_EP_act) pour le variant GWAS considere  
              ! ==================================================================================
              interm_1=0.0d0
              call dsymm('R','L',neff_GWAS,neq_ep_2,1.0d0,inf_LEFT_EP_act,neq_ep_2,ZpXZpM(1:neff_GWAS,:),neff_GWAS,1.0d0,interm_1(1:neff_GWAS,:),neff_GWAS)


              ! ==================================================================================
              ! On calcule interm_2 = interm1 t[Z'X Z'M] pour le variant GWAS considere  
              ! ==================================================================================
              interm_2=0.0d0
              call dgemm('N','T',neff_GWAS,neff_GWAS,neq_ep_2,1.0d0,interm_1(1:neff_GWAS,:),neff_GWAS,ZpXZpM(1:neff_GWAS,:),neff_GWAS,0.0d0,interm_2(1:neff_GWAS,1:neff_GWAS),neff_GWAS)


              !print*,' '
              !print*,'interm2 = '
              !do i=1,neff_GWAS
              !   print*,interm_2(i,:)
              !enddo

              ! ==================================================================================
              ! On calcule interm_3 = interm1 [X'Y M'Y] pour le variant GWAS considere  
              ! ==================================================================================

              ! ! ! ! ! ! ! ! LE PROBLEME SEMBLE SE PASSER DANS LES 2 LIGNES CI DESSOUS

              interm_3=0.0d0
              !call dgemm('N','N',neff_GWAS,1,neq_ep_2,1.0d0,interm_1(1:neff_GWAS,1:neq_ep_2),neff_GWAS,RIGHT_EP_act(1:neq_ep_2,1),neq_ep_2,0.0d0,interm_3(1:neff_GWAS,1),neff_GWAS)
              call dgemm('N','N',neff_GWAS,1,neq_ep_2,1.0d0,interm_1(1:neff_GWAS,1:neq_ep_2),neff_GWAS,RIGHT_EP_act(1:neq_ep_2,1:1),neq_ep_2,0.0d0,interm_3(1:neff_GWAS,1:1),neff_GWAS) ! modif JVDP




              ! ==================================================================================
              ! On calcule interm_4 = ZpZ - interm_2 pour le variant GWAS considere  
              ! ==================================================================================
              interm_4 = ZpZ - interm_2

              ! ==================================================================================
              ! On calcule interm_5 = ZpY - interm_3 pour le variant GWAS considere  
              ! ==================================================================================
              interm_5=0.0d0
              interm_5 = ZpY - interm_3

              ! ============================================================================================================
              ! On calcule les solutions pour l effet additif et eventuellement de dominance pour le variant GWAS en cours  
              ! ============================================================================================================    

              call invmatlapack(interm_4(1:neff_GWAS,1:neff_GWAS),neff_GWAS)

              !call dgemm('N','N',neff_GWAS,1,neff_GWAS,1.0d0,interm_4(1:neff_GWAS,1:neff_GWAS),neff_GWAS,interm_5(1:neff_GWAS,1),neff_GWAS,0.0d0,sol_VarGWAS(1:neff_GWAS,VarGWAS),neff_GWAS)
              call dgemm('N','N',neff_GWAS,1,neff_GWAS,1.0d0,interm_4(1:neff_GWAS,1:neff_GWAS),neff_GWAS,interm_5(1:neff_GWAS,1:1),neff_GWAS,0.0d0,sol_VarGWAS(1:neff_GWAS,VarGWAS:VarGWAS),neff_GWAS) ! modifJVDP


              ! modif du 05 janvier 2023
              ! on calcule les solutions pour les effets environnementaux et SNPparente pour calculer correctement la variance residuelle



              !===========================================================================!
              !                                                                           !
              ! CALCUL DE LA VARIANCE D ESTIMATION DES EFFETS DU VARIANT GWAS CONSIDERE   !
              ! COMPUTATION OF VARIANCE OF ESTIMATION FOR CURRENT GWAS VARIANT EFFECTS    !
              !                                                                           !
              !===========================================================================!

              ! _______________________________________________________________________________________________________________________________________
              !
              ! methode 1 : on considere que les effets environnementaux (beta) et des SNPparente ne dependent pas des effets du variant GWAS considere

              restemp=0.0d0

              !call dgemm('T','N',neff_GWAS,1,nbind,1.0d0,vecDosesPond(:,1:neff_GWAS),nbind,vecPERFpond,nbind,0.0d0,restemp(1:neff_GWAS,1),neff_GWAS)
              call dgemm('T','N',neff_GWAS,1,nbind,1.0d0,vecDosesPond(:,1:neff_GWAS),nbind,reshape(vecPERFpond, [size(vecPERFpond),1]),nbind,0.0d0,restemp(1:neff_GWAS,1:1),neff_GWAS) ! modif JVDP


              contrib_VarGWAS = dot_product(sol_VarGWAS(1:neff_GWAS,VarGWAS),restemp(1:neff_GWAS,1))
              !!print*,'contrib_VarGWAS = ',contrib_VarGWAS

              if((CalcVarRes.eq.'approx').or.(CalcVarRes.eq.'optim')) then

                 exp4_1=0.0d0

                 exp4_1 = (BASE1 - contrib_VarGWAS)*invdenom ! inutile a terme car BASE2 calcule juste apres mais utile pour comparaison

                 vec_VarRes(1,VarGWAS)=exp4_1

                 ! On utilise matC l inverse de interm_4 = inv(ZpZ - (ZpX ZpM) inv(LEFT_EP_act) (XpZ MpZ)) = matrice (neff_GWAS x neff_GWAS) qu on a deja inverse pour calculer les solutions add et dom
                 Test_effVarGWAS_1=0.0d0
                 do i=1,neff_GWAS
                    !if((interm_4(i,i).ne.0.0d0).and.(exp4_1.gt.0.0d0)) then  ! modif 04/09/2025 car quelques variants avec STD(DOSES)tres proche de 0 ont un interm4 < 0 du a limite precision machine
                    if((interm_4(i,i).gt.0.0d0).and.(exp4_1.gt.0.0d0)) then
                       Test_effVarGWAS_1(i) = abs(sol_VarGWAS(i,VarGWAS)) / (sqrt(interm_4(i,i)*exp4_1))   !  formule a utiliser mais plante sur exemple de 4 individus car interm_4*exp4 negatif
                    endif
                    mat_TestEff_1(i,VarGWAS)=Test_effVarGWAS_1(i)
                    mat_pval_1(i,VarGWAS) = normal_ttr(Test_effVarGWAS_1(i),dx)
                 enddo

              endif


              ! **********************************************************************************************************************************
              ! ETAPE 2 : on calcule les solutions add_hat et dom_hat de maniere exacte pour les effets environnementaux sans SNPparente a exclure
              ! **********************************************************************************************************************************

              if((CalcVarRes.eq.'exact').or.(exp4_1.le.0.0d0).or.((CalcVarRes.eq.'optim').and.(maxval(Test_effVarGWAS_1(1:neff_GWAS)).gt.seuil_test1))) then ! on recalcule exactement egalement si VarResapprox < 0

                 BetaSnew=0.0d0

                 exp4_2=0.0d0

                 interm_6=RIGHT_EP_act ! on initialise interm_6

                 !call dgemm('T','N',neq_ep_2 ,1 ,neff_GWAS ,-1.0d0,ZpXZpM(1:neff_GWAS,:),neff_GWAS ,sol_VarGWAS(1:neff_GWAS,VarGWAS),neff_GWAS ,1.0d0,interm_6(1:neq_ep_2,1),neq_ep_2 )
                 call dgemm('T','N',neq_ep_2 ,1 ,neff_GWAS,-1.0d0,ZpXZpM(1:neff_GWAS,:),neff_GWAS,sol_VarGWAS(1:neff_GWAS,VarGWAS:VarGWAS),neff_GWAS ,1.0d0,interm_6(1:neq_ep_2,1:1),neq_ep_2 ) ! modif JVDP

                 ! on calcule inverse(membre de gauche sans variant GWAS actualise) %*% interm_6 = solutions pour effets environnementaux et SNPparente avec variant GWAS dans modele
                 !call dsymm('L','L',neq_ep_2,1,1.0d0,inf_LEFT_EP_act,neq_ep_2,interm_6(1:neq_ep_2,1),neq_ep_2,0.0d0,BetaSnew,neq_ep_2)
                 call dsymm('L','L',neq_ep_2,1,1.0d0,inf_LEFT_EP_act,neq_ep_2,interm_6(1:neq_ep_2,1:1),neq_ep_2,0.0d0,BetaSnew,neq_ep_2) ! modif JVDP

                 BASE2 = dot_product(vecPERFpond,vecPERFpond) - (dot_product( BetaSnew(:,1) , RIGHT_EP_act(1:neq_ep_2,1)))
                 !!print*,'BASE2=',BASE2

                 exp4_2 = (BASE2 - contrib_VarGWAS)*invdenom

                 vec_VarRes(2,VarGWAS)=exp4_2

                 ! On utilise matC l inverse de interm_4 = inv(ZpZ - (ZpX ZpM) inv(LEFT_EP_act) (XpZ MpZ)) = matrice (neff_GWAS x neff_GWAS) qu on a deja inversee pour calculer les solutions add et dom
                 Test_effVarGWAS_2=0.0d0
                 do i=1,neff_GWAS
                    !if ((interm_4(i,i).ne.0.0d0).and.(exp4_2.gt.0.0d0)) then  ! some variants with STD(DOSE) very close to 0 have interm_4 <0 because of machine accuracy limits
                    if ((interm_4(i,i).gt.0.0d0).and.(exp4_2.gt.0.0d0)) then
                       Test_effVarGWAS_2(i) = abs(sol_VarGWAS(i,VarGWAS)) / (sqrt(interm_4(i,i)*exp4_2))
                    endif
                    mat_TestEff_2(i,VarGWAS)=Test_effVarGWAS_2(i)
                    mat_pval_2(i,VarGWAS)=normal_ttr(Test_effVarGWAS_2(i),dx)
                 enddo

              endif

           endif  ! fin du test  le VarGWAS est il monomorphe ?

        endif  ! fin du test if(VarDiffCojo.eq.1) 



        if(mod(VarGWAS,100).eq.0) then
           call fdate(jour)
           print*,'GWAS for Variant ',VarGWAS,' completed ',jour
        endif


        if(mod(VarGWAS,StepFichOut).eq.0) then

           do k=(CurrLine+1),VarGwas



              if(CalcVarRes.eq.'optim') then
                 write(io_s,fmt=format_res) k,GWASuttoMAP_p3(k),vec_stat_doses(:,k),vecTemDom(k),sol_VarGWAS(:,k),mat_TestEff_1(:,k),mat_pval_1(:,k),mat_TestEff_2(:,k),mat_pval_2(:,k),vec_VarRes(1,k),vec_VarRes(2,k),detZpZ(k) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
              endif

              if(CalcVarRes.eq.'approx') then
                 write(io_s,fmt=format_res) k,GWASuttoMAP_p3(k),vec_stat_doses(:,k),vecTemDom(k),sol_VarGWAS(:,k),mat_TestEff_1(:,k),mat_pval_1(:,k),vec_VarRes(1,k),detZpZ(k) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
              endif

              if(CalcVarRes.eq.'exact') then
                 write(io_s,fmt=format_res) k,GWASuttoMAP_p3(k),vec_stat_doses(:,k),vecTemDom(k),sol_VarGWAS(:,k),mat_TestEff_2(:,k),mat_pval_2(:,k),vec_VarRes(2,k),detZpZ(k) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
              endif


           enddo

           CurrLine=CurrLine+StepFichOut

        endif

     enddo   ! fin de la boucle sur le variant VarGWAS

     close(io_dosesbin0)


     print*,' '
     call fdate(jour)
     print*,'Last GWAS Variant has been treated'
     print*,jour
     print*,' '



     !deallocate(inf_LEFT_EP_act)
     !deallocate(RIGHT_EP_act)
     !deallocate(ZpXZpM)
     !deallocate(vecDosesPond)
     !deallocate(vecDoses)
     !deallocate(vecGENO1)
     !deallocate(ZpY)
     !deallocate(ZpZ)
     !deallocate(interm_1)
     !deallocate(interm_2)
     !print*,'Probleme quand on tente desallocation de interm_3 et interm_4'
     !deallocate(interm_3)
     !deallocate(interm_4)
     !deallocate(interm_5)

     ! on ecrit les resultats des derniers variants GWAS au dela de la derniere tranche de StepFichOut

     if(Currline.lt.nbGWASexec) then

        do i=(Currline+1),nbGWASexec

           if(CalcVarRes.eq.'optim') then
              write(io_s,fmt=format_res) i,GWASuttoMAP_p3(i),vec_stat_doses(:,i),vecTemDom(i),sol_VarGWAS(:,i),mat_TestEff_1(:,i),mat_pval_1(:,i),mat_TestEff_2(:,i),mat_pval_2(:,i),vec_VarRes(1,i),vec_VarRes(2,i),detZPZ(i) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
           endif

           if(CalcVarRes.eq.'approx') then
              write(io_s,fmt=format_res) i,GWASuttoMAP_p3(i),vec_stat_doses(:,i),vecTemDom(i),sol_VarGWAS(:,i),mat_TestEff_1(:,i),mat_pval_1(:,i),vec_VarRes(1,i),detZPZ(i) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
           endif

           if(CalcVarRes.eq.'exact') then
              write(io_s,fmt=format_res) i,GWASuttoMAP_p3(i),vec_stat_doses(:,i),vecTemDom(i),sol_VarGWAS(:,i),mat_TestEff_2(:,i),mat_pval_2(:,i),vec_VarRes(2,i),detZPZ(i) ! A TESTER : si effets_estimes Tests et varRes dans meme table I/O plus rapide ?
           endif

        enddo

     endif


     close(io_s)

     print*,'Step_3 completed'

  endif  ! End of paft executed if steptodo = 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  STOP 0

  call fdate(jour)
  print*,'End of program'
  print*,jour
  print*,'*******************************'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains


  !========================================================================
  subroutine read_parameters_hm()
    !========================================================================
    ! reads and decodes parameter file
    implicit none

    integer n,ef,stat,curr_eff,i,j,k
    real(kind=8) x(100)
    !real(kind=8) xreal(100)
    character xc(100)*150

    call getarg(1,parfile)
    parfile=adjustl(parfile)               ! ignore leading spaces 

    print*,'Parameter file read:'
    write(*,'(5X,a)') trim(parfile)
    open(io_p,file=parfile)

    call chkfmtnewline(io_p)

    call readline(io_p,x,xc,n)
    if (n == -1) then
       print*,'INPUT FILE EMPTY'
       stop 21
    endif


    if (xc(1) /= 'STEPS') then
       print*,'STEPS not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    steptodo=x(1)

    print*,' '
    print*,'Step that will be performed = ',steptodo
    if((steptodo.ne.1).and.(steptodo.ne.2).and.(steptodo.ne.3)) then
       print*,'STEPS different from 1 2 or 3'
       STOP 21
    endif

    if(steptodo.eq.3) then
       call readline(io_p,x,xc,n)
       if (xc(1) /= 'Chrom_numGR') then
          print*,'STEP 3 but Chrom_numGR not found'
          stop 21
       endif
       call readline(io_p,x,xc,n)
       if(n.ne.3) then
          print*,'3 integer values expected after keyword Chrom_numGR'
          stop 21
       endif
       chrom_p3=x(1)
       numgr_p3=x(2)
       numSgr_p3=x(3)

       print*,' '
       print*,'chrom_p3 read =',chrom_p3
       print*,'numgr_p3 read =',numgr_p3
       print*,'numSgr_p3 read =',numSgr_p3

    endif


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'EXEC_PATH') then
       print*,'EXEC_PATH not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    execpath=xc(1)
    print*,' '
    !if(detailed_log.eq.3) print*,'EXEC_PATH = ',execpath


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'BIN_DIR') then
       print*,'BIN_DIR not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    location=xc(1)
    !if(detailed_log.eq.3) print*,'BIN_DIR = ',location

    ! verification de l existence du repertoire des fichiers binaires
    INQUIRE(DIRECTORY =trim(location), EXIST=testficdir)    !!! ATTENTION DIRECTORY ne fonctionne pas en gfortran (a remplacer par file)
    if (testficdir.eq..false.) then
       print *
       print "(a)", " WORNING!!"
       print "(a)", " -----------"
       print "(a)", " Directory "//trim(location)//" does not exist"
       print "(a)", "  => CHECK"
       print*
       STOP 11
    end if

    call readline(io_p,x,xc,n) 
    if (xc(1) /= 'DATAFILE') then
       print*,'DATAFILE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    datafile=xc(1)

    !if(detailed_log.eq.3) print*,'datafile = ',datafile

    !call readline(io_p,x,xc,n)
    !if (xc(1) /= 'NUMBER_OF_TRAITS') then
    !   print*,'NUMBER_OF_TRAITS not found'
    !   stop 21
    !endif
    !call readline(io_p,x,xc,n) 
    !ntrait=x(1)
    !print*,'ntrait lu =',ntrait

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'NUMBER_OF_EFFECTS') then
       print*,'NUMBER_OF_EFFECTS not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    neff=x(1)
    !if(detailed_log.eq.3) print*,'neff read =',neff

    !allocate(rep(0:neff))
    allocate(pos_y(ntrait),pos_weight(ntrait),pos_eff(neff,ntrait),&
         nlev(neff),effecttype(neff),nestedcov(neff,ntrait),&
         randomtype(neff),randomnumb(neff),&
         r(ntrait,ntrait),rinv(ntrait,ntrait),g(ntrait,ntrait),ginv(ntrait,ntrait),rand(neff,maxcorr,maxcorr),randinv(neff,maxcorr,maxcorr),RhoDiag(neff,maxcorr,maxcorr))
    randomtype=g_fixed
    randomnumb=0
    nestedcov=0
    r=0
    g=0
    rand=0
    typparfile=' '
    fmttyppar=' '
    mapparfile=' '
    typgwasfile=' '
    fmttypgwas=' '
    mapgwasfile=' '

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'POS_ID_CHAR') then
       print*,' POS_ID_CHAR not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    posIdchar=x(1)
    !if(detailed_log.eq.3) print*,'posIdchar read =',posIdchar

    !call readline(io_p,x,xc,n)
    !if (xc(1) /= 'NB_IND') then
    !   print*,' NB_IND not found'
    !   stop 21
    !endif
    !call readline(io_p,x,xc,n) 
    !nbind=x(1)
    !print*,'nbind lu =',nbind



    call readline(io_p,x,xc,n)
    if (xc(1) /= 'OBSERVATION(S)') then
       print*,' OBSERVATION(S) not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    if (n == ntrait) then
       pos_y(1:ntrait)=x(1:ntrait)
    else
       print*,ntrait,' numbers expected after OBSERVATION(S)'
       stop 21
    endif

    !if(detailed_log.eq.3) print*,'Position of phenotype in DATA file = ',pos_y(1:ntrait)

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'WEIGHT(S)') then
       print*,' WEIGHT(S) not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    if (n == ntrait) then
       pos_weight(1:ntrait)=x(1:ntrait)
    elseif (n == 0) then
       pos_weight=0
    else
       print*,ntrait,' numbers expected after WEIGHT(S) (1 weight per trait)'
       stop 21
    endif

    !if(detailed_log.eq.3) print*,'Position of weights of phenotypes in DATA file = ',pos_weight(1:ntrait)

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'EFFECTS:') then
       print*,' EFFECTS: not found'
       stop 21
    endif
    do ef=1,neff
       call readline(io_p,x,xc,n)
       !print*,'print de x() pour verif covnested) : '
       !print*,'effet ',ef,' : ',' / n = ',n,' / x=',x(:)
       if (n == ntrait+2 .or. n == 2*ntrait+2) then
          pos_eff(ef,:)=x(1:ntrait)
          nlev(ef)=x(ntrait+1)
          select case (xc(ntrait+2)(1:3))
          case ('cro')         !crossclassified
             effecttype(ef)=effcross
          case ('cov')         !covariate
             effecttype(ef)=effcov
             if (n==2*ntrait+2) then
                nestedcov(ef,:)=x(ntrait+3:n)
             endif
          case default
             print*,'Unknown type of effect: ',xc(ntrait+2)
          end select
       else
          print*,'too few or too many numbers for effect ',ef,'. Was ',n
          stop 21
       endif
    enddo

    !if(detailed_log.eq.3) then 
    !   print*,'Non genetic effects of the model : position / nb of levels / effect type / random_type / parent effect if nested regr coef'
    !   do ef=1,neff
    !      print*,pos_eff(ef,1),'  /  ',nlev(ef),'  /  ',effecttype(ef),'  /  ',randomtype(ef),'  /  ',nestedcov(ef,1)
    !   enddo
    !endif

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'RESIDUAL_VARIANCE') then
       print*,' RESIDUAL_VARIANCE not found'
       stop 21
    endif
    read(io_p,*,iostat=stat)((r(i,j),i=1,ntrait),j=1,ntrait)
    if (stat.ne.0) then
       print*,'error reading RESIDUAL_VARIANCE'
       stop 21
    endif
    !if(detailed_log.eq.3) print*,'RESIDUAL VARIANCE read =',r(1,1)


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'GENETIC_VARIANCE') then
       print*,' GENETIC_VARIANCE not found'
       stop 21
    endif
    read(io_p,*,iostat=stat)((g(i,j),i=1,ntrait),j=1,ntrait)
    if (stat.ne.0) then
       print*,'error reading GENETIC_VARIANCE'
       stop 21
    endif
    !if(detailed_log.eq.3) print*,'GENETIC VARIANCE read =',g(1,1)


    ! code for random diagonal effects was not removed, but not tested so far --> don't use without testing!!!
    do
       call readline(io_p,x,xc,n)
       if (n == -1) exit
       if (xc(1) /='RANDOM_GROUP') exit

       call readline(io_p,x,xc,n)
       if (n<1 .or. n>neff) then
          print*,'line after RANDOM_GROUP has too few or too many numbers'
          stop 21
       endif
       if (x(1)<1 .or. x(n)> neff) then
          print*,'number ', x(1:n),' after RANDOM_GROUP out of range'
          stop 21
       endif
       do i=1,n
          if (i+x(1)-1 /= x(i)) then
             print*,' correlated effects:',x(1:n),' should be consecutive'
             stop 21
          endif
       enddo
       curr_eff=x(1)
       randomnumb(curr_eff)=n
       ! bridage temporaire pour limiter a des effets diagonaux independants : 1 seul effet aleatoire
       if(randomnumb(curr_eff).gt.1) then
          print*,'WARNING : Currently impossible to have several random diagonal effects correlated together'
          STOP 10
       endif

       call readline(io_p,x,xc,n)
       if (xc(1) /= 'RANDOM_TYPE') then
          print*,' RANDOM_TYPE not found'
          stop 21
       endif

       call readline(io_p,x,xc,n) 
       select case (xc(1)(1:15))
       case ('diag')
          randomtype(curr_eff)=g_diag
          temRand=1
          !case('add')
          !   randomtype(curr_eff)=g_AD
          !case('addmat')
          !   randomtype(curr_eff)=g_AD
          !   randomtype(curr_eff+1)=g_AM
       case default
          print*,'unknown RANDOM_TYPE: ', xc(1)
          print*,' Should be diag, other random effects currently not supported'
          stop 21
       end select

       call readline(io_p,x,xc,n)
       if (xc(1) /= 'VARIANCE') then
          print*,' VARIANCE not found'
          stop 21
       endif

       k=ntrait*randomnumb(curr_eff)
       if (k <= maxcorr) then
          read(io_p,*,iostat=stat)((rand(curr_eff,i,j),i=1,k),j=1,k)
          if (stat /= 0) then
             print*,'error reading variances for effect ',curr_eff
             stop 21
          endif
       else
          print*,'maxcorr should be increased to at least ',k
          stop 21
       endif
    enddo

    !if(temRand.eq.1) then
    !   print*,' '
    !   print*,'Non genetic effects : location in DATA file / nb of levels / effect type / random_type / parent effect if nested regr coef / associated variance'
    !   do ef=1,neff
    !      print*,pos_eff(ef,1),'  /  ',nlev(ef),'  /  ',effecttype(ef),'  /  ',randomtype(ef),'  /  ',nestedcov(ef,1),'  /  ',rand(ef,1,1)
    !   enddo
    !   print*,' '
    !endif


    !call readline(io_p,x,xc,n)
    if (xc(1) /= 'TYP_PAR_FILE') then
       print*,'keyword TYP_PAR_FILE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    typparfile=xc(1)

    !if(detailed_log.eq.3) print*,'Genotype file for Markers used to model polygenic BV = ',typparfile


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'FORMAT_TYP_PAR') then
       print*,'keyword FORMAT_TYP_PAR not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    fmttyppar=xc(1)

    if((fmttyppar.ne.'typ_eval').and.(fmttyppar.ne.'plink')) then
       print*,'unknown : FORMAT_TYP_PAR', xc(1)
       print*,' Should be typ_eval or plink'
       stop 21
    endif

    !if(detailed_log.eq.3) print*,'Format for Marker Genotype file = ',fmttyppar


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'MAP_PAR_FILE') then
       print*,'keyword MAP_PAR_FILE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) 
    mapparfile=xc(1)

    !if(detailed_log.eq.3) print*,'MAP file for Markers used to model poygenic BV = ',mapparfile

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'TYP_GWAS_FILE') then
       print*,'keyword TYP_GWAS_FILE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! not obligatory in STEP_1, obligatory in STEP_2 and STEP_3
    if(n.ne.0) then
       typgwasfile=xc(1)
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'Variant genotype/allelic dosage file is missing'
          stop 21
       elseif(steptodo.eq.1) then
          typgwasfile=' '
       endif
    endif
    !if(steptodo.ne.1) print*,'Variant genotype/allelic dosage file = ',typgwasfile

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'FORMAT_TYP_GWAS') then
       print*,'keyword FORMAT_TYP_GWAS not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! non obligatoire a renseigner en STEP_1, obligatoire en STEP_2 et STEP_3
    if(n.ne.0) then
       fmttypgwas=xc(1)
       if((fmttypgwas.ne.'minimac').and.(fmttypgwas.ne.'typ_eval').and.(fmttypgwas.ne.'plink')) then
          print*,'unknown : FORMAT_TYP_GWAS', xc(1)
          print*,' Should be minimac or typ_eval or plink'
          stop 21
       endif
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'Format for variant genotype/allelic dosage file is missing'
          stop 21
       elseif(steptodo.eq.1) then
          fmttypgwas=' '
       endif
    endif
    !if(steptodo.ne.1) print*,'Format for variant genotype/allelic dosage file = ',fmttypgwas

    !if(fmttypgwas.eq.'plink') then
    call readline(io_p,x,xc,n) 
    if (xc(1) /= 'NUMCHR_GWAS') then
       !print*,'FORMAT_TYP_GWAS = plink but keyword NUMCHR_GWAS not found'
       print*,'keyword NUMCHR_GWAS not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! non obligatoire a renseigner en STEP_1, obligatoire en STEP_2 et STEP_3
    if(n.ne.0) then
       NUMCHR_GWAS=x(1)
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'Id of Chromosome carrying the variants for GWAS is missing'
          stop 21
       elseif(steptodo.eq.1) then
          NUMCHR_GWAS=0
       endif
    endif
    !if(steptodo.ne.1) print*,'Id of chromosome carrying the variants for GWAS = ',NUMCHR_GWAS


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'MAP_GWAS_FILE') then
       print*,'keyword MAP_GWAS_FILE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! non obligatoire a renseigner en STEP_1, obligatoire en STEP_2 et STEP_3
    if(n.ne.0) then
       mapgwasfile=xc(1)
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'Map file of the variants for GWAS is missing'
          stop 21
       elseif(steptodo.eq.1) then
          mapgwasfile=' '
       endif
    endif
    !if(steptodo.ne.1) print*,'Map file of the variants for GWAS = ',mapgwasfile


    call readline(io_p,x,xc,n)
    if (xc(1) /= 'GWAS_TYPE') then
       print*,' GWAS_TYPE not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! non obligatoire a renseigner en STEP_1, obligatoire en STEP_2 et STEP_3
    if(n.ne.0) then
       gwastype=xc(1)
       if ((gwastype.ne.'add').and.(gwastype.ne.'add_dom')) then
          print*,'unknown GWAS_TYPE: ', gwastype
          print*,' Should be add or add_dom'
          stop 21
       endif
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'GWAS type is missing / should be among add add_dom'
          stop 21
       elseif(steptodo.eq.1) then
          gwastype=' '
       endif
    endif
    if(steptodo.ne.1) print*,'GWAS_TYPE read = ',gwastype

    call readline(io_p,x,xc,n)
    if (xc(1) /= 'SEGM_REM') then
       print*,' SEGM_REM not found'
       stop 21
    endif
    call readline(io_p,x,xc,n) ! non obligatoire a renseigner en STEP_1, obligatoire en STEP_2 et STEP_3
    if(n.ne.0) then
       segmrem=x(1)
    elseif(n.eq.0) then
       if(steptodo.ne.1) then
          print*,'Choice about markers surrounding GWAS variant to remove is missing'
          stop 21
       elseif(steptodo.eq.1) then
          segmrem=0
       endif
    endif
    !if(steptodo.ne.1) print*,'SEGM_REM read = ',segmrem


    call readline(io_p,x,xc,n)

    do
       if (xc(1).eq.'COJO') then  ! l utilisateur peut indiquer qu il veut (1) ou pas (0) l option COJO
          ! si il n indique pas COJO=0 c est qu il ne veut pas COJO
          cojofile=' '
          read(io_p,*,iostat=stat) indic_cojo
          !call readline(io_p,x,xc,n) 
          !indic_cojo=xc(1)
          if((indic_cojo.ne.0).and.(indic_cojo.ne.1)) then
             print*,'The indicator provided by user for COJO option is different from 0 and 1'
             STOP 21
          endif

          if(indic_cojo.eq.1) then
             call readline(io_p,x,xc,n)
             if (xc(1) /= 'COJO_FILE') then
                print*,'COJO_FILE not found'
                stop 21
             endif
             call readline(io_p,x,xc,n) 
             if(n.ne.0) then
                cojofile=xc(1)
             elseif(n.eq.0) then
                print*,'Choice about markers surrounding GWAS Variant to remove is missing'
                stop 21
             endif
          endif
          ! if((indic_cojo.eq.1).and.(segmrem.ne.9999999999)) then
          !    print*,' '
          !    print*,'PROBLEM : COJO only possible if ALL Markers on the chromosome carrying GWAS variants are excluded <-> segmrem=9999999999'
          !    STOP 21
          ! endif
          ! if((indic_cojo.eq.1).and.(gwastype.ne.'add')) then
          !    print*,' '
          !    print*,'PROBLEM : COJO only possible if GWAS_TYPE for Variants = ADDITIVE EFFECT <-> gwastype=add'
          !    STOP 21
          ! endif



          !if(detailed_log.eq.3) then
          !   print*,'Option COJO read = ',indic_cojo
          !   if(indic_cojo.eq.1) print*,'File for COJO = ',cojofile
          !endif

       endif
       call readline(io_p,x,xc,n) 
       exit
    enddo

    if(detailed_log.eq.3) then
       print*,' '
       print*,'Execution of subroutine read_parameters_hm before reading OPTIONS is completed'
       print*,' '
    endif

    do
       !call readline(io_p,x,xc,n)
       if (xc(1).eq.'OPTION') exit
    enddo

  end subroutine read_parameters_hm
  !========================================================================




  !========================================================================
  subroutine crea_fichpar_p3(chr,GR,SGr)
    !========================================================================
    ! reads and copy parameter file in sub-folders for PART 3
    implicit none

    integer:: chr,GR,SGr
    integer,parameter:: io_parout=68
    integer :: ioparam
    !real(kind=8) xreal(100)
    character*15 :: fx3,fx4,fx42
    character*300 :: ligne_param

    ioparam=0

    write(fx3,'(i0)') chr
    write(fx4,'(i0)') GR
    write(fx42,'(i0)') SGr

    if(indic_cojo.eq.0) then
       open(io_parout,file=trim("chr")//trim(fx3)//trim("_gr")//trim(fx4)//trim("/SG")//trim(fx42)//trim("/fichpar_part3_GrVar.par"),form="formatted")
    endif
    if(indic_cojo.eq.1) then
       open(io_parout,file=trim("cojo_")//trim(fxcojo)//trim("/chr")//trim(fx3)//trim("_gr")//trim(fx4)//trim("/SG")//trim(fx42)//trim("/fichpar_part3_GrVar.par"),form="formatted")
    endif

    rewind(io_p)

    do

       read(io_p,"(a300)",iostat=ioparam) ligne_param
       if(ioparam.ne.0) exit
       write(io_parout,"(a)") adjustl(trim(ligne_param))

       if(adjustl(trim(ligne_PARAM)).eq.'STEPS') then
          read(io_p,*,iostat=ioparam) ligne_param
          ligne_param='3'
          write(io_parout,"(a)") adjustl(trim(ligne_param))
          ligne_param='Chrom_numGR'
          write(io_parout,"(a)") adjustl(trim(ligne_param))
          ligne_param=trim(fx3)//"  "//trim(fx4)//"  "//trim(fx42)
          write(io_parout,"(a)") adjustl(trim(ligne_param))
       endif

       ! si option COJO demandee on remplace le fichier donnant l identite du variant COJO par le fichier contenant les doses pour ce variant cree en partie 2
       if(adjustl(trim(ligne_PARAM)).eq.'COJO_FILE') then
          read(io_p,*,iostat=ioparam) ligne_param
          ligne_param=trim("cojo_")//trim(fxcojo)//trim("/Doses_VarCojo_")//trim(fxcojo)//trim(".bin")
          write(io_parout,"(a)") adjustl(trim(ligne_param))
       endif

    enddo

    close(io_parout)

  end subroutine crea_fichpar_p3
  !========================================================================




  !========================================================================
  subroutine defaults_and_checks()
    !========================================================================
    ! check optional parameters and check for positive definitness
    !integer::i
    integer::n

    do

       call getoption('missing',n,x,xc)
       if(n > 0) then
          read(xc(1),*) mis
          ! print'('' * missing observation (default=-9999.0):'',F8.1)',mis
       endif

       call getoption('print_effects',n,x,xc)
       if (n > 0) then
          read(xc(1),'(i10)') iprint_eff
          print*,'! characteristics of first and last ',iprint_eff,' effects are printed'
          if(iprint_eff.lt.1) then
             print*,'Number of levels per effect asked by user to be printed < 1'
             stop 20
          endif
       endif

       ! modif du 13 avril 2018 pour fixer le seuil de MAF minimum en deca duquel on considere un SNP comme monomorphe
       call getoption('MAF_minimum',n,x,xc)
       if (n > 0) then
          read(xc(1),*) MAFmin
          ! print*,'seuil de MAF en deca duquel on considere un SNP comme monomorphe = ',MAFmin
          if((MAFmin.lt.0.0d0).or.(MAFmin.gt.0.5d0)) then
             print*,'The value for MAFmin specified by user is <0 or >0.5'
             stop 20
          endif
       endif


       ! modif du 22 mai 2019 pour choisir niveau de details dans la LOG d execution 
       call getoption('detailed_log',n,x,xc)
       if (n > 0) then
          read(xc(1),*) detailed_log
          ! print*,'Niveau de details contenus dans la LOG (1=basique / 2= qlq details / 3=maxi) :',detailed_log
          if((detailed_log.ne.1).and.(detailed_log.ne.2).and.(detailed_log.ne.3)) then
             print*,'Value for detailed_log specified by user different from 1 2 3'
             stop 20
          endif
       endif

       ! modif du 18 jan 2023 pour lire l option choisi par utilisateur pour segments de SNPpar a exclure des GWAS (strict ou window)
       call getoption('StratElim',n,x,xc)
       if (n > 0) then
          read(xc(1),*) StratElim
          ! print*,'Strategie d eimination des SNPpar autour de VarGWAS demandee (strict / window) :',StratElim
          !if((StratElim.ne.'strict').and.(StratElim.ne.'window')) then ! 14/11/2024 : option strict not used for a long time, not sure it still works after program evolutions
          if(StratElim.ne.'window') then
             print*,'Only current possible value for StratElim is window'
             stop 20
          endif
       endif

       ! modif du 25 jan 2023 pour lire le nombre max de variants GWAS par job en partie 3 choisi par utilisateur 
       call getoption('NbMaxVarGWAS',n,x,xc)
       if (n > 0) then
          read(xc(1),*) NbMaxVarGWAS
          print*,'Maximum number of GWAS variant per job in Step_3 fixed by user :',NbMaxVarGWAS
          if((NbMaxVarGWAS.le.0)) then
             print*,'Maximum number of GWAS variant per job in Step_3 fixed by user <= 0'
             stop 20
          endif
       endif

       ! modif du 26 jan 2023 pour lire l option de calcul de la variance residuelle par Variant GWAS
       call getoption('CalcVarRes',n,x,xc)
       if (n > 0) then
          read(xc(1),*) CalcVarRes
          print*,'Strategy for computing Residuel Variance for GWAS Variants in Step_3 chosen by user: ',CalcVarRes
          if((CalcVarRes.ne.'approx').and.(CalcVarRes.ne.'exact').and.(CalcVarRes.ne.'optim')) then
             print*,'Strategy for computing Residuel Variance for GWAS Variants in Step_3 chosen by user different from approx exact optim'
             stop 20
          endif
       endif


       ! modif du 26 jan 2023 pour lire la valeur mini du Test_1 d effet d 1 Variant GWAS pour qu on lui calcule VarRes exacte quand CalcVarRes=optim
       call getoption('seuil_test1',n,x,xc)
       if (n > 0) then
          read(xc(1),*) seuil_test1
          print*,'Threshold value for Test_approx to calculate exact residual variance in Step_3 when CalcVarRes=optim :',seuil_test1
          if(seuil_test1.le.0) then
             print*,'Threshold value for Test_approx (seuil_test1) negative or zero'
             stop 20
          endif
       endif


       ! modif du 27 jan 2023 pour lire la valeur de StepFichOut = taille du paquet de VarGWAS qu on sort regulierement dans fichier Resultats.txt en partie 3
       call getoption('StepFichOut',n,x,xc)
       if (n > 0) then
          read(xc(1),*) StepFichOut
          print*,'Writing in result file will be done every ',StepFichOut,' GWAS variants'
          if((StepFichOut.le.0).or.(StepFichOut.gt.NbMaxVarGWAS)) then
             print*,'Step for writing in result file (StepFichOut) is negative or larger than NbMaxVarGWAS specified by user'
             stop 20
          endif
       endif

       ! modif du 05 sep 2025 pour ecrire dan la LOG tous les individus avec genotype inconnu a au moins 1 maruquer parente
       call getoption('PrintMissGeno',n,x,xc)
       if (n > 0) then
          read(xc(1),*) PrintMissGeno
          print*,'Printing total list of animals with unknown genotype at kinship Markers: PrintMissGeno=',PrintMissGeno
          if((PrintMissGeno.ne.0).or.(PrintMissGeno.ne.1)) then
             print*,'Value for OPTION PrintMissGeno specified by user has to be 0 (default) or 1'
             stop 20
          endif
       endif

       ! modif du 27 jan 2023 pour lire la valeur de nbGWAS_test si l utilisateur veut arreter l execution des GWAS de la cellule Group x Batch
       ! apres seulement nbGWAStest variants = utile pour developper et tester des modifs de programme sans faire les GWAS sur tous les variants
       ! pour vérifier que le programme fonctionne jusqu'au bout sans traiter tous les variants du Group x Batch
       call getoption('nbGWAS_limite',n,x,xc)
       if (n > 0) then
          read(xc(1),*) nbGWAStest
          nbGWAS_limite=1
          print*,'nbGWAS_limite indicator=',nbGWAS_limite
          print*,'Nunmber of variants for which the program performs GWAS before stopping (useful for tests) with OPTION nbGWAS_limite :',nbGWAStest
          if((nbGWAStest.le.0).or.(nbGWAStest.gt.NbMaxVarGWAS)) then
             print*,'Value given by user for OPTION nbGWAS_limite is negative or bigger than the total number of Variants in the Group x Batch cell'
             stop 20
          endif
       endif

       ! modif du 17 fev 2023 pour lire la valeur de nbSNPparAuto si l utilisateur veut que le programme choisisse automatiquement les SNPparente a utiliser
       call getoption('SelAutoSNPpar',n,x,xc)
       if (n > 0) then
          read(xc(1),*) nbSNPparAuto
          SelAutoSNPpar=1
          print*,'Upper bound of Marker number to be automatically chosen via OPTION SelAutoSNPpar :',nbSNPparAuto
          if(nbSNPparAuto.le.0) then
             print*,'nbSNPparAuto value given by user is negative: ',nbSNPparAuto
             stop 20
          endif
       endif

       ! modif du 17 mars 2023 pour lire la valeur de MaxMem si l utilisateur veut modifier la valeur par defaut (10Go)
       call getoption('MaxMem',n,x,xc)
       if (n > 0) then
          read(xc(1),*) MaxMem
          print*,' '
          print*,'Upper bound of RAM usable for reading Variant Genotype / allelic dosage per pass in Step_2:',MaxMem
          print*,'NB: this parameter only limits the Table allocation size during Step_2, not a System parameter'
          print*,' To ajust system memory usable by the program, change -hvmem value when launching the program'
          if (MaxMem.le.0) then
             print*,'MaxMem value given by user is negative or zero: ',MaxMem
             stop 20
          endif
       endif

       ! modif du 20 mars 2023 pour lire la valeur de meth_mpm si l utilisateur ne veut pas (=1) ou veut (=2) stocker les genotypes centres aux SNPpar en memoire
       call getoption('meth_MpM',n,x,xc)
       if (n > 0) then
          read(xc(1),*) meth_mpm
          print*,' '
          print*,'Storage method chosen for Marker genotypes: non centered integers (=1) or centered reals (=2) :',meth_mpm
          print*,'NB : 1 = lower memory usage by slower because genotypes must be centered when genotypes are needed'
          if ((meth_mpm.ne.1).and.(meth_mpm.ne.2)) then
             print*,'meth_MPM value is different from 1 and 2: ',meth_mpm
             stop 20
          endif
       endif

       ! modif du 28 juil 2023 pour convertir les DOSES des variants GWAS en genotypes discrets = nb alleles 2 portes par individus
       call getoption('doses_to_geno',n,x,xc)
       if (n > 0) then
          read(xc(1),*) doses_to_geno
          print*,'Conversion from Variant allelic dosages to discrete genotypes 0 1 2 (0 = no conversion = default / 1 = conversion) :',doses_to_geno
          if ((doses_to_geno.ne.0).and.(doses_to_geno.ne.1)) then
             print*,'doses_to_geno value is different from 0 and 1 : ',doses_to_geno
             stop 20
          endif
       endif

       call getoption('end',n,x,xc)
       if(n <= 0) exit
    enddo

  end subroutine defaults_and_checks
  !========================================================================



  !========================================================================
  subroutine chkfmtnewline(un)
    !========================================================================
    character(50000) :: a,i_name ! big number to allow check genotypes files
    integer :: i,un
    !
    inquire(un,name=i_name)
    read(un,'(a)') a 
    do i=1,len(a) 
       if (ichar(a(i:i)) == 13) then ! Mac files
          print '(/a,a,a)', ' Mac old newline format detected in file: "',trim(i_name(index(i_name,'/',back=.true.)+1:)),'"'
          print '(a)'  , ' Convert it to Unix newline format, e.g. "flip -u <file>"'
          stop 21
       endif
    enddo

    if  (scan(achar(9),a)/=0) then 
       print '(/a,a,a)',' TAB was found, Use spaces to separeate columns in file: "',trim(i_name(index(i_name,'/',back=.true.)+1:)),'"'
       print '(a)'     ,' Convert to TAB to spaces, e.g. "expand file > newfile"'
       stop 21
    endif
    rewind(un)
  endsubroutine chkfmtnewline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !========================================================================




  !======================================================================== 
  subroutine readline(unit,x,xc,n,nohash)
    !========================================================================
    ! Reads a line from unit and decomposes it into numeric fields in x and
    ! alphanumeric fields in xc; the number of fields is stored in n.
    ! If optional variable nohash is absent, 
    !  all characters after character # are ignored.

    integer unit,n,stat,i
    real(kind=8)::x(:)
    integer,optional::nohash
    character xc(:)*(*)
    character(len=max_string_readline) :: a   

    do
       read(unit,'(a)',iostat=stat)a
       if (stat /= 0) then
          n=-1
          return
       endif
       a=adjustl(a)               ! ignore leading spaces

       if (present(nohash)) then 
          exit
       else	!ignore characters after #
          i=scan(a,'#')
          if (i == 0) exit
          if (i>1) then
             a=a(1:i-1)
             exit
          endif
       endif
    enddo

    !call nums(a,x,xc,size(x),n)
    call nums2(a,n,x,xc)
    return
  end subroutine readline
  !======================================================================== 



  !======================================================================== 
  subroutine nums2(a,n,x,xc)
    !======================================================================== 
    ! separates array a into items delimited by blanks. character elements are
    ! put into optional character vector xc, decoded numeric values 
    ! into optional real vector x, and n contains the number of items. The 
    ! dimension of x and xc can be lower than n.
    ! A modification of nums() from f77 to f90
    ! Now accepts real numbers
    ! 2/23/2000

    character (*)::a
    character (*),optional::xc(:)
    real(kind=8),optional::x(:)
    integer::n,curr,first,last,lena,stat,i

    curr=1;lena=len(a);n=0
    if (present(xc)) xc=' '
    if (present(x)) x=0

    do 
       ! search for first nonspace
       first=0
       do i=curr,lena
          if (a(i:i) /= ' ') then
             first=i
             exit
          endif
       enddo
       if (first == 0) exit


       ! search for first space
       curr=first+1
       last=0
       do i=curr,lena
          if (a(i:i) == ' ') then
             last=i
             exit
          endif
       enddo

       if (last == 0) last=lena

       n=n+1
       if (present(xc)) then
          if (size(xc) >= n) then
             xc(n)=a(first:last)
          else
             print*, "Error in nums2 splitting string starting: "
             print*,a(1:80)
             print*,'into m:',size(xc),'items'
             stop 21
          endif
       endif
       if (present(x)) then
          if (size(x) >= n) then
             read(a(first:last),'(f12.0)',iostat=stat)x(n)
             if (stat /=0) x(n)=0
          else
             print*, "Error in nums2 splitting string starting: "
             print*,a(1:80)
             print*,'into m:',size(xc),'items'            
             stop 21
          endif
       endif

       curr=last+1
    enddo
  end subroutine nums2
  !======================================================================== 






  !======================================================================== 
  subroutine nb_chps_ttr(unit,n)
    !======================================================================== 
    ! separates array a into items delimited by blanks. character elements are
    ! put into optional character vector xc, decoded numeric values 
    ! into optional real vector x, and n contains the number of items. The 
    ! dimension of x and xc can be lower than n.
    ! A modification of nums() from f77 to f90
    ! Now accepts real numbers
    ! 2/23/2000

    integer unit,n,i
    !character(len=max_string_readline) :: a   
    character*250::a
    integer::curr,first,last,lena

    read(unit,"(a30)")a
    !a=adjustl(a)               ! ignore leading spaces
    a=adjustl(adjustr(a))

    curr=1
    lena=len(a)
    n=0

    do 
       ! search for first nonspace and first non tab
       first=0
       do i=curr,lena
          if ((a(i:i) /= ' ').and.(a(i:i).ne.tab)) then
             first=i
             exit
          endif
       enddo
       if (first == 0) exit


       ! search for first space or first tab
       curr=first+1
       last=0
       do i=curr,lena
          if ((a(i:i) == ' ').or.(a(i:i) .eq. tab)) then
             last=i
             exit
          endif
       enddo

       if (last == 0) last=lena

       n=n+1

       curr=last+1
    enddo

  end subroutine nb_chps_ttr
  !======================================================================== 




  !========================================================================  
  subroutine getoption(typnam,n,x,xc,norewind)
    !========================================================================  
    ! In unit io_p=40, which reads parameter line in BLUPF90, locates line:
    ! OPTION typnam str1 str2
    ! where str1, str2 are strings separated by spaces.
    ! Then, it assigns: xc(1)=str1,  xc(2)=str2,...
    ! and attempts to decode strings into real values: x1=value(str1),....
    !
    ! n contains the number of strings.  x and xc are optional and their 
    ! dimensions may be smaller than n in which case some strings/values are
    ! not stored.
    !
    ! Upon exit, unit io_p=40 points to line next to the one located.
    !
    ! If the line cannot be located, n=-1
    ! 
    ! if present norewind continue reading from last line in file
    !

    character (*)::typnam
    integer::n
    real(kind=8),optional::x(:)
    character (*),optional::xc(:)
    integer,optional :: norewind
    real(kind=8)::x1(1000)
    integer::stat,m
    character (400)::xc1(1000)
    character (400)::a

    if (.not. present(norewind)) rewind io_p
    n=-1

    do
       read (io_p,'(a)',iostat=stat)a
       if (stat /= 0) exit
       call nums2(a,m,x1,xc1)       
       if ( m>1 .and. xc1(1) == 'OPTION' .and. xc1(2) == typnam) then
          n=m-2
          if (present(xc)) xc=xc1(3:size(xc)+2)
          if (present(x)) x=x1(3:size(x)+2)
          exit
       endif
    enddo
  end subroutine getoption
  !========================================================================  



  !========================================================================  
  subroutine chkgetoption(xa)
    !========================================================================  
    ! check all arguments to be in the list of valid ones. 
    character (*) :: xa(:)
    real(kind=8)::x1(1000)
    integer::stat,m,i
    character (50)::xc1(1000)
    character (400)::a
    logical :: found
    rewind io_p
    do
       read (io_p,'(a)',iostat=stat)a
       if (stat /= 0) exit
       call nums2(a,m,x1,xc1)       
       if ( m>1 .and. xc1(1) == 'OPTION') then
          found=.false.
          do i=1,size(xa)
             if (trim(xa(i))==trim(xc1(2)) ) then
                found=.true. 
                exit
             endif
          enddo
          if (.not. found) then 
             print '(/a)'  ,' ***** OPTION NOT valid for this program !!!! '
             print '(a,a,a)' ,'       "',trim(adjustl(xc1(2))),'"'
             print '(a,a)' ,'       Check for typing errors (e.g. upper/lower cases) in the key'
             stop 21
          endif
       endif

    enddo
  end subroutine chkgetoption
  !========================================================================  





  !========================================================================  
  subroutine print_parameters_hm
    !========================================================================  
    ! Prints parameters of the model
    !use denseop
    implicit none
    integer :: i,j,k,l
    character (20):: typnam
    character (40):: fmt
    !logical::not_pd

    testpresfic=.true.

    !on verifie que le fichier performances existe bien
    INQUIRE (FILE=trim(datafile),EXIST=testpresfic)
    if(testpresfic.eq..false.) then
       print*,'DATA file does not exist'
       STOP 21
    endif


    print*,' '
    print '('' EXEC_PATH:'',t90,a)',trim(execpath)
    print '('' Folder containing binary files that will be created:'',t90,a)',trim(location)

    print*,' '
    print '('' Data file:'',t90,a)',trim(datafile)
    !on verifie que le fichier performances existe bien
    INQUIRE (FILE=trim(datafile),EXIST=testpresfic)
    if(testpresfic.eq..false.) then
       print*,'DATA file does not exist'
       STOP 21
    endif

    print '('' Marker genotype file used to model breeding values: '',t90,a)',trim(typparfile)
    !on verifie que le fichier TYPAGES pour la parente existe bien
    INQUIRE (FILE=trim(typparfile),EXIST=testpresfic)
    if(testpresfic.eq..false.) then
       print*,'Marker Genotype file does not exist'
       STOP 21
    endif

    print '('' Format of Marker genotype file used to model breeding values: '',t90,a)',trim(fmttyppar)


    print '('' Map file for markers for breeding values: '',t90,a)',trim(mapparfile)
    !on verifie que le fichier CARTE pour la parente existe bien
    INQUIRE (FILE=trim(mapparfile),EXIST=testpresfic)
    if(testpresfic.eq..false.) then
       print*,'Marker MAP file does not exist'
       STOP 21
    endif

    ! Verifications realisees uniquement pour STEP_2 ou STEP3 car information facultative en STEP_1
    if(steptodo.ne.1) then

       print '('' Allelic dosage / Genotype file for GWAS variants: '',t90,a)',trim(typgwasfile)
       !on verifie que le fichier TYPAGES pour les GWAS existe bien
       INQUIRE (FILE=trim(typgwasfile),EXIST=testpresfic)
       if(testpresfic.eq..false.) then
          print*,'Variant genotype/allelic dosage file does not exist'
          STOP 21
       endif

       print '('' Format of genotype file for GWAS: '',t90,a)',trim(fmttypgwas)
       if((trim(fmttypgwas).ne.'minimac').and.(trim(fmttypgwas).ne.'typ_eval').and.(trim(fmttypgwas).ne.'plink')) then
          print*,'Unknown format for GWAS variant genotype/allelic dosage file - has to be minimac or typ_eval or plink'
          STOP 21
       endif


       print '('' Map file for Variants considered for GWAS: '',t90,a)',trim(mapgwasfile)
       !on verifie que le fichier CARTE pour la parente existe bien
       INQUIRE (FILE=trim(mapgwasfile),EXIST=testpresfic)
       if(testpresfic.eq..false.) then
          print*,'GWAS Variant MAP file does not exist'
          STOP 21
       endif

       print'('' Id of chromosome carrying the variants for GWAS: '',t90,i3)',NUMCHR_GWAS
       if(NUMCHR_GWAS.lt.1) then
          print*,'Id of chromosome carrying the variants for GWAS is incorrect'
          stop 21
       endif


       print*,' '
       print '('' GWAStype:'',t90,a)',trim(gwastype)

       if(indic_cojo.eq.1) then

          print '('' Indicator for OPTION COJO: indic_cojo='',t90,i3)',indic_cojo

          print '('' File giving the ID of COJO variant: '',t90,a)',trim(cojofile)
          !on verifie que le fichier COJO existe bien
          INQUIRE (FILE=trim(cojofile),EXIST=testpresfic)
          if(testpresfic.eq..false.) then
             print*,'File giving the ID of COJO variant does not exist'
             STOP 21
          endif
       endif

    endif
    ! Fin des verifications realisees uniquement pour STEP_2 ou STEP3 car information facultative en STEP_1

    print*,' '
    if(detailed_log.eq.3) print '('' Number of Traits'',t90,i3)',ntrait

    print '('' Number of Effects'',t90,i8)',neff

    print '('' Position of alphanum Id of individuals in DATA file'',t90,i4)',posIdchar

    print '('' Position of Phenotypes'',t90,20i3)',(pos_y(i),i=1,ntrait)

    print '('' Position of Weights (if any)'',t90,20i3)',(pos_weight(i),i=1,ntrait)

    print '('' Value of Missing Trait/Observation'',t90,F8.1)',mis

    print '('' MAF under wich a marker is considered monomorphous and is excluded (default=1d-3)= '',t110,D14.6)',MAFmin

    print '('' Storage method chosen for Marker genotypes: non centered integers (=1) or centered reals (=2) : '',t110,i3)',meth_mpm


    !print '('' Nombre de coeurs pour factorisation Ann par MKL = '',t90,i3)',nb_threads_mkl

    print*,' '
    print '(/,''NON GENETIC EFFECTS DECLARED IN PARAMETER FILE'')'
    print '(a,t25,a,t45,a)',' #  type','position (2)','levels   [positions for nested]'
    !print*,'iprint_eff=',iprint_eff
    do i=1,neff
       select case (effecttype(i))
       case (effcross)
          typnam='cross-classified'
       case (effcov)
          typnam='covariable'
       case default
          typnam='???'
       end select

       ! writing with number of fields dependent on ntrait and presence of cov
       write(fmt,'(''(i4,2x,a19,'',i4,'',i4,t40,i10,t65,20i3)'')') ntrait
       if (i.le.iprint_eff .or. i.ge.neff+1-iprint_eff) then
          if (nestedcov(i,1) /= 0) then
             !a l origine VDU print fmt,i,typnam,pos_eff(i,1),nlev(i),nestedcov(i,1)
             !print fmt,i,typnam,pos_eff(i,:),nlev(i),nestedcov(i,:)
             print *,i,typnam,pos_eff(i,:),nlev(i),nestedcov(i,:)
          else
             !  print '(''(i4,5x,a19,5x,i6,t40,i6)'')',i,typnam,pos_eff(i,1),nlev(i)
             !a l origine VDU print *,i,typnam,pos_eff(i,1),nlev(i)
             print *,i,typnam,pos_eff(i,:),nlev(i)
          endif
       else if (i==iprint_eff+1) then
          print'(''  ..................'',/)'   
       endif
    enddo

    do i=1,neff
       if (i>iprint_eff +1 .and. i<neff+1-iprint_eff) cycle
       if (i==iprint_eff+1) print'(''  ..................'')'  
       if (randomnumb(i) == 0)then
          cycle
       else if (randomnumb(i) == 1) then
          print '(/,'' Random Effect(s)'',t20,20i6)',(/(j,j=i,i+randomnumb(i)-1)/)
       else if(randomnumb(i) >1) then
          print '(/,'' correlated random effects'',t30,20i3)',&
               (/(j,j=i,i-1+randomnumb(i)) /)
       endif
       select case (randomtype(i))
       case (g_diag)
          !if (rep(i)) then
          !    print '('' Type of Random Effect:'',t30,a)','diagonal rep'
          !else
          print '('' Type of Random Effect:'',t30,a)','diagonal'
          !endif

          !case(g_AD)
          !   print '('' Type of Random Effect:'',t30,a,a)','additive animal direct'
          !   if (randomnumb(i) .eq. 2) then
          !      if (randomtype(i+1) .eq. g_AM) then  
          !         print '('' Type of Random Effect:'',t30,a,a)','additive animal maternal'
          !      endif
          !   endif
       case default
          print*,'unknown RANDOM-TYPE'
          stop 21
       end select

       j=ntrait*randomnumb(i)   
       print '('' trait   effect    (CO)VARIANCES'')'

       do k=0,randomnumb(i)-1
          do j=1,ntrait
             l=j+k*ntrait
             print '(i3,i8,2x,50g12.4)',j,k+i,rand(i,l,1:randomnumb(i)*ntrait)
          enddo
       enddo
    enddo

    print '('' Residual Variance read:'',t90,f15.6)',r(1,1)
    print '('' Genetic Variance read:'',t90,f15.6)',g(1,1)


  end subroutine print_parameters_hm
  !========================================================================  



  !========================================================================  
  function checksym8(x) result (a) 
    !========================================================================  
    ! check that matrix x is symmetric
    real(kind=8) :: x(:,:) 
    integer :: i,j
    logical :: a 
    !
    a=.true. 
    do i=1,size(x,1)
       do j=i,size(x,1)
          if (x(i,j) /= x(j,i) ) then 
             a = .false.
             exit 
          endif
       enddo
    enddo
  end function checksym8
  !========================================================================  




  !========================================================================  
  function normal_ttr(x, dx)  result (pval_ttr)
  !========================================================================  

    ! the function normal_ttr calculates the P-value associated to the null hypothesis "x diff 0"
    ! x is the test staistic
    ! dx is the step for approximate calculation of surface under the standard centered normal distribution (0.0001 by default)
    ! Pval_ref is a vector containing precalculated bilateral P-value for integers from 1 to 37
    ! The function approximates the surface Sx under the normal distribution curve from abs(x) to ceiling(abs(x)) with trapeze approxiamtion
    ! then the bilateral Pvalue for x is (2*Sx)+Pval_ref(ceiling(abs(x)))
    ! using a smaller dx increases approximation but increases the number of "slices" and therfore the computing time

    real(kind=8), intent(in) :: x, dx
    real(kind=8) :: t, sum, f, invden
    integer :: n, i, ref
    real(kind=8) :: abs_x,pval_ttr,y
    real(kind=8)::pval_ref(37)  ! vecteur des valeurs de Pvalue dans une loi normale pour les entiers de 1 a 37
    save pval_ref
    save invden

    data pval_ref/0.317310507862914, 4.550026389635837d-02, 2.699796063260190d-03, 6.334248366623979d-05, 5.733031437583881d-07, 1.973175290075398d-09, &
    2.559625087771674d-12, 1.244192114854365d-15, 2.257176811907689d-19, 1.523970604832098d-23, 3.821319148997350d-28, 3.552964224155365d-33, &
    1.223432879909973d-38, 1.558707363838535d-44, 7.341932398625426d-51, 1.277750880107608d-57, 8.211992404197855d-65, 1.948189783787418d-72, &
    1.705444790526215d-80, 5.507248237212310d-89, 6.558556037957993d-98, 2.879784870290108d-107, 4.661274012441260d-117, 2.780784237099379d-127, &
    6.113393412765075d-138, 4.952126631006798d-149, 1.477896201376999d-160, 1.624773893931854d-172, 6.579570533408875d-185, 9.813427854296267d-198, &
    5.390500162400810d-211, 1.090416120702428d-224, 8.122371241832289d-239, 2.227797571148807d-253, 2.249821412944903d-268, 8.365248131594677d-284, &
    1.145114244504795d-299 /

    data invden /0.3989422804014326/ ! equivalent a invden= 1.0d0 / (sqrt(2.0d0 * 3.141592653589793d0))

    if(x.eq.0.d0) then 
       pval_ttr=1.0d0
       return
    endif

    ! if x is negative we work on abs(x)
    abs_x=abs(x)

    if(abs_x.gt.37) then ! already very close to the accuracy limit of double precision reals
       pval_ttr= Pval_ref(37)
       return
    endif

    pval_ttr=0.0d0

    sum=0.0d0

    ref=ceiling(abs_x)

    y=0.0d0

    n = int( (dfloat(ref) - abs_x) / dx) ! number of slices for approximation of complement of surface under the normal curve from abs_x to ceiling(abs_x)

    if(n.eq.0) then  ! pour le cas tres peu probable que la statistique de test soit un entier entre 1 et 37
       pval_ttr = Pval_ref(ref)  
       return
    endif
       
    if (n.ne.0) then

       sum = 0.5d0 * exp(-0.5d0 * (abs_x*abs_x)) ! initialisation of the complement of surface

       y=abs_x

       do i=1,n-1
          y=y+dx
          sum = sum + exp(-0.5d0 * y * y)
       enddo

       y=y+dx    ! y = abs_x + n*dx

       sum = sum + 0.5d0* exp(-0.5d0 * y * y) ! closing the right side of trapeze from abs_x+(n-1)*dx to abs_x+n*dx

       sum = sum * dx ! base of the n trapezes

       ! add surface of ultimate trapeze from abs_x+n*dx to ref if necessary
       if(y.ne.dfloat(ref)) sum = sum + 0.5d0 * (  exp(-0.5d0 * dfloat(ref) * dfloat(ref)) + exp(-0.5d0 * y * y) ) * ( dfloat(ref) - y )

       ! divide sum by sqrt(2 * PI)
       sum = sum * invden

       pval_ttr = 2.0d0*sum + Pval_ref(ref) ! 2*complement of surface bicause bilateral test ; Pval_ref is already bilateral

    endif

    !if(abs_x.gt.37) pval_ttr= 0.5d0* ( Pval_ref(37) + (exp(-0.5d0 * (abs_x**2))*invden) ) * (abs_x - 37)  ! au dela de la derniere valeur de reference on suppose que la fonction est linéaire

  end function normal_ttr
  !========================================================================  





  !========================================================================  
function normal_ttr2(x, dx) result(pval_ttr)
!========================================================================  

    ! This function computes the bilateral P-value for a test statistic x
    ! under the standard normal distribution using trapezoidal approximation.

    real(kind=8), intent(in) :: x, dx
    real(kind=8) :: abs_x, pval_ttr, y, sum
    integer :: n, i, ref
    real(kind=8) :: pval_ref(37), invden
    save pval_ref, invden

    ! Precalculated P-values for integers 1 to 37 under standard normal distribution
    data pval_ref /0.317310507862914, 4.550026389635837d-02, 2.699796063260190d-03, 6.334248366623979d-05, 5.733031437583881d-07, 1.973175290075398d-09, &
                   2.559625087771674d-12, 1.244192114854365d-15, 2.257176811907689d-19, 1.523970604832098d-23, 3.821319148997350d-28, 3.552964224155365d-33, &
                   1.223432879909973d-38, 1.558707363838535d-44, 7.341932398625426d-51, 1.277750880107608d-57, 8.211992404197855d-65, 1.948189783787418d-72, &
                   1.705444790526215d-80, 5.507248237212310d-89, 6.558556037957993d-98, 2.879784870290108d-107, 4.661274012441260d-117, 2.780784237099379d-127, &
                   6.113393412765075d-138, 4.952126631006798d-149, 1.477896201376999d-160, 1.624773893931854d-172, 6.579570533408875d-185, 9.813427854296267d-198, &
                   5.390500162400810d-211, 1.090416120702428d-224, 8.122371241832289d-239, 2.227797571148807d-253, 2.249821412944903d-268, 8.365248131594677d-284, &
                   1.145114244504795d-299 /

    data invden /0.3989422804014326/ ! 1 / sqrt(2 * PI)

    ! Special case: x = 0
    if (x == 0.d0) then
        pval_ttr = 1.0d0
        return
    endif

    ! Compute absolute value of x
    abs_x = abs(x)

    ! Special case: x > 37 (beyond precomputed P-values)
    if (abs_x > 37.d0) then
        pval_ttr = pval_ref(37)
        return
    endif

    ! Initialize variables
    ref = ceiling(abs_x)
    n = int((dfloat(ref) - abs_x) / dx)

    ! Special case: abs_x is an integer
    if (n == 0) then
        pval_ttr = pval_ref(ref)
        return
    endif

    ! Compute the complement of the surface under the curve
    sum = 0.5d0 * exp(-0.5d0 * abs_x * abs_x) ! Start with the first half-trapezoid
    y = abs_x

    do i = 1, n-1
        y = y + dx
        sum = sum + exp(-0.5d0 * y * y)
    enddo

    ! Add the final slice up to y = abs_x + n * dx
    y = y + dx
    sum = sum + 0.5d0 * exp(-0.5d0 * y * y)

    ! Scale sum by the step size and normalize
    sum = sum * dx * invden

    ! Handle the last segment from y to ref, if necessary
    if (y /= dfloat(ref)) then
        sum = sum + 0.5d0 * (exp(-0.5d0 * ref*ref) + exp(-0.5d0 * y*y)) * (dfloat(ref) - y) * invden
    endif

    ! Compute final bilateral P-value
    pval_ttr = 2.0d0 * sum + pval_ref(ref)

end function normal_ttr2
!========================================================================






  !========================================================================  
  function MV(i,j,NL,NC)
    !========================================================================  
    ! donne en resultat la position MV de l element i,j de la matrice M de dim (NL*NC) dans le vecteur m en ColMajor pour BLAS/LAPACK
    integer :: i,j,NL,NC
    integer(kind=8)::MV
    !
    if((i.le.0).or.(i.gt.NL).or.(j.le.0).or.(j.gt.NC)) then
       print*,'ERROR : one matrix element has coordinates out of the matrix: ',i,j,NL,NC
    endif
    MV = int( dfloat(i) +   dfloat(NL)*dfloat((j-1))  )
    !
  end function MV
  !========================================================================  



  !========================================================================  
  function TI(i,j,NL)
    !========================================================================  
    ! donne en resultat la position TI de l element iti,jti de la matrice symetrique M de dim (NL*NL) dans le vecteur m en Condense Triangulaire inferieur pour BLAS/LAPACK
    integer :: i,j,iti,jti,kti,NL
    integer(kind=8)::TI
    !
    iti=0
    jti=0
    kti=0
    if((i.le.0).or.(i.gt.NL).or.(j.le.0).or.(j.gt.NL)) then
       print*,'ERROR in fonction TI: one matrix element has coordinates out of the matrix: ',i,j,NL
    endif
    iti=i
    jti=j
    if(jti.gt.iti) then
       kti=iti
       iti=jti
       jti=kti
    endif
    TI = int( dfloat((jti-1)) * dfloat(((2*NL)-jti)) * 0.5d0) + iti
    !
  end function TI
  !========================================================================  






  !========================================================================  
  subroutine ecriMat_inf(A,NL,NE)
    !========================================================================  

    ! subroutine permettant d ecrire dans fichier log une matrice carree symetrique de NL lignes stockee sous forme Triangulaire inferieure

    integer, intent(in) :: NL ! modif JVDP
    real(kind=8),intent(in) :: A(NE)
    real(kind=8)::vecprinttemp(min(NL,20))
    integer::i,j,NE ! modif JVDP

    if(NL.gt.20) then
       print*,'The squared matrix is bigger than 20 x 20 --> orint only 1:20 x 1:20 block'
    endif

    print*,' '
    print*,'********************************'
    do i=1,min(NL,20)
       do j=1,min(NL,20)
          vecprinttemp(j)=A(TI(i,j,NL))
       enddo
       print*,vecprinttemp(:)
    enddo
    print*,'********************************'
    print*,' '

  end subroutine ecriMat_inf
  !========================================================================  





  !========================================================================  
  subroutine ecriBlock_inf(A,NL,NE8,NK)
    !========================================================================  

    ! subroutine permettant d ecrire dans fichier log le bloc (1:NK*1:NK) d une matrice carree symetrique M de NL lignes stockee sous forme Triangulaire inferieure dans A(1:NE8)

    ! NE8 = nombre d elements dans le vecteur contenant la triangulaire inferieure = NL * (NL+1) *0.5
    ! NL = nombre de ligne de la matrice carree stockee sous forme de triangulaire inferieure dans le vecteur A
    ! NK = nombre de 1eres lignes et de 1eres colonnes de M qu on veut ecrire

    !real(kind=8),intent(in) :: A(NE8)
    integer, intent(in) :: NL, NK   ! modif JVDP
    integer(kind=8)::NE8
    real(kind=8) :: A(NE8)
    real(kind=8)::vecprinttemp(min(NK,20))
    integer::i,j   ! modif JVDP

    if(NK.gt.20) then
       print*,'The squared Block-mmatrix is bigger than 20 x 20 --> print only 1:20 x 1:20 block'
    endif

    print*,' '
    print*,'********************************'
    do i=1,min(NK,20)
       do j=1,min(NK,20)
          vecprinttemp(j)=A(TI(i,j,NL))
       enddo
       print*,vecprinttemp(:)
    enddo
    print*,'********************************'
    print*,' '

  end subroutine ecriBlock_inf
  !========================================================================  












  !========================================================================  
  subroutine invmatlapack(A,N)
    !========================================================================  

    real(kind=8),intent(inout) :: A(N,N)
    !real(kind=8) :: work(N,N)
    real(kind=8) :: work(N) !work should be a 1D array ! modif JVDP
    integer:: IPIV(N)
    integer::N,info

    ! factorisation LU de A
    info=0
    call dgetrf(N,N,A,N,IPIV,info)

    !print*,'apres dgetrf info=',info

    ! inversion de la matrice apres factorisation LU
    info=0
    call dgetri(N,A,N,IPIV,WORK,N,info)

    !print*,'apres dgetri info=',info

  end subroutine invmatlapack
  !========================================================================  



  !========================================================================  
  subroutine invmatlapack_ti(A,N)
    !========================================================================  

    real(kind=8),intent(inout) :: A(N,N)
    !real(kind=8) :: work(N,N)
    real(kind=8) :: work(N) !work should be a 1D array  ! modif JVDP
    integer:: IPIV(N)
    integer::N,info

    ! factorisation LU de A
    info=0
    call dsytrf('L',N,A,N,IPIV,WORK,N,info)

    !print*,'apres dgetrf info=',info

    ! inversion de la matrice apres factorisation LU
    info=0
    call dsytri('L',N,A,N,IPIV,WORK,info)

    !print*,'apres dgetri info=',info

  end subroutine invmatlapack_ti
  !========================================================================  





  !========================================================================  
  subroutine invmatlapack_inf(A,N,dimA)
    !========================================================================  

    real(kind=8),intent(inout) :: A(dimvec)
    real(kind=8) :: work(2*N)
    integer:: IPIV(N)
    integer::N,info
    integer(kind=8)::dimA

    work=0.0d0

    print*,' '
    print*,'invmatlapack_inf : dimA=',dimA

    ! factorisation LU de A
    info=0
    !call dsptrf('L',N,A,IPIV,work,info)
    call dsptrf('L',N,A,IPIV,info) !work not in the standard API ! modif JVDP

    if(detailed_log.eq.3) print*,'After dsptrf info=',info

    !print*,'resultat factorisation'
    !print*,A(:)
    !print*,' '


    ! inversion de la matrice apres factorisation LU
    info=0
    work=0.0d0
    call dsptri('L',N,A,IPIV,work,info)

    if(detailed_log.eq.3) print*,'After dsptri info=',info

  end subroutine invmatlapack_inf
  !========================================================================  





  !========================================================================  
  subroutine decode_record_effix
    !========================================================================  
    ! decodes data record into y's and weight and NivAniv = levels of fixed effects for Animal
    ! and weight_cov

    integer::i,j

    if (pos_weight(1) /= 0) then
       do j=1,ntrait
          weight_y(j) = indata(pos_weight(j))
       enddo
    else
       weight_y=1.0d0
    endif

    do j=1,ntrait
       y(j)=indata(pos_y(j))
    enddo

    do i=1,neff
       ! effcross
       if ((effecttype(i).eq.0) .and. (randomtype(i).eq.1)) then
          do j=1,ntrait
             if(pos_eff(i,j).ne.0) then ! ajout test sur position differente de 0 car plantage avec option -CB car indata(0)
                Nivanim(i,j)=int(indata(pos_eff(i,j)))
                if(Nivanim(i,j).ne.0) weight_cov(i,j)=1.0d0
                if(Nivanim(i,j).eq.0) weight_cov(i,j)=0.0d0
             endif
             if(pos_eff(i,j).eq.0) then
                Nivanim(i,j)=0
                weight_cov(i,j)=0.0d0 ! pour mettre a 0 cov pour effet absent du modele du car
             endif
             if (Nivanim(i,j).gt.nlev(i)) then
                print*,'Some levels of non genetic effects in DATA file are higher than the number of levels specified in parameter file'
                STOP 20
             endif
          enddo
       endif
       ! effcov 
       if ((effecttype(i).eq.1) .and. (randomtype(i).eq.1)) then
          do j=1,ntrait
             if ((nestedcov(i,j).ne.0).and.(pos_eff(i,j).ne.0)) then
                Nivanim(i,j)=int(indata(nestedcov(i,j)))
                weight_cov(i,j)=indata(pos_eff(i,j))
             endif

             if ((nestedcov(i,j) .eq. 0).and.(pos_eff(i,j).ne.(0))) then  ! effcov
                Nivanim(i,j)=1
                weight_cov(i,j)=indata(pos_eff(i,j))
             endif
             if(pos_eff(i,j).eq.0) then
                Nivanim(i,j)=0
                weight_cov(i,j)=0.0d0
             endif
          enddo
       endif
       ! random_diagonal_effect
       if ((effecttype(i).eq.0) .and. ((randomtype(i).eq.2).or.(randomtype(i).eq.3).or.(randomtype(i).eq.4))) then
          do j=1,ntrait
             if(pos_eff(i,j).ne.0) then
                Nivanim(i,j)=int(indata(pos_eff(i,j)))
                if(Nivanim(i,j).ne.0) weight_cov(i,j)=1.0d0
                if(Nivanim(i,j).eq.0) weight_cov(i,j)=0.0d0
             endif
             if(pos_eff(i,j).eq.0) then
                Nivanim(i,j)=0
                weight_cov(i,j)=0.0d0 
             endif
             if (Nivanim(i,j).gt.nlev(i)) then 
                print*,'Some levels of non genetic effects in DATA file are higher than the number of levels specified in parameter file'
                STOP 20
             endif
          enddo
       endif
    enddo

  end subroutine decode_record_effix
  !========================================================================  



  !========================================================================  
  function address1(e,l,t)
    !========================================================================  
    ! returns address for level l of effect e and trait t
    integer :: e,l,t, address1,i
    logical,save::first=.true.
    integer,allocatable,save::offset(:)
    !
    if (first) then
       first=.false.
       allocate(offset(neff))
       do i=1,neff
          offset(i)=sum(nlev(1:i-1))*ntrait
       enddo
    endif
    !
    address1= offset(e)+(l-1)*ntrait+t
  end function address1
  !========================================================================  




  !========================================================================  
  subroutine STAT_NUM (vec, n, moy_doses, std_doses)
    !========================================================================  
    integer i,n,n_notmiss
    real(kind=8) vec(n)
    real*8 som, somcar, moy_doses, var, std_doses

    som=0.0d0
    somcar=0.0d0
    moy_doses=0.0d0
    var=0.0d0
    std_doses=0.0d0

    if(n.eq.0) return         
    n_notmiss=0
    do i=1,n
       som=som+vec(i)
       somcar=somcar+vec(i)*vec(i)
       n_notmiss=n_notmiss+1
    end do
    moy_doses=som/dfloat(n_notmiss)
    var=somcar/dfloat(n_notmiss) - (moy_doses*moy_doses)
    if(var.ge.0.d0) then
       std_doses=dsqrt(var)
    else
       std_doses=0.d0
    end if

    return
  end subroutine STAT_NUM

  !========================================================================  






  !========================================================================
  subroutine vversion(cdate,anlim,molim,jlim)
    !========================================================================

    implicit none
    character(len=128) :: cdate
    integer :: annee,imo,anlim,i,molim,jlim,diff,jour
    character(len=10) :: cal,mois,heure
    character(len=3) :: mm(12)
    !logical :: francais,english

    data mm/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

    read (cdate,*) cal,mois,jour,heure,annee 
    imo=0
    do i=1,12
       if (mois.eq.mm(i)) imo=i
    end do
    if (imo.eq.0) then
       imo=1
       print *,'probleme dans la date, Janvier suppose'
    end if

    call ecdate(diff,anlim,molim,jlim,annee,imo,jour)

    print *
    if (diff.gt.0) then
       print *,'******************************************************'
       print *,'***                                                ***'
       print *,'***            Votre licence est perimee           ***'
       print *,'***                                                ***'
       print *,'***   Contactez INRAE G2B      pour la renouveler  ***'
       print *,'***                                                ***'
       print *,'******************************************************'
       STOP 20

    else if (diff.ge.-60) then
       print *,'******************************************************'
       print *,'***                                                ***'
       print '(a36,i3,1x,a3,i5,a7)',' ***  Votre licence se terminera le ',jlim,mm(molim),anlim,'    ***'
       print *,'***                                                ***'
       print *,'***  Contactez INRAE G2B       pour la renouveler  ***'
       print *,'***                                                ***'
       print *,'******************************************************'
    end if

    print *
    print *
    return

    !========================================================================
  end subroutine vversion
  !========================================================================



  !========================================================================
  subroutine ecdate(diff,alim,mlim,jlim,annee,mois,jour)
    !========================================================================

    integer mc(12)
    integer jour,mois,annee,alim,mlim,jlim
    integer j1,m1,a1,j2,m2,a2
    integer i,diff
    logical test
    data mc/0,31,59,90,120,151,181,212,243,273,304,334/
    j1=jour ; m1=mois ; a1=annee
    j2=jlim ; m2=mlim ; a2=alim

    test=.false.
    if (a1.gt.a2) test=.true.
    if (a1.eq.a2.and.m1.gt.m2) test=.true.
    if (a1.eq.a2.and.m1.eq.m2.and.j1.gt.j2) test=.true.

    if (test) then
       i=j1 ; j1=j2; j2=i
       i=m1 ; m1=m2; m2=i
       i=a1 ; a1=a2; a2=i
    end if


    diff=365*a1 + mc(m1) + j1
    diff=diff - (365*a2 + mc(m2) + j2)

    do i=a2, a1
       if (mod(i,4).eq.0) then 
          diff=diff+1
          if (i.eq.a2.and.m2.gt.2) diff=diff-1
          if (i.eq.a1.and.m1.le.2) diff=diff-1
       end if
    end do

    if (test) diff=-diff

    return

    !========================================================================
  end subroutine ecdate
  !========================================================================




  !========================================================================
  !DECK HPSORT
  SUBROUTINE HPSORT (HX, N, STRBEG, STREND, IPERM, KFLAG, WORK, IER)
    !========================================================================
    !***BEGIN PROLOGUE  HPSORT
    !***PURPOSE  Return the permutation vector generated by sorting a
    !  substring within a character array and, optionally,
    !  rearrange the elements of the array.  The array may be
    !  sorted in forward or reverse lexicographical order.  A
    !  slightly modified quicksort algorithm is used.
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A1C, N6A2C
    !***TYPE CHARACTER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
    !***KEYWORDS  PASSIVE SORTING, SINGLETON QUICKSORT, SORT, STRING SORTING
    !***AUTHOR  Jones, R. E., (SNLA)
    !      Rhoads, G. S., (NBS)
    !      Sullivan, F. E., (NBS)
    !      Wisniewski, J. A., (SNLA)
    !***DESCRIPTION
    !
    !   HPSORT returns the permutation vector IPERM generated by sorting
    !   the substrings beginning with the character STRBEG and ending with
    !   the character STREND within the strings in array HX and, optionally,
    !   rearranges the strings in HX.   HX may be sorted in increasing or
    !   decreasing lexicographical order.  A slightly modified quicksort
    !   algorithm is used.
    !
    !   IPERM is such that HX(IPERM(I)) is the Ith value in the
    !   rearrangement of HX.  IPERM may be applied to another array by
    !   calling IPPERM, SPPERM, DPPERM or HPPERM.
    !
    !   An active sort of numerical data is expected to execute somewhat
    !   more quickly than a passive sort because there is no need to use
    !   indirect references. But for the character data in HPSORT, integers
    !   in the IPERM vector are manipulated rather than the strings in HX.
    !   Moving integers may be enough faster than moving character strings
    !   to more than offset the penalty of indirect referencing.
    !
    !   Description of Parameters
    ! HX - input/output -- array of type character to be sorted.
    !      For example, to sort a 80 element array of names,
    !      each of length 6, declare HX as character HX(100)*6.
    !      If ABS(KFLAG) = 2, then the values in HX will be
    !      rearranged on output; otherwise, they are unchanged.
    ! N  - input -- number of values in array HX to be sorted.
    ! STRBEG - input -- the index of the initial character in
    !     the string HX that is to be sorted.
    ! STREND - input -- the index of the final character in
    !     the string HX that is to be sorted.
    ! IPERM - output -- permutation array such that IPERM(I) is the
    !    index of the string in the original order of the
    !    HX array that is in the Ith location in the sorted
    !    order.
    ! KFLAG - input -- control parameter:
    !  =  2  means return the permutation vector resulting from
    !   sorting HX in lexicographical order and sort HX also.
    !  =  1  means return the permutation vector resulting from
    !   sorting HX in lexicographical order and do not sort
    !   HX.
    !  = -1  means return the permutation vector resulting from
    !   sorting HX in reverse lexicographical order and do
    !   not sort HX.
    !  = -2  means return the permutation vector resulting from
    !   sorting HX in reverse lexicographical order and sort
    !   HX also.
    ! WORK - character variable which must have a length specification
    !   at least as great as that of HX.
    ! IER - output -- error indicator:
    !     =  0  if no error,
    !     =  1  if N is zero or negative,
    !     =  2  if KFLAG is not 2, 1, -1, or -2,
    !     =  3  if work array is not long enough,
    !     =  4  if string beginning is beyond its end,
    !     =  5  if string beginning is out-of-range,
    !     =  6  if string end is out-of-range.
    !
    !     E X A M P L E  O F  U S E
    !
    ! CHARACTER*2 HX, W
    ! INTEGER STRBEG, STREND
    ! DIMENSION HX(10), IPERM(10)
    ! DATA (HX(I),I=1,10)/ '05','I ',' I','  ','Rs','9R','R9','89',
    !     1     ',*','N"'/
    ! DATA STRBEG, STREND / 1, 2 /
    ! CALL HPSORT (HX,10,STRBEG,STREND,IPERM,1,W)
    ! PRINT 100, (HX(IPERM(I)),I=1,10)
    ! 100 FORMAT (2X, A2)
    ! STOP
    ! END
    !
    !***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
    !       for sorting with minimal storage, Communications of
    !       the ACM, 12, 3 (1969), pp. 185-187.
    !***ROUTINES CALLED  XERMSG
    !***REVISION HISTORY  (YYMMDD)
    !   761101  DATE WRITTEN
    !   761118  Modified by John A. Wisniewski to use the Singleton
    !      quicksort algorithm.
    !   811001  Modified by Francis Sullivan for string data.
    !   850326  Documentation slightly modified by D. Kahaner.
    !   870423  Modified by Gregory S. Rhoads for passive sorting with the
    !      option for the rearrangement of the original data.
    !   890620  Algorithm for rearranging the data vector corrected by R.
    !      Boisvert.
    !   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
    !   920507  Modified by M. McClain to revise prologue text.
    !   920818  Declarations section rebuilt and code restructured to use
    !      IF-THEN-ELSE-ENDIF.  (SMR, WRB)
    !***END PROLOGUE  HPSORT
    !     .. Scalar Arguments ..
    INTEGER IER, KFLAG, N, STRBEG, STREND
    CHARACTER * (*) WORK
    !     .. Array Arguments ..
    INTEGER IPERM(*)
    CHARACTER * (*) HX(*)
    !     .. Local Scalars ..
    REAL R
    INTEGER I, IJ, INDX, INDX0, IR, ISTRT, J, K, KK, L, LM, LMT, M,NN, NN2
    !     .. Local Arrays ..
    INTEGER IL(30), IU(30)
    character*10 pgm
    character*80 msg
    !     .. External Subroutines ..
    !EXTERNAL XERMSG
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, INT, LEN
    !***FIRST EXECUTABLE STATEMENT  HPSORT
    pgm='HPSORT'
    IER = 0
    NN = N
    IF (NN .LT. 1) THEN
       IER = 1
       msg='The number of values to be sorted, N, is not positive.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF
    KK = ABS(KFLAG)
    IF (KK.NE.1 .AND. KK.NE.2) THEN
       IER = 2
       msg='The sort control parameter, KFLAG, is not 2, 1, -1, or -2.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF

    IF(LEN(WORK) .LT. LEN(HX(1))) THEN
       IER = 3
       msg='The length of the work variable, WORK, is too short.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF
    IF (STRBEG .GT. STREND) THEN
       IER = 4
       msg='The string beginning, STRBEG, is beyond its end, STREND.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF
    IF (STRBEG .LT. 1 .OR. STRBEG .GT. LEN(HX(1))) THEN
       IER = 5
       msg='The string beginning, STRBEG, is out-of-range.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF
    IF (STREND .LT. 1 .OR. STREND .GT. LEN(HX(1))) THEN
       IER = 6
       msg='The string end, STREND, is out-of-range.'
       CALL XERMSG (pgm,msg,IER, 1)
       RETURN
    ENDIF
    !
    !     Initialize permutation vector
    !


    DO 10 I=1,NN
       !print*,'I=',I
       IPERM(I) = I
10     CONTINUE
       !
       !     Return if only one value is to be sorted
       !
       IF (NN .EQ. 1) RETURN
       !
       !     Sort HX only
       !
       M = 1
       I = 1
       J = NN
       R = .375E0

20     IF (I .EQ. J) GO TO 70
       IF (R .LE. 0.5898437E0) THEN
          R = R+3.90625E-2
       ELSE
          R = R-0.21875E0
       ENDIF

30     K = I
       !
       !     Select a central element of the array and save it in location L
       !
       IJ = I + INT((J-I)*R)
       LM = IPERM(IJ)
       !
       !     If first element of array is greater than LM, interchange with LM
       !
       !print*,'IPERM(I)=',IPERM(I)
       !print*,'LM=',LM

       IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND)) THEN
          IPERM(IJ) = IPERM(I)
          IPERM(I) = LM
          LM = IPERM(IJ)
       ENDIF
       L = J
       !
       !     If last element of array is less than LM, interchange with LM
       !

       !print*,'IPERM(J)=',IPERM(J)
       !print*,'LM=',LM

       IF (HX(IPERM(J))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND)) THEN
          IPERM(IJ) = IPERM(J)
          IPERM(J) = LM
          LM = IPERM(IJ)
          !
          !   If first element of array is greater than LM, interchange
          !   with LM
          !
          IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND))THEN
             IPERM(IJ) = IPERM(I)
             IPERM(I) = LM
             LM = IPERM(IJ)
          ENDIF
       ENDIF
       GO TO 50
40     LMT = IPERM(L)
       IPERM(L) = IPERM(K)
       IPERM(K) = LMT
       !
       !     Find an element in the second half of the array which is smaller
       !     than LM
       !
50     L = L-1
       !PRINT*,'IPERM(L)=',IPERM(L)
       !print*,'LM=',LM
       IF (HX(IPERM(L))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND)) GO TO 50
       !
       !     Find an element in the first half of the array which is greater
       !     than LM
       !
60     K = K+1
       !print*,'IPERM(K)=',IPERM(K)
       !print*,'LM=',LM
       IF (HX(IPERM(K))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND)) GO TO 60
       !
       !
       !    Interchange these elements

       IF (K .LE. L) GO TO 40
       !
       !     Save upper and lower subscripts of the array yet to be sorted
       !
       IF (L-I .GT. J-K) THEN
          IL(M) = I
          IU(M) = L
          I = K
          M = M+1
       ELSE
          IL(M) = K
          IU(M) = J
          J = L
          M = M+1
       ENDIF
       GO TO 80
       !
       !     Begin again on another portion of the unsorted array
       !
70     M = M-1
       IF (M .EQ. 0) GO TO 110
       I = IL(M)
       J = IU(M)
       !
80     IF (J-I .GE. 1) GO TO 30
       IF (I .EQ. 1) GO TO 20
       I = I-1

90     I = I+1
       IF (I .EQ. J) GO TO 70
       LM = IPERM(I+1)
       !print*,'IPERM(I)=',IPERM(I)
       !print*,'LM=',LM
       IF (HX(IPERM(I))(STRBEG:STREND) .LE. HX(LM)(STRBEG:STREND)) GO TO 90
       K = I

100    IPERM(K+1) = IPERM(K)
       K = K-1

       !print*,'LM=',LM
       !print*,'IPERM(K)=',IPERM(K)

       IF (HX(LM)(STRBEG:STREND) .LT. HX(IPERM(K))(STRBEG:STREND)) GO TO 100
       IPERM(K+1) = LM
       GO TO 90

       !     Clean up
       !
110    IF (KFLAG .LE. -1) THEN
          !
          !   Alter array to get reverse order, if necessary
          !
          NN2 = NN/2
          DO 120 I=1,NN2
             IR = NN-I+1
             LM = IPERM(I)
             IPERM(I) = IPERM(IR)
             IPERM(IR) = LM
120          CONTINUE
          ENDIF
          !
          !     Rearrange the values of HX if desired
          !
          IF (KK .EQ. 2) THEN
             !
             !   Use the IPERM vector as a flag.
             !   If IPERM(I) < 0, then the I-th value is in correct location
             !
             DO 140 ISTRT=1,NN
                IF (IPERM(ISTRT) .GE. 0) THEN
                   INDX = ISTRT
                   INDX0 = INDX
                   !print*,'ISTRT=',ISTRT
                   WORK = HX(ISTRT)
                   !print*,'INDX=',INDX
                   !print*,'IPERM(INDX)=',IPERM(INDX)
130                IF (IPERM(INDX) .GT. 0) THEN
                      HX(INDX) = HX(IPERM(INDX))
                      INDX0 = INDX
                      IPERM(INDX) = -IPERM(INDX)
                      INDX = ABS(IPERM(INDX))
                      GO TO 130
                   ENDIF
                   !print*,'INDX0=',INDX0
                   HX(INDX0) = WORK
                ENDIF
140             CONTINUE
                !
                !   Revert the signs of the IPERM values
                !
                DO 150 I=1,NN
                   IPERM(I) = -IPERM(I)
150                CONTINUE

                ENDIF

                RETURN
              END subroutine HPSORT
              !========================================================================


              !========================================================================
              !*DECK HPPERM
              SUBROUTINE HPPERM (HX, N, IPERM, WORK, IER)
                !========================================================================

                INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
                CHARACTER*(*) HX(*), WORK
                character*10 pgm
                character*80 msg
                pgm='HPPERM'

                IER=0
                IF(N.LT.1)THEN
                   IER=164276
                   msg='The number of values to be rearranged, N, is not positive.'
                   CALL XERMSG (pgm,msg,IER, 1)
                   RETURN
                ENDIF
                IF(LEN(WORK).LT.LEN(HX(1)))THEN
                   IER=2
                   msg='The length of the work variable, WORK, is too short.'
                   CALL XERMSG (pgm,msg,IER, 1)
                   RETURN
                ENDIF

                !     CHECK WHETHER IPERM IS A VALID PERMUTATION

                DO 100 I=1,N
                   INDX=ABS(IPERM(I))
                   IF((INDX.GE.1).AND.(INDX.LE.N))THEN
                      IF(IPERM(INDX).GT.0)THEN
                         IPERM(INDX)=-IPERM(INDX)
                         GOTO 100
                      ENDIF
                   ENDIF
                   IER=3
                   msg='The permutation vector, IPERM, is not valid.'
                   CALL XERMSG (pgm,msg,IER, 1)
                   RETURN
100                CONTINUE

                   !     REARRANGE THE VALUES OF HX

                   !     USE THE IPERM VECTOR AS A FLAG.
                   !     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION

                   DO 330 ISTRT = 1 , N
                      IF (IPERM(ISTRT) .GT. 0) GOTO 330
                      INDX = ISTRT
                      INDX0 = INDX
                      WORK = HX(ISTRT)
320                   CONTINUE
                      IF (IPERM(INDX) .GE. 0) GOTO 325
                      HX(INDX) = HX(-IPERM(INDX))
                      INDX0 = INDX
                      IPERM(INDX) = -IPERM(INDX)
                      INDX = IPERM(INDX)
                      GOTO 320
325                   CONTINUE
                      HX(INDX0) = WORK
330                   CONTINUE

                      RETURN
                    END subroutine HPPERM
                    !========================================================================


                    !========================================================================
                    subroutine XERMSG (c2,c3,j,i)
                      !========================================================================
                      character*10 c2
                      character*80 c3
                      integer i,j
                      print *,'Message d erreur'
                      print *,'****************'
                      print *,'Subroutine : ',c2
                      print *,c3
                      if (i.eq.1) stop 20
                      return
                    end subroutine XERMSG
                    !========================================================================




                    !========================================================================
                    SUBROUTINE IPSORT (IX, N, IPERM, KFLAG, WORK, IER)
                      !========================================================================

                      !     .. Scalar Arguments ..
                      INTEGER IER, KFLAG, N
                      integer(kind=8) WORK
                      !     .. Array Arguments ..
                      INTEGER IPERM(*)
                      integer(kind=8) IX(*)
                      !     .. Local Scalars ..
                      REAL R
                      INTEGER I, IJ, INDX, INDX0, IR, ISTRT, J, K, KK, L, LM, LMT, M,NN, NN2
                      !     .. Local Arrays ..
                      INTEGER IL(30), IU(30)
                      character*10 pgm
                      character*80 msg
                      !     .. External Subroutines ..
                      !EXTERNAL XERMSG
                      !     .. Intrinsic Functions ..
                      INTRINSIC ABS, INT, LEN
                      !***FIRST EXECUTABLE STATEMENT  HPSORT
                      pgm='HPSORT'
                      IER = 0
                      NN = N
                      IF (NN .LT. 1) THEN
                         IER = 1
                         msg='The number of values to be sorted, N, is not positive.'
                         CALL XERMSG (pgm,msg,IER, 1)
                         RETURN
                      ENDIF
                      KK = ABS(KFLAG)
                      IF (KK.NE.1 .AND. KK.NE.2) THEN
                         IER = 2
                         msg='The sort control parameter, KFLAG, is not 2, 1, -1, or -2.'
                         CALL XERMSG (pgm,msg,IER, 1)
                         RETURN
                      ENDIF

                      !     Initialize permutation vector

                      DO 10 I=1,NN
                         IPERM(I) = I
10                       CONTINUE

                         !     Return if only one value is to be sorted

                         IF (NN .EQ. 1) RETURN

                         !     Sort HX only

                         M = 1
                         I = 1
                         J = NN
                         R = .375E0

20                       IF (I .EQ. J) GO TO 70
                         IF (R .LE. 0.5898437E0) THEN
                            R = R+3.90625E-2
                         ELSE
                            R = R-0.21875E0
                         ENDIF

30                       K = I

                         !     Select a central element of the array and save it in location L

                         IJ = I + INT((J-I)*R)
                         LM = IPERM(IJ)

                         !     If first element of array is greater than LM, interchange with LM

                         IF (IX(IPERM(I)) .GT. IX(LM)) THEN
                            IPERM(IJ) = IPERM(I)
                            IPERM(I) = LM
                            LM = IPERM(IJ)
                         ENDIF
                         L = J

                         !     If last element of array is less than LM, interchange with LM

                         IF (IX(IPERM(J)).LT. IX(LM)) THEN
                            IPERM(IJ) = IPERM(J)
                            IPERM(J) = LM
                            LM = IPERM(IJ)

                            !   If first element of array is greater than LM, interchange
                            !   with LM

                            IF (IX(IPERM(I)) .GT. IX(LM)) THEN
                               IPERM(IJ) = IPERM(I)
                               IPERM(I) = LM
                               LM = IPERM(IJ)
                            ENDIF
                         ENDIF
                         GO TO 50
40                       LMT = IPERM(L)
                         IPERM(L) = IPERM(K)
                         IPERM(K) = LMT

                         !     Find an element in the second half of the array which is smaller
                         !     than LM

50                       L = L-1
                         IF (IX(IPERM(L)).GT. IX(LM)) GO TO 50

                         !     Find an element in the first half of the array which is greater
                         !     than LM

60                       K = K+1
                         IF (IX(IPERM(K)).LT. IX(LM)) GO TO 60

                         !     Interchange these elements

                         IF (K .LE. L) GO TO 40

                         !     Save upper and lower subscripts of the array yet to be sorted

                         IF (L-I .GT. J-K) THEN
                            IL(M) = I
                            IU(M) = L
                            I = K
                            M = M+1
                         ELSE
                            IL(M) = K
                            IU(M) = J
                            J = L
                            M = M+1
                         ENDIF
                         GO TO 80

                         !     Begin again on another portion of the unsorted array

70                       M = M-1
                         IF (M .EQ. 0) GO TO 110
                         I = IL(M)
                         J = IU(M)

80                       IF (J-I .GE. 1) GO TO 30
                         IF (I .EQ. 1) GO TO 20
                         I = I-1

90                       I = I+1
                         IF (I .EQ. J) GO TO 70
                         LM = IPERM(I+1)
                         IF (IX(IPERM(I)) .LE. IX(LM)) GO TO 90
                         K = I

100                      IPERM(K+1) = IPERM(K)
                         K = K-1

                         IF (IX(LM).LT. IX(IPERM(K))) GO TO 100
                         IPERM(K+1) = LM
                         GO TO 90

                         !     Clean up

110                      IF (KFLAG .LE. -1) THEN

                            !   Alter array to get reverse order, if necessary

                            NN2 = NN/2
                            DO 120 I=1,NN2
                               IR = NN-I+1
                               LM = IPERM(I)
                               IPERM(I) = IPERM(IR)
                               IPERM(IR) = LM
120                            CONTINUE
                            ENDIF

                            !     Rearrange the values of HX if desired

                            IF (KK .EQ. 2) THEN

                               !   Use the IPERM vector as a flag.
                               !   If IPERM(I) < 0, then the I-th value is in correct location

                               DO 140 ISTRT=1,NN
                                  IF (IPERM(ISTRT) .GE. 0) THEN
                                     INDX = ISTRT
                                     INDX0 = INDX
                                     WORK = IX(ISTRT)
130                                  IF (IPERM(INDX) .GT. 0) THEN
                                        IX(INDX) = IX(IPERM(INDX))
                                        INDX0 = INDX
                                        IPERM(INDX) = -IPERM(INDX)
                                        INDX = ABS(IPERM(INDX))
                                        GO TO 130
                                     ENDIF
                                     IX(INDX0) = WORK
                                  ENDIF
140                               CONTINUE
                                  !
                                  !  Revert the signs of the IPERM values

                                  DO 150 I=1,NN
                                     IPERM(I) = -IPERM(I)
150                                  CONTINUE

                                  ENDIF

                                  RETURN
                                END subroutine IPSORT
                                !========================================================================



                                !========================================================================
                                !*DECK IPPERM
                                SUBROUTINE IPPERM (IX, N, IPERM, IER)
                                  !========================================================================
                                  INTEGER(kind=8) IX(*)
                                  INTEGER  N, IPERM(*), I, IER, INDX, INDX0, ITEMP, ISTRT
                                  character*10 pgm
                                  character*80 msg
                                  !***FIRST EXECUTABLE STATEMENT  IPPERM
                                  IER=0
                                  IF(N.LT.1)THEN
                                     IER=1
                                     msg='The number of values to be rearranged, N, is not positive.'
                                     CALL XERMSG (pgm,msg,IER, 1)
                                     RETURN
                                  ENDIF

                                  !     CHECK WHETHER IPERM IS A VALID PERMUTATION

                                  DO 100 I=1,N
                                     INDX=ABS(IPERM(I))
                                     IF((INDX.GE.1).AND.(INDX.LE.N))THEN
                                        IF(IPERM(INDX).GT.0)THEN
                                           IPERM(INDX)=-IPERM(INDX)
                                           GOTO 100
                                        ENDIF
                                     ENDIF
                                     IER=2
                                     msg='The permutation vector, IPERM, is not valid.'
                                     CALL XERMSG (pgm,msg,IER, 1)
                                     RETURN
100                                  CONTINUE

                                     !     REARRANGE THE VALUES OF IX

                                     !     USE THE IPERM VECTOR AS A FLAG.
                                     !     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION

                                     DO 330 ISTRT = 1 , N
                                        IF (IPERM(ISTRT) .GT. 0) GOTO 330
                                        INDX = ISTRT
                                        INDX0 = INDX
                                        ITEMP = IX(ISTRT)
320                                     CONTINUE
                                        IF (IPERM(INDX) .GE. 0) GOTO 325
                                        IX(INDX) = IX(-IPERM(INDX))
                                        INDX0 = INDX
                                        IPERM(INDX) = -IPERM(INDX)
                                        INDX = IPERM(INDX)
                                        GOTO 320
325                                     CONTINUE
                                        IX(INDX0) = ITEMP
330                                     CONTINUE

                                        RETURN
                                      END subroutine IPPERM
                                      !========================================================================




                                      ! ************** SUBROUTINE CHERCH DE D. BOICHARD ********************
                                      ! recherche l element iout egal a a dans la liste triee vecteur de longueur n
                                      ! iout=0 si pas trouve

                                      SUBROUTINE cherch(a,vecteur,nl,iout,lc)          
                                        IMPLICIT none 
                                        integer nl,iout,l,nd,jj,kk,lc       
                                        character(lc) vecteur(nl),a         

                                        !print *, "a : ", a

                                        iout=0           
                                        if (nl.eq.0) return
                                        if (nl.eq.1) then
                                           if (trim(a).eq.trim(vecteur(1))) iout=1
                                           return
                                        end if

                                        jj=1            
                                        nd=-1           
2                                       if (jj.gt.nl) GOTO 3           
                                        jj=jj*2              
                                        nd=nd+1              
                                        goto 2          
3                                       kk=jj/2              
                                        jj=kk/2              
                                        do l=1,nd   
                                           !print *, "kk, vecteur(kk) :", kk, vecteur(kk)	       
                                           if (trim(a).eq.trim(vecteur(kk))) goto 10           
                                           if (trim(a).lt.trim(vecteur(kk))) then          
                                              kk=kk-jj         
                                           else          
                                              kk=kk+jj         
                                              if (kk.gt.nl) kk=kk-jj          
                                           end if
                                           jj=jj/2            
                                        end do
                                        if (trim(a).gt.trim(vecteur(kk))) THEN       
                                           kk=kk+1            
                                        end if
10                                      continue

                                        if (trim(a).ne.trim(vecteur(kk))) return
                                        iout=kk

                                        RETURN          
                                      END subroutine cherch
                                      !*****************************************************************




                                      !========================================================================
                                    end program SNP_GWAS_publi
                                       !========================================================================





