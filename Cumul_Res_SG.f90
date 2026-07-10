!========================================================================
program cumul_res_SG
  !========================================================================
  
  ! Cumul_Res_SG
  ! Copyright (C) 2026 Thierry Tribout and Didier Boichard

  !  This program is free software: you can redistribute it and/or modify
  !  it under the terms of the GNU General Public License as published by
  !  the Free Software Foundation, either version 3 of the License, or
  !  (at your option) any later version.
  !
  !  This program is distributed in the hope that it will be useful,
  !  but WITHOUT ANY WARRANTY; without even the implied warranty of
  !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  !
  !  See the LICENSE file for more details.


! 21/11/2024 :  program reads GWAS result files  for each Group x Batch cell for 1 single Chromosom and gather them in 1 single Result file
  !               it also calculates _LOG10(P-value) for each estimated effect and each variant

  implicit none

  logical::testpresfic

  character(len=150) :: RepTrav0,RepTrav,FichInfos0,FichInfos,FicOut0,FicOut,methode0,methode,FicMapGwas0,FicMapGwas,FormatMapGwas0,FormatMapGwas
  character(len=20) :: modele0,modele ! indique si on a uniquement un effet additif ou un effet additif et un effet de dominance pour les VarGWAS
  character (250)    ::   parfile    !name of parameter file
  character*15 fx1, fx2, fx22
  character(320)::outformatA='', outformat='',outformatOPT='',outformatR=''

  integer::io,i,nblinFichInfos,chrom_anal, prem, dern, nlres, iores, iomapgwas, nlmapgwas,posgwasPrec,numgwasPrec

  real(kind=8)::moy_doses, std_doses, VarRes_1,VarRes_2,VarRes

  real(kind=8),allocatable::sol(:),test_1(:),pval_1(:),test_2(:),pval_2(:),test_F(:),pval_F(:),mlogP(:)

  integer::test_opt=0

  integer::posAdd=1,posDom

  integer,allocatable::mat_Infos(:,:)

  integer::numintra, numori

  integer,allocatable::MAPGWAS(:,:)
  integer::mapgwastemp(4)

  integer::nbeff  ! number of effects in the model for GWAS Variants: 1 if add , 2 if add_dom
  integer::TemDom,modDom=0

  character(len=20)::nomvar
  integer::pospb,poscm
  integer::iprec,nltot
  integer::firstVarGWAS,lastVarGWAS ! Id (numeric) of the first and last GWAS variants on the chromosom studied in Step_2 and Step_3 if the genotype format file is plink
  integer::nbGWAS
  

  character(len=128)::jour

  ! end of declarations 


  call fdate(jour)

  print*,'****************************************************'
  print*,'*                                                  *'
  print*,'*                Cumul_Res_SG V1.1                 *'
  print*,'*               Version 08 jul 2026                *'
  print*,'*                                                  *'
  print*,'   ',jour
  print*,'*                                                  *'
  print*,'****************************************************'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! reading parameter file

  PRINT*,'Reading parameter file'

  call getarg(1,parfile)
  parfile=adjustl(parfile)               ! ignore leading spaces
  write(*,'(5X,a)') trim(parfile)
  open(40,file=parfile)



  READ(40,*)
  READ(40,"(a)")RepTrav0
  RepTrav=trim(adjustl(RepTrav0))
  PRINT*,'Workig directory: ',RepTrav

  READ(40,*)
  READ(40,"(a)")FichInfos0
  FichInfos=trim(adjustl(FichInfos0))
  PRINT*,'File containing Group and Batch cells information for the chromosom : ',FichInfos

  READ(40,*)
  READ(40,"(a)")modele0
  modele=trim(adjustl(modele0))
  PRINT*,'Effects considered in the GWAS for the Variants (add or add_dom) : ',modele

  if((modele.ne.'add').and.(modele.ne.'add_dom')) then
     print*,'Incorrect Model for GWAS Variants: should be add or add_dom'
  endif

  if(modele.eq.'add') nbeff=1
  if(modele.eq.'add_dom') nbeff=2
  print*,' '
  print*,'Number of effect(s) considered for Variants in GWAS model =',nbeff
  print*,' '


  if(modele.eq.'add') modDom=0
  if(modele.eq.'add_dom') modDom=1
  print*,'modDom=',modDom
  print*,' '  

  READ(40,*)
  READ(40,"(a)")methode0
  methode=trim(adjustl(methode0))
  PRINT*,'Strategy for computing the Residual Variance (approx or optim or exact) : ',methode

  if((methode.ne.'approx').and.(methode.ne.'optim').and.(methode.ne.'exact')) then
     print*,'Incorrect value for Residual Variance computing strategy: must be among approx optim exact'
     STOP 20
  endif

  READ(40,*)
  READ(40,"(a)")FicMapGwas0
  FicMapGwas=trim(adjustl(FicMapGwas0))
  PRINT*,'MAP file for the GWAS Variants that was used during GWAS: ',FicMapGwas

  READ(40,*)
  READ(40,"(a)")FormatMapGwas0
  FormatMapGwas=trim(adjustl(FormatMapGwas0))
  PRINT*,'MAP file format for the GWAS Variants: ',FormatMapGwas

  if((FormatMapGwas.ne.'typ_eval').and.(FormatMapGwas.ne.'plink')) then
     print*,'Incorrect MAP file format for the GWAS Variants: must be among typ_eval plink'
     STOP 20
  endif
     

  READ(40,*)
  READ(40,"(a)")FicOut0
  FicOut=trim(adjustl(FicOut0))
  PRINT*,'Name for final result file in Working Directory, gathering all Greoup x Batch results for the considered chromosom: ',FicOut

  ! Checking if File with Group x Batch information and GWAS Variant MAP file exist

   INQUIRE (FILE=trim(RepTrav)//trim("/")//trim(FichInfos),EXIST=testpresfic)
   if(testpresfic.eq..false.) then
      print*,'File with Group x Batch information does not exist'
      STOP 30
   endif

   INQUIRE (FILE=trim(FicMapGwas),EXIST=testpresfic)
   if(testpresfic.eq..false.) then
      print*,'Variant MAP file does not exist'
      STOP 30
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! opening input and output files
   open(unit=10,file=trim(RepTrav)//trim("/")//trim(FichInfos),status='old')

   ! Result file for ADDITIVE effects
   if((modele.eq.'add').or.(modele.eq.'add_dom')) then
      open(unit=13,file=trim(RepTrav)//trim("/ADD_")//trim(FicOut))
      write(13,*) 'Chrom    Pos     NumOriVar numGR  numSG   moy(doses)       std(doses)         sol_add           test_add          p-val_add      -LOG(P-val_add)    VarRES'
   endif

   ! Result file for dominance effects
   if(modele.eq.'add_dom') then
      open(unit=14,file=trim(RepTrav)//trim("/DOM_")//trim(FicOut))
      write(14,*) 'Chrom    Pos     NumOriVar numGR  numSG   moy(doses)       std(doses)         sol_dom           test_dom          p-val_dom      -LOG(P-val_dom)    varRES'
   endif



   ! on lit le fichier FichInfos une premiere fois pour connaitre le nb de fichiers resultats elementaires qu on va devoir lire

   ! Checking file exists and is not empty
   read(10,*,iostat=io)
   if (io.lt.0) then
      print*,'Group x Batch cells information file is empty'
      stop 30
   endif
   rewind(10)

   ! first reading to count nb of lines in the file
   i=0
   do
      read(10,*,iostat=io)
      if (io/=0) exit
      i=i+1
   enddo
   nblinFichInfos = i
   rewind(10)

   allocate(mat_Infos(3,nblinFichInfos))
   mat_Infos=0

   ! second reading to fill table mat_Infos
   do i=1,nblinFichInfos
      read(10,*) chrom_anal,mat_Infos(1:2,i),prem,dern,mat_Infos(3,i)
   enddo

   close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ! opening GWAS Variants MAP file, storing in memory position of each Variant
   open(20,file=FicMapGwas,status='old')

   if(FormatMapGwas.eq.'typ_eval') then

      iomapgwas=0
      nlmapgwas=0
      do
         read(20,*,iostat=iomapgwas)
         if (iomapgwas .ne. 0) exit
         nlmapgwas=nlmapgwas+1
      enddo

      rewind(20)

      print*,'Number of lines in GWAS Variant MAP file: ',nlmapgwas

      allocate(MAPGWAS(3,nlmapGWAS))  ! 1st line = num_chromos / 2nd line = position on chromosom / 3rd line = inclusion/exclusion of Variant
      MAPGWAS=0
      MAPGWAStemp=0
      posgwasPrec=0
      numgwasPrec=0

      ! 2nd reading to fill table MAPGWAS
      nlmapgwas=0
      do
         read(20,*,iostat=iomapgwas) MAPGWAStemp(:)
         if (iomapgwas .ne. 0) exit
         nlmapgwas=nlmapgwas+1
         if((MAPGWAStemp(2).le.posgwasPrec).or.(MAPGWAStemp(3).le.numgwasPrec).or.(MAPGWAStemp(3).ne.nlmapgwas)) then
            print*,'GWAS Variant MAP file is not sorted by ascending position / Variant order'
            STOP 30
         endif
         if((MAPGWAStemp(4).ne.0).and.(MAPGWAStemp(4).ne.1)) then
            print*,'INCL/EXCL indicator in MAP file is different from 0 and 1 for Variant: ',nlmapgwas
            STOP 30
         endif
         MAPGWAS(1:2,nlmapgwas)=MAPGWAStemp(1:2)
         MAPGWAS(3,nlmapgwas)=MAPGWAStemp(4)
         posgwasPrec=MAPGWAStemp(2)
         numgwasPrec=MAPGWAStemp(3)
      enddo

   endif

   ! if FORMAT for MAP file=plink : the file contains positions of all variants on all chromosomes
   ! we must know the number of the chromosom carrying the Variants BEFORE reading MAP file to keep only Information of the Variants on the considered Chromosom

   if(FormatMapGwas.eq.'plink') then

        ! on lit une 1ere fois le fichier CARTE_GWAS plink pour connaitre le nb de variants pour GWAS sur le chromosome considere pour dimensionner vecteurs et tableaux
        ! a priori pas d indicateur inclusion/exclusion dans fichier carte variants GWAS plink pour l instant

        iomapgwas=0
        nlmapgwas=0
        do
           read(20,*,iostat=iomapgwas) i,nomvar,poscM,pospb
           if (iomapgwas .ne. 0) exit
           if(i.eq.chrom_anal) nlmapgwas=nlmapgwas+1
        enddo

        rewind(20)

        print*,'Number of lines read in GWAS Variant MAP file for chromosom ',chrom_anal,' is: ',nlmapgwas

        allocate(MAPGWAS(3,nlmapGWAS))  ! 1st line = num_chromos / 2nd line = position on chromosom / 3rd line = inclusion/exclusion of Variant
        MAPGWAS=0
        posgwasPrec=0
        nbGWAS=0
        i=0
        iprec=0
        pospb=0
        firstVarGWAS=0
        lastVarGWAS=0
        nltot=0
        
        ! 2nd reading of GWAS Variant MAP file to fill table MAPGWAS
        nlmapgwas=0
        do
           read(20,*,iostat=iomapgwas) i,nomvar,poscM,pospb
           nltot=nltot+1
           if (iomapgwas .ne. 0) exit
           if(i.eq.chrom_anal) then ! GWAS variant is on the chromosom studied in Step_2 and Step_3
              nlmapgwas=nlmapgwas+1
              if(firstVarGWAS.eq.0) firstVarGWAS=nltot ! overall order of first GWAS Variant on the chromosom considered
              lastVarGWAS=nltot                        ! overall order of last Variant on the chromosom considered 
              if(pospb.le.posgwasPrec) then
                 print*,'GWAS Variant MAP file for Chromosom ',chrom_anal,' is not sorted by ascending position on chromosom'
                 STOP 30
              endif
              MAPGWAS(1,nlmapgwas)=chrom_anal
              MAPGWAS(2,nlmapgwas)=pospb
              MAPGWAS(3,nlmapgwas)=1
              posgwasPrec=pospb
           endif
           iprec=i
           if(i.lt.iprec) then
              print*,'GWAS Variant MAP file for Chromosom ',chrom_anal,' is not sorted by ascending chromosom'
              STOP 30
           endif
        enddo
        nbGWAS=nlmapgwas

   endif


   close(20)


  ! writing format for final result file
  !write(outformat,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,d17.9,1x,d17.9,1x,d17.9,1x,d17.9,1x,d17.9,1x,d17.9)'
  write(outformatA,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3)'
  write(outformat,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3)'
  write(outformatOPT,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3)'

  write(outformatR,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,i7,1x,i7,1x,i7,1x,i7,1x,i7,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3,1x,E24.12E3)'

  !write(outformat,'(a,i0,a)') '(i3,1x,i9,1x,i9,1x,i5,1x,i6,1x,*,1x,*,1x,*,1x,*,1x,*,1x,*)'

  ! on alloue les vecteurs pour lire les solutions, tests, pvalues en fonction du nombre d effets dans la partie GWAS du modele en partie 3
  allocate(sol(nbeff),test_1(nbeff),pval_1(nbeff),test_2(nbeff),pval_2(nbeff),test_F(nbeff),pval_F(nbeff),mlogP(nbeff))
  sol=0.0d0
  test_1=0.0d0
  pval_1=0.0d0
  test_2=0.0d0
  pval_2=0.0d0
  test_F=0.0d0
  pval_F=0.0d0
  mlogP=0.0d0


  ! Reading of each result file 1 by 1 and adding of necessary information in the final file

  do i=1,nblinFichInfos

     write(fx1,'(i0)') chrom_anal
     write(fx2,'(i0)') mat_Infos(1,i)
     write(fx22,'(i0)') mat_Infos(2,i)

     INQUIRE (FILE=trim(RepTrav)//trim("/")//trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Resultats_chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim(".txt"),EXIST=testpresfic)
     if(testpresfic.eq..false.) then
        print*,'Result file for  chrom ',trim(fx1),' Group ',trim(fx2),' Batch ',trim(fx22),'is missing'
        STOP 30
     endif

     open(unit=11,file=trim(RepTrav)//trim("/")//trim("chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim("/SG")//trim(fx22)//trim("/Resultats_chr")//trim(fx1)//trim("_gr")//trim(fx2)//trim(".txt"),status='old')

     nlres=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  DEBUT NOUVELLE VERSION DU CODE  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if((modele.eq.'add').or.(modele.eq.'add_dom')) then

        posAdd=1

        do

           if(methode.eq.'optim')then

              posDom=2
              sol=0.0d0
              test_1=0.0d0
              pval_1=0.0d0
              test_2=0.0d0
              pval_2=0.0d0
              VarRes_1=0.0d0
              VarRes_2=0.0d0
              VarRes=0.0d0
              mlogP=0.0d0
              test_F=0.0d0
              pval_F=0.0d0

              read(11,*,iostat=iores) numintra, numori, moy_doses, std_doses, TemDOM, sol(:), test_1(:), pval_1(:), test_2(:), pval_2(:), VarRes_1, VarRes_2
              if (iores .ne. 0) exit
              nlres=nlres+1
              if((TemDom+1).gt.nbeff) then
                 print*,' '
                 print*,'ERROR for GWAS Variant ',numori,' : more effects found in result file than in model: TemDom=',TemDom,' nbeff=',nbeff
                 STOP 30
              endif

              if(nbeff.eq.2)then
                 if(temDom.eq.1) posDom=2
              endif

              ! solution for ADDITIVE EFFECT of GWAS Variant
              if(test_2(posAdd).ne.0.0d0) then     ! exact test statistic for ADDITIVE effect has been calculated for the Variant 
                 test_F(posAdd)=test_2(posAdd)
                 pval_F(posAdd)=pval_2(posAdd)
                 VarRes=VarRes_2

                 !if((nbeff.eq.2).or.(nbeff.eq.2)) then
                 if(temDom.eq.1) then
                    test_F(posDom)=test_2(posDom)
                    pval_F(posDom)=pval_2(posDom)
                 endif
              else
                 test_F(posAdd)=test_1(posAdd)        ! exact test statistics for ADDITIVE effect have NOT been calculated for the Variant --> use approx test statistics
                 pval_F(posAdd)=pval_1(posAdd)
                 VarRes=VarRes_1

                 if(temDom.eq.1) then
                    test_F(posDom)=test_1(posDom)
                    pval_F(posDom)=pval_1(posDom)
                 endif
              endif

           endif

           if((methode.eq.'approx').or.(methode.eq.'exact'))then

              posDom=0
              sol=0.0d0
              test_1=0.0d0
              pval_1=0.0d0
              test_2=0.0d0
              pval_2=0.0d0
              VarRes_1=0.0d0
              VarRes_2=0.0d0
              VarRes=0.0d0
              mlogP=0.0d0
              test_F=0.0d0
              pval_F=0.0d0

              read(11,*,iostat=iores) numintra, numori, moy_doses, std_doses, TemDOM, sol(:), test_1(:), pval_1(:), VarRes_1
              if (iores .ne. 0) exit
              nlres=nlres+1
              if((TemDom+1).gt.nbeff) then
                 print*,' '
                 print*,'ERROR for GWAS Variant ',numori,' : more effects found in result file than in model: TemDom=',TemDom,' nbeff=',nbeff
                 STOP 30
              endif

              if(nbeff.eq.2)then
                 if(temDom.eq.1) posDom=2
              endif

              test_F(posAdd)=test_1(posAdd)          ! exact test statistic for effects has NOT been calculated for the Variant --> use approx test statistic
              pval_F(posAdd)=pval_1(posAdd)
              VarRes=VarRes_1

              if(temDom.eq.1) then
                 test_F(posDom)=test_1(posDom)
                 pval_F(posDom)=pval_1(posDom)
              endif

           endif


           if(pval_F(posAdd).ne.0.0d0)  mlogP(posAdd)=-1.0d0*LOG10(pval_F(posAdd))

           if(temDom.eq.1) then
              if(pval_F(posDom).ne.0.0d0) mlogP(posDom)=-1.0d0*LOG10(pval_F(posDom))
           endif

           ! Writing in result file for ADDITIVE EFFECTS
           write(13,outformatA) chrom_anal,MAPGWAS(2,numori),numori,mat_Infos(1,i),mat_Infos(2,i), moy_doses, std_doses, sol(posAdd), test_F(posAdd), pval_F(posAdd), mlogP(posAdd), VarRes

           ! Writing in result file for DOMINANCE EFFECTS
           if(modDom.eq.1) then
              write(14,outformatA) chrom_anal,MAPGWAS(2,numori),numori,mat_Infos(1,i),mat_Infos(2,i), moy_doses, std_doses, sol(posDom), test_F(posDom), pval_F(posDom), mlogP(posDom), VarRes
           endif

        enddo

     endif  ! end of test if model = add or add_dom



     print*,' '
     print*,'Number of lines read in Result file for chr ',trim(fx1),' Group ',trim(fx2),' Batch ',trim(fx22),' is :',nlres
     print*,' '

     if(nlres.ne.mat_Infos(3,i)) then
        print*,'Number of lines read in Result file for chr ',trim(fx1),' Group ',trim(fx2),' Batch ',trim(fx22),' is ',nlres,' whereas expecting ',mat_Infos(3,i),' lines'
        !STOP 33
     endif

     close(unit=11)

  enddo ! end of do loop i=1,nblinFichInfos

  close(unit=13)
  if(modDom.eq.1) close(unit=14)

  if((modele.eq.'add_dom').and.(test_opt.eq.1)) close(unit=21)

  print*,' '
  print*,'********************'
  print*,' '
  print*,'Program completed'

  STOP 0

  ! Program completed


  !========================================================================
end program cumul_res_SG
!========================================================================





