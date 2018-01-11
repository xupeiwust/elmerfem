!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) asny later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Denis Cohen, Juha Hartikainen
! *  Email:  thomas Zwinger [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date:  January 2017                
! * 
! *****************************************************************************
!>  Module containing solver for enhanced permafrost problem and material
!>  settings
MODULE PermafrostMaterials
  
  USE Types
  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE
  !---------------------------------
  ! type for solvent (water and ice)
  !---------------------------------
  TYPE SolventMaterial_t
     REAL(KIND=dp) :: &
          Mw,rhow0,rhoi0,hw0,hi0,vi0,bccw0,&
          kw0th,ki0th,bw,bi, &
          cw0,acw(0:5),bcw(0:5), &
          ci0,aci(0:5),&
          aw0,kw0,zw0,aaw(0:5),bzw(0:5),ckw(0:5), &
          ai0,ki0,aai(0:5),cki(0:5),&
          muw0,nu10,anw(0:5),bnw(0:5)
     INTEGER :: &
          acwl,bcwl,aawl,bzwl,ckwl,&
          acil,aail,ckil,anwl,bnwl 
  END type SolventMaterial_t
  
  !---------------------------------
  ! type for solute (ions)
  !---------------------------------
  TYPE SoluteMaterial_t
     REAL(KIND=dp) ::  Mc,vc0,kc0th,&
          d1,d2,dc0,dc1,bc,&
          cc0,acc(0:5),bcc(0:5),&
          rhoc0,ac0,kc0,zc0,aac(0:5),ckc(0:5),bzc(0:5),&
          nu20,anc(0:5),bnc(0:5)
     INTEGER :: accl, bccl,aacl,ckcl,bzcl,ancl,bncl
     CHARACTER(LEN=MAX_NAME_LEN) :: SoluteName
  END type SoluteMaterial_t
  
  !---------------------------------
  ! type for rock material
  !---------------------------------
  TYPE RockMaterial_t
     INTEGER :: NumerOfRockRecords
     REAL(KIND=dp), ALLOCATABLE :: ks0th(:),e1(:),bs(:),rhos0(:),&
          Xi0(:),eta0(:),hs0(:),Kgwh0(:,:,:),qexp(:),alphaL(:),alphaT(:),RadGen(:),&
          cs0(:),acs(:,:),as0(:),aas(:,:),ks0(:),cks(:,:)
     INTEGER, ALLOCATABLE :: acsl(:),aasl(:),cksl(:)
     CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  END TYPE RockMaterial_t
  
CONTAINS
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------

  !-------------------------------------------------
  ! I/O related functions
  !-------------------------------------------------

  SUBROUTINE SetPermafrostSolventMaterial( CurrentSolventMaterial)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params, Constants
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    ! ----------- local
    TYPE(SolventMaterial_t), TARGET :: LocalSolventMaterial
    LOGICAL :: FirstTime=.TRUE.
    CHARACTER(LEN=MAX_NAME_LEN) :: SubroutineName='SetPermafrostSolventMaterial'
    SAVE LocalSolventMaterial, FirstTime

    IF (FirstTime) THEN
      !------------------------------------------------------------------------------
      ! set constants for water and ice 
      ! Mw,rhow0,rhoi0,hw0,hi0,vi0,cw0,ci0,acw(3),bcw(0:2),aci(0:1),kw0th,ki0th, bi, bw
      !------------------------------------------------------------------------------
      LocalSolventMaterial % Mw =    0.018_dp 
      LocalSolventMaterial % rhow0 = 999.9_dp 
      LocalSolventMaterial % rhoi0 = 916.8_dp 
      LocalSolventMaterial % hw0 =   0.0_dp 
      LocalSolventMaterial % hi0 =  -333360.0_dp 
      LocalSolventMaterial % vi0 =   0.0010908_dp
      LocalSolventMaterial % kw0th = 0.56_dp 
      LocalSolventMaterial % ki0th = 2.24_dp
      LocalSolventMaterial % muw0=    1.792d-03
      LocalSolventMaterial % bi = 0.0_dp ! ?????????????
      LocalSolventMaterial % bw = 0.0_dp ! ?????????????
      ! --------------------- polynomials
      LocalSolventMaterial % cw0 = 4207.7_dp
      LocalSolventMaterial % kw0 = 4.4534d-10
      LocalSolventMaterial % aw0 = 5.3358d-05
      LocalSolventMaterial % zw0 = -6.6237d-01
      LocalSolventMaterial % ci0 = 2088.8_dp
      LocalSolventMaterial % ki0 = 1.1417d-10
      LocalSolventMaterial % ai0 = 1.6781d-04
      LocalSolventMaterial % nu10 = -0.031950_dp
      LocalSolventMaterial % acwl=2
      LocalSolventMaterial % bcwl=2
      LocalSolventMaterial % ckwl=0
      LocalSolventMaterial % aawl=5
      LocalSolventMaterial % bzwl=1
      LocalSolventMaterial % acil=1
      LocalSolventMaterial % ckil=0
      LocalSolventMaterial % aail=1
      LocalSolventMaterial % anwl=3
      LocalSolventMaterial % bnwl=0
      LocalSolventMaterial % acw(0:5) = &
           RESHAPE([1.0,-0.088,0.2859,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % acw))
      LocalSolventMaterial % bcw(0:5) = &
           RESHAPE([1.0,3.0201,14.7607,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % bcw))
      LocalSolventMaterial % ckw(0:5) = &
           RESHAPE([1.0,0.0,0.0,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % ckw))
      LocalSolventMaterial % aaw(0:5) = &
           RESHAPE([1.0,-79.1305,207.4836,-403.8270,395.5347,-166.1466], SHAPE(LocalSolventMaterial % aaw))
      LocalSolventMaterial % bzw(0:5) = &
           RESHAPE([1.0,5.1806,0.0,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % bzw))
      LocalSolventMaterial % aci(0:5) = &
           RESHAPE([1.0,0.9557,0.0,0.0,0.0,0.0],SHAPE(LocalSolventMaterial % aci))
      LocalSolventMaterial % cki(0:5) = &
           RESHAPE([1.0,0.0,0.0,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % cki))
      LocalSolventMaterial % aai(0:5) = &
           RESHAPE([1.0,1.1923,0.0,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % aai))
      LocalSolventMaterial % anw(0:5) = &
           RESHAPE([1.0,-1.5771,1.3920,-0.4784,0.0,0.0], SHAPE(LocalSolventMaterial % anw))
      LocalSolventMaterial % bnw(0:5) = &
           RESHAPE([1.0,0.0,0.0,0.0,0.0,0.0], SHAPE(LocalSolventMaterial % anw))
      CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
      CALL INFO(SubroutineName,"Solvent related constants",Level=9)
      WRITE(Message,*) "Mw",LocalSolventMaterial % Mw,"rhow0",LocalSolventMaterial % rhow0,"rhoi0",LocalSolventMaterial % rhoi0,&
           "hw0",LocalSolventMaterial % hw0,"hi0",LocalSolventMaterial % hi0,"vi0",LocalSolventMaterial % vi0,&
           "cw0",LocalSolventMaterial % cw0,"ci0",LocalSolventMaterial % ci0,"acw(3)",LocalSolventMaterial % acw(0:2)
      CALL INFO(SubroutineName,Message,Level=9)
      WRITE(Message,*) "bcw(0:2)",LocalSolventMaterial % bcw(0:2),"aci(0:1)",LocalSolventMaterial % aci(0:1),&
           "kw0th",LocalSolventMaterial % kw0th,"ki0th",LocalSolventMaterial % ki0th," bi",LocalSolventMaterial % bi,&
           "bw",LocalSolventMaterial % bw
      CALL INFO(SubroutineName,Message,Level=9)
      CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
      FirstTime = .FALSE.
    END IF
    CurrentSolventMaterial => LocalSolventMaterial
  END SUBROUTINE SetPermafrostSolventMaterial
  
  SUBROUTINE ReadPermafrostSoluteMaterial( Params,Constants,CurrentSoluteMaterial )
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params, Constants
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    ! ----------- local
    TYPE(SoluteMaterial_t), TARGET :: LocalSoluteMaterial
    INTEGER :: i,j,k,l, n,t, active, DIM, ok,InitialNumerOfSoluteRecords, EntryNumber
    INTEGER,parameter :: io=20
    LOGICAL :: Found, DataRead=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  SoluteFileName, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SubroutineName='ReadPermafrostSoluteMaterial'

    SAVE DataRead,SoluteFileName, LocalSoluteMaterial

    IF (DataRead) THEN
      CurrentSoluteMaterial => LocalSoluteMaterial   
      RETURN
    ELSE
      DIM = CoordinateSystemDimension()
      !------------------------------------------------------------------------------
      ! Inquire and open file
      !------------------------------------------------------------------------------
      ! give preference to a defined material database
      SoluteFileName = GetString( Params, 'Solute Material File', Found )
      IF (.NOT.Found) THEN
        CALL INFO(SubroutineName," 'Solute Material File' keyword not found - assigning default values!",Level=1)
        DataRead=.TRUE.
        LocalSoluteMaterial % SoluteName = TRIM('Sea Salt')
        LocalSoluteMaterial % Mc=        0.0323_dp
        LocalSoluteMaterial % vc0=       0.00057372_dp
        LocalSoluteMaterial % kc0th=     0.56_dp      
        LocalSoluteMaterial % d1=       0.87_dp		       
        LocalSoluteMaterial % d2=       1.30_dp		       
        LocalSoluteMaterial % dc0=       0.81_dp
        LocalSoluteMaterial % dc1=       2.60_dp
        LocalSoluteMaterial % bc=        0.0_dp
        ! heat capacity
        LocalSoluteMaterial % cc0=       2413.8_dp
        LocalSoluteMaterial % acc(0:5) = &
             RESHAPE([1.0_dp,-0.0887_dp,0.2859_dp,0.0_dp,0.0_dp,0.0_dp], SHAPE(LocalSoluteMaterial % acc))
        LocalSoluteMaterial % bcc(0:5) = &
             RESHAPE([1.0_dp,-3.0201_dp,14.7607_dp,0.0_dp,0.0_dp,0.0_dp], SHAPE(LocalSoluteMaterial % bcc))
        LocalSoluteMaterial % accl = 2
        LocalSoluteMaterial % bccl = 2
        ! density
        LocalSoluteMaterial % rhoc0 = 1743.0_dp
        LocalSoluteMaterial % ac0 =   5.3358d-05 
        LocalSoluteMaterial % kc0 =   4.4534d-10 
        LocalSoluteMaterial % zc0 =   0.66237_dp
        LocalSoluteMaterial % aac = &
             RESHAPE([1.0_dp, -0.0877_dp, 2.2859_dp, 0.0_dp, 0.0_dp, 0.0_dp], SHAPE(LocalSoluteMaterial % aac))
        LocalSoluteMaterial % aacl = 2
        LocalSoluteMaterial % ckc = &
             RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], SHAPE(LocalSoluteMaterial % ckc))
        LocalSoluteMaterial % ckcl = 0 
        LocalSoluteMaterial % bzc = &
             RESHAPE([1.0_dp, -5.1806_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], SHAPE(LocalSoluteMaterial % bzc))
        LocalSoluteMaterial % bzcl = 1
        LocalSoluteMaterial % nu20 = 1.9175_dp
        LocalSoluteMaterial % anc = &
             RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], SHAPE(LocalSoluteMaterial % anc))
        LocalSoluteMaterial % bnc= &
             RESHAPE([1.0_dp, 5.4188_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], SHAPE(LocalSoluteMaterial % bnc))
        LocalSoluteMaterial % ancl = 0
        LocalSoluteMaterial % bncl = 1
        DataRead=.TRUE.
      ELSE      
        OPEN(unit = io, file = TRIM(SoluteFileName), status = 'old',iostat = ok)
        IF (ok /= 0) THEN
          WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(SoluteFileName)
          CALL FATAL(Trim(SubroutineName),Trim(message))
        ELSE
          !------------------------------------------------------------------------------
          ! Read in the number of records in file (first line integer)
          !------------------------------------------------------------------------------
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % SoluteName
          WRITE(Message,'(A,A)') 'Reading entry', TRIM(LocalSoluteMaterial % SoluteName)
          CALL INFO(Trim(SubroutineName),Trim(Message),Level=3)
          !------------------------------------------------------------------------------
          ! Read in information for (currently fixed) solute (= salts)
          ! 
          !     Mc, rhoc0
          !     vc0,cc0,acc(0:2),bcc(0:2)
          !     kc0th,muw0,d1,d2,dc0,dc1,bc
          !------------------------------------------------------------------------------
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % Mc, Comment 
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % vc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % kc0th, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % d1, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % d2, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % dc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % dc1, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bc, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % cc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % acc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bcc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % accl, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bccl, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % rhoc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ac0, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % kc0, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % zc0 , Comment 
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % aac(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % aacl, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ckc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ckcl, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bzc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bzcl , Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % nu20, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % anc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ancl , Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bnc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bncl, Comment
          DataRead=.TRUE.
10        CLOSE(io)
          IF (.NOT.DataRead) THEN
            WRITE(Message,'(A,A,A)')  'Not all entries in "Solute material File" ',TRIM(SoluteFileName),' found.'
            CALL FATAL(Trim(SubroutineName),Trim(message))
          END IF
        END IF
      END IF
      CurrentSoluteMaterial => LocalSoluteMaterial
    END IF
    CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(SubroutineName,"Solute related constants",Level=9)
    WRITE(Message,*) "Mc",CurrentSoluteMaterial % Mc,"vc0",CurrentSoluteMaterial % vc0,&
         "kc0th", CurrentSoluteMaterial %kc0th
    CALL INFO(SubroutineName,Message,Level=9)    
    WRITE(Message,*) "d1",CurrentSoluteMaterial % d1,"d2",&
         CurrentSoluteMaterial % d2,"dc0",CurrentSoluteMaterial % dc0,"dc1",&
         CurrentSoluteMaterial % dc1,"bc",CurrentSoluteMaterial % bc
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*) "cc0",CurrentSoluteMaterial % cc0,"acc(0:5)",CurrentSoluteMaterial % acc(0:5),&
         "bcc(0:5)",CurrentSoluteMaterial % bcc(0:5)
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*) "rhoc0",CurrentSoluteMaterial % rhoc0,"ac0",CurrentSoluteMaterial % ac0,"kc0",&
         CurrentSoluteMaterial % kc0,"zc0",CurrentSoluteMaterial % zc0
    WRITE(Message,*)  "aac(0:5)",CurrentSoluteMaterial % aac(0:5),"ckc(0:5)",CurrentSoluteMaterial % ckc(0:5),&
         "bzc(0:5)",CurrentSoluteMaterial % bzc(0:5)
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "bzc(0:5)",CurrentSoluteMaterial % bzc(0:5),"bnc(0:5)",CurrentSoluteMaterial % bnc(0:5),&
         "nu20",CurrentSoluteMaterial % nu20    
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "aacl",LocalSoluteMaterial % aacl,"ckcl", LocalSoluteMaterial % ckcl    
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "bzcl",LocalSoluteMaterial % bzcl,"ancl", LocalSoluteMaterial % ancl,&
         "bncl",LocalSoluteMaterial % bncl
    CALL INFO(SubroutineName,Message,Level=9)
    CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
    RETURN
20  WRITE(Message,'(A,A,A)')  'Not all entries in "Solute material File" ',TRIM(SoluteFileName),' found.'
    CLOSE(io)
    CALL FATAL(Trim(SubroutineName),Trim(message))
  END SUBROUTINE ReadPermafrostSoluteMaterial

  !---------------------------------------------------------------------------------------------
  FUNCTION ReadPermafrostRockMaterial( Params,Constants,CurrentRockMaterial ) RESULT(NumerOfRockRecords)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params, Constants
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(RockMaterial_t), TARGET :: LocalRockMaterial
    Integer :: NumerOfRockRecords

    INTEGER :: i,j,k,l, n,t, active, DIM, ok,InitialNumerOfRockRecords, EntryNumber
    INTEGER,parameter :: io=21
    LOGICAL :: Found, fexist, FirstTime=.TRUE., AllocationsDone=.FALSE., DataRead=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  MaterialFileName, NewMaterialFileName, str, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='ReadPermafrostRockMaterial'

    SAVE AllocationsDone,DataRead,InitialNumerOfRockRecords,LocalRockMaterial,MaterialFileName

    IF (DataRead) THEN
      NumerOfRockRecords = InitialNumerOfRockRecords
      CurrentRockMaterial => LocalRockMaterial
      RETURN
    ELSE ! we read Data from file database
      DIM = CoordinateSystemDimension()
      !------------------------------------------------------------------------------
      ! Inquire and open file
      !------------------------------------------------------------------------------
      ! give preference to a defined material database
      MaterialFileName = GetString( Params, 'Rock Material File', Found )
      IF (.NOT.Found) THEN
        CALL INFO(FunctionName," 'Rock Material File' keyword not found - looking for default DB!")
        fexist = .FALSE.
#ifdef USE_ISO_C_BINDINGS
        str = 'ELMER_LIB'
#else
        str = 'ELMER_LIB'//CHAR(0)
#endif
        CALL envir( str,MaterialFileName,k ) 
        IF ( k > 0  ) THEN
          MaterialFileName = MaterialFileName(1:k) // '/permafrostmaterialdb.dat'
          INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
        END IF
        IF (.NOT. fexist) THEN
#ifdef USE_ISO_C_BINDINGS
          str = 'ELMER_HOME'
#else
          str = 'ELMER_HOME'//CHAR(0)
#endif
          CALL envir( str,MaterialFileName,k ) 
          IF ( k > 0 ) THEN
            MaterialFileName = MaterialFileName(1:k) // '/share/elmersolver/lib/' // 'permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
          END IF
          IF ((.NOT. fexist) .AND. k>0) THEN
            MaterialFileName = MaterialFileName(1:k) // '/permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
          END IF
        END IF
        IF (.NOT. fexist) THEN
          CALL Fatal('CheckKeyWord', 'permafrostmaterialdb.dat not found')
        END IF
      END IF

      ! if we are still here, we open the file (what ever it may be)
      OPEN(unit = io, file = TRIM(MaterialFileName), status = 'old',iostat = ok)
      IF (ok /= 0) THEN
        WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
        CALL FATAL(Trim(FunctionName),Trim(message))
      ELSE
        !------------------------------------------------------------------------------
        ! Read in the number of records in file (first line integer)
        !------------------------------------------------------------------------------
        READ (io, *, END=30, IOSTAT=OK, ERR=40) NumerOfRockRecords, Comment
        WRITE (Message,*) "Attempting to read ",NumerOfRockRecords," ",&
             TRIM(Comment)," from data file ",TRIM(MaterialFileName)
        CALL INFO(FunctionName,Message,level=3)
        InitialNumerOfRockRecords = NumerOfRockRecords
      END IF
      !------------------------------------------------------------------------------
      ! Allocate and read stuff
      !------------------------------------------------------------------------------
      !M = Model % Mesh % NumberOfNodes
      IF (AllocationsDone) THEN
        DEALLOCATE(&
             LocalRockMaterial % ks0th,&
             LocalRockMaterial % e1,&
             LocalRockMaterial % bs,&
             LocalRockMaterial % rhos0,&
             LocalRockMaterial % cs0,&
             LocalRockMaterial % Xi0,&
             LocalRockMaterial % eta0,&
             LocalRockMaterial % hs0,&
             LocalRockMaterial % Kgwh0, &
             LocalRockMaterial % qexp, &
             LocalRockMaterial % alphaL, &
             LocalRockMaterial % alphaT, &
             LocalRockMaterial % RadGen, &
             LocalRockMaterial % acs, &
             LocalRockMaterial % as0, &
             LocalRockMaterial % aas, &
             LocalRockMaterial % ks0, &
             LocalRockMaterial % cks, &
             LocalRockMaterial % acsl, &
             LocalRockMaterial % aasl, &
             LocalRockMaterial % cksl, &
             LocalRockMaterial % VariableBaseName)
      END IF
      ALLOCATE(&
           LocalRockMaterial % ks0th(NumerOfRockRecords),&
           LocalRockMaterial % e1(NumerOfRockRecords),&
           LocalRockMaterial % bs(NumerOfRockRecords),&
           LocalRockMaterial % rhos0(NumerOfRockRecords),&
           LocalRockMaterial % cs0(NumerOfRockRecords),&
           LocalRockMaterial % Xi0(NumerOfRockRecords),&
           LocalRockMaterial % eta0(NumerOfRockRecords),&
           LocalRockMaterial % hs0(NumerOfRockRecords),&
           LocalRockMaterial % Kgwh0(3,3,NumerOfRockRecords),&
           LocalRockMaterial % qexp(NumerOfRockRecords), &
           LocalRockMaterial % alphaL(NumerOfRockRecords), &
           LocalRockMaterial % alphaT(NumerOfRockRecords), &
           LocalRockMaterial % RadGen(NumerOfRockRecords), &
           LocalRockMaterial % acs(0:5,NumerOfRockRecords), &
           LocalRockMaterial % as0(NumerOfRockRecords), &
           LocalRockMaterial % aas(0:5,NumerOfRockRecords), &
           LocalRockMaterial % ks0(NumerOfRockRecords), &
           LocalRockMaterial % cks(0:5,NumerOfRockRecords), &
           LocalRockMaterial % acsl(NumerOfRockRecords), &     
           LocalRockMaterial % aasl(NumerOfRockRecords), &
           LocalRockMaterial % cksl(NumerOfRockRecords), &
           LocalRockMaterial % VariableBaseName(NumerOfRockRecords),&
           STAT=OK)
      AllocationsDone = .TRUE.
      DataRead = .TRUE.
      IF (OK /= 0) THEN
        CLOSE(io)
        CALL FATAL(FunctionName, 'Allocation Error of input data array')
      END IF
      
      DO I=1,NumerOfRockRecords
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % VariableBaseName(I), EntryNumber
        IF (EntryNumber /= I) THEN
          WRITE(Message,'(A,I3,A,I3)') &
               "Entry number", EntryNumber, "does not match expected number ",I
          CLOSE(io)
          CALL FATAL(FunctionName,Message)
        ELSE
          WRITE(Message,'(A,A,A,I3,A)')&
               "Material ", TRIM(LocalRockMaterial % VariableBaseName(I)),&
               " entry number ", EntryNumber, " will be read in"
          CALL INFO(FunctionName,Message,Level=3)
        END IF
        WRITE(Message,'(A,I2,A,A)') "Input for Variable No.",I,": ", LocalRockMaterial % VariableBaseName(I)
        CALL INFO(FunctionName,Message,Level=9)
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % Xi0(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % eta0(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % ks0th(I), Comment          
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % e1(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % bs(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % rhos0(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % cs0(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % hs0(I), Comment
        DO J=1,3
          DO K=1,3
            READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % Kgwh0(J,K,I), Comment
          END DO
        END DO
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % qexp(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % alphaL(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % alphaT(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % RadGen(I), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % acs(0:5,I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % as0(I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % aas(0:5,I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % ks0(I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % cks(0:5,I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % acsl(I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % aasl(I),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) LocalRockMaterial % cksl(I),  Comment
      END DO
      WRITE(Message,'(A,I2,A,A)') "Read ",NumerOfRockRecords," rock material records from file ", TRIM(MaterialFileName)
      CALL INFO(FunctionName,Message,Level=1)
30    CLOSE(io)
      IF (I < NumerOfRockRecords) THEN
        WRITE(Message,'(I3,A,I3)') I,"records read, which is smaller than given number ", NumerOfRockRecords
        CALL FATAL(FunctionName,Message)
      ELSE
        CurrentRockMaterial => LocalRockMaterial
        WRITE(Message,'(A,I2,A,A)') "Read ",NumerOfRockRecords," rock material records from file ", TRIM(MaterialFileName)
        CALL INFO(FunctionName,Message,Level=1)
      END IF
      RETURN
    END IF

40  CALL WARN(FunctionName,"I/O error! Last successfully read variable:")
    CALL WARN(FunctionName,Comment)
    CALL FATAL(FunctionName,"Stopping simulation")    
  END FUNCTION ReadPermafrostRockMaterial
  
  !---------------------------------------------------------------------------------------------  
  FUNCTION ReadPermafrostElementRockMaterial(CurrentRockMaterial,MaterialFileName,NoElements) RESULT(NumberOfRockRecords)
    IMPLICIT NONE
    CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: MaterialFileName
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER, INTENT(IN) :: NoElements
    INTEGER :: NumberOfRockRecords
    !-----------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: SubroutineName="ReadPermafrostElementRockMaterial"
    LOGICAL :: FirstTime=.TRUE.
    TYPE(RockMaterial_t), TARGET :: LocalRockMaterial
    INTEGER :: OK, CurrentNo, I, io
    REAL(KIND=dp) :: ReceivingArray(50)

    SAVE LocalRockMaterial, FirstTime

    IF (FirstTime) THEN
      ALLOCATE(&
           LocalRockMaterial % ks0th(NoElements),&
           LocalRockMaterial % e1(NoElements),&
           LocalRockMaterial % bs(NoElements),&
           LocalRockMaterial % rhos0(NoElements),&
           LocalRockMaterial % cs0(NoElements),&
           LocalRockMaterial % Xi0(NoElements),&
           LocalRockMaterial % eta0(NoElements),&
           LocalRockMaterial % hs0(NoElements),&
           LocalRockMaterial % Kgwh0(3,3,NoElements),&
           LocalRockMaterial % qexp(NoElements), &
           LocalRockMaterial % alphaL(NoElements), &
           LocalRockMaterial % alphaT(NoElements), &
           LocalRockMaterial % RadGen(NoElements), &           
           LocalRockMaterial % acs(NoElements,0:5), &
           LocalRockMaterial % as0(NoElements), &
           LocalRockMaterial % aas(0:5,NoElements), &
           LocalRockMaterial % ks0(NoElements), &
           LocalRockMaterial % cks(0:5,NoElements), &
           LocalRockMaterial % acsl(NoElements), &
           LocalRockMaterial % aasl(NoElements), &
           LocalRockMaterial % cksl(NoElements), &
           LocalRockMaterial % VariableBaseName(NoElements),&
           STAT=OK)
      OPEN(unit = io, file = TRIM(MaterialFileName), status = 'old',iostat = ok)
      IF (ok /= 0) THEN
        WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
        CALL FATAL(Trim(SubroutineName),Trim(message))        
      ELSE        
        !------------------------------------------------------------------------------
        ! Read in the number of records in file (first line integer)
        !------------------------------------------------------------------------------
        WRITE (Message,'(A,I2,A,A,A,A)') "Attempting read ",NoElements,&
             " from data file ",TRIM(MaterialFileName)
        CALL INFO(SubroutineName,Message,level=3)
        DO I=1,NoElements
          READ (io, *, END=50, ERR=60, IOSTAT=OK) CurrentNo, ReceivingArray(1:50)
          LocalRockMaterial % ks0th(I) = ReceivingArray(12) ! shall be changed to tensor
          LocalRockMaterial % e1(I) = ReceivingArray(34) ! e1 (mail from Juha 11.10.)
          LocalRockMaterial % bs(I) = ReceivingArray(24) ! b11,1 (mail from Juha 11.10.)
          LocalRockMaterial % rhos0(I) = ReceivingArray(2)
          LocalRockMaterial % cs0(I) = ReceivingArray(15)
          LocalRockMaterial % Xi0(I) = ReceivingArray(33)
          LocalRockMaterial % eta0(I) = ReceivingArray(31) ! eta_t (mail from Juha 11.10.)
          LocalRockMaterial % hs0(I) = 0.0_dp! will be removed
          LocalRockMaterial % Kgwh0(1,1,I) = ReceivingArray(36)
          LocalRockMaterial % Kgwh0(2,2,I) = ReceivingArray(37)
          LocalRockMaterial % Kgwh0(3,3,I) = ReceivingArray(38)
          LocalRockMaterial % Kgwh0(1,2,I) = ReceivingArray(39)
          LocalRockMaterial % Kgwh0(1,3,I) = ReceivingArray(40)
          LocalRockMaterial % Kgwh0(2,3,I) = ReceivingArray(41)
          LocalRockMaterial % Kgwh0(2,1,I) = LocalRockMaterial % Kgwh0(1,2,I)
          LocalRockMaterial % Kgwh0(3,1,I) = LocalRockMaterial % Kgwh0(1,3,I)
          LocalRockMaterial % Kgwh0(3,2,I) = LocalRockMaterial % Kgwh0(2,3,I)
          LocalRockMaterial % qexp(I) = ReceivingArray(42)
          LocalRockMaterial % alphaL(I) = ReceivingArray(48)
          LocalRockMaterial % alphaT(I) = ReceivingArray(49)
          LocalRockMaterial % RadGen(I) = ReceivingArray(30)
          LocalRockMaterial % acs(NoElements,0) =  ReceivingArray(16)
          LocalRockMaterial % acs(NoElements,1) =  ReceivingArray(17)
          LocalRockMaterial % acs(NoElements,2) =  ReceivingArray(18)
          LocalRockMaterial % acs(NoElements,2:5) = 0.0_dp
          !LocalRockMaterial % bcs(NoElements,0) =  ReceivingArray(19)
          !LocalRockMaterial % bcs(NoElements,1) =  ReceivingArray(20)
          !LocalRockMaterial % bcs(NoElements,2) =  ReceivingArray(21)
          LocalRockMaterial % as0(NoElements)= ReceivingArray(3)
          LocalRockMaterial % aas(NoElements,0) =  ReceivingArray(4)
          LocalRockMaterial % aas(NoElements,1) =  ReceivingArray(5)
          LocalRockMaterial % aas(NoElements,2) =  ReceivingArray(6)
          LocalRockMaterial % aas(NoElements,3) =  ReceivingArray(7)
          LocalRockMaterial % aas(NoElements,4) =  ReceivingArray(8)
          LocalRockMaterial % aas(NoElements,5) =  ReceivingArray(9)
          LocalRockMaterial % ks0(NoElements)= ReceivingArray(13)
          LocalRockMaterial % cks(0,NoElements)= ReceivingArray(14)
          LocalRockMaterial % cks(1:5,NoElements)= 0.0_dp
          LocalRockMaterial % acsl(NoElements)= 2
          LocalRockMaterial % aasl(NoElements)= 5
          LocalRockMaterial % cksl(NoElements)= 0
          WRITE(Message,*) 'Element',I
          LocalRockMaterial % VariableBaseName(I) = TRIM(Message)
        END DO
        LocalRockMaterial % NumerOfRockRecords = NoElements
        NumberOfRockRecords = NoElements
50      CLOSE(io)
        IF (CurrentNo < NoElements) THEN
          WRITE (Message,*) 'Found only ',CurrentNo,' entries in file ',TRIM(MaterialFileName),&
               ' for ', NoElements, ' elements in mesh.'
          CALL FATAL(TRIM(SubroutineName),Message)
        ELSE
          WRITE (Message,*) 'Read ',CurrentNo,' entries in file ',TRIM(MaterialFileName)
          CALL INFO(TRIM(SubroutineName),Message,Level=3)
        END IF
      END IF
      CurrentRockMaterial => LocalRockMaterial
      FirstTime = .FALSE.
    ELSE
      CurrentRockMaterial => LocalRockMaterial
      NumberOfRockRecords = NoElements
    END IF
    RETURN
60  WRITE (Message,*) 'I/O error at entry ',CurrentNo,' of file ',TRIM(MaterialFileName)
    CALL FATAL(TRIM(SubroutineName),Message)
  END FUNCTION ReadPermafrostElementRockMaterial
  
  !---------------------------------------------------------------------------------------------
  FUNCTION ReadPermafrostConstants(Model, FunctionName,&
       DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity) RESULT(Constantsread)
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: DIM
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps, Gravity(3) ! constants read only once 
    LOGICAL :: Constantsread
    !------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: gWork(:,:)
    LOGICAL :: Found
    INTEGER :: I
    !------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    gWork => ListGetConstRealArray( Model % Constants,"Gravity",Found)
    IF (.NOT.Found) THEN
      Gravity = 0.0
      CALL WARN(FunctionName,'Gravity not found in Constants section. Setting to zero')
    ELSE
      Gravity = gWork(1:3,1)*gWork(4,1)
    END IF
    !------------------------------------------------------------------------------
    ! Constants
    ! GasConstant, N0, T0, p0, DeltaT, eps
    !------------------------------------------------------------------------------
    !IF (.NOT.ASSOCIATED(Model % Constants)) STOP
    GasConstant = GetConstReal(Model % Constants, "Gas Constant", Found)
    IF (.NOT.Found) THEN
      GasConstant = 8.3145_dp
      CALL INFO(FunctionName, ' "Gas Constant" not found in Constants and set to default value 8.3145',Level=3)
    END IF
    N0 = GetConstReal(Model % Constants, "Avogadro Number", Found)
    IF (.NOT.Found) THEN
      N0 = 6.022140857d23
      CALL INFO(FunctionName, ' "Avogadro Number" not found in Constants and set to default value 6.022140857E23',Level=3)
    END IF
    T0 = GetConstReal(Model % Constants, 'Reference Temperature', Found)
    IF (.NOT.Found) THEN
      T0 = 273.15_dp
      CALL INFO(FunctionName, ' "Reference Temperature" not found in Constants and set to default value T0=273.15',Level=3)
    END IF
    p0 = GetConstReal(Model % Constants, 'Reference Pressure', Found)
    IF (.NOT.Found) THEN
      p0 = 100132.0_dp
      CALL INFO(FunctionName, ' "Reference Pressure not found in Constants and set to default value p0=100132.0',Level=3)
    END IF
    DeltaT = GetConstReal(Model % Constants,"Permafrost DeltaT",Found)
    IF (.NOT.Found) THEN
      DeltaT = 1.0_dp
      CALL INFO(FunctionName, ' "Permafrost DeltaT" not found in Constants and set to default value DeltaT=1.0',Level=3)
    END IF
    Eps = GetConstReal(Model % Constants,"Permafrost eps",Found)
    IF (.NOT.Found) THEN
      eps = 0.99_dp
      CALL INFO(FunctionName, ' "Permafrost eps" not found in Constants and set to default value eps=0.99',Level=3)
    END IF
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(FunctionName,"Model Constants:", Level=9)
    WRITE(Message,*) "GasConstant, T0, p0, DeltaT, eps:"
    CALL INFO(FunctionName,Message, Level=9)
    WRITE(Message,*) GasConstant, T0, p0, DeltaT, eps
    CALL INFO(FunctionName,Message, Level=9)
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
  END FUNCTION ReadPermafrostConstants
  
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------
  ! model parameters functions
  !---------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------------------
  ! general functions 
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION GeneralPolynomial(Variable,ReferenceValue,Normation,coeff,pdeg)
    IMPLICIT NONE
    !-------
    REAL(KIND=dp), INTENT(IN) :: Variable,ReferenceValue,Normation,coeff(0:5)
    INTEGER, INTENT(IN) :: pdeg
    REAL(KIND=dp) outval
    ! ------
    REAL(KIND=dp) currpot
    INTEGER :: i

    outval = 0.0_dp
    currpot = 1.0_dp
    DO i=0,pdeg
      !IF (pdeg == 0) PRINT*,i,"of",pdeg, currpot, Variable, ReferenceValue, Normation
      outval = outval + coeff(i) * currpot
      currpot = currpot * (Variable - ReferenceValue)/Normation
    END DO
    GeneralPolynomial = outval
  END FUNCTION GeneralPolynomial
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION GeneralIntegral(Variable,ReferenceValue,Normation,coeff0,coeff,pdeg)
    IMPLICIT NONE
    !-------
    REAL(KIND=dp), INTENT(IN) :: Variable,ReferenceValue,Normation,coeff0,coeff(0:5)
    INTEGER, INTENT(IN) :: pdeg
    REAL(KIND=dp) outval
    ! ------
    REAL(KIND=dp) currpot
    INTEGER :: currdeg
    !PRINT *,"Integral"
    outval = coeff0*(Variable - ReferenceValue)
    currpot = 1.0_dp
    DO currdeg=0,pdeg
      outval = outval * coeff(currdeg) * currpot/(DBLE(currdeg) + 1.0_dp)
      currpot = currpot * (Variable - ReferenceValue)/Normation
    END DO
    GeneralIntegral = outval
  END FUNCTION GeneralIntegral
  !---------------------------------------------------------------------------------------------
  ! functions specific to heat transfer and phase change
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXiAnderson(A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) RESULT(XiAnderson)
    REAL(KIND=dp), INTENT(IN) :: A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar, XiAnderson
    IF (Porosity <= 0.0) &
         CALL FATAL("Permafrost(GetXiAnderson)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    XiAnderson =  MAX(MIN((rhos0/rhow)*(A*(Tstar**B)/Porosity),1.0_dp),0.0_dp)
  END FUNCTION GetXiAnderson
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonT(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity)
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0) &
             CALL FATAL("Permafrost(GetXiAndersonT)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1.0_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonT = 0.0_dp
    ELSE
      XiAndersonT = -(rhos0/rhow)*(A*B*(Tstar**(B - 1.0_dp)))/Porosity
    END IF
  END FUNCTION XiAndersonT
    !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonP(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity)
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0_dp) &
             CALL FATAL("Permafrost(GetXiAndersonT)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonP = 0.0_dp
    ELSE
      XiAndersonP = -Beta*(rhos0/rhow)*(A*B*(Tstar**(B - 1.0_dp)))/Porosity
    END IF
  END FUNCTION XiAndersonP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonEta(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) 
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0) &
         CALL FATAL("Permafrost(GetXiAndersonEta)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1.0_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonEta = 0.0_dp
    ELSE
      XiAndersonEta =  -(rhos0/rhow)*(A*(Tstar**B))/(Porosity*Porosity)
    END IF
  END FUNCTION XiAndersonEta
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION delta(CurrentSolventMaterial,&
       eps,DeltaT,T0,GasConstant)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: eps,DeltaT,T0,GasConstant
    REAL(KIND=dp) :: aux,Mw,hi0,cw0,ci0
    LOGICAL :: FirstTime=.TRUE.
    SAVE FirstTime,Mw,hi0,ci0,cw0
    IF(FirstTime) THEN
      !assign needed variables
      Mw = CurrentSolventMaterial % Mw
      hi0 = CurrentSolventMaterial % hi0
      ci0 = CurrentSolventMaterial % ci0
      cw0 = CurrentSolventMaterial % cw0
      FirstTime = .FALSE.
    END IF
    aux = -0.5_dp*hi0*DeltaT/T0 &
         + (cw0 - ci0)*((T0 + 0.5_dp*DeltaT)*LOG(1.0_dp + 0.5_dp*DeltaT/T0) - 0.5_dp*DeltaT)
    delta = aux*(eps*(1.0_dp - eps)/(2.0_dp*eps - 1.0_dp))* Mw/(GasConstant*(T0 + 0.5_dp*DeltaT))
    IF (delta < 1.0d-10) PRINT *, "delta=", delta
  END FUNCTION delta
  !---------------------------------------------------------------------------------------------
  FUNCTION GetAcwtilde(CurrentSolventMaterial) RESULT(acwtilde)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: acwtilde(0:5)
    REAL(KIND=dp) :: cw0,kw0,acw(0:5)
    INTEGER :: acwl,I
 
    !assign needed properties
    acwl = CurrentSolventMaterial % acwl
    acw(0:5) =  CurrentSolventMaterial % acw(0:5)
    acwtilde(0:5) = 0.0_dp
    !compute
    I=acwl
    acwtilde(0)= (-(1.0_dp)**(1.0_dp*I))*acw(i)
    acwtilde(i)= (-(1.0_dp)**(1.0_dp*I))*(1.0_dp/(1.0_dp*I) - (1.0_dp/(1.0_dp*I + 1.0_dp)))*acw(i)
    DO I=acwl-1,1,-1
      acwtilde(0)= acwtilde(0) + (-1.0_dp**(1.0_dp*I))*acw(i)
      acwtilde(i) = acwtilde(i+1) + (-1.0_dp**(1.0_dp*I))*(1.0_dp/(1.0_dp*I) &
           - (1.0_dp/(1.0_dp*I + 1.0_dp)))*acw(i)
    END DO
    acwtilde(0)= acwtilde(0) + acw(0)
  END FUNCTION GetAcwtilde
   !---------------------------------------------------------------------------------------------
  FUNCTION GetAcitilde(CurrentSolventMaterial) RESULT(acitilde)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: acitilde(0:5)
    REAL(KIND=dp) :: aci(0:5)
    INTEGER :: acil,I
    !assign needed properties
    acitilde(0:5) = 0.0_dp
    acil = CurrentSolventMaterial % acil
    aci(0:5) = CurrentSolventMaterial % aci(0:5)
    ! compute
    I=acil    
    acitilde(0)= (-(1.0_dp)**(1.0_dp*I))*aci(I)
    acitilde(I)= (-(1.0_dp)**(1.0_dp*I))*(1.0_dp/(1.0_dp*I) - (1.0_dp/(1.0_dp*I + 1.0_dp)))*aci(i)
    DO I=acil-1,1,-1
      acitilde(0)= acitilde(0) + (-1.0_dp**(1.0_dp*I))*aci(i)
      acitilde(i) = acitilde(i+1) + (-1.0_dp**(1.0_dp*I))*(1.0_dp/(1.0_dp*I)&
           - (1.0_dp/(1.0_dp*I + 1.0_dp)))*aci(i)
    END DO
    acitilde(0)= acitilde(0) + aci(0)
  END FUNCTION GetAcitilde
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gwa(CurrentSolventMaterial,&
       p0,T0,rhow,Xi,Temperature,Pressure,Salinity)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhow,Xi,Temperature,Pressure,Salinity
    REAL(KIND=dp) :: cw0,kw0,acw(0:5),bcw(0:5)
    INTEGER :: I, acwl,bcwl
    REAL(KIND=dp) :: watercont,acwtilde(0:5),aux,aux1
    LOGICAL :: FirstTime=.TRUE.

    SAVE FirstTime,cw0,kw0,acw,acwl,bcw,bcwl,acwtilde

    IF (FirstTime) THEN
      !assign needed properties
      cw0 = CurrentSolventMaterial % cw0
      kw0 = CurrentSolventMaterial %  kw0
      acw(0:5) = CurrentSolventMaterial %  acw(0:5)
      bcw(0:5) = CurrentSolventMaterial %  bcw(0:5)
      acwl = CurrentSolventMaterial %  acwl
      bcwl = CurrentSolventMaterial %  bcwl
      !acwtilde = GetAcwtilde(CurrentSolventMaterial)
      FirstTime = .FALSE.
    END IF
    
    watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
    !simplified version
 !   aux = -cw0*( Temperature * LOG(Temperature/T0) &
    !        - (Temperature - T0) )
    !original version
 !   aux = -cw0*(acwtilde(0) * Temperature * LOG(Temperature/T0) &
 !        - (Temperature - T0) * GeneralPolynomial(Temperature,T0,T0,acwtilde,acwl))
    !   aux = aux * GeneralPolynomial(watercont,1.0_dp,1.0_dp,bcw,bcwl)
    !    gwa = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*kw0*(Pressure - p0))/rhow)
    ! written-out version
    aux1=(Temperature/T0 - 1.0_dp)
    aux = -cw0*((acw(0) - acw(1) + acw(2))* Temperature * LOG(Temperature/T0) &
         - (acw(0) - acw(1) + acw(2) - (0.5_dp* acw(1) -&
         acw(2)/6.0_dp)*aux1 + acw(2)*(aux1**2.0_dp)/6.0_dp )*(Temperature - T0))
    aux = aux * GeneralPolynomial(watercont,1.0_dp,1.0_dp,bcw,bcwl)
    gwa = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*kw0*(Pressure - p0))/rhow)
    IF (gwa .NE. gwa) THEN
      PRINT *, "gwa"
      STOP
    END IF
    !PRINT *,"gwa:",gwa,aux,((Pressure - p0)*(1.0_dp + 0.5_dp*kw0*(Pressure - p0))/rhow),&
    !     GeneralPolynomial(watercont,1.0_dp,1.0_dp,bcw,bcwl)
  END FUNCTION gwa
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gia(CurrentSolventMaterial,&
       p0,T0,rhoi,Temperature,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhoi,Temperature,Pressure
    !---------------------------------
    REAL(KIND=dp) :: ci0,ki0,hi0,aci(0:5)
    INTEGER :: acil,I
    REAL(KIND=dp) ::acitilde(0:5),aux,aux1
    LOGICAL :: FirstTime=.TRUE.
    !-----------------------------------
    SAVE FirstTime,ci0,ki0,hi0,aci,acil,acitilde

    IF (FirstTime) THEN
      ci0 = CurrentSolventMaterial %  ci0 
      ki0 = CurrentSolventMaterial %  ki0
      hi0 = CurrentSolventMaterial %  hi0
      acil = CurrentSolventMaterial % acil
      aci(0:5) = CurrentSolventMaterial % aci(0:5)
      !acitilde = GetAcitilde(CurrentSolventMaterial)
      FirstTime = .FALSE.
    END IF
    ! simple version
    !aux = -hi0*((Temperature - T0)/T0)-ci0*(Temperature * LOG(Temperature/T0) &
    !     - (Temperature - T0) )
    !original version
    !aux = -hi0*((Temperature - T0)/T0)-ci0*(acitilde(0) * Temperature * LOG(Temperature/T0) &
    !      - (Temperature - T0) * GeneralPolynomial(Temperature,T0,T0,acitilde,acil))
    !gia = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*ki0*(Pressure - p0))/rhoi)
    ! written-out version
    aux1 = (Temperature/T0 - 1.0_dp)
    aux = -hi0*aux1 &
         - ci0*((aci(0) - aci(1))* Temperature * LOG(Temperature/T0)&
         - (aci(0) - aci(1) - 0.5_dp*aci(1)*aux1)*(Temperature - T0))
    gia = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*ki0*(Pressure - p0))/rhoi)
  END FUNCTION gia
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gwaT(CurrentSolventMaterial,&
       p0,T0,rhow,Xi,Temperature,Salinity)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhow,Xi,Temperature,Salinity
    !----------------------------
    INTEGER :: I
    REAL(KIND=dp) :: watercont,acwtilde(0:5),aux, aux1
    LOGICAL :: FirstTime=.TRUE.

    !SAVE FirstTime, acwtilde

    !IF (FirstTime) THEN
 
      !acwtilde = GetAcwtilde(CurrentSolventMaterial)
      !DO I=1,acwl
      !  acwtilde(I) = (1.0_dp + (1.0_dp*I)) * acwtilde(I)
      !END DO
    !  FirstTime = .FALSE.
    !END IF
    

    watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
    !original version
!     aux = -cw0*(acwtilde(0) * (1.0_dp + LOG(Temperature/T0)) &
!         -  GeneralPolynomial(Temperature,T0,T0,acwtilde,acwl)) 
    !gwaT = aux !* GeneralPolynomial(watercont,1.0_dp,1.0_dp,bcw,bcwl)
    ! written-out version
    aux1= (Temperature/T0 - 1.0_dp)
    aux = -(CurrentSolventMaterial %cw0)*&
         ((CurrentSolventMaterial % acw(0) - CurrentSolventMaterial % acw(1) &
         + CurrentSolventMaterial % acw(2)) * LOG(Temperature/T0) &
         + (CurrentSolventMaterial % acw(1) - CurrentSolventMaterial % acw(2)/3.0_dp)*aux1 &
         - 0.5_dp*(CurrentSolventMaterial % acw(2))*(aux1**2.0_dp))
    !PRINT *,"gwaT:",aux, GeneralPolynomial(watercont,1.0_dp,1.0_dp,bcw,bcwl),acw,acwl,T0,Temperature,cw0
    gwaT = aux * GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
         CurrentSolventMaterial % bcw(0:5),&
         CurrentSolventMaterial % bcwl)
    IF (gwaT .NE. gwaT) THEN
      PRINT *, "gwaT"
      STOP
    END IF
  END FUNCTION gwaT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION giaT(CurrentSolventMaterial,&
       p0,T0,rhoi,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhoi,Temperature
    !---------------------------------
    !INTEGER :: I
    !REAL(KIND=dp) ::acitilde(0:5)
    !LOGICAL :: FirstTime=.TRUE.
    !-----------------------------------
    !SAVE FirstTime,acitilde

   ! IF (FirstTime) THEN
      !acitilde = GetAcitilde(CurrentSolventMaterial)
      !DO I=1,acil
      !  acitilde(I) = (1.0_dp + (1.0_dp*I)) * acitilde(I)
      !END DO
      !FirstTime = .FALSE.
    !END IF
    ! original version
 !   giaT = -hi0/T0 - ci0*(acitilde(0) * (1.0_dp + LOG(Temperature/T0)) &
    !         -  GeneralPolynomial(Temperature,T0,T0,acitilde,acil))
    ! written-out version
    giaT = -(CurrentSolventMaterial % hi0)/T0 &
         - (CurrentSolventMaterial %  ci0)*&
         ((CurrentSolventMaterial % aci(0) - CurrentSolventMaterial % aci(1)) * ( LOG(Temperature/T0)) &
         + ( (CurrentSolventMaterial % aci(1)) *(Temperature/T0 - 1.0_dp)))
  END FUNCTION giaT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION deltaG(gwa,gia)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: gwa,gia
    deltaG = gwa - gia
  END FUNCTION deltaG  
  !---------------------------------------------------------------------------------------------
  FUNCTION GetBij(CurrentSoluteMaterial,&
       Xi0tilde,Salinity) RESULT(bij)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,Salinity
    REAL(KIND=dp):: bij(2,2)
    !----------
    REAL(KIND=dp)::  aux,d1,d2

    d1 = CurrentSoluteMaterial % d1
    d2 = CurrentSoluteMaterial % d2
    
    aux = Salinity/(Xi0tilde - Salinity)
    bij(1,1) = aux*(d1 + 0.5_dp*d2*aux)
    bij(1,2) = aux*(d1 + d2*aux)*Xi0tilde/(Xi0tilde - Salinity)
    aux = Salinity/(1.0_dp - Salinity)
    bij(2,1) = aux*(d1 + 0.5_dp*d2*aux)
    bij(2,2) = aux*((d1 + d2*aux))/(1.0_dp - Salinity)     
  END FUNCTION GetBij
    !---------------------------------------------------------------------------------------------
  FUNCTION GetBijYc(CurrentSoluteMaterial,&
       Xi0tilde,bij,Salinity) RESULT(bijYc)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,bij(2,2),Salinity
    REAL(KIND=dp):: bijYc(2,2)
    !----------
    REAL(KIND=dp)::  aux,d1,d2

    d1 = CurrentSoluteMaterial % d1
    d2 = CurrentSoluteMaterial % d2

    aux = 1.0_dp/(Xi0tilde - Salinity)
    bijYc(1,1) = (d1 + d2*Salinity*aux)*Xi0tilde*aux*aux
    bijYc(1,2) = aux*((d1 + 2.0_dp*d2*Salinity*aux)*(Xi0tilde**2.0_dp)*aux*aux + bij(1,2))
    aux = 1.0_dp/(1.0_dp - Salinity)
    bijYc(2,1) = aux*(d1 + d2*Salinity*aux)
    bijYc(2,2) = aux*((d1 + 2.0_dp*d2*Salinity*aux)*aux*aux + bij(2,2))        
  END FUNCTION GetBijYc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
       delta,deltaG,GasConstant,bij,Temperature) RESULT(B1)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: delta,deltaG,GasConstant,bij(2,2),Temperature
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: B1
    REAL(KIND=dp) :: e1,Mw
    Mw = CurrentSolventMaterial % Mw
    e1 = CurrentRockMaterial % e1(RockMaterialID)
    B1 =(Mw*deltaG/(GasConstant*Temperature) - bij(1,1))/(e1 + delta + bij(1,2))
    IF (B1 .NE. B1) THEN
      PRINT *, "B1"
      STOP
    END IF
  END FUNCTION GetB1
  !---------------------------------------------------------------------------------------------
  FUNCTION GetB2(CurrentSolventMaterial,delta,deltaG,GasConstant,bij,Temperature) RESULT(B2)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: delta,deltaG,GasConstant,bij(2,2),Temperature
    REAL(KIND=dp) :: B2
    REAL(KIND=dp) :: Mw
    IF ((delta + bij(2,2)) == 0.0_dp) CALL FATAL("GetB2","Singular expression")
    Mw = CurrentSolventMaterial % Mw
    B2 = (Mw*deltaG/(GasConstant*Temperature) - bij(2,1))/(delta + bij(2,2))
    IF (B2 .NE. B2) THEN
      PRINT *, "B2"
      STOP
    END IF
    !PRINT *,"B2:",Mw,deltaG,GasConstant,Temperature,bij(2,1),bij(2,2),delta
  END FUNCTION GetB2
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION D1(CurrentRockMaterial,RockMaterialID,delta,bij)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: delta,bij(2,2)
    REAL(KIND=dp) :: e1
    e1 = CurrentRockMaterial % e1(RockMaterialID)
    ! local
    D1 = delta/(e1 + delta + bij(1,2))
    IF (D1 .NE. D1) THEN
      PRINT *, "D1"
      STOP
    END IF
  END FUNCTION D1
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION D2(delta,bij)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,bij(2,2)
    ! local
    D2 = delta/(delta + bij(2,2))
    IF (D2 .NE. D2) THEN
      PRINT *, "D2"
      STOP
    END IF
  END FUNCTION D2
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXi0Tilde(CurrentRockMaterial,RockMaterialID,Porosity) RESULT(Xi0tilde)
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: Porosity
    REAL(KIND=dp) :: Xi0tilde
    REAL(KIND=dp) :: Xi0,eta0
    LOGICAL :: FirstTime = .TRUE.
    SAVE FirstTime
    
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    IF (Porosity <= 0.0_dp) THEN
      IF (Xi0 == 0.0_dp) THEN
        Xi0tilde = 1.0_dp
      ELSE
        CALL FATAL("Permafrost(GetXi)","Zero or negative porosity detected")
      END IF
    ELSE
      Xi0tilde = MIN(Xi0 * (eta0/Porosity) * (1.0_dp - Porosity)/(1.0_dp - eta0),1.0_dp)
    END IF
  END FUNCTION GetXi0Tilde
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
       Xi0tilde,rhow,Xi,GasConstant,Temperature)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,rhow,Xi,GasConstant,Temperature
    !-----------------
    REAL(KIND=dp) :: Mw, e1
    !----------------
 !   IF (Xi > Xi0tilde) THEN
 !     fw = 0.0_dp
 !   ELSE
 !     e1 = CurrentRockMaterial % e1(RockMaterialID)
 !     Mw = CurrentSolventMaterial % Mw
 !     fw = rhow*GasConstant*Temperature*e1*Xi0tilde/(Mw*Xi)
 !   END IF
    fw = 0.0_dp !! CHANGE BACK WHEN JUHA TELLS US TO DO SO !!
  END FUNCTION fw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXi(B1,B2,D1,D2,Xi0tilde) RESULT(Xi)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0tilde
    REAL(KIND=dp) :: Xi
    Xi= Xi0tilde/(1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1)) &
         + (1.0_dp - Xi0tilde)/(1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))
    IF (Xi < 0.0_dp) Xi = 0.0_dp
    IF (Xi > 1.0_dp) Xi = 1.0_dp
  END FUNCTION GetXi
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
       B1,B2,D1,D2,Xi0Tilde,bij,p0,delta,deltaG,T0,gwaT,giaT,GasConstant,Temperature, Pressure)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0Tilde,bij(2,2),p0,delta,deltaG,&
         T0,gwaT,giaT,GasConstant,Temperature, Pressure    
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3, Mw,e1, hi0,hw0,rhow0,rhoi0,cw0,ci0
    LOGICAL :: FirstTime=.TRUE.

    SAVE Mw,hi0,hw0,rhow0,rhoi0,cw0,ci0, FirstTime
    
    IF (FirstTime) THEN
      Mw = CurrentSolventMaterial % Mw
      hi0   = CurrentSolventMaterial % hi0  
      hw0   = CurrentSolventMaterial % hw0  
      rhow0 = CurrentSolventMaterial % rhow0
      rhoi0 = CurrentSolventMaterial % rhoi0
      cw0   = CurrentSolventMaterial % cw0  
      ci0   = CurrentSolventMaterial % ci0  
      FirstTime=.FALSE.
    END IF
    e1 = CurrentRockMaterial % e1(RockMaterialID)
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
!    aux3 = ((hw0 - hi0) + (cw0 - ci0)*(Temperature - T0) & 
    !         + (-(1.0_dp/rhoi0) + (1.0_dp/rhow0))*(Pressure - p0))/Temperature
    aux3 = (deltaG/Temperature -(gwaT - giaT))
    !if (aux2 .NE. aux2) THEN
    !  PRINT *,"Problem in XiT:B2=",B2,"D2=",D2     
    !END if
    !    IF ((Temperature < 273.65) .AND. (Temperature > 272.65))&
    !IF (Temperature > 274.0)&
    !     PRINT *,"XiT:",0.5_dp*Xi0Tilde*aux1/(e1 + delta + bij(1,2)),&
    !     0.5_dp*(1.0_dp - Xi0Tilde)*aux2/(delta + bij(2,2)),&
    !     Mw * aux3/(GasConstant*Temperature)
    XiT = (0.5_dp*Xi0Tilde*aux1/(e1 + delta + bij(1,2)) &
         + 0.5_dp*(1.0_dp - Xi0Tilde)*aux2/(delta + bij(2,2)))&
         *Mw * aux3/(GasConstant*Temperature)
  END FUNCTION XiT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
       B1,B2,D1,D2,bij,gwap,giap,Xi0Tilde,delta,GasConstant,Temperature)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,bij(2,2),gwap,giap,Xi0Tilde,delta,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2, rhow0,rhoi0, Mw,e1
    LOGICAL :: FirstTime=.TRUE.
    
    SAVE Mw,rhow0,rhoi0,FirstTime
    IF (FirstTime) THEN
      Mw = CurrentSolventMaterial % Mw
      rhow0 = CurrentSolventMaterial % rhow0
      rhoi0 = CurrentSolventMaterial % rhoi0
      FirstTime=.FALSE.
    END IF
    e1 = CurrentRockMaterial % e1(RockMaterialID)
    IF (Temperature <= 0.0_dp) CALL FATAL("Permafrost (XiP)","(sub-)Zero Temperature detected")
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
    XiP = (0.5_dp*Xi0Tilde*aux1/(e1 + delta + bij(1,2)) + 0.5_dp*(1.0_dp - Xi0Tilde)*aux2/(delta + bij(2,2))) &
         *(giap - gwap)* Mw/(GasConstant*Temperature)
  END FUNCTION XiP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiYc(CurrentRockMaterial,RockMaterialID,&
       B1,B2,D1,D2,bij,bijYc,Xi0Tilde,Xi,delta,Salinity)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: B1,B2,bij(2,2),bijYc(2,2),D1,D2,&
         Xi0Tilde,Xi,delta,Salinity
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3, e1

    e1 = CurrentRockMaterial % e1(RockMaterialID)
    IF ((Salinity < 1.0_dp) .AND. (Salinity >= 0.0_dp)) THEN
      aux1 = 0.5_dp*((1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))*(bijYc(1,1) + B1*bijYc(1,2)) + bijYc(1,2)/SQRT(B1*B1 + 4.0_dp*D1))
      aux1 = aux1*(0.5_dp*Xi0Tilde*aux1/(e1 + delta + bij(1,2)))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)      
      aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))*(bijYc(2,1) + B2*bijYc(2,2)) + bijYc(1,2)/SQRT(B2*B2 + 4.0_dp*D2)
      aux2 = 0.5_dp*((1.0_dp - Xi0Tilde)/delta + bij(2,2))*aux2/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
      XiYc = aux1 + aux2
      IF (XiYc .NE. XiYc) THEN
        PRINT *, "XiYc"
        STOP
      END IF
    ELSE
      CALL WARN("Permafrost(XiXc)","Salinity out of physical range - returning zero")
      XiYc = 0.0_dp
    END IF
    
  END FUNCTION XiYc
  
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiEta(CurrentRockMaterial,RockMaterialID,&
       B1,B2,D1,D2,Xi0,Xi0Tilde,eta0,Porosity)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0Tilde,Porosity
    !local
    REAL(KIND=dp) :: aux1, aux2, Xi0,eta0

    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    
    IF (Porosity >= 0.0_dp) THEN
      IF (Xi0Tilde < 1.0_dp) THEN
        aux1 = 1.0_dp/(1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))
        aux2 = 1.0_dp/(1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))
        XiEta = (aux1 + aux2) * (Xi0*eta0/(1.0_dp - eta0))*(1.0_dp/(Porosity**2.0_dp))
      ELSE
        XiEta = 0.0_dp
      END IF
    ELSE
      CALL WARN("Permafrost(XiEta)","Porosity out of physical range - returning zero")
      XiEta = 0.0_dp
    END IF
  END FUNCTION XiEta
  !---------------------------------------------------------------------------------------------
  ! Densities and their derivatives, thermal expansion, isothermal chemical compaction and
  !     compressibility coefficients
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhos(CurrentRockMaterial,RockMaterialID,&
       T0,p0,Temperature,Pressure,ConstVal)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Temperature,Pressure
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1,aux2
    !----------------------
    rhos = CurrentRockMaterial % rhos0(RockMaterialID)
    IF (.NOT.ConstVal) THEN
      aux1 = GeneralIntegral(Pressure,p0,p0,&
           CurrentRockMaterial % ks0(RockMaterialID),&
           CurrentRockMaterial % cks(0:5,RockMaterialID),&
           CurrentRockMaterial % cksl(RockMaterialID))
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentRockMaterial % as0(RockMaterialID),&
           CurrentRockMaterial % aas(0:5,RockMaterialID),&
           CurrentRockMaterial % aasl(RockMaterialID))
      rhos = rhos * EXP(aux1 - aux2)
    END IF
  END FUNCTION rhos
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhosT(CurrentRockMaterial,RockMaterialID,rhos,T0,Temperature)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    REAL(KIND=dp), INTENT(IN) :: rhos,T0,Temperature
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: alphaS
    
    alphaS = CurrentRockMaterial % as0(RockMaterialID) *&
	 GeneralPolynomial(Temperature,T0,T0,&
         CurrentRockMaterial % aas(0:5,RockMaterialID),&
         CurrentRockMaterial % aasl(RockMaterialID))
    rhosT = rhos * alphaS
  END FUNCTION rhosT
  
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhow(CurrentSolventMaterial,&
       T0,p0,Xi,Temperature,Pressure,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Xi,Temperature,Pressure,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, aux3, watercont
    !----------------------
    rhow = CurrentSolventMaterial % rhow0
    IF (.NOT.ConstVal) THEN
      aux1 = GeneralIntegral(Pressure,p0,p0,&
           CurrentSolventMaterial % kw0,&
           CurrentSolventMaterial % ckw(0:5),&
           CurrentSolventMaterial % ckwl)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % aw0,&
           CurrentSolventMaterial % aaw(0:5),&
           CurrentSolventMaterial % aawl)
      watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
      aux3 = GeneralIntegral(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % zw0,&
           CurrentSolventMaterial % bzw(0:5),&
           CurrentSolventMaterial % bzwl)
      rhow = rhow * EXP(aux1 - aux2 + aux3)
      IF (rhow .NE. rhow) THEN
        PRINT *, "rhow"
        STOP
      END IF
    END IF
  END FUNCTION rhow
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhowT(CurrentSolventMaterial,rhow,T0,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhow,T0,Temperature
    !--------------------
    REAL(KIND=dp) :: alphaW
    !--------------------
    alphaW = (CurrentSolventMaterial % aw0) *&
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSolventMaterial % aaw(0:5),&
         CurrentSolventMaterial % aawl)
    rhowT = rhow * alphaW
  END FUNCTION rhowT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhowP(CurrentSolventMaterial,rhow,p0,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhow,p0,Pressure
    !--------------------
    REAL(KIND=dp) ::  kappaW
    !--------------------
    kappaW = (CurrentSolventMaterial % kw0) *&
         GeneralPolynomial(Pressure,p0,p0,&
         CurrentSolventMaterial % ckw(0:5),&
         CurrentSolventMaterial % ckwl)
    rhowP = rhow * kappaW
  END FUNCTION rhowP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhoi(CurrentSolventMaterial,T0,p0,Temperature,Pressure,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Temperature,Pressure
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2
    !----------------------
    rhoi = CurrentSolventMaterial % rhoi0
    IF (.NOT.ConstVal) THEN
      aux1 = GeneralIntegral(Pressure,p0,p0,&
           CurrentSolventMaterial % ki0 ,&
           CurrentSolventMaterial % cki(0:5),&
           CurrentSolventMaterial % ckil)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % ai0,&
           CurrentSolventMaterial % aai(0:5),&
           CurrentSolventMaterial % aail)
      rhoi = rhoi * EXP(aux1 - aux2)
    END IF
  END FUNCTION rhoi
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhoiT(CurrentSolventMaterial,rhoi,T0,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoi,T0,Temperature
    !------------------------
    REAL(KIND=dp) :: alphaI
    !------------------------
    alphaI = (CurrentSolventMaterial % ai0) *&
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSolventMaterial % aai(0:5),&
         CurrentSolventMaterial % aail) 
    rhoiT = rhoi * alphaI
  END FUNCTION rhoiT
    !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhoiP(CurrentSolventMaterial,rhoi,T0,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoi,T0,Temperature
    !------------------------
    REAL(KIND=dp):: kappaI
    !------------------------
    kappaI = (CurrentSolventMaterial % ai0) *&
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSolventMaterial % aai(0:5),&
         CurrentSolventMaterial % aail) 
    rhoiP = rhoi * kappaI
  END FUNCTION rhoiP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhoc(CurrentSoluteMaterial,T0,p0,Xi,Temperature,Pressure,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Xi,Temperature,Pressure,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, aux3, xc
    !----------------------
    rhoc = CurrentSoluteMaterial % rhoc0
    IF (.NOT.ConstVal) THEN
      xc = Salinity/Xi
      aux1 = GeneralIntegral(Pressure,p0,p0,&
           CurrentSoluteMaterial % kc0,&
           CurrentSoluteMaterial % ckc(0:5),&
           CurrentSoluteMaterial % ckcl)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSoluteMaterial % ac0,&
           CurrentSoluteMaterial % aac(0:5),&
           CurrentSoluteMaterial % aacl)
      aux3 = GeneralIntegral(xc,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % zc0,&
           CurrentSoluteMaterial % bzc(0:5),&
           CurrentSoluteMaterial % bzcl)
      rhoc = rhoc * EXP(aux1 - aux2 + aux3)
    END IF
  END FUNCTION rhoc
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocT(CurrentSoluteMaterial,rhoc,T0,Temperature)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc,T0,Temperature
    !-------------------------
    REAL(KIND=dp):: alphaC
    !-------------------------
    alphaC = (CurrentSoluteMaterial % ac0) * &
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSoluteMaterial % aac(0:5),&
         CurrentSoluteMaterial % aacl) 
    rhocT = rhoc * alphaC
  END FUNCTION rhocT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocP(CurrentSoluteMaterial,rhoc)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc
    !---------------
    rhocP = rhoc * (CurrentSoluteMaterial % kc0)
  END FUNCTION rhocP
    !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocYc(CurrentSoluteMaterial,rhoc,Xi,Salinity)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc, Xi, Salinity
    !---------------
    REAL(KIND=dp):: xc, zc
    !---------------
    xc = Salinity/Xi
    zc = (CurrentSoluteMaterial % zc0) * &
         GeneralPolynomial(xc,0.0_dp,1.0_dp,&
         CurrentSoluteMaterial % bzc(0:5),&
         CurrentSoluteMaterial % bzcl)
    rhocYc = rhoc * zc /Xi
  END FUNCTION rhocYc
 !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogw(rhow,rhoc,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow,rhoc,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc    
    !------------
    xc = Salinity/Xi
    rhogw = rhow + xc*(rhoc - rhow)
  END FUNCTION rhogw  
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cs(CurrentRockMaterial,RockMaterialID,T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux
    !----------------------
    cs = CurrentRockMaterial % cs0(RockMaterialID)
      IF (.NOT.ConstVal) &
         cs = cs * GeneralPolynomial(Temperature,T0,T0,&
         CurrentRockMaterial % acs(0:5,RockMaterialID),&
         CurrentRockMaterial % acsl(RockMaterialID))
  END FUNCTION cs
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cw(CurrentSolventMaterial,T0,Xi,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Xi,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, watercont
    !----------------------
    cw = CurrentSolventMaterial % cw0
    IF (.NOT.ConstVal) THEN
      watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
      aux1 = GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % acw(0:5),&
           CurrentSolventMaterial % acwl)
      aux2 = GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % bcw(0:5),&
           CurrentSolventMaterial % bcwl)
      cw = cw * aux1 * aux2
    END IF
  END FUNCTION cw
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION ci(CurrentSolventMaterial,&
       T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    REAL(KIND=dp) :: ci0
    REAL(KIND=dp), DIMENSION(0:5) :: aci
    INTEGER :: acil
    LOGICAL :: ConstVal
    !----------------------
    ci = CurrentSolventMaterial % ci0     
    IF (.NOT.ConstVal) THEN
      !PRINT *, "ci:",ci0,aci(0:5),acil
      ci = ci *&
           GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % aci(0:5),&
           CurrentSolventMaterial % acil)
    END IF
  END FUNCTION ci
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cc(CurrentSoluteMaterial,&
       T0,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2
    !----------------------
    cc = CurrentSoluteMaterial % cc0
    IF (.NOT.ConstVal) THEN
      aux1 = GeneralPolynomial(Temperature,T0,T0,&
           CurrentSoluteMaterial % acc(0:5),&
           CurrentSoluteMaterial % accl)
      aux2 = GeneralPolynomial(Salinity,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % bcc(0:5),&
           CurrentSoluteMaterial % bccl)
      cc = cc*aux1*aux2
    END IF
  END FUNCTION cc
  !---------------------------------------------------------------------------------------------
  ! latent heat of water  
  REAL (KIND=dp) FUNCTION hw(CurrentSolventMaterial,&
       T0,Xi,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Xi,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, watercont
    !----------------------
    hw = CurrentSolventMaterial % hw0 
    IF (.NOT.ConstVal) THEN
      watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
      aux1 = GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % bcw(0:5),&
           CurrentSolventMaterial % bcwl)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % cw0,&
           CurrentSolventMaterial % acw(0:5),&
           CurrentSolventMaterial % acwl)
      hw = hw + aux1 * aux2
    END IF
  END FUNCTION hw
  !---------------------------------------------------------------------------------------------
  ! latent heat of ice  
  REAL (KIND=dp) FUNCTION hi(CurrentSolventMaterial,&
       T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=DP) :: hi0,ci0
    REAL(KIND=DP), DIMENSION(0:5) :: aci
    INTEGER :: acil
    LOGICAL :: FirstTime = .TRUE.
    SAVE FirstTime,hi0,ci0,aci,acil
    !----------------------
    hi0 = CurrentSolventMaterial % hi0
    IF (ConstVal) THEN
      hi = hi0
    ELSE
      IF (FirstTime) THEN
        ci0 = CurrentSolventMaterial % ci0 
        aci(0:5) = CurrentSolventMaterial % aci(0:5)
        acil = CurrentSolventMaterial % acil
        FirstTime = .FALSE.
      END IF
      hi = hi0 + GeneralIntegral(Temperature,T0,T0,ci0,aci,acil)
    END IF
  END FUNCTION hi
  !---------------------------------------------------------------------------------------------
  ! General constituent thermal conductivity: kalpha0th and balpha have to be directly transferred
  FUNCTION GetKAlphaTh(kalpha0th,balpha,T0,Temperature)RESULT(kalphath)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: kalpha0th,balpha,T0,Temperature
    REAL(KIND=dp) :: kalphath
    !-------------------------
    kalphath = kalpha0th/( 1.0_dp + balpha*(Temperature - T0)/T0)
    !kalphath = kalpha0th
  END FUNCTION GetKAlphaTh
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTT(Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,hi,hw,&
       Porosity,Salinity)RESULT(CGTT)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,&
         hi,hw,Porosity,Salinity
    REAL(KIND=dp) :: CGTT
    !-------------------------
    REAL(KIND=dp) :: xc
    !-------------------------
    xc = Salinity/Xi
    CGTT = (1.0_dp - Porosity)*rhos*cs &
         + (Xi - Salinity) * Porosity * rhow * cw & ! mind xc * Xi = Salinity
         + Salinity * Porosity * rhoc * cc & ! mind xc * Xi = Salinity
         + (1.0_dp - Xi)*Porosity*rhoi*ci &
         + rhoi*(hw - hi)*Porosity*XiT
  END FUNCTION GetCGTT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTp(rhoi,hi,hw,XiP,Porosity)RESULT(CGTp)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoi,hi,hw,XiP,Porosity
    REAL(KIND=dp) :: CGTp
    !-------------------------
    CGTp = Porosity*rhoi*(hw - hi)*XiP
  END FUNCTION GetCGTp  
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTyc(rhoi,hi,hw,XiYc,Porosity)RESULT(CGTyc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoi,hi,hw,XiYc,Porosity
    REAL(KIND=dp) :: CGTyc
    !-------------------------
    CGTyc = Porosity*rhoi*(hw - hi)*XiYc
  END FUNCTION GetCGTyc
  !---------------------------------------------------------------------------------------------
  ! functions specific to groundwater flow
  !---------------------------------------------------------------------------------------------
  FUNCTION GetJgwD(Kgwpp,KgwpT,Kgw,gradp,gradT,Gravity,rhogw,DIM) RESULT(JgwD)
    IMPLICIT NONE
    REAL (KIND=dp), INTENT(IN) :: Kgwpp(3,3),KgwpT(3,3),Kgw(3,3),gradp(3),gradT(3),Gravity(3),&
         rhogw
    REAL (KIND=dp)  :: JgwD(3)
    !-------------------------
    INTEGER, INTENT(IN) :: DIM
    INTEGER :: i
    REAL (KIND=dp) :: fluxp(3),fluxT(3),fluxg(3)
    !-------------------------
    !PRINT *,"Flux",fluxpAtIP,"gradp",gradpAtIP(1:DIM),"K",KgwppAtIP(1:DIM,1:DIM)
    DO i=1,dim
      fluxp(i) = -1.0_dp * SUM(Kgwpp(i,1:DIM)*gradp(1:DIM))
      fluxT(i) = -1.0_dp * SUM(KgwpT(i,1:DIM)*gradT(1:DIM))
      fluxg(i) =   rhogw * SUM(Kgw(i,1:DIM)*Gravity(1:DIM))
    END DO
    JgwD(1:3) = 0.0
    JgwD(1:DIM) = fluxp(1:DIM) + fluxT(1:DIM) + fluxg(1:DIM)    
  END FUNCTION GetJgwD
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKGTT(ksth,kwth,kith,kcth,Xi,&
       Salinity,Porosity,meanfactor)RESULT(KGTT) ! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ksth,kwth,kith,kcth,Xi,&
         Salinity,Porosity,meanfactor
    !-------------------------
    REAL(KIND=dp) :: KGTT(3,3), KGaTT, KghTT, unittensor(3,3),xc
    !-------------------------
    xc = Salinity/Xi
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    KGhTT = 1.0_dp/((1.0_dp - Porosity)/ksth + (1.0_dp - xc)*Xi*Porosity/kwth &
         + xc*Porosity/kcth + (1.0_dp - Xi)*Porosity/kith)
    KGaTT = (1.0_dp - Porosity)*ksth + (1.0_dp - xc)*Xi*Porosity*kwth &
         + xc*Porosity*kcth + (1.0_dp - Xi)*Porosity*kith
    KGTT = unittensor*((1.0_dp - meanfactor)*KGhTT + meanfactor * KGaTT)
  END FUNCTION GetKGTT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwTT(rhow,rhoc,cw,cc,Xi,Salinity)RESULT(CgwTT)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow,rhoc,cw,cc,Xi,Salinity
    REAL(KIND=dp) :: CgwTT
    !-------------------------
    REAL(KIND=dp) :: xc
    !-------------------------
    xc = Salinity/Xi
    CgwTT = (1.0_dp - xc)*rhow*cw + xc*rhoc*cc
  END FUNCTION GetCgwTT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
       Xi,T0,Salinity,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi,T0,Salinity,Temperature
    LOGICAL :: ConstVal
    !-------------------------
    REAL(KIND=dp) :: nu1, nu2, xc
    !-------------------------
    mugw = CurrentSolventMaterial % muw0
    IF (.NOT.ConstVal) THEN
      xc = Salinity/Xi
      nu1 = (CurrentSolventMaterial % nu10) *&
           GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % anw(0:5),&
           CurrentSolventMaterial % anwl)
!      PRINT *,"mugw:anw,anwl,nu1", anw(0:5),anwl, nu1
!      PRINT *,"mugw:bnc", bnc(0:5),bncl
!      PRINT *,"bnc=",bnc(0:5),"bncl=",bncl,"xc=",xc
!      STOP
      nu2 = (CurrentSoluteMaterial % nu20) *&
           GeneralPolynomial(xc,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % bnc(0:5),&
           CurrentSoluteMaterial % bncl)
!      PRINT *,"mugw:", nu1,nu2,Temperature, xc
      mugw = mugw * EXP(nu1 * (Temperature - T0) + nu2 * (xc - 0.0))
    END IF    
  END FUNCTION mugw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
       mugw,Xi,MinKgw)RESULT(Kgw)
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: Xi,MinKgw,mugw
    REAL(KIND=dp) :: Kgw(3,3)
    !--------------------------
    REAL(KIND=dp) :: muw0,rhow0,qexp,Kgwh0(3,3),kGpe(3,3),kG0pe(3,3),factor
    REAL(KIND=dp), PARAMETER :: gval=9.81_dp !hard coded, so match Kgwh0 with this value
    INTEGER :: I, J
    !-------------------------
    IF (mugw <= 0.0_dp) &
         CALL FATAL("Permafrost(GetKgw)","Unphysical viscosity detected")
    muw0 = CurrentSolventMaterial % muw0
    rhow0 = CurrentSolventMaterial % rhow0
    qexp = CurrentRockMaterial % qexp(RockMaterialID)
    Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
    ! transformation factor from hydr. conductivity to hydr. conductivity tensor
    factor = (muw0/mugw)*(Xi**qexp)/(rhow0*gval)
    !PRINT *,"Kgw:",muw0,mugw,rhow0,Kgwh0,Xi,factor
    Kgw = 0.0_dp
    DO I=1,3
      DO J=1,3
        Kgw(i,j) = Kgwh0(i,j)*factor
      END DO
    END DO
    DO I=1,3
      Kgw(i,i) = MAX(Kgw(i,i),MinKgw)
    END DO
  END FUNCTION GetKgw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgwpT(fw,XiT,Kgw)RESULT(KgwpT) ! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: fw,XiT,Kgw(3,3)
    REAL(KIND=dp) :: KgwpT(3,3)
    !-------------------------
    KgwpT(1:3,1:3) = fw*XiT*Kgw(1:3,1:3)
  END FUNCTION GetKgwpT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgwpp(fw,XiP,Kgw)RESULT(Kgwpp)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: fw,XiP,Kgw(3,3)
    REAL(KIND=dp) :: Kgwpp(3,3)
    !-------------------------
    Kgwpp(1:3,1:3) = (1.0_dp + fw*XiP)*Kgw(1:3,1:3)
  END FUNCTION GetKgwpp
  !---------------------------------------------------------------------------------------------
  ! functions specific to solute transport
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKc(CurrentRockMaterial,RockMaterialID,Dm,Xi,JgwD,Porosity)RESULT(Kc) 
    IMPLICIT NONE
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    REAL(KIND=dp), INTENT(IN) :: Dm,Xi,JgwD(3),Porosity
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: alphaL,alphaT,Kc(3,3), unittensor(3,3), aux, eL(3),absJgwD
    INTEGER :: I,J
    !-------------------------
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    IF (Porosity <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero Porosity detected")
    IF (Xi <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero water content detected")
    Kc =  Dm * unittensor
    absJgwD = SQRT(SUM(JgwD(1:3) * JgwD(1:3)))
    IF (absJgwD > 0.0_dp) THEN
      alphaL = CurrentRockMaterial % alphaL(RockMaterialID)
      alphaT = CurrentRockMaterial % alphaT(RockMaterialID)
      eL = JgwD/absJgwD
      aux = absJgwD/(Porosity * Xi)   
      DO I=1,3
        DO J=1,3
          Kc(I,J) = Kc(I,J)  &
               + aux*((alphaL - alphaT)*eL(I)*eL(J)  + alphaT * unittensor(I,J))
        END DO
      END DO
    END IF
  END FUNCTION GetKc  
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION Dm(CurrentSoluteMaterial,N0,GasConstant,rhoc,mugw,Temperature)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: N0,GasConstant,rhoc,mugw,Temperature
    !-------------------------
    REAL(KIND=dp), PARAMETER :: bconst = 3.0 * PI
    REAL(KIND=dp) :: Mc, lcbar
    !-------------------------
    Mc = CurrentSoluteMaterial % Mc
    lcbar = (Mc/(rhoc * N0))**(1.0_dp/3.0_dp)
    Dm = GasConstant * Temperature / (bconst * mugw * lcbar * N0)
  END FUNCTION Dm
  !---------------------------------------------------------------------------------------------
  FUNCTION GetR(CurrentSoluteMaterial,GasConstant,rhoc,Xi,Temperature,Salinity) RESULT(r12)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: GasConstant,rhoc,Xi,Salinity,Temperature
    REAL(KIND=dp) :: r12(2)
    REAL(KIND=dp) :: d1, d2, Mc, aux
    !-------------------------
    d1 = CurrentSoluteMaterial % d1
    d2 = CurrentSoluteMaterial % d2
    Mc = CurrentSoluteMaterial % Mc
    aux = Salinity/(Xi - Salinity)
    r12(1) = Mc*(1.0_dp - Salinity/Xi)/(rhoc * GasConstant * Temperature)
    r12(2) = d1 + (d1 + d2)*aux + d2*aux*aux    
  END FUNCTION GetR
  !---------------------------------------------------------------------------------------------
  FUNCTION  GetKcYcYc(Kc,r12) RESULT(KcYcYc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Kc(3,3),r12(2)
    REAL(KIND=dp) :: KcYcYc(3,3)
    !-------------------------
    KcYcYc(1:3,1:3) = r12(2) * Kc(1:3,1:3) 
  END FUNCTION GetKcYcYc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetFc(rhoc,rhow,Gravity,r12,XiT,XiP,Xi,gradP,gradT) RESULT(fc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoc,rhow,Gravity(3),r12(2),XiT,XiP,Xi,gradP(3),gradT(3)
    REAL(KIND=dp) :: fc(3)
    !-------------------------
    fc(1:3) = r12(1)*(rhoc - rhow)*Gravity(1:3) + r12(2)*(XiT*gradT(1:3) + XiP*gradP(1:3))/Xi
  END FUNCTION GetFc
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcT(rhocT,Porosity,Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhocT,Porosity, Salinity
    CcYcT = Porosity*Salinity*rhocT
  END FUNCTION CcYcT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcP(rhocP,Porosity, Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhocP,Porosity, Salinity
    !-------------------------
    CcYcP = Porosity*Salinity*rhocp
  END FUNCTION CcYcP
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcYc(rhoc,rhocYc,Porosity, Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoc,rhocYc,Porosity, Salinity
    !-------------------------
    CcYcYc = Porosity*(rhoc + Salinity*rhocYc)
  END FUNCTION CcYcYc

  !---------------------------------------------------------------------------------------------
END MODULE PermafrostMaterials
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
! Solvers for problems
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> Solver for groundwater flow of the enhanced permafrost model
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlow( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: TemperatureVar,PorosityVar,SalinityVar!,GWfluxVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, Active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=22
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:)!,GWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:)!,GWflux(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
       NodalPressure(:) !NodalGWflux(:,:),
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.,ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, VarName, PhaseChangeModel,&
       ElementRockMaterialName 

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure,ElementWiseRockMaterial
  !------------------------------------------------------------------------------
  CALL DefaultStart()

  Params => GetSolverParams()

  Pressure => Solver % Variable % Values
  PressurePerm => Solver % Variable % Perm
  VarName = Solver % Variable % Name

  IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
    N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    IF (AllocationsDone) &
         DEALLOCATE(NodalPressure,NodalPorosity,NodalTemperature,NodalSalinity)!,NodalGWflux)
    ALLOCATE(NodalPressure(N),NodalPorosity(N),NodalTemperature(N),&
         NodalSalinity(N),STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL FATAL(SolverName,"Allocation error")
    END IF
    AllocationsDone = .TRUE.
  END IF

  maxiter = ListGetInteger( Params, &
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! find variables for dependencies
  !--------------------------------
  TemperatureName = ListGetString(Params, &
       'Temperature Variable', Found )
  IF (.NOT.Found) THEN
    CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
    WRITE(TemperatureName,'(A)') 'Temperature'
  ELSE
    WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
    CALL INFO(SolverName,Message,Level=9)
  END IF
  TemperatureVar => VariableGet(Solver % Mesh % Variables,TemperatureName)
  IF (.NOT.ASSOCIATED(TemperatureVar)) THEN
    CALL FATAL(SolverName,'Temperature Variable not found.')
  ELSE
    Temperature => TemperatureVar % Values
    TemperaturePerm => TemperatureVar % Perm
  END IF

  PorosityName = ListGetString(Params, &
       'Porosity Variable', Found )
  IF (.NOT.Found) THEN
    CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
    WRITE(PorosityName,'(A)') 'Porosity'
  ELSE
    WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
    CALL INFO(SolverName,Message,Level=9)
  END IF
  ConstantPorosity= GetLogical(Params,'Constant Porosity', Found)
  IF ((.NOT.Found) .OR. (.NOT.ConstantPorosity)) THEN
    PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
    IF (.NOT.ASSOCIATED(PorosityVar)) THEN
      CALL FATAL(SolverName,'Porosity Variable not found')
    ELSE
      Porosity => PorosityVar % Values
      PorosityPerm => PorosityVar % Perm
    END IF
  ELSE
    NULLIFY(PorosityVar)
  END IF

  SalinityName = ListGetString(Params, &
       'Salinity Variable', Found )
  IF (.NOT.Found) THEN
    CALL WARN(SolverName," 'Salinity Variable' not found. Using default 'Salinity' ")
    WRITE(SalinityName,'(A)') 'Salinity'
  ELSE
    WRITE(Message,'(A,A)') "'Salinity Variable' found and set to: ", SalinityName
    CALL INFO(SolverName,Message,Level=9)
  END IF
  SalinityVar => VariableGet(Solver % Mesh % Variables,SalinityName)
  IF (.NOT.ASSOCIATED(SalinityVar)) THEN
    CALL WARN(SolverName,'Salinity Variable not found. Switching Salinity off')
    NoSalinity = .TRUE.
  ELSE
    Salinity => SalinityVar % Values
    SalinityPerm => SalinityVar % Perm
  END IF


  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL DefaultInitialize()
    !------------------------------------------------------------------------------
    Active = Solver % NumberOfActiveElements
    DO t=1,Active
      Element => GetActiveElement(t)
      IF (.NOT.ASSOCIATED(Element)) CYCLE
      ! cycle halo elements
      !-------------------
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      Material => GetMaterial(Element)
      
      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF
      
      IF (FirstTime) THEN
        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Active)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF      
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
      END IF
      IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
        CALL FATAL(SolverName,Message)
      END IF
      N  = GetElementNOFNodes()
      ND = GetElementNOFDOFs()
      NB = GetElementNOFBDOFs()

      ! Nodal variable dependencies
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
      meanfactor = GetConstReal(Material,"Conductivity Aritmetic Mean Weight",Found)
      IF (.NOT.Found) THEN
        meanfactor = 1.0_dp
      END IF
      IF (ConstantPorosity) THEN
        NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
        IF (.NOT.Found) THEN
          WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
          CALL FATAL(SolverName,Message)
        END IF
      ELSE
        NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
      END IF
      NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
      IF (NoSalinity) THEN
        NodalSalinity(1:N) = 0.0_dp
      ELSE
        NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
      END IF
      
      ! compose element-wise contributions to matrix and R.H.S
      CALL LocalMatrixDarcy( Model, Element, t, N, ND+NB, Active, NodalPressure, NodalTemperature, &
           NodalPorosity, NodalSalinity, CurrentRockMaterial,CurrentSoluteMaterial,&
           CurrentSolventMaterial,PhaseChangeModel,ElementWiseRockMaterial)
    END DO
    CALL DefaultFinishBulkAssembly()
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCDarcy(  Element, n, nd+nb )
      END IF
    END DO
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    !Solve the system:
    !--------------------
    Norm = DefaultSolve()
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT


  END DO

  CALL DefaultFinish()

CONTAINS
  ! PermafrostGroundWaterFlow : Assembly of the matrix entries arising from the bulk elements 

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixDarcy( Model, Element, ElementID, n, nd, NoElements, NodalPressure, NodalTemperature, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, CurrentSoluteMaterial,&
       CurrentSolventMaterial, PhaseChangeModel, ElementWiseRockMaterial )
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    LOGICAL :: ElementWiseRockMaterial
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),meanfactor,MinKgw,gradTAtIP(3),&
         gradPAtIP(3),fluxTAtIP(3),fluxgAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP,bijAtIP(2,2),bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: ks0th,e1,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,RadGen,acs(0:5),&
         as0,aas(0:5),ks0,cks(0:5)  ! stuff comming from RockMaterial
    INTEGER :: acsl,aasl,cksl       ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps,Gravity(3)! real constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp) , POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow', &
         FunctionName='Permafrost (LocalMatrixDarcy)'

    SAVE Nodes, ConstantsRead, ConstVal, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)      
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      LOAD(1:n) = GetReal( BodyForce,'Groundwater source', Found )   
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)


    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:N) * LOAD(1:N) )

      ! from coordinate system
      Weight = IP % s(t) * DetJ

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N) )
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N) )
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N) )
      
      !Materialproperties needed at IP
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW
      !PRINT *,"Densities",rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
        bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
        gwaAtIP = gwa(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP) ! NEW
        gwaTAtIP =  gwaT(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)!
        gwapAtIP = 1.0_dp/rhowAtIP
        giaAtIP = gia(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP) !NEW        
        giaTAtIP = giaT(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP) ! NEW
        giapAtIP = 1.0_dp/rhoiAtIP
        deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW
        D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
        D2AtIP= D2(deltaInElement,bijAtIP)! NEW
        B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
             GasConstant,bijAtIP,TemperatureAtIP)
        B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
        XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,gwapAtIP,giapAtIP,Xi0Tilde,&
             deltaInElement,GasConstant,TemperatureAtIP)         
      END SELECT
      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      
      rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP,MinKgw)
      !PRINT *,"KgwppAt(1,1)",KgwAtIP(1,1)
      fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
           Xi0Tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
      !PRINT *,"mugw,fw", mugwAtIP, fwAtIP
      KgwpTAtIP = 0.0_dp
      KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
      KgwppAtIP = 0.0_dp
      KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
      !PRINT *,"KgwppAtIP(1,1)",KgwppAtIP(1,1)
      !PRINT *, "KgwpTAtIP=",KgwpTAtIP(i,1:DIM),"rhowAtIP=",rhowAtIP,&
      !         "SalinityAtIP=",SalinityAtIP,"rhocAtIP=",rhocAtIP

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      DO p=1,nd
        DO q=1,nd
          StiffPQ = 0.0

          ! diffusion term ( grad(u),grad(v))
          ! div(J_gwp) = d(Kgwpp_i,j dp/dx_j)/dx_i
          ! -----------------------------------
          DO i=1,DIM
            DO j=1,DIM
              StiffPQ = StiffPQ +  rhogwAtIP * KgwppAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)              
            END DO
          END DO
          !PRINT *,"KgwppAtIP,KgwAtIP",KgwppAtIP(1,1),KgwAtIP(1,1)
          !StiffPQ = 1.0d-04 *  SUM(dBasisdx(p,1:DIM)* dBasisdx(q,1:DIM))
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
        END DO
      END DO
      !IF (XiAtIP > 0.0_dp) THEN
      ! body forces
      DO i=1,DIM
        gradTAtIP(i) =  SUM(NodalTemperature(1:n)*dBasisdx(1:n,i))
      END DO
      DO i=1,DIM
        fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
        fluxgAtIP(i) = rhogwAtIP * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))   !!! NEW
        IF ((fluxgAtIP(i) .NE. fluxgAtIP(i)) .OR. (fluxTAtIP(i) .NE. fluxTAtIP(i))) THEN
          PRINT *, "NaN in r.h.s. of Darcy fluxes"
          PRINT *, "flux(",i,")= Jgwg",fluxgAtIP(i),"+ JgwpT", fluxTAtIP(i)
          PRINT *, "KgwpTAtIP=",KgwpTAtIP(i,1:DIM),"rhowAtIP=",rhowAtIP,&
               "SalinityAtIP=",SalinityAtIP,"rhocAtIP=",rhocAtIP
          STOP
        END IF
      END DO
      DO p=1,nd     
        FORCE(p) = FORCE(p) - Weight * SUM(fluxTAtIP(1:DIM)*dBasisdx(p,1:DIM)) !!! NEW
        FORCE(p) = FORCE(p) + Weight * rhogwAtIP * SUM(fluxgAtIP(1:DIM)*dBasisdx(p,1:DIM))        
      END DO
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
      !END IF
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixDarcy
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCDarcy( Element, n, nd)

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Pressure(n), F,Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,PressureAtIP
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    LOGICAL :: Stat,Found,FluxCondition,WeakDirichletCond
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlow (LocalMatrixBCDarcy)'

    SAVE Nodes,DIM
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BoundaryCondition,'Groundwater Flux', FluxCondition )
    ! Check, whether we have a weakly imposed Dirichlet condition
    Pressure(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)

    ! Numerical integration:
    !-----------------------
    IF (FluxCondition .OR. WeakDirichletCond) THEN ! spare us, if natural BC
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ
        ! Given flux:
        ! -----------
        IF (Fluxcondition) THEN

          F = SUM(Basis(1:n)*flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
          ! Given pressure, weakly imposed
          !----------------------------------------------------------------------
        ELSE IF (WeakDirichletCond) THEN
          PressureAtIP = SUM(Pressure(1:n)*Basis(1:n))
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * PressureAtIP * Basis(1:nd)
        END IF
      END DO
      CALL DefaultUpdateEquations(STIFF,FORCE)
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCDarcy
  !------------------------------------------------------------------------------

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------  
END SUBROUTINE PermafrostGroundwaterFlow

!-----------------------------------------------------------------------------
!> Solver for groundwater flow of the enhanced permafrost model
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlux( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName, PotName
  INTEGER :: i,j,k,dim,DOFs,firstmag
  LOGICAL :: GotIt
  REAL(KIND=dp) :: Unorm, Totnorm, val
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER :: SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: FluxSol
  TYPE FieldTable_t
     REAL(KIND=dp), POINTER :: Values(:)
     INTEGER, POINTER :: Perm(:)
  END TYPE FieldTable_t
  TYPE(FieldTable_t) :: Fields(3)

  TYPE(Variable_t), POINTER :: PressureVar,TemperatureVar,PorosityVar,SalinityVar
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:)
  INTEGER :: NumberOfRockRecords
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:)
  LOGICAL :: ConstantPorosity, NoSalinity, UnfoundFatal=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, PressureName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux"

  SAVE SaveRHS

  CALL Info( SolverName, '-------------------------------------',Level=3 )
  CALL Info( SolverName, 'Computing the groundwater flux       ',Level=3 )
  CALL Info( SolverName, '-------------------------------------',Level=3 )

  dim = CoordinateSystemDimension()
  !------------------------------------------------------------------------------
  !  Check what needs to be computed
  !------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN

  SolverParams => GetSolverParams()
  Dofs = Dim

  ! Read Variables
  CALL AssignVarsGWFlux()

  !-------------------------------------------------------------------------------
  ! If only one component is used use the scalar equation, otherwise use an
  ! auxiliary variable to store all the dimensions
  !-------------------------------------------------------------------------------
  Varname = TRIM('Groundwater')

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 1',UnFoundFatal=UnFoundFatal )
  Fields(1) % Values => FluxSol % Values
  Fields(1) % Perm => FluxSol % Perm

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 2',UnFoundFatal=UnFoundFatal )
  Fields(2) % Values => FluxSol % Values
  Fields(2) % Perm => FluxSol % Perm

  IF( dim == 3 ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 3',UnFoundFatal=UnFoundFatal )
    Fields(3) % Values => FluxSol % Values
    Fields(3) % Perm => FluxSol % Perm
  END IF
  at0 = RealTime()

  CALL DefaultInitialize()

  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),DOFs))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS

  CALL BulkAssembly()
  CALL DefaultFinishAssembly()

  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( SolverName, Message, Level=5 )
  !        
  !------------------------------------------------------------------------------     


  TotNorm = 0.0_dp
  DO i=1,Dofs
    WRITE(Message,'(A,I1,A,I1)') "Working on DOF ",i," out of ",Dofs
    CALL INFO(SolverName,Message,Level=3)
    Solver % Matrix % RHS => ForceVector(:,i)
    UNorm = DefaultSolve()
    WRITE( Message, * ) 'Norm of DOF: ',i,'=',UNorm ** 2.0_dp
    CALL INFO(SolverName,Message,Level=3)
    TotNorm = TotNorm + Unorm ** 2.0_dp
    Fields(i) % Values = Solver % Variable % Values
    !Fields(i) % Values = 1.0_dp * i
  END DO

  DEALLOCATE( ForceVector )  
  Solver % Matrix % RHS => SaveRHS
  TotNorm = SQRT(TotNorm)
  Solver % Variable % Norm = Totnorm


  !------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( SolverName, Message, Level=4 )

  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( SolverName, Message, Level=4 )

  CALL Info( SolverName, 'All done',Level=4 )
  CALL Info( SolverName, '-------------------------------------',Level=6 )



CONTAINS


  !------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: elem,t,i,j,k,p,q,n,nd, DIM,Rank, RockMaterialID, Active
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,detJ,GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    LOGICAL :: Found, ConstVal=.FALSE.,ConstantsRead=.FALSE., FirstTime=.TRUE., ElementWiseRockMaterial
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: GasConstant, N0, meanfactor,DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         JgwDAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, &
         bijAtIP(2,2), bijYcAtIP(2,2),gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuf
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP  
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:), NodalSalinity(:), NodalPressure(:),NodalPorosity(:)
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlux (BulkAssembly)'
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel,ElementRockMaterialName
    ! -------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, meanfactor,GasConstant, N0,DeltaT, eps,T0, p0,Gravity,&
         CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
         FirstTime,ElementWiseRockMaterial
    ! -------------------------------------------------------------
          
    n = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )

    ALLOCATE( STIFF(n,n), FORCE(dofs,n) )
    ALLOCATE( Basis(n), dBasisdx(n,3) )
    ALLOCATE( NodalPressure(N),NodalPorosity(N),NodalTemperature(N),NodalSalinity(N) )

    Active = Solver % NumberOFActiveElements
    DO elem = 1,Active

      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      Material => GetMaterial(Element)
      
      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF
      IF (FirstTime) THEN
        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Active)
          PRINT *, "NumberOfRockRecords", NumberOfRockRecords
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF      
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )        
      END IF

      IF(.NOT.ConstantsRead) THEN
        ConstantsRead = &
             ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      END IF

      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

      ! Nodal Variable Values
      NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
      IF (ConstantPorosity) THEN
        NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
        IF (.NOT.Found) THEN
          WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
          CALL FATAL(SolverName,Message)
        END IF
      ELSE
        NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
      END IF
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))      
      IF (NoSalinity) THEN
        NodalSalinity(1:N) = 0.0_dp
      ELSE
        NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
      END IF

      ! Get stuff from SIF Material section
      Material => GetMaterial(Element)     
      meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
      IF (.NOT.Found) THEN
        CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
        meanfactor = 1.0_dp
      END IF
      MinKgw = GetConstReal( Material, &
           'Hydraulic Conductivity Limit', Found)
      IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
           MinKgw = 1.0D-14
      
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      IF (ElementWiseRockMaterial) THEN
        RockMaterialID = elem  ! each element has it's own set of parameters
      ELSE
        RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
      END IF
      
      ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
      IF (.NOT.Found) THEN
        ConstVal = .FALSE.
      ELSE
        IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
      END IF

      deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
             IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

        ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
        TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
        PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
        PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
        SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))        
        gradTAtIP = 0.0_dp
        gradPAtIP = 0.0_dp
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradpAtIP(i) =  SUM(NodalPressure(1:N)*dBasisdx(1:N,i))
        END DO
        
        !Materialproperties needed at IP
        rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('Anderson')
          XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
        CASE DEFAULT ! Hartikainen model
          Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
          bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
          bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
!          PRINT *, "1_Flux:",SalinityAtIP,Xi0tilde,bijAtIP(1:2,1:2)
          gwaAtIP = gwa(CurrentSolventMaterial,&
               p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP) ! NEW
          gwaTAtIP =  gwaT(CurrentSolventMaterial,&
               p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)!
          gwapAtIP = 1.0_dp/rhowAtIP
          giaAtIP = gia(CurrentSolventMaterial,&
               p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP) !NEW
          giaTAtIP = giaT(CurrentSolventMaterial,&
               p0,T0,rhoiAtIP,TemperatureAtIP)
          giapAtIP = 1.0_dp/rhoiAtIP
          deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW
          D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
          D2AtIP= D2(deltaInElement,bijAtIP)! NEW
          !PRINT *,"Flux:",deltaInElement,bijAtIP,SalinityAtIP
          B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
               GasConstant,bijAtIP,TemperatureAtIP)
          B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
          !IF (B2AtIP .NE. B2AtIP) THEN
          !  PRINT *,"Flux:deltaInElement,deltaGAtIP,",deltaInElement,deltaGAtIP,GasConstant,bijAtIP,Xi0tilde,SalinityAtIP
          !  STOP
          !END IF
!          PRINT *,"2_Flux:",D1AtIP,D2AtIP,B1AtIP,B2AtIP,deltaInElement
          XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
          XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
          XiPAtIP= XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,gwapAtIP,giapAtIP,Xi0Tilde,&
             deltaInElement,GasConstant,TemperatureAtIP)
        END SELECT
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = 0.0_dp
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = 0.0_dp
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        !PRINT *,"Flux:", rhogwAtIP,rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP
        Weight = IntegStuff % s(t) * detJ

        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
          END DO
        END DO

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM)
        DO i=1,dim
          FORCE(i,1:nd) = FORCE(i,1:nd) + Weight *  JgwDAtIP(i) * Basis(1:nd)
        END DO
      END DO

      !------------------------------------------------------------------------------
      !      Update global matrices from local matrices 
      !------------------------------------------------------------------------------
      Solver % Matrix % Rhs => SaveRhs
      CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      !      END IF

      DO i=1,Dofs
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO

    END DO

    DEALLOCATE( STIFF, FORCE, Basis, dBasisdx,&
         NodalPressure,NodalPorosity,NodalTemperature,NodalSalinity )

    !------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
  !------------------------------------------------------------------------------


  SUBROUTINE AssignVarsGWFlux()
    ! find variables for dependencies
    !--------------------------------
    TemperatureName = ListGetString(SolverParams, &
         'Temperature Variable', GotIt )
    IF (.NOT.GotIt) THEN
      CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
      WRITE(TemperatureName,'(A)') 'Temperature'
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    TemperatureVar => VariableGet(Solver % Mesh % Variables,TemperatureName)
    IF (.NOT.ASSOCIATED(TemperatureVar)) THEN
      CALL FATAL(SolverName,'Temperature Variable not found.')
    ELSE
      Temperature => TemperatureVar % Values
      TemperaturePerm => TemperatureVar % Perm
    END IF

    PressureName = ListGetString(SolverParams, &
         'Pressure Variable', GotIt )
    IF (.NOT.GotIt) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      PRINT *,PressureName
      CALL FATAL(SolverName,'Pressure Variable not found.')
    ELSE
      Pressure => PressureVar % Values
      PressurePerm => PressureVar % Perm
    END IF

    PorosityName = ListGetString(SolverParams, &
         'Porosity Variable', GotIt )
    IF (.NOT.GotIt) THEN
      CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
      WRITE(PorosityName,'(A)') 'Porosity'
    ELSE
      WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    ConstantPorosity= GetLogical(SolverParams,'Constant Porosity', GotIt)
    IF ((.NOT.GotIt) .OR. (.NOT.ConstantPorosity)) THEN
      PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
      IF (.NOT.ASSOCIATED(PorosityVar)) THEN
        CALL FATAL(SolverName,'Porosity Variable not found')
      ELSE
        Porosity => PorosityVar % Values
        PorosityPerm => PorosityVar % Perm
      END IF
    ELSE
      NULLIFY(PorosityVar)
    END IF

    SalinityName = ListGetString(SolverParams, &
         'Salinity Variable', GotIt )
    IF (.NOT.GotIt) THEN
      CALL WARN(SolverName," 'Salinity Variable' not found. Using default 'Salinity' ")
      WRITE(SalinityName,'(A)') 'Salinity'
    ELSE
      WRITE(Message,'(A,A)') "'Salinity Variable' found and set to: ", SalinityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    SalinityVar => VariableGet(Solver % Mesh % Variables,SalinityName)
    IF (.NOT.ASSOCIATED(SalinityVar)) THEN
      CALL WARN(SolverName,'Salinity Variable not found. Switching Salinity off')
      NoSalinity = .TRUE.
    ELSE
      NoSalinity = .FALSE.
      Salinity => SalinityVar % Values
      SalinityPerm => SalinityVar % Perm
    END IF
  END SUBROUTINE AssignVarsGWFlux
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostGroundwaterFlux
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Initialization for the primary solver, PermafrostGroundwaterFlux. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlux_Init( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: DT
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: dim
  CHARACTER(LEN=MAX_NAME_LEN) :: EqName, VarName, FluxName, GradName
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux_Init", &
       FluxVariableName
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  dim = CoordinateSystemDimension()

  IF( dim < 2 .OR. dim > 3 ) THEN
    CALL Fatal('PermafrostGroundwaterFlux_init','Flux computation makes sense only in 2D and 3D')
  END IF


  VarName = TRIM('GroundWater')
  !VarName = TRIM('GW')

  IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    EqName = ListGetString( SolverParams,'Equation')
    CALL ListAddString( SolverParams, 'Variable','-nooutput '//TRIM(EqName)//'_temp' )
  END IF

  FluxName = TRIM(VarName)//' Flux'
  CALL Info('PermafrostGroundwaterFlux_init','Saving flux to: '//TRIM(FluxName), Level=3) 
  IF(dim == 2) THEN
    FluxVariableName=TRIM(FluxName)//'['//TRIM(FluxName)//':2]'
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':2]')
  ELSE IF(dim == 3) THEN
    FluxVariableName=TRIM(FluxName)//'['//TRIM(FluxName)//':3]'
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':3]')
  ELSE
    CALL FATAL('PermafrostGroundwaterFlux_init','Wrong dimension of problem')
  END IF
  CALL ListAddString( SolverParams,&
       NextFreeKeyword('Exported Variable',SolverParams),&
       FluxVariableName)
  WRITE(Message,*) 'Added ',TRIM(FluxVariableName),' as variable'
  CALL INFO('PermafrostGroundwaterFlux_init',Message,Level=3)
  IF( GetLogical( SolverParams,'Calculate Flux Abs',GotIt) ) THEN
    FluxName = TRIM(VarName)//' Flux_abs'
    CALL Info('PermafrostGroundwaterFlux_init','Saving flux abs to: '//FluxName) 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),TRIM(FluxName))
  END IF

  CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

  CALL ListAddLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)

  ! Add linear system defaults: cg+diagonal
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
       CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
       CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
       CALL ListAddString(SolverParams,'Linear System Preconditioning','diagonal')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
       CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-10_dp)

  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostGroundwaterFlux_Init
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> solute transport equation for enhanced permafrost model
!-----------------------------------------------------------------------------
SUBROUTINE PermafrostSoluteTransport( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,TemperatureVar,GWfluxVar1,GWfluxVar2,GWfluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostSoluteTransportEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, VarName, TemperatureName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,ElementWiseRockMaterial!,NodalGWflux,NoGWflux
  !------------------------------------------------------------------------------

  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()

  CALL AssignVarsSolute()

  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      Material => GetMaterial()
      IF (FirstTime) THEN

        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Active)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      END IF

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF

      CALL ReadVarsSolute(N)

      CALL LocalMatrixSolute(  Element, t, Active, n, nd+nb,&
           NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
           NodalGWflux,GivenGWflux,&
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()

    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCSolute(  Element, n, nd+nb )
        !PRINT *,t,"of",Active,":",n, nb
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()

CONTAINS
  
  ! find needed variables and assign pointers
  SUBROUTINE AssignVarsSolute()

    IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
      DIM = CoordinateSystemDimension()
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (AllocationsDone) &
           DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,&
           NodalSalinity,NodalGWflux)
      ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),&
           NodalSalinity(N),NodalGWflux(3,N),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        AllocationsDone = .TRUE.
        CALL INFO(SolverName,"Allocations Done",Level=1)
      END IF

    END IF
    Salinity => Solver % Variable % Values
    SalinityPerm => Solver % Variable % Perm

    PressureName = ListGetString(Params, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      NoPressure = .TRUE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " not associated"
      CALL WARN(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      PressurePerm => PressureVar % Perm
      NoPressure = .FALSE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PorosityName = ListGetString(Params, &
         'Porosity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
      WRITE(PorosityName,'(A)') 'Porosity'
    ELSE
      WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    ConstantPorosity= GetLogical(Params,'Constant Porosity', Found)
    IF ((.NOT.Found) .OR. (.NOT.ConstantPorosity)) THEN
      PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
      IF (.NOT.ASSOCIATED(PorosityVar)) THEN
        CALL FATAL(SolverName,'Porosity Variable not found')
      ELSE
        Porosity => PorosityVar % Values
        PorosityPerm => PorosityVar % Perm
      END IF
    ELSE
      NULLIFY(PorosityVar)
    END IF

    TemperatureName = ListGetString(Params, &
         'Temperature Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
      WRITE(TemperatureName,'(A)') 'Temperature'
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    TemperatureVar => VariableGet(Solver % Mesh % Variables,TemperatureName)
    IF (.NOT.ASSOCIATED(TemperatureVar)) THEN
      CALL FATAL(SolverName,'Temperature Variable not found')
    ELSE
      Temperature => TemperatureVar % Values
      TemperaturePerm => TemperatureVar % Perm
    END IF


    GWfluxName = ListGetString(Params, &
         'Groundwater Flux Variable', GivenGWFlux )
    IF (GivenGWFlux) THEN
      WRITE(Message,'(A,A)') "'Groundwater flux Variable' found and set to: ", GWfluxName
      CALL INFO(SolverName,Message,Level=9)
      GWfluxVar1 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 1")
      IF (.NOT.ASSOCIATED(GWfluxVar1)) THEN
        PRINT *, TRIM(GWfluxName) // " 1", " not found"
        GivenGWflux = .FALSE.
      END IF
      IF (DIM > 1) THEN
        GWfluxVar2 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 2")
        IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
          PRINT *, TRIM(GWfluxName) // " 2", " not found"
          GivenGWflux = .FALSE.
        END IF
        IF (DIM > 2) THEN
          GWfluxVar3 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 3")
          IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
            PRINT *, TRIM(GWfluxName) // " 3", " not found"
            GivenGWflux = .FALSE.
          END IF
        END IF
      END IF
      GWflux1 => GWfluxVar1 % Values
      GWfluxPerm1 => GWfluxVar1 % Perm
      IF (DIM > 1) THEN
        GWflux2 => GWfluxVar2 % Values
        GWfluxPerm2 => GWfluxVar2 % Perm
        IF (DIM > 2) THEN
          GWflux3 => GWfluxVar3 % Values
          GWfluxPerm3 => GWfluxVar3 % Perm
        END IF
      END IF
      CALL INFO(SolverName,'Groundwater flux Variable found. Using this as prescribed groundwater flux',Level=9)
    END IF
  END SUBROUTINE AssignVarsSolute

  
  ! compute element-wise nodal variables
  SUBROUTINE ReadVarsSolute(N)
    INTEGER :: N,I
    REAL(KIND=dp) :: p0

    NodalPressure(1:N) = 0.0_dp
    NodalSalinity(1:N) = 0.0_dp
    NodalGWflux(1:3,1:N) = 0.0_dp
    NodalPorosity(1:N) = 0.0_dp
    ! Nodal variable dependencies
    NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
    IF (ConstantPorosity) THEN
      NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
      IF (.NOT.Found) THEN
        WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
        CALL FATAL(SolverName,Message)
      END IF
    ELSE
      NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoPressure) THEN
      CALL INFO(SolverName,'No Pressure variable found - setting to "Reference Pressure"',Level=9)
      p0 = GetConstReal(Model % Constants, 'Reference Pressure', Found)
      IF (.NOT.Found) THEN
        p0 = 100132.0_dp
        CALL INFO(SolverName, ' "Reference Pressure not found in Constants and set to default value p0=100132.0',Level=9)
      END IF
      NodalPressure(1:N) = p0
    ELSE
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoSalinity) THEN
      NodalSalinity(1:N) = 0.0_dp
    ELSE
      NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (GivenGWflux) THEN
      NodalGWflux(1,1:N) = &
           GWflux1(GWfluxPerm1(Element % NodeIndexes(1:N)))
      IF (DIM > 1) THEN
        NodalGWflux(2,1:N) = &
             GWflux2(GWfluxPerm2(Element % NodeIndexes(1:N)))
        IF (DIM > 2) THEN
          NodalGWflux(3,1:N) = &
               GWflux3(GWfluxPerm3(Element % NodeIndexes(1:N)))
        END IF
      END IF
    ELSE
      NodalGWflux(1:3,1:N) = 0.0_dp
    END IF
  END SUBROUTINE ReadVarsSolute
  !------------------------------------------------------------------------------
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSolute(  Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, GivenGWflux,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementNo, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: vstarAtIP(3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,KgwppAtIP(3,3),fwAtIP,mugwAtIP!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP, & ! material properties at IP
         CcYcTAtIP, CcYcPAtIP, CcYcYcAtIP, rhocPAtIP, rhocYcAtIP, rhocTAtIP,& ! material properties at IP
         DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3), fcAtIP(3), extforceFlux(3)! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixSolute)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Solute Source', Found )

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementNo  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))
      vstarAtIP = 0.0_dp ! CHANGE to SUM(  Basis(1:N) * NodalRockVelocity(1:N) )
      
      ! Gradients of Variables at IP - moved from flux computation to here
      DO i=1,DIM
        gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
        gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
      END DO
      
      !Materialproperties needed at IP

      ! water/ice densitities

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      
      !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
        bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
        gwaAtIP = gwa(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP)
        gwaTAtIP =  gwaT(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)!        
        gwapAtIP = 1.0_dp/rhowAtIP
        giaAtIP = gia(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP)
        giaTAtIP = giaT(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP)
        giapAtIP = 1.0_dp/rhoiAtIP
        deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW        
        D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
        D2AtIP= D2(deltaInElement,bijAtIP)! NEW
        B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
             GasConstant,bijAtIP,TemperatureAtIP)
        B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
        XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,gwapAtIP,giapAtIP,Xi0Tilde,&
             deltaInElement,GasConstant,TemperatureAtIP)
        !PRINT *,"Solute: gwaTAtIP,giaTAtIP,gwaAtIP,giaAtIP,TemperatureAtIP",&
        !     gwaTAtIP,giaTAtIP,gwaAtIP,giaAtIP,TemperatureAtIP
        XiYcAtIP = XiYc(CurrentRockMaterial,RockMaterialID,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,bijYcAtIP,&
             Xi0Tilde,XiAtIP,deltaInElement,SalinityAtIP)
      END SELECT    
      ! solute and rock densities and derivatives
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      rhocPAtIP = rhocP(CurrentSoluteMaterial,rhocAtIP)
      rhocYcAtIP = rhocYc(CurrentSoluteMaterial,rhocAtIP,XiAtIP,SalinityAtIP)
      rhocTAtIP = rhocT(CurrentSoluteMaterial,rhocAtIP,T0,TemperatureAtIP)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW

      ! capacities of solutes      
      CcYcTAtIP = CcYcT(rhocTAtIP,PorosityAtIP,SalinityAtIP)
      CcYcPAtIP = CcYcP(rhocPAtIP,PorosityAtIP, SalinityAtIP)
      CcYcYcAtIP = CcYcYc(rhocAtIP,rhocYcAtIP,PorosityAtIP, SalinityAtIP)

      ! groundwater viscosity is pulled outtside flux computation, as needed anyway
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)

      JgwDAtIP = 0.0_dp
      IF (GivenGWFlux) THEN
        !PRINT *, "Solute: Interpolate Flux"
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N)) 
        END DO
      ELSE        
        !PRINT *, "Solute: Compute Flux"
        
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        !PRINT *,"Solute: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        
        ! gradT and gradP have been moved upwards, as needed elsewhere
        
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM)

      END IF
      
      ! parameters for diffusion-dispersion flow
      r12AtIP = GetR(CurrentSoluteMaterial,GasConstant,rhocAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)
      DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)
      KcAtIP = GetKc(CurrentRockMaterial,RockMaterialID,DmAtIP,XiAtIP,JgwDAtIP,PorosityAtIP)
      KcYcYcAtIP = GetKcYcYc(KcAtIP,r12AtIP)
      fcAtIP = GetFc(rhocAtIP,rhowAtIP,Gravity,r12AtIP,XiTAtIP,XiPAtIP,XiAtIP,gradPAtIP,gradTAtIP) 
      
      Weight = IP % s(t) * DetJ

      DO p=1,nd
        DO q=1,nd
          ! diffusion term (Porosity rhoc KcYcYc.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + &
                   Weight * PorosityAtIP * rhocAtIP * KcYcYcAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
              !PRINT *,"cond", Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          ! advection term (CcYcYc * (v* .grad(u)),v) ! V* not implemented yet and set to zero
          ! -----------------------------------
          STIFF (p,q) = STIFF(p,q) &
               + Weight * CcYcYcAtIP * SUM(vstarAtIP(1:DIM)*dBasisdx(q,1:dim)) * Basis(p)

          ! left overs from partial integration of fluxes
          ! (rhoc/Xi) (u, Jgw.grad(v)) 
          STIFF (p,q) = STIFF(p,q) &
               + Weight * rhocAtIP * Basis(q) * SUM(JgwDAtIP(1:DIM) * dBasisdx(p,1:DIM))/ XiAtIP

          ! porosity rhoc  (u,(Kc.fc).grad(v))
          DO i=1,DIM
            extforceFlux =  SUM(KcAtIP(i,1:DIM)*fcAtIP(1:DIM))
          END DO          
          STIFF (p,q) = STIFF(p,q) &
               + Weight * PorosityAtIP * rhocAtIP * Basis(q) * SUM(extforceFlux(1:DIM) * dBasisdx(p,1:DIM))
          ! time derivative (CcYcYc*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * CcYcYcAtIP  * Basis(q) * Basis(p) 
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSolute
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSolute( Element, n, nd )
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Salinity(n), F,Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,SalinityAtIP
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    LOGICAL :: Stat,Found,FluxCondition,WeakDirichletCond
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostSoluteTransport (LocalMatrixBCSolute)'

    SAVE Nodes,DIM
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
    
    Flux(1:n)  = GetReal( BoundaryCondition,'Solute Flow', FluxCondition )
    ! Check, whether we have a weakly imposed Dirichlet condition
    Salinity(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)

    ! Numerical integration:
    !-----------------------
    IF (FluxCondition .OR. WeakDirichletCond) THEN ! spare us, if natural BC
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ
        ! Given flux:
        ! -----------
        IF (Fluxcondition) THEN

          F = SUM(Basis(1:n)*flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
          ! Given pressure, weakly imposed
          !----------------------------------------------------------------------
        ELSE IF (WeakDirichletCond) THEN
          SalinityAtIP = SUM(Salinity(1:n)*Basis(1:n))
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * SalinityAtIP * Basis(1:nd)
        END IF
      END DO
      CALL DefaultUpdateEquations(STIFF,FORCE)
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCSolute
  !------------------------------------------------------------------------------

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostSoluteTransport


!-----------------------------------------------------------------------------
!> heat transfer equation for enhanced permafrost model
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE PermafrostHeatTransfer( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,GWfluxVar1,GWfluxVar2,GWfluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,ElementWiseRockMaterial!,NodalGWflux,NoGWflux
  !------------------------------------------------------------------------------

  CALL DefaultStart()


  Params => GetSolverParams()

  CALL AssignVarsHTEQ()

  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      Material => GetMaterial()
      IF (FirstTime) THEN

        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Active)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      END IF

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF

      CALL ReadVarsHTEQ(N)

      CALL LocalMatrixHTEQ(  Element, t, Active, n, nd+nb,&
           NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
           NodalGWflux,GivenGWflux,&
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()

    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCHTEQ(  Element, n, nd+nb )
        !PRINT *,t,"of",Active,":",n, nb
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()

CONTAINS

  ! compute element-wise nodal variables
  SUBROUTINE ReadVarsHTEQ(N)
    INTEGER :: N,I
    REAL(KIND=dp) :: p0

    NodalPressure(1:N) = 0.0_dp
    NodalSalinity(1:N) = 0.0_dp
    NodalGWflux(1:3,1:N) = 0.0_dp
    NodalPorosity(1:N) = 0.0_dp
    ! Nodal variable dependencies
    NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
    IF (ConstantPorosity) THEN
      NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
      IF (.NOT.Found) THEN
        WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
        CALL FATAL(SolverName,Message)
      END IF
    ELSE
      NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoPressure) THEN
      CALL INFO(SolverName,'No Pressure variable found - setting to "Reference Pressure"',Level=9)
      p0 = GetConstReal(Model % Constants, 'Reference Pressure', Found)
      IF (.NOT.Found) THEN
        p0 = 100132.0_dp
        CALL INFO(SolverName, ' "Reference Pressure not found in Constants and set to default value p0=100132.0',Level=9)
      END IF
      NodalPressure(1:N) = p0
    ELSE
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoSalinity) THEN
      NodalSalinity(1:N) = 0.0_dp
    ELSE
      NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (GivenGWflux) THEN
      NodalGWflux(1,1:N) = &
           GWflux1(GWfluxPerm1(Element % NodeIndexes(1:N)))
      IF (DIM > 1) THEN
        NodalGWflux(2,1:N) = &
             GWflux2(GWfluxPerm2(Element % NodeIndexes(1:N)))
        IF (DIM > 2) THEN
          NodalGWflux(3,1:N) = &
               GWflux3(GWfluxPerm3(Element % NodeIndexes(1:N)))
        END IF
      END IF
    ELSE
      NodalGWflux(1:3,1:N) = 0.0_dp
    END IF
  END SUBROUTINE ReadVarsHTEQ

  ! find needed variables and assign pointers
  SUBROUTINE AssignVarsHTEQ()

    IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
      DIM = CoordinateSystemDimension()
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (AllocationsDone) &
           DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,&
           NodalSalinity,NodalGWflux)
      ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),&
           NodalSalinity(N),NodalGWflux(3,N),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        AllocationsDone = .TRUE.
        CALL INFO(SolverName,"Allocations Done",Level=1)
      END IF

    END IF
    Temperature => Solver % Variable % Values
    TemperaturePerm => Solver % Variable % Perm

    PressureName = ListGetString(Params, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      NoPressure = .TRUE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " not associated"
      CALL WARN(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      PressurePerm => PressureVar % Perm
      NoPressure = .FALSE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PorosityName = ListGetString(Params, &
         'Porosity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
      WRITE(PorosityName,'(A)') 'Porosity'
    ELSE
      WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    ConstantPorosity= GetLogical(Params,'Constant Porosity', Found)
    IF ((.NOT.Found) .OR. (.NOT.ConstantPorosity)) THEN
      PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
      IF (.NOT.ASSOCIATED(PorosityVar)) THEN
        CALL FATAL(SolverName,'Porosity Variable not found')
      ELSE
        Porosity => PorosityVar % Values
        PorosityPerm => PorosityVar % Perm
      END IF
    ELSE
      NULLIFY(PorosityVar)
    END IF

    SalinityName = ListGetString(Params, &
         'Salinity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Salinity Variable' not found. Using default 'Salinity' ")
      WRITE(SalinityName,'(A)') 'Salinity'
    ELSE
      WRITE(Message,'(A,A)') "'Salinity Variable' found and set to: ", SalinityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    SalinityVar => VariableGet(Solver % Mesh % Variables,SalinityName)
    IF (.NOT.ASSOCIATED(SalinityVar)) THEN
      CALL WARN(SolverName,'Salinity Variable not found. Switching Salinity off')
      NoSalinity = .TRUE.
    ELSE
      Salinity => SalinityVar % Values
      SalinityPerm => SalinityVar % Perm
      NoSalinity=.FALSE.
    END IF


    GWfluxName = ListGetString(Params, &
         'Groundwater Flux Variable', GivenGWFlux )
    IF (GivenGWFlux) THEN
      WRITE(Message,'(A,A)') "'Groundwater flux Variable' found and set to: ", GWfluxName
      CALL INFO(SolverName,Message,Level=9)
      GWfluxVar1 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 1")
      IF (.NOT.ASSOCIATED(GWfluxVar1)) THEN
        PRINT *, TRIM(GWfluxName) // " 1", " not found"
        GivenGWflux = .FALSE.
      END IF
      IF (DIM > 1) THEN
        GWfluxVar2 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 2")
        IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
          PRINT *, TRIM(GWfluxName) // " 2", " not found"
          GivenGWflux = .FALSE.
        END IF
        IF (DIM > 2) THEN
          GWfluxVar3 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 3")
          IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
            PRINT *, TRIM(GWfluxName) // " 3", " not found"
            GivenGWflux = .FALSE.
          END IF
        END IF
      END IF
      GWflux1 => GWfluxVar1 % Values
      GWfluxPerm1 => GWfluxVar1 % Perm
      IF (DIM > 1) THEN
        GWflux2 => GWfluxVar2 % Values
        GWfluxPerm2 => GWfluxVar2 % Perm
        IF (DIM > 2) THEN
          GWflux3 => GWfluxVar3 % Values
          GWfluxPerm3 => GWfluxVar3 % Perm
        END IF
      END IF
      CALL INFO(SolverName,'Groundwater flux Variable found. Using this as prescribed groundwater flux',Level=9)
    END IF
  END SUBROUTINE AssignVarsHTEQ

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHTEQ(  Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, GivenGWflux,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementNo, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,KgwppAtIP(3,3),fwAtIP,mugwAtIP!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Heat Source', Found )

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementNo  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))
      
      !Materialproperties needed at IP

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW

      !PRINT *,"HTEQ: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
        bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
        gwaAtIP = gwa(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP)
        gwaTAtIP =  gwaT(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)!        
        gwapAtIP = 1.0_dp/rhowAtIP
        giaAtIP = gia(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP)
        giaTAtIP = giaT(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP)
        giapAtIP = 1.0_dp/rhoiAtIP
        deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW        
        D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
        D2AtIP= D2(deltaInElement,bijAtIP)! NEW
        B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
             GasConstant,bijAtIP,TemperatureAtIP)
        B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
        XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
        !PRINT *,"HTEQ: gwaTAtIP,giaTAtIP,gwaAtIP,giaAtIP,TemperatureAtIP",&
        !     gwaTAtIP,giaTAtIP,gwaAtIP,giaAtIP,TemperatureAtIP
        XiYcAtIP = XiYc(CurrentRockMaterial,RockMaterialID,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,bijYcAtIP,&
             Xi0Tilde,XiAtIP,deltaInElement,SalinityAtIP)
      END SELECT

      !Materialproperties needed at IP:
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      
      ! heat capacities
      csAtIP   = cs(CurrentRockMaterial,RockMaterialID,&
           T0,TemperatureAtIP,ConstVal)
      cwAtIP   = cw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"cw",T0,TemperatureAtIP,SalinityAtIP,cw0,&
      !     acw,bcw,acwl,bcwl
      !PRINT *, "cwAtIP", cwAtIP, "cw0",cw0,"acw",acw,"bcw",bcw,"T0",T0,SalinityAtIP,TemperatureAtIP,PressureAtIP
      ciAtIP   = ci(CurrentSolventMaterial,&
        T0,TemperatureAtIP,ConstVal)
      !ci(ci0,aci,T0,TemperatureAtIP,PressureAtIP)  !!! NEW
      ccAtIP   = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"HTEQ: cw,ci,cs,cc",cwAtIP,ciAtIP,csAtIP,ccAtIP

      ! latent heat
      hiAtIP = hi(CurrentSolventMaterial,&
        T0,TemperatureAtIP,ConstVal)
      hwAtIP = hw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
      !IF ((TemperatureAtIP < 273.65) .AND. (TemperatureAtIP > 272.65)) PRINT *,"hw/hi/XiT/Xi",hwAtIP,hiAtIP,XiTAtIP,XiAtIP

      ! heat conductivity at IP
      ksthAtIP = GetKalphath(CurrentRockMaterial % ks0th(RockMaterialID),&
           CurrentRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
      kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)      
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP,&
           SalinityATIP,PorosityAtIP,meanfactor)
      !IF ((TemperatureAtIP < 273.25))&! .AND. (TemperatureAtIP > 273.05))&
      !     PRINT *, "HTEQ: KGTTAtIP",KGTTAtIP(1,1),TemperatureAtIP,XiAtIP,XiTAtIP

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      !IF ((TemperatureAtIP < 273.65) .AND. (TemperatureAtIP > 272.65))&
      !     PRINT *, "HTEQ: CGTTAtIP", CGTTAtIP, XiAtIP,XiTAtIP,hiAtIP,hwAtIP
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
      !PRINT *,"CgwTTAtIP",CgwTTAtIP,rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,SalinityAtIP
      ! compute groundwater flux for advection term


      CGTpAtIP = GetCGTp(rhoiAtIP,hiAtIP,hwAtIP,XiPAtIP,PorosityAtIP) !NEW
      CGTycAtIP = GetCGTyc(rhoiAtIP,hiAtIP,hwAtIP,XiYcAtIP,PorosityAtIP) !NEW
      ! groundwater flux
      !PRINT *, "KGTTAtIP", KGTTAtIP,"CgwTTAtIP",CgwTTAtIP
      JgwDAtIP = 0.0_dp
      IF (GivenGWFlux) THEN
        !PRINT *, "HTEQ: Interpolate Flux"
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N)) 
        END DO
      ELSE        
        !PRINT *, "HTEQ: Compute Flux"
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        !PRINT *,"HTEQ: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO
        
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM)

      END IF

      Weight = IP % s(t) * DetJ

      DO p=1,nd
        DO q=1,nd
          ! diffusion term (KGTTAtIP.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
              !PRINT *,"cond", Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          ! advection term (CgwTT * (Jgw.grad(u)),v)
          ! -----------------------------------
          !IF (.NOT.NoGWFlux .OR. ComputeGWFlux) THEN
          STIFF (p,q) = STIFF(p,q) + Weight * &
               CgwTTAtIP * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p) !
          ! PRINT *,JgwDAtIP(1:dim),CgwTTAtIP,Weight,"adv",Weight * CgwTTAtIP * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)
          !END IF
          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * (CGTTAtIP + 0.0_dp*CGTycAtIP) * Basis(q) * Basis(p) ! check CGTycAtIP
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixHTEQ
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCHTEQ( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Fluxcondition,Robincondition
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BoundaryCondition

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BoundaryCondition,'Heat Flux', FluxCondition )
    Coeff(1:n) = GetReal( BoundaryCondition,'Heat Transfer Coefficient', RobinCondition )
    Ext_t(1:n) = GetReal( BoundaryCondition,'External Temperature', RobinCondition )

    IF (FluxCondition .OR. RobinCondition)  THEN
      ! Numerical integration:
      !-----------------------
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ

        ! Evaluate terms at the integration point:
        !------------------------------------------

        ! Given flux:
        ! -----------
        F = SUM(Basis(1:n)*flux(1:n))

        ! Robin condition (C*(u-u_0)):
        ! ---------------------------
        C = SUM(Basis(1:n)*coeff(1:n))
        Ext = SUM(Basis(1:n)*ext_t(1:n))

        IF (Robincondition) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C*Ext * Basis(1:nd)
        ELSE IF (Fluxcondition) THEN
          FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
        END IF
      END DO
    END IF
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCHTEQ
  !------------------------------------------------------------------------------

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostHeatTransfer
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!> output of unfrozen water content as variable (post-processing)
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE PermafrostUnfrozenWaterContent( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),WaterContentPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),WaterContent(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeXiT=.FALSE.
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostUnfrozenWaterContent'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName, PhaseChangeModel

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NumberOfRockRecords,NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  ComputeXiT = GetLogical(Params,"Compute XiT",Found)
  IF (.NOT.Found) ComputeXiT=.FALSE.
  
  CALL DefaultInitialize()
  
  ! Assign output variables
  WaterContent => Solver % Variable % Values
  WaterContentPerm => Solver % Variable % Perm

  ! Read Variables
  CALL AssignVarsXi()
  Active = GetNOFActive()
  

  DO t=1,Active
    Element => GetActiveElement(t)      
    n  = GetElementNOFNodes(Element)
    Material => GetMaterial(Element)
    PhaseChangeModel = ListGetString(Material, &
         'Permafrost Phase Change Model', Found )
    IF (Found) THEN
      WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
      CALL INFO(SolverName,Message,Level=9)
    END IF
    IF (FirstTime) THEN
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants, CurrentRockMaterial )
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
        FirstTime = .FALSE.
      END IF
      CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
      CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      dim = CoordinateSystemDimension()
    END IF
    CALL ReadVarsXi(N)
    CALL LocalMatrixXi(  Element, n, NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
          CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
          PhaseChangeModel,ComputeXiT)
  END DO
  
  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()

  ! Trim values into [0,1]
  DO I=1,Model % NumberOfNodes
    WaterContent(WaterContentPerm(I)) =  MAX(WaterContent(WaterContentPerm(I)),0.0_dp)
    WaterContent(WaterContentPerm(I)) =  MIN(WaterContent(WaterContentPerm(I)),1.0_dp)
  END DO

  CALL INFO("SolverName","Computation of unfrozen water content (Xi) for post-processing done",Level=1)

CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixXi( Element, n, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, CurrentSoluteMaterial,&
       CurrentSolventMaterial, PhaseChangeModel, ComputeXit)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    LOGICAL :: ComputeXiT
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP, giaAtIP, gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fwAtIp, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0,DeltaT, T0, p0,eps,Gravity(3)
    REAL(KIND=dp) :: rhowAtIP, rhoiAtIP, rhosAtIP, rhocAtIP
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=DP) :: gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),fluxPAtIP(3),fluxgAtIP(3)
    REAL(KIND=dp) :: MASS(n,n), STIFF(n,n), FORCE(n), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixXi)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, GasConstant, N0,DeltaT, T0, p0,eps,Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      dim = CoordinateSystemDimension()
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF
    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      ! Material properties at IP
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiTAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
     

      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
        bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
        gwaAtIP = gwa(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP) ! NEW
        gwaTAtIP =  gwaT(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)!
        gwapAtIP = 1.0_dp/rhowAtIP
        giaAtIP = gia(CurrentSolventMaterial,& 
             p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP) !NEW
        giaTAtIP = giaT(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP)
        giapAtIP = 1.0_dp/rhoiAtIP
        deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW
        D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
        D2AtIP= D2(deltaInElement,bijAtIP)! NEW
        B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
             GasConstant,bijAtIP,TemperatureAtIP)
        B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
        XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,gwapAtIP,giapAtIP,Xi0Tilde,&
             deltaInElement,GasConstant,TemperatureAtIP)
      END SELECT
      
      Weight = IP % s(t) * DetJ

      DO p=1,n
        DO q=1,n
          Stiff(p,q) = Stiff(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO
      IF (ComputeXiT) THEN
        FORCE(1:n) = FORCE(1:n) + Weight * XiTAtIP * Basis(1:n)
      ELSE
        FORCE(1:n) = FORCE(1:n) + Weight * XiAtIP * Basis(1:n)
      END IF
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixXi

  !------------------------------------------------------------------------------
  SUBROUTINE ReadVarsXi(N)
    INTEGER :: N
    ! Nodal variable dependencies
    NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
    !PRINT *, NodalTemperature(1:N), TemperaturePerm(Element % NodeIndexes(1:N)), Element % NodeIndexes(1:N)
    IF (ConstantPorosity) THEN
      NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
      IF (.NOT.Found) THEN
        WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
        CALL FATAL(SolverName,Message)
      END IF
    ELSE
      NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoPressure) THEN
      NodalPressure(1:N) = 0.0_dp
    ELSE
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoSalinity) THEN
      NodalSalinity(1:N) = 0.0_dp
    ELSE
      NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    END IF
  END SUBROUTINE ReadVarsXi

  SUBROUTINE AssignVarsXi()

    IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (AllocationsDone) &
           DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,&
           NodalSalinity)
      ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),&
           NodalSalinity(N),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        AllocationsDone = .TRUE.
        CALL INFO(SolverName,"Allocations Done",Level=1)
      END IF

    END IF

    TemperatureName = ListGetString(Params, &
         'Temperature Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
      WRITE(TemperatureName,'(A)') 'Temperature'
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    TemperatureVar => VariableGet(Solver % Mesh % Variables,TemperatureName)
    IF (.NOT.ASSOCIATED(TemperatureVar)) THEN
      WRITE(Message,'(A,A,A)') "'Temperature Variable ", TRIM(TemperatureName), " not associated"
      CALL FATAL(SolverName,Message)
    ELSE
      Temperature => TemperatureVar % Values
      TemperaturePerm => TemperatureVar % Perm
      WRITE(Message,'(A,A,A)') "'Temperature Variable ", TRIM(TemperatureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PressureName = ListGetString(Params, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      NoPressure = .TRUE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " not associated"
      CALL WARN(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      PressurePerm => PressureVar % Perm
      NoPressure = .FALSE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PorosityName = ListGetString(Params, &
         'Porosity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
      WRITE(PorosityName,'(A)') 'Porosity'
    ELSE
      WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    ConstantPorosity= GetLogical(Params,'Constant Porosity', Found)
    IF ((.NOT.Found) .OR. (.NOT.ConstantPorosity)) THEN
      PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
      IF (.NOT.ASSOCIATED(PorosityVar)) THEN
        CALL FATAL(SolverName,'Porosity Variable not found')
      ELSE
        Porosity => PorosityVar % Values
        PorosityPerm => PorosityVar % Perm
      END IF
    ELSE
      NULLIFY(PorosityVar)
    END IF

    SalinityName = ListGetString(Params, &
         'Salinity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Salinity Variable' not found. Using default 'Salinity' ")
      WRITE(SalinityName,'(A)') 'Salinity'
    ELSE
      WRITE(Message,'(A,A)') "'Salinity Variable' found and set to: ", SalinityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    SalinityVar => VariableGet(Solver % Mesh % Variables,SalinityName)
    IF (.NOT.ASSOCIATED(SalinityVar)) THEN
      CALL WARN(SolverName,'Salinity Variable not found. Switching Salinity off')
      NoSalinity = .TRUE.
    ELSE
      Salinity => SalinityVar % Values
      SalinityPerm => SalinityVar % Perm
      NoSalinity=.FALSE.
    END IF


  END SUBROUTINE AssignVarsXi
END SUBROUTINE PermafrostUnfrozenWaterContent


!-----------------------------------------------------------------------------
!>  salinity transport equation for enhanced permafrost model
!----------------------------------------------------------------------------

!!! TEMPORARILY REMOVED 


!-----------------------------------------------------------------------------
!>  initialization of Porosity to given reference value in material
!----------------------------------------------------------------------------
!==============================================================================
SUBROUTINE PorosityInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
  INTEGER, POINTER :: PorosityPerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:)
  REAL(KIND=dp), ALLOCATABLE :: NodalHits(:)
  INTEGER :: DIM, i, j, k, NumberOfRockRecords,RockMaterialID,CurrentNode,Active
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PorosityInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,ElementRockMaterialName
  LOGICAL :: Visited = .False., Found, GotIt,ElementWiseRockMaterial

  SAVE Visited,ElementWiseRockMaterial,NodalHits
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------
  
  ! Execute solver only once at beginning
  if (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing porosity to reference ', Level=1)
  CALL Info(SolverName, 'levels in material file            ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  ! Get variables
  DIM = CoordinateSystemDimension()

  ! Get info
  SolverParams => GetSolverParams()

  PorosityName = ListGetString(SolverParams, &
       'Porosity Variable', GotIt )
  IF (.NOT.GotIt) THEN
    PorosityName = "Porosity"
    CALL WARN(SolverName, ' "Porosity Variable" not found - trying default "Porosity"')
  END IF
  PorosityVariable => VariableGet( Solver % Mesh % Variables, PorosityName,GotIt )

  IF ( ASSOCIATED( PorosityVariable ) ) THEN
    PorosityPerm    => PorosityVariable % Perm
    PorosityValues  => PorosityVariable % Values
  ELSE
    CALL FATAL(SolverName, 'Could not find "Porosity Variable"')
  END IF
  PorosityValues = 0.0_dp
  
  ! Loop over elements
  Active = Solver % NumberOFActiveElements
  IF (.NOT.Visited) THEN
    ALLOCATE(NodalHits(Solver % Mesh % NumberOfNodes))
    NodalHits = 0.0_dp
  ELSE
    NodalHits = 0.0_dp
  END IF
  DO i = 1, Active
    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial()
    IF (.NOT.Visited) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Active)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
      END IF
      dim = CoordinateSystemDimension()
      Visited=.True.
    END IF
    
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF
    
    ! Loop over nodes of element
    DO k = 1, GetElementNOFNodes(CurrentElement)
      CurrentNode = PorosityPerm(CurrentElement % NodeIndexes(k))
      IF (PorosityValues(CurrentNode) >= 0.0_dp) THEN
        PorosityValues(CurrentNode) = &
             0.5_dp*(PorosityValues(CurrentNode) + CurrentRockMaterial % eta0(RockMaterialID))
      ELSE
        PorosityValues(CurrentNode) = CurrentRockMaterial % eta0(RockMaterialID)
      END IF
    END DO
  END DO

  !==============================================================================
END SUBROUTINE PorosityInit
!==============================================================================

!-----------------------------------------------------------------------------
!>  initialization of arbitrary scalar given by nodal file
!----------------------------------------------------------------------------
!==============================================================================
SUBROUTINE NodalVariableInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: NodalVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  INTEGER, POINTER :: NodalVariablePerm(:)
  INTEGER,PARAMETER :: io=26
  REAL(KIND=dp), POINTER :: NodalVariableValues(:)
  REAL(KIND=dp) :: InputField
  INTEGER :: DIM, i, CurrentNode, NumberOfNodes, OK,  counter
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="NodalVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: NodalVariableName,NodalVariableFileName
  LOGICAL :: Visited = .False., Found, GotIt,ElementWiseRockMaterial

  SAVE Visited,ElementWiseRockMaterial
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------
  
  ! Execute solver only once at beginning
  if (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing variable to reference ', Level=1)
  CALL Info(SolverName, 'levels in material file            ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  ! Get variables
  DIM = CoordinateSystemDimension()

  ! Get variable to fill in
  SolverParams => GetSolverParams()

  NodalVariableName = ListGetString(SolverParams, &
       'Nodal Variable', GotIt )
  IF (.NOT.GotIt) THEN
    CALL FATAL(SolverName, ' "Nodal Variable" not found')
  END IF
  NodalVariable => VariableGet( Solver % Mesh % Variables, NodalVariableName,GotIt )

  IF ( ASSOCIATED( NodalVariable ) ) THEN
    NodalVariablePerm    => NodalVariable % Perm
    NodalVariableValues  => NodalVariable % Values
    WRITE (Message,*) 'Reading variable ',TRIM(NodalVariableName)
    CALL INFO(SolverName,Message,Level=3)
  ELSE
    WRITE (Message,*) 'Could not find ',TRIM(NodalVariableName)
    CALL FATAL(SolverName, Message)
  END IF
  NodalVariableValues = 0.0_dp

  NodalVariableFileName = ListGetString(SolverParams, &
       'Nodal Variable File', GotIt, Unfoundfatal=.TRUE. )
  
  NumberOfNodes = Model % Mesh % NumberOfNodes
  
  OPEN(unit = io, file = TRIM(NodalVariableFileName), status = 'old',iostat = ok)
  IF (ok /= 0) THEN
    WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(NodalVariableFileName)
    CALL FATAL(TRIM(SolverName),TRIM(message))
  ELSE
    !------------------------------------------------------------------------------
    ! Read in the number of records in file (first line integer)
    !------------------------------------------------------------------------------
    DO i=1,NumberOfNodes
      READ (io, *, END=70, IOSTAT=OK, ERR=80) counter, InputField
      IF (counter .NE. i) CALL FATAL(SolverName,'No concecutive numbering in file')
      IF (NodalVariablePerm(i)==0) CALL FATAL(SolverName,'No corresponding entry of target variable')
      NodalVariableValues(NodalVariablePerm(i)) = InputField
     ! PRINT *,i,counter
    END DO
    !PRINT *, "END", i,counter
70  IF (i-1 .NE. NumberOfNodes) THEN
      WRITE (Message,*) 'Number of records ',i,' in file ',&
           TRIM(NodalVariableFileName),' does not match number of nodes ', NumberOfNodes, ' in mesh'
      CALL FATAL(SolverName,Message)
    END IF
    CLOSE (io)
    RETURN
80  CALL FATAL(SolverName,"I/O error")
  END IF
END SUBROUTINE NodalVariableInit

!-----------------------
!
!--------------------
!-----------------------------------------------------------------------------
!> output of material parameter
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE PermafrostMaterialOutput( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,TensorComponent(2)
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),OutputPropertyPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),OutputProperty(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeXiT=.FALSE.
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostUnfrozenWaterContent'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName, PhaseChangeModel, VariableName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  VariableName = ListGetString(Params,"Output Property",Found)

  IF (.NOT.Found) &
    CALL INFO(SolverName, ' "Output Property" not found')
    
  TensorComponent(1:2) = ListGetInteger(Params,"Output Property Component",Found)
  IF (.NOT.Found) TensorComponent = 1
  
  CALL DefaultInitialize()
  
  ! Assign output variables
  OutputProperty => Solver % Variable % Values
  OutputPropertyPerm => Solver % Variable % Perm

  ! Read Variables
  CALL AssignVarsMaterialOutput()
  Active = GetNOFActive()
  

  DO t=1,Active
    Element => GetActiveElement(t)      
    n  = GetElementNOFNodes(Element)
    Material => GetMaterial(Element)
    PhaseChangeModel = ListGetString(Material, &
         'Permafrost Phase Change Model', Found )
    IF (Found) THEN
      WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
      CALL INFO(SolverName,Message,Level=9)
    END IF
    IF (FirstTime) THEN
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants, CurrentRockMaterial )
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
        FirstTime = .FALSE.
      END IF
      CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
      CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      dim = CoordinateSystemDimension()
    END IF
    CALL ReadVarsMaterialOutput(N)
    CALL LocalMatrixMaterialOutput(  Element, n, NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
          CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
          PhaseChangeModel,VariableName,TensorComponent)
  END DO
  
  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()
  IF (TRIM(VariableName) == 'kgw') THEN
    DO I=1,Solver % Mesh % NumberOfNodes
      OutputProperty(OutputPropertyPerm(i)) = EXP(OutputProperty(OutputPropertyPerm(i)))
    END DO
  END IF
  CALL INFO("SolverName","Computation of " // TRIM(VariableName) // "done",Level=1)

CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixMaterialOutput( Element, n, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, CurrentSoluteMaterial,&
       CurrentSolventMaterial, PhaseChangeModel,VariableName,TensorComponent)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,TensorComponent(2)
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP,bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP, giaAtIP, gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), gwAtIp, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0,DeltaT, T0, p0,eps,Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhowAtIP, rhoiAtIP, rhosAtIP, rhocAtIP
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=DP) :: PropertyAtIP
    REAL(KIND=dp) :: MASS(n,n), STIFF(n,n), FORCE(n), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixMaterialOutput)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, GasConstant, N0,DeltaT, T0, p0,eps,Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))


      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiTAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal) !!! NEW
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        bijAtIP = GetBij(CurrentSoluteMaterial,Xi0tilde,SalinityAtIP)! NEW
        bijYcAtIP = GetBijYc(CurrentSoluteMaterial,Xi0tilde,bijAtIP,SalinityAtIP)! NEW
        gwaAtIP = gwa(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP) ! NEW
        gwaTAtIP =  gwaT(CurrentSolventMaterial,&
             p0,T0,rhowAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)! NEW
        gwapAtIP = 1.0_dp/rhowAtIP ! NEW
        giaAtIP = gia(CurrentSolventMaterial,& 
             p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP) !NEW
        giaTAtIP = giaT(CurrentSolventMaterial,&
             p0,T0,rhoiAtIP,TemperatureAtIP)
        giapAtIP = 1.0_dp/rhoiAtIP
        deltaGAtIP = deltaG(gwaAtIP,giaAtIP) ! NEW
        D1AtIP= D1(CurrentRockMaterial,RockMaterialID,deltaInElement,bijAtIP)! Changed Argument
        D2AtIP= D2(deltaInElement,bijAtIP)! NEW
        B1AtIP = GetB1(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,deltaInElement,deltaGAtIP,&
             GasConstant,bijAtIP,TemperatureAtIP)
        B2AtIP = GetB2(CurrentSolventMaterial,deltaInElement,deltaGAtIP,GasConstant,bijAtIP,TemperatureAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde)
        XiTAtIP= XiT(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,Xi0Tilde,bijAtIP,p0,&
             deltaInElement,deltaGAtIP,T0,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             B1AtIP,B2AtIP,D1AtIP,D2AtIP,bijAtIP,gwapAtIP,giapAtIP,Xi0Tilde,&
             deltaInElement,GasConstant,TemperatureAtIP)
      END SELECT
      ! those etwo are not really needed at the moment
      !rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!! NEW
      !rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      
      SELECT CASE(VariableName)
      CASE('kgw')

        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP =  GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        IF (KgwAtIP(TensorComponent(1),TensorComponent(2)) <= 0.0_dp) STOP
        PropertyAtIP = LOG(KgwAtIP(TensorComponent(1),TensorComponent(2)))
        !IF (PropertyAtIP < 0.0_dp) THEN
        !  PRINT *, "KgwAtIP < 0 =", KgwAtIP, muw0,muw0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw!
         ! STOP
       ! ELSE
       !   PRINT *, "KgwAtIP=", PropertyAtIP
       ! END IF
      !CASE('rhow')
      !  CALL FATAL(SolverName,' "rhow" not implemented yet")
      CASE DEFAULT
        WRITE(Message,*) ' Variable "', TRIM(VariableName), '" not implemented.'
         CALL FATAL(SolverName,Message)
      END SELECT
      
      Weight = IP % s(t) * DetJ

      DO p=1,n
        DO q=1,n
          Stiff(p,q) = Stiff(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO
      FORCE(1:n) = FORCE(1:n) + Weight * PropertyAtIP * Basis(1:n)
      
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixMaterialOutput

  !------------------------------------------------------------------------------
  SUBROUTINE ReadVarsMaterialOutput(N)
    INTEGER :: N
    ! Nodal variable dependencies
    NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
    !PRINT *, NodalTemperature(1:N), TemperaturePerm(Element % NodeIndexes(1:N)), Element % NodeIndexes(1:N)
    IF (ConstantPorosity) THEN
      NodalPorosity(1:N) = ListGetReal(Material,PorosityName,N,Element % NodeIndexes, Found)
      IF (.NOT.Found) THEN
        WRITE (Message,'(A,A,A)') "No '",TRIM(PorosityName) ,"'found in Material"
        CALL FATAL(SolverName,Message)
      END IF
    ELSE
      NodalPorosity(1:N) = Porosity(PorosityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoPressure) THEN
      NodalPressure(1:N) = 0.0_dp
    ELSE
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
    END IF
    IF (NoSalinity) THEN
      NodalSalinity(1:N) = 0.0_dp
    ELSE
      NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    END IF
  END SUBROUTINE ReadVarsMaterialOutput

  SUBROUTINE AssignVarsMaterialOutput()

    IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (AllocationsDone) &
           DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,&
           NodalSalinity)
      ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),&
           NodalSalinity(N),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        AllocationsDone = .TRUE.
        CALL INFO(SolverName,"Allocations Done",Level=1)
      END IF

    END IF

    TemperatureName = ListGetString(Params, &
         'Temperature Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
      WRITE(TemperatureName,'(A)') 'Temperature'
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    TemperatureVar => VariableGet(Solver % Mesh % Variables,TemperatureName)
    IF (.NOT.ASSOCIATED(TemperatureVar)) THEN
      WRITE(Message,'(A,A,A)') "'Temperature Variable ", TRIM(TemperatureName), " not associated"
      CALL FATAL(SolverName,Message)
    ELSE
      Temperature => TemperatureVar % Values
      TemperaturePerm => TemperatureVar % Perm
      WRITE(Message,'(A,A,A)') "'Temperature Variable ", TRIM(TemperatureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PressureName = ListGetString(Params, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      NoPressure = .TRUE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " not associated"
      CALL WARN(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      PressurePerm => PressureVar % Perm
      NoPressure = .FALSE.
      WRITE(Message,'(A,A,A)') "'Pressure Variable ", TRIM(PressureName), " associated"
      CALL INFO(SolverName,Message,Level=9)
    END IF

    PorosityName = ListGetString(Params, &
         'Porosity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
      WRITE(PorosityName,'(A)') 'Porosity'
    ELSE
      WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    ConstantPorosity= GetLogical(Params,'Constant Porosity', Found)
    IF ((.NOT.Found) .OR. (.NOT.ConstantPorosity)) THEN
      PorosityVar => VariableGet(Solver % Mesh % Variables,PorosityName)
      IF (.NOT.ASSOCIATED(PorosityVar)) THEN
        CALL FATAL(SolverName,'Porosity Variable not found')
      ELSE
        Porosity => PorosityVar % Values
        PorosityPerm => PorosityVar % Perm
      END IF
    ELSE
      NULLIFY(PorosityVar)
    END IF

    SalinityName = ListGetString(Params, &
         'Salinity Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Salinity Variable' not found. Using default 'Salinity' ")
      WRITE(SalinityName,'(A)') 'Salinity'
    ELSE
      WRITE(Message,'(A,A)') "'Salinity Variable' found and set to: ", SalinityName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    SalinityVar => VariableGet(Solver % Mesh % Variables,SalinityName)
    IF (.NOT.ASSOCIATED(SalinityVar)) THEN
      CALL WARN(SolverName,'Salinity Variable not found. Switching Salinity off')
      NoSalinity = .TRUE.
    ELSE
      Salinity => SalinityVar % Values
      SalinityPerm => SalinityVar % Perm
      NoSalinity=.FALSE.
    END IF


  END SUBROUTINE AssignVarsMaterialOutput
END SUBROUTINE PermafrostMaterialOutput














