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
! *  Authors: Thomas Zwinger
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

  TYPE RockMaterial_t
     INTEGER :: NumberOfRecords
     REAL(KIND=dp), ALLOCATABLE :: ks0th(:),ew(:),bs(:),rhos0(:),&
          cs0(:),Xi0(:),eta0(:),hs0(:),Kgwh0(:,:,:),qexp(:),alphaL(:),alphaT(:),As0(:)
     REAL(KIND=dp) :: GasConstant, Mw,Mc,DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
          hw0,hi0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc
     CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  END TYPE RockMaterial_t

CONTAINS

  FUNCTION ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
       NumberOfRecords,GasConstant,Mw,Mc,DeltaT,T0,p0,rhow0,rhoi0,rhoc0,&
       l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,&
       Gravity) RESULT(Constantsread)
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER :: DIM, NumberOfRecords
    REAL(KIND=dp) :: GasConstant, Mw, Mc,DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,&
         Gravity(3)
    LOGICAL :: Constantsread
    !------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: gWork(:,:)
    LOGICAL :: Found
    INTEGER :: I
    !------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    gWork => ListGetConstRealArray( Model % Constants,'Gravity',Found)
    IF (.NOT.Found) THEN
      Gravity = 0.0
      CALL WARN(FunctionName,'Gravity not found in Constants section. Setting to zero')
    ELSE
      Gravity = gWork(1:3,1)*gWork(4,1)
    END IF
    GasConstant= CurrentRockMaterial % GasConstant
    Mw= CurrentRockMaterial % Mw
    Mc= CurrentRockMaterial % Mc
    DeltaT= CurrentRockMaterial % DeltaT
    T0= CurrentRockMaterial % T0
    p0= CurrentRockMaterial % p0
    rhow0= CurrentRockMaterial % rhow0
    rhoi0= CurrentRockMaterial % rhoi0
    rhoc0= CurrentRockMaterial % rhoc0
    cw0= CurrentRockMaterial % cw0
    ci0= CurrentRockMaterial % ci0
    cc0= CurrentRockMaterial % cc0
    eps=CurrentRockMaterial % eps
    kw0th= CurrentRockMaterial % kw0th
    ki0th= CurrentRockMaterial % ki0th
    kc0th= CurrentRockMaterial % kc0th
    mu0 = CurrentRockMaterial % mu0
    Dm0 = CurrentRockMaterial % Dm0
    dw1 = CurrentRockMaterial % dw1
    dw2 = CurrentRockMaterial % dw2
    dc0 = CurrentRockMaterial % dc0  
    dc1 = CurrentRockMaterial % dc1
    bw = CurrentRockMaterial % bw
    bi = CurrentRockMaterial % bi
    bc = CurrentRockMaterial % bc
    l0= (CurrentRockMaterial % hw0) - (CurrentRockMaterial % hi0)
    ConstantsRead=.TRUE.

    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(FunctionName,"General Constants:", Level=9)
    WRITE(Message,'(A)') "GasConstant,Mw,Mc,DeltaT,T0,p0,rhow0,rhoi0,rhoc0,hw0,hi0,cw0,ci0"
    CALL INFO(FunctionName,Message,Level=9)
    WRITE(Message,'(A)') "eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1:"
    CALL INFO(FunctionName,Message,Level=9)
    WRITE(Message,'(24E12.5)') CurrentRockMaterial % GasConstant, &
         CurrentRockMaterial % Mw, CurrentRockMaterial % Mc, &
         CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0, CurrentRockMaterial % p0,&
         CurrentRockMaterial % rhow0, CurrentRockMaterial % rhoi0, CurrentRockMaterial % rhoc0,&
         CurrentRockMaterial % hw0, CurrentRockMaterial % hi0,&
         CurrentRockMaterial % cw0,CurrentRockMaterial % ci0, CurrentRockMaterial % cc0,&
         CurrentRockMaterial % eps, CurrentRockMaterial % kw0th,&
         CurrentRockMaterial % ki0th, CurrentRockMaterial % kc0th, CurrentRockMaterial % mu0,&
         CurrentRockMaterial % Dm0, CurrentRockMaterial % dw1, CurrentRockMaterial % dw2, &
         CurrentRockMaterial % dc0, CurrentRockMaterial % dc1
    CALL INFO(FunctionName,Message,Level=9)
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(FunctionName,"Material Constants:", Level=9)
    ! Read in material specific values
    DO I=1,NumberOfRecords
      WRITE(Message,'(I2,A,A,A)') I,": ", CurrentRockMaterial % VariableBaseName(I),":"
      WRITE(Message,'(A)') "Xi0,eta0,ks0th,Xi0,ew,b,rhos0,cs0:"
      CALL INFO(FunctionName,Message,Level=9)
      WRITE(Message,'(E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5)') CurrentRockMaterial % Xi0(I),&
           ",",CurrentRockMaterial % eta0(I), ",", CurrentRockMaterial % Ks0th(I), "," ,&
           CurrentRockMaterial % ew(I), ",", CurrentRockMaterial % bs(I),",", CurrentRockMaterial % rhos0(I),&
           ",",CurrentRockMaterial % cs0(I)
      CALL INFO(FunctionName,Message,Level=9)
    END DO
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9) 

  END FUNCTION ReadPermafrostRockMaterialConstants

  SUBROUTINE ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
       ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,GasConstant, Mw,Mc,&
       DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)

    IMPLICIT NONE
    TYPE(Element_t) :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    REAL(KIND=dp), INTENT(IN) :: GasConstant,Mw, Mc,DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity(3) ! constants that need to have been set
    INTEGER :: RockMaterialID,DIM
    REAL(KIND=dp), INTENT(OUT) :: meanfactor,MinKgw,ks0th,&
         ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(1:3,1:3),qexp,alphaL,alphaT,As0
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    CHARACTER(LEN=MAX_NAME_LEN) :: SubroutineName="ReadPermafrostRockMaterialVariables"

    Material => GetMaterial(Element)
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(SubroutineName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    ew = CurrentRockMaterial % ew(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
    IF (DIM == 2) THEN
      Kgwh0(1,2) =  Kgwh0(1,3)
      Kgwh0(2,2) =  Kgwh0(3,3)
      Kgwh0(2,1) =  Kgwh0(3,1)
      Kgwh0(3,1:3) = 0.0_dp
      Kgwh0(1:3,3) = 0.0_dp
    END IF
    qexp = CurrentRockMaterial % qexp(RockMaterialID)
    alphaL = CurrentRockMaterial % alphaL(RockMaterialID)
    alphaT = CurrentRockMaterial % alphaT(RockMaterialID)
    As0 = CurrentRockMaterial % As0(RockMaterialID)
    ! derive element rock material specific parameters
    deltaInElement = delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    D1InElement = D1(deltaInElement,ew)
    D2InElement = 1.0_dp
  END SUBROUTINE ReadPermafrostRockMaterialVariables

  FUNCTION ReadPermafrostRockMaterial(Params,CurrentRockMaterial ) RESULT(NumberOfRecords)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(RockMaterial_t), TARGET :: LocalRockMaterial
    Integer :: NumberOfRecords

    INTEGER :: i,j,k,l, n,t, active, DIM, ok,InitialNumberOfRecords, EntryNumber
    INTEGER,parameter :: io=20
    LOGICAL :: Found, fexist, FirstTime=.TRUE., AllocationsDone=.FALSE., DataRead=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  MaterialFileName, NewMaterialFileName, str, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='ReadPermafrostRockMaterial'

    SAVE AllocationsDone,DataRead,InitialNumberOfRecords,LocalRockMaterial,MaterialFileName

    IF (DataRead) THEN
      NewMaterialFileName = GetString( Params, 'Rock Material File', Found )
      IF (.NOT.Found) THEN
        CALL INFO(FunctionName," 'Rock Material File' keyword not found - looking for default DB!")
        fexist = .FALSE.
#ifdef USE_ISO_C_BINDINGS
        str = 'ELMER_LIB'
#else
        str = 'ELMER_LIB'//CHAR(0)
#endif
        CALL envir( str,NewMaterialFileName,k ) 
        IF ( k > 0  ) THEN
          NewMaterialFileName = NewMaterialFileName(1:k) // '/permafrostmaterialdb.dat'
          INQUIRE(FILE=TRIM(NewMaterialFileName), EXIST=fexist)
        END IF
        IF (.NOT. fexist) THEN
#ifdef USE_ISO_C_BINDINGS
          str = 'ELMER_HOME'
#else
          str = 'ELMER_HOME'//CHAR(0)
#endif
          CALL envir( str,NewMaterialFileName,k ) 
          IF ( k > 0 ) THEN
            NewMaterialFileName = NewMaterialFileName(1:k) // '/share/elmersolver/lib/' // 'permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(NewMaterialFileName), EXIST=fexist)
          END IF
          IF ((.NOT. fexist) .AND. k>0) THEN
            NewMaterialFileName = NewMaterialFileName(1:k) // '/permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(NewMaterialFileName), EXIST=fexist)
          END IF
        END IF
        IF (.NOT. fexist) THEN
          CALL Fatal('CheckKeyWord', 'permafrostmaterialdb.dat not found')
        END IF
      END IF
      IF (NewMaterialFileName /= MaterialFileName) THEN
        WRITE (Message, '(A,A,A,A)') NewMaterialFileName,' does not match existing datafile ', MaterialFileName,'. Exiting!'
        CALL FATAL(FunctionName,Message)
      END IF
      NumberOfRecords = InitialNumberOfRecords
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
        READ (io, *, END=10, IOSTAT=OK, ERR=30) NumberOfRecords, Comment
        WRITE (Message,'(A,I2,A,A,A,A)') "Attempting read ",NumberOfRecords," ",&
             TRIM(Comment)," from data file ",TRIM(MaterialFileName)
        CALL INFO(FunctionName,Message,level=3)
        InitialNumberOfRecords = NumberOfRecords
      END IF
      !------------------------------------------------------------------------------
      ! Allocate and read stuff
      !------------------------------------------------------------------------------
      !M = Model % Mesh % NumberOfNodes
      IF (AllocationsDone) THEN
        DEALLOCATE(&
             LocalRockMaterial % ks0th,&
             LocalRockMaterial % ew,&
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
             LocalRockMaterial % As0, &
             LocalRockMaterial % VariableBaseName)
      END IF
      ALLOCATE(&
           LocalRockMaterial % ks0th(NumberOfRecords),&
           LocalRockMaterial % ew(NumberOfRecords),&
           LocalRockMaterial % bs(NumberOfRecords),&
           LocalRockMaterial % rhos0(NumberOfRecords),&
           LocalRockMaterial % cs0(NumberOfRecords),&
           LocalRockMaterial % Xi0(NumberOfRecords),&
           LocalRockMaterial % eta0(NumberOfRecords),&
           LocalRockMaterial % hs0(NumberOfRecords),&
           LocalRockMaterial % Kgwh0(3,3,NumberOfRecords),&
           LocalRockMaterial % qexp(NumberOfRecords), &
           LocalRockMaterial % alphaL(NumberOfRecords), &
           LocalRockMaterial % alphaT(NumberOfRecords), &
           LocalRockMaterial % As0(NumberOfRecords), &
           LocalRockMaterial % VariableBaseName(NumberOfRecords),&
           STAT=OK)
      AllocationsDone = .TRUE.
      DataRead = .TRUE.

      IF (OK /= 0) THEN
        CLOSE(io)
        CALL FATAL(FunctionName, 'Allocation Error of input data array')
      END IF
      !------------------------------------------------------------------------------
      ! Read in information from material file
      ! General constants
      !  GasConstant, Mw, Mc,DeltaT, T0, p0, rhow0,rhoi0,rhoc0,rhos0
      !     hw0,hi0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc
      !------------------------------------------------------------------------------
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % GasConstant, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Mw, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Mc, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % DeltaT,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % T0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % p0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhow0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhoi0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhoc0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % hw0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % hi0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % cw0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ci0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % cc0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % eps, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % kw0th, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ki0th, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % kc0th, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % mu0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Dm0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % dw1, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % dw2, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % dc0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % dc1, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % bw, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % bi, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % bc, Comment
      ! for each material
      !       ks0th(:),ew(:),b(:),rhos0(:),cs0(:)
      DO I=1,NumberOfRecords
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % VariableBaseName(I), EntryNumber
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
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Xi0(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % eta0(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ks0th(I), Comment          
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ew(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % bs(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhos0(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % cs0(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % hs0(I), Comment
        DO J=1,3
          DO K=1,3
            READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Kgwh0(J,K,I), Comment
          END DO
        END DO
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % qexp(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % alphaL(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % alphaT(I), Comment
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % As0(I), Comment
      END DO
      WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRecords," rock material records from file ", TRIM(MaterialFileName)
      CALL INFO(FunctionName,Message,Level=1)
10    CLOSE(io)
      IF (I < NumberOfRecords) THEN
        WRITE(Message,'(I3,A,I3)') I,"records read, which is smaller than given number ", NumberOfRecords
        CALL FATAL(FunctionName,Message)
      ELSE
        CurrentRockMaterial => LocalRockMaterial
        WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRecords," rock material records from file ", TRIM(MaterialFileName)
        CALL INFO(FunctionName,Message,Level=1)
        CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
        CALL INFO(FunctionName,"General Constants:", Level=9)
        WRITE(Message,'(A)') "GasConstant,Mw,Mc,DeltaT,T0,p0,rhow0,rhoi0,hw0,hi0,cw0,ci0,eps,kw0th,ki0th,mu0:"
        CALL INFO(FunctionName,Message,Level=9)
        WRITE(Message,'(16E12.5)') CurrentRockMaterial % GasConstant, &
             CurrentRockMaterial % Mw, CurrentRockMaterial % Mc, CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0,&
             CurrentRockMaterial % p0, CurrentRockMaterial % rhow0, CurrentRockMaterial % rhoi0,&           
             CurrentRockMaterial % hw0, CurrentRockMaterial % hi0, CurrentRockMaterial % cw0,&
             CurrentRockMaterial % ci0, CurrentRockMaterial % eps, CurrentRockMaterial % kw0th,&
             CurrentRockMaterial % ki0th, CurrentRockMaterial % mu0
        CALL INFO(FunctionName,Message,Level=9)
        CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
        CALL INFO(FunctionName,"Material Constants:", Level=9)
        DO I=1,NumberOfRecords
          WRITE(Message,'(I2,A,A,A)') I,": ", CurrentRockMaterial % VariableBaseName(I),":"
          WRITE(Message,'(A)') "Xi0,eta0,Ks0th,Xi0,ew,b,rhos0,cs0:"
          CALL INFO(FunctionName,Message,Level=9)
          WRITE(Message,'(E10.5,E10.5,E10.5,E10.5,E10.5,E10.5,E10.5)') CurrentRockMaterial % Xi0(I),&
               CurrentRockMaterial % eta0(I), CurrentRockMaterial % Ks0th(I), &
               CurrentRockMaterial % ew(I),CurrentRockMaterial % bs(I),CurrentRockMaterial % rhos0(I),&
               CurrentRockMaterial % cs0(I)
          CALL INFO(FunctionName,Message,Level=9)
        END DO
        CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)  
      END IF
      RETURN
    END IF

30  CALL WARN(FunctionName,"I/O error! Last successfully read variable:")
    CALL WARN(FunctionName,Comment)
    CALL FATAL(FunctionName,"Stopping simulation")    
  END FUNCTION ReadPermafrostRockMaterial

  FUNCTION groundwaterflux(Salinity,Kgwpp,KgwpT,Kgw,gradp,gradT,Gravity,rhow,rhoc,DIM) RESULT(JgwD)
    IMPLICIT NONE
    REAL (KIND=dp), INTENT(IN) :: Salinity,Kgwpp(3,3),KgwpT(3,3),Kgw(3,3),gradp(3),gradT(3),Gravity(3),&
         rhow,rhoc
    REAL (KIND=dp)  :: JgwD(3)
    INTEGER, INTENT(IN) :: DIM
    INTEGER :: i
    REAL (KIND=dp) :: Fluxp(3),fluxT(3),gFlux(3)

    JgwD(1:3) = 0.0
    DO i=1,DIM
      fluxp(i) = -1.0_dp * SUM(Kgwpp(i,1:DIM)*gradp(1:DIM))          
      !fluxT(i) =  -1.0_dp * SUM(KgwpT(i,1:DIM)*gradT(1:DIM))
      gFlux(i) = ((1.0_dp - Salinity) * rhow + Salinity * rhoc) * SUM(Kgw(i,1:DIM)*Gravity(1:DIM))
    END DO
    JgwD(1:DIM) = fluxp(1:DIM) + fluxT(1:DIM) + gFlux(1:DIM)    
  END FUNCTION groundwaterflux

  REAL (KIND=dp) FUNCTION delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant
    delta = 0.5_dp*l0*DeltaT/T0 &
         + (cw0 - ci0)*((T0 + 0.5_dp*DeltaT)*LOG(1.0_dp + 0.5_dp*DeltaT/T0) - 0.5_dp*DeltaT)
    delta = delta*(eps*(1.0_dp - eps)/(2.0_dp*eps -1.0_dp))* Mw/(GasConstant*(T0 + 0.5_dp*DeltaT))
  END FUNCTION delta

  REAL (KIND=dp) FUNCTION deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,dw1,dw2,&
       GasConstant,Temperature,Pressure,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,dw1,dw2,&
         GasConstant,Temperature,Pressure,Salinity
    REAL(KIND=dp) :: relSalinity
    IF ((Salinity < 1.0_dp) .AND. (Salinity >= 0.0_dp)) THEN
      relSalinity = Salinity/(1.0_dp - Salinity)
    ELSE
      CALL WARN("deltaG", "Salinity either too small or too large - removing salinity effect")
      relSalinity = 0.0_dp
    END IF
    deltaG = -l0*(Temperature - T0)/T0 &  ! the first one is L-Zero, do not add a _dp to it!
         + ((1.0_dp/rhow0) + (1.0_dp/rhoi0))*(Pressure - p0) &
         - (cw0 - ci0)*(Temperature*LOG(Temperature/T0) - (Temperature - T0)) &
         - GasConstant * Temperature *(dw1 * relSalinity + dw2 * (relSalinity**2.0_dp))/Mc
    !deltaG = 1.0_dp ! REMOVE THIS LINE !
  END FUNCTION deltaG

  REAL (KIND=dp) FUNCTION B1(delta,deltaG,ew,Mw,GasConstant,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,deltaG,ew,Mw,GasConstant,Temperature
    B1 = Mw*deltaG/(GasConstant*Temperature*(ew + delta))
  END FUNCTION B1

  REAL  (KIND=dp) FUNCTION B2(delta,deltaG,GasConstant,Mw,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,deltaG,Mw,GasConstant,Temperature
    B2 = Mw*deltaG/(GasConstant*Temperature*delta)
  END FUNCTION B2

  REAL (KIND=dp) FUNCTION D1(delta,ew)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,ew
    ! local
    D1 = delta/(ew + delta)
  END FUNCTION D1

  FUNCTION GetXiAnderson(A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) RESULT(XiAnderson)
    REAL(KIND=dp), INTENT(IN) :: A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar, XiAnderson
    IF (Porosity <= 0.0) &
         CALL FATAL("Permafrost(GetXiAnderson)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    XiAnderson =  MAX(MIN((rhos0/rhow)*(A*(Tstar**B)/Porosity),1.0_dp),0.0_dp)
  END FUNCTION GetXiAnderson

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

!  REAL (KIND=dp) FUNCTION XiAndersonEta(A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) RESULT(XiAnderson)
!    REAL(KIND=dp), INTENT(IN) :: A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
!    REAL(KIND=dp) :: Tstar
!    IF (Porosity <= 0.0) &
!         CALL FATAL("Permafrost(GetXiAndersonEta)","Zero or negative porosity detected")
!    Tstar = T0 - Beta * Pressure - Temperature
!    IF (Xi == 1_dp .OR. Xi == 0.0_dp) THEN
!      XiAndersonEta = 0.0_dp
!    ELSE
!      XiAndersonEta =  -(rhos0/rhow)*(A*(Tstar**B)/(Porosity*Porosity)
!    END IF
!  END FUNCTION GetXiAnderson
  
  FUNCTION GetXi0Tilde(Xi0,mu0,Porosity) RESULT(Xi0tilde)
    REAL(KIND=dp), INTENT(IN) :: Xi0,mu0,Porosity
    REAL(KIND=dp) Xi0tilde    
    IF (Porosity <= 0.0_dp) &
         CALL FATAL("Permafrost(GetXi)","Zero or negative porosity detected")
    Xi0tilde = MIN(Xi0 * (mu0/Porosity) * (1.0_dp - Porosity)/(1 - mu0),1.0_dp)
  END FUNCTION GetXi0Tilde
  
  FUNCTION GetXi(B1,B2,D1,D2,Xi0tilde) RESULT(Xi)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0tilde
    REAL(KIND=dp) :: Xi
    Xi= Xi0tilde/(1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1)) &
         + (1.0_dp - Xi0tilde)/(1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))
    IF (Xi < 0.0_dp) Xi = 0.0_dp
    IF (Xi > 1.0_dp) Xi = 1.0_dp
  END FUNCTION GetXi

  REAL (KIND=dp) FUNCTION XiT(B1,B2,D1,D2,Xi0,p0,Mw,ew,delta,rhow0,rhoi0,cw0,ci0,&
       l0,T0,GasConstant,Temperature, Pressure)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0,p0,Mw,ew,delta,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,Temperature, Pressure
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
    aux3 = (l0 + (cw0 - ci0)*(Temperature - T0) & ! this l0 is small L and 0 - don't change
         + (-(1.0_dp/rhoi0) + (1.0_dp/rhow0))*(Pressure - p0))/Temperature
    XiT = (0.5_dp*Xi0*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0)*aux2/delta) *Mw*aux3/(T0*GasConstant*Temperature)
  END FUNCTION XiT

  REAL (KIND=dp) FUNCTION XiP(B1,B2,D1,D2,Xi0Tilde,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0Tilde,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2
    IF (Temperature <= 0.0_dp) CALL FATAL("Permafrost (XiP)","(sub-)Zero Temperature detected")
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
    XiP = (0.5_dp*Xi0Tilde*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0Tilde)*aux2/delta) &
         *((1.0_dp/rhoi0) - (1.0_dp/rhow0))* Mw/(GasConstant*Temperature)
  END FUNCTION XiP

  REAL (KIND=dp) FUNCTION XiXc(B1,B2,D1,D2,Xi0Tilde,Mw,Mc,ew,dw1,dw2,delta,GasConstant,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0Tilde,Mw,Mc,ew,dw1,dw2,delta,GasConstant,Salinity
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3
    IF ((Salinity < 1.0_dp) .AND. (Salinity >= 0.0_dp)) THEN
      aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
      aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
      aux3 = 1.0_dp - Salinity
      XiXc = (0.5_dp*Xi0Tilde*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0Tilde)*aux2/delta) &
           *(dw1/(aux3**2.0_dp) + dw2*Salinity/(aux3**3.0_dp)) * Mw/Mc
    ELSE
      CALL WARN("Permafrost(XiXc)","Salinity out of physical range - returning zero")
      XiXc = 0.0_dp
    END IF
  END FUNCTION XiXc

  REAL (KIND=dp) FUNCTION XiEta(B1,B2,D1,D2,Xi0,Xi0Tilde,eta0,Porosity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0,Xi0Tilde,eta0,Porosity
    !local
    REAL(KIND=dp) :: aux1, aux2
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

  REAL (KIND=dp) FUNCTION rhos(rhos0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhos0,TemperatureAtIP,PressureAtIP
    rhos = rhos0
  END FUNCTION rhos
  
  REAL (KIND=dp) FUNCTION rhow(rhow0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,TemperatureAtIP,PressureAtIP
    rhow = rhow0
  END FUNCTION rhow

    REAL (KIND=dp) FUNCTION rhoi(rhoi0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoi0,TemperatureAtIP,PressureAtIP
    rhoi = rhoi0
  END FUNCTION rhoi
  
  REAL (KIND=dp) FUNCTION rhoc(rhoc0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoc0,TemperatureAtIP,PressureAtIP
    rhoc = rhoc0
  END FUNCTION rhoc
  
  REAL (KIND=dp) FUNCTION cs(cs0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: cs0,TemperatureAtIP,PressureAtIP
    cs = cs0
  END FUNCTION cs
  
  REAL (KIND=dp) FUNCTION cw(cw0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: cw0,TemperatureAtIP,PressureAtIP
    cw = cw0
  END FUNCTION cw

  REAL (KIND=dp) FUNCTION ci(ci0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ci0,TemperatureAtIP,PressureAtIP
    ci = ci0
  END FUNCTION ci
  
  REAL (KIND=dp) FUNCTION cc(cc0,TemperatureAtIP,PressureAtIP)  !!! Replace with function
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: cc0,TemperatureAtIP,PressureAtIP
    cc = cc0
  END FUNCTION cc
  
  FUNCTION GetKAlphaTh(kalpha0th,balpha,T0,Temperature)RESULT(kalphath)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: kalpha0th,balpha,T0,Temperature
    REAL(KIND=dp) :: kalphath
    kalphath = kalpha0th/( 1.0_dp + balpha*(Temperature - T0)/T0)
  END FUNCTION GetKAlphaTh

  FUNCTION GetKGTT(ksth,kwth,kith,kcth,Xi,&
       Salinity,Porosity,meanfactor)RESULT(KGTT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ksth,kwth,kith,kcth,Xi,&
         Salinity,Porosity,meanfactor
    ! local
    REAL(KIND=dp) :: KGTT(3,3), KGaTT, KghTT, unittensor(3,3)
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    KGhTT = 1.0_dp/((1.0_dp - Porosity)/ksth + (1.0_dp - Salinity)*Xi*Porosity/kwth &
         + Salinity*Xi*Porosity/kcth + (1.0_dp - Xi)*Porosity/kith)
    KGaTT = (1.0_dp - Porosity)*ksth + (1.0_dp - Salinity)*Xi*Porosity*kwth &
         + Salinity*Xi*Porosity*kcth + (1.0_dp - Xi)*Porosity*kith
    !KGhTT =1.0_dp/(Xi*Porosity/kwth + (1.0_dp - Xi)*Porosity/kith) ! REMOVE
    !KGaTT = Xi*Porosity*kwth + (1.0_dp - Xi)*Porosity*kith ! REMOVE
    KGTT = unittensor*((1.0_dp - meanfactor)*KGhTT + meanfactor * KGaTT)
    
  END FUNCTION GetKGTT

  FUNCTION GetCGTT(Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,l0,&
       Porosity,Salinity)RESULT(CGTT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,&
         l0,Porosity,Salinity
    REAL(KIND=dp) :: CGTT
    CGTT = (1.0_dp - Porosity)*rhos*cs &
         + (1.0_dp - Salinity) * Xi * Porosity * rhow * cw &
         + Salinity * Xi * Porosity * rhoc * cc &
         + (1.0_dp - Xi)*Porosity*rhoi*ci &
         + rhoi*l0*Porosity*XiT
  END FUNCTION GetCGTT

  FUNCTION GetCgwTT(rhow0,rhoc0,cw0,cc0,Salinity)RESULT(CgwTT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,rhoc0,cw0,cc0,Salinity
    REAL(KIND=dp) :: CgwTT
    CgwTT = (1.0_dp - Salinity)*rhow0*cw0 + Salinity*rhoc0*cc0
  END FUNCTION GetCgwTT

  REAL (KIND=dp) FUNCTION fTildewT(B1,Temperature,D1,delta,ew,l0,cw0,ci0,T0,Xi,Xi0)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,Temperature,D1,delta,ew,l0,cw0,ci0,T0,Xi,Xi0
    REAL(KIND=dp) :: aux1, aux2, aux3
    IF (Xi > Xi0) THEN
      fTildewT = 0.0_dp
    ELSE
      aux1 = (l0/T0) + ((cw0  - ci0)*(Temperature  - T0)/T0)
      aux2 = 1.0_dp + B1/(SQRT(B1*B1 + 4.0_dp*D1))
      aux3 = ew/(ew + delta)
      fTildewT = 0.5_dp * aux3 * aux2 * aux1
    END IF
  END FUNCTION fTildewT

  REAL (KIND=dp) FUNCTION fTildewp(B1,D1,delta,ew,rhow0,rhoi0,Xi,Xi0)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,D1,delta,ew,rhow0,rhoi0,Xi,Xi0
    REAL(KIND=dp) :: aux1, aux2, aux3
    IF (Xi > Xi0) THEN
      fTildewp = 0.0_dp
    ELSE
      aux1 = (1.0_dp/rhoi0) - (1.0_dp/rhow0)
      aux2 = 1.0_dp + B1/(SQRT(B1*B1 + 4.0_dp*D1))
      aux3 = ew/(ew + delta)
      fTildewp = 0.5_dp * aux3 * aux2 * aux1
    END IF
  END FUNCTION fTildewp

  FUNCTION GetKgw(mu0,mu,Xi,rhow0,qexp,Kgwh0,MinKgw)RESULT(Kgw)
    IMPLICIT NONE
    !    REAL(KIND=dp), INTENT(IN) :: mu0,Xi,rhow0,qexp,Kgwh0(3,3),MinKgw,mu
    REAL(KIND=dp) :: mu0,Xi,rhow0,qexp,Kgwh0(3,3),MinKgw,mu
    REAL(KIND=dp) :: Kgw(3,3), factor
    REAL(KIND=dp), PARAMETER :: gval=9.81_dp !hard coded, so match Kgwh0 with this value
    INTEGER :: I, J
    IF (mu <= 0.0_dp) &
         CALL FATAL("Permafrost(GetKgw)","Unphysical viscosity detected")
    factor = (mu0/mu)*(Xi**qexp)/(rhow0*gval)
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

  FUNCTION GetKgwpT(rhow0,fTildewT,Kgw)RESULT(KgwpT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,fTildewT,Kgw(3,3)
    REAL(KIND=dp) :: KgwpT(3,3)
    KgwpT = Kgw*rhow0*fTildewT
  END FUNCTION GetKgwpT

  FUNCTION GetKgwpp(rhow0,fTildewp,Kgw)RESULT(Kgwpp)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,fTildewp,Kgw(3,3)
    REAL(KIND=dp) :: Kgwpp(3,3)
    Kgwpp = Kgw*(1.0_dp + rhow0*fTildewp)
  END FUNCTION GetKgwpp

  FUNCTION GetKc(alphaL,alphaT,Dm0,Xi,absJgwD,eL,Porosity)RESULT(Kc)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: alphaL,alphaT,Dm0,Xi,absJgwD,eL(3),Porosity
    REAL(KIND=dp) :: Kc(3,3), unittensor(3,3), aux
    INTEGER :: I,J
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    IF (Porosity <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero Porosity detected")
    IF (Xi <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero water content detected")
    aux = absJgwD/(Porosity * Xi)
    Kc = 0.0_dp
    DO I=1,3
      DO J=1,3
        Kc(I,J) = Kc(I,J) + Dm0 * unittensor(I,J) &
             + aux*((alphaL - alphaT)*eL(I)*eL(J)  + alphaT * unittensor(I,J))
      END DO
    END DO
  END FUNCTION GetKc

  FUNCTION GetKcXcXc(T0,rhoc0,dw1,dw2,dc0,dc1,Kc,Temperature,Salinity,Pressure)RESULT(KcXcXc)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: T0,rhoc0,dw1,dw2,dc0,dc1,Kc(3,3),Temperature,Salinity,Pressure
    REAL(KIND=dp) :: KcXcXc(3,3), aux, relSalinity

    aux = rhoc0 * Temperature/(T0 * rhoc0)
    KcXcXc(1:3,1:3) = aux * Kc(1:3,1:3)
    IF (Salinity >= 1.0_dp) &
         CALL FATAL("GetKcXcXc","(Larger than) unity Salinity detected")
    relSalinity = Salinity/(1.0_dp - Salinity)
    aux = 1.0_dp + ((dw1 + dc1)/dc0) * relSalinity + 2.0_dp * (dw2/dc0) * (relSalinity**2.0_dp)
    KcXcXc(1:3,1:3) = aux * KcXcXc(1:3,1:3)    
  END FUNCTION GetKcXcXc
END MODULE PermafrostMaterials


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
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, Active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:)!,GWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:)!,GWflux(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
       NodalPressure(:) !NodalGWflux(:,:),
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, VarName, PhaseChangeModel 

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,&
       NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure
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
        NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )
        IF (NumberOfRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
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
      CALL LocalMatrixDarcy(  Element, N, ND+NB, NodalPressure, NodalTemperature, &
           NodalPorosity, NodalSalinity, CurrentRockMaterial,PhaseChangeModel)
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
  SUBROUTINE LocalMatrixDarcy( Element, n, nd,  NodalPressure, NodalTemperature, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, PhaseChangeModel )
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),&
         gradPAtIP(3),fluxTAtIP(3),fluxgAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: fTildewTAtIP, fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,As0  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, Mc, DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity(3)    ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP

    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp) , POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow', &
         FunctionName='Permafrost (LocalMatrixDarcy)'

    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,dw1,dw2,dc0,dc1,bw,bi,bc, MinKgw, Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
           l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      LOAD(1:n) = GetReal( BodyForce,'Groundwater source', Found )   
    END IF

    ! read variable material parameters from CurrentRockMaterial
    CALL ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
         ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,&
         GasConstant, Mw, Mc,DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)    

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
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,GasConstant,dw1,dw2,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP)
        B1AtIP = B1(deltaInElement,deltaGAtIP,ew,Mw,GasConstant,TemperatureAtIP)
        B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
        Xi0Tilde = GetXi0Tilde(Xi0,mu0,PorosityAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0Tilde)
        XiTAtIP= XiT(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,p0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,GasConstant,TemperatureAtIP)
      END SELECT

      !Materialproperties needed at IP
      rhosAtIP = rhos(rhos0,TemperatureAtIP,PressureAtIP)  !!! NEW
      rhowAtIP = rhow(rhow0,TemperatureAtIP,PressureAtIP)  !!! NEW
      rhoiAtIP = rhoi(rhoi0,TemperatureAtIP,PressureAtIP)  !!! NEW
      !rhocAtIP = rhoc(rhoc0,TemperatureAtIP,PressureAtIP)  !!! NEW
      !csAtIP   = cs(cs0,TemperatureAtIP,PressureAtIP)  !!! NEW
      !cwAtIP   = cw(cw0,TemperatureAtIP,PressureAtIP)  !!! NEW
      !ciAtIP   = ci(ci0,TemperatureAtIP,PressureAtIP)  !!! NEW
      !ccAtIP   = cc(cc0,TemperatureAtIP,PressureAtIP)  !!! NEW
      
      fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
      fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
      
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw) ! NB it is meant to be rhow0 and kghh0
      KgwpTAtIP = 0.0_dp
      KgwpTAtIP = GetKgwpT(rhowAtIP,fTildewTATIP,KgwAtIP)
      KgwppAtIP = 0.0_dp
      KgwppAtIP = GetKgwpp(rhowAtIP,fTildewpATIP,KgwAtIP)

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
              StiffPQ = StiffPQ + KgwppAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)              
            END DO
          END DO
          !PRINT *,"KgwppAtIP(1,1)",KgwppAtIP(1,1)
          !StiffPQ = 1.0d-04 *  SUM(dBasisdx(p,1:DIM)* dBasisdx(q,1:DIM))
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
        END DO
      END DO
      !IF (XiAtIP > 0.0_dp) THEN
        ! body forcexs
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:n)*dBasisdx(1:n,i))
        END DO
        DO i=1,DIM
          !fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
          fluxgAtIP(i) = ( (1.0_dp - SalinityAtIP) * rhowAtIP  + SalinityAtIP * rhocAtIP)& 
               * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))   !!! NEW
          !fluxgAtIP(i) =rhow0*SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
        END DO
        DO p=1,nd     
          FORCE(p) = FORCE(p) - Weight * SUM(fluxTAtIP(1:DIM)*dBasisdx(p,1:DIM)) !!! NEW
          FORCE(p) = FORCE(p) + Weight * SUM(fluxgAtIP(1:DIM)*dBasisdx(p,1:DIM))        
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
    LOGICAL :: Stat,Found,ConstantsRead=.FALSE.,FluxCondition,WeakDirichletCond
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlow (LocalMatrixBCDarcy)'

    SAVE Nodes, ConstantsRead, DIM
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
  INTEGER :: NumberOfRecords
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
    !------------------------------------------------------------------------------
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER :: elem,t,i,j,k,p,q,n,nd, DIM,Rank, RockMaterialID,NumberOfRecords
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,coeff,detJ,GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    LOGICAL :: Found, ConstantsRead=.FALSE., FirstTime=.TRUE.
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), POINTER :: Conductivity(:,:,:)=>NULL()
    REAL(KIND=dp) :: GasConstant, Mw, Mc, DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity(3)    ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         fluxTAtIP(3),fluxpAtIP(3),fluxgAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: fTildewTAtIP, fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: meanfactor,ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,As0  ! stuff comming from RockMaterial
    REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:), NodalSalinity(:), NodalPressure(:),NodalPorosity(:)
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlux (BulkAssembly)'
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    
    SAVE Conductivity, Nodes, ConstantsRead,DIM,&
         GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,&
         bw,bi,bc,Gravity,NumberOfRecords,CurrentRockMaterial,FirstTime

    n = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )

    ALLOCATE( STIFF(n,n), FORCE(dofs,n) )
    ALLOCATE( Basis(n), dBasisdx(n,3) )
    ALLOCATE( NodalPressure(N),NodalPorosity(N),NodalTemperature(N),NodalSalinity(N) )

    DO elem = 1,Solver % NumberOFActiveElements

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
        NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )
        IF (NumberOfRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
      END IF
      IF(.NOT.ConstantsRead) THEN
        ConstantsRead = ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
             NumberOfRecords,GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
             l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity)
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

      ! read variable material parameters from CurrentRockMaterial
      CALL ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
           ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,&
           GasConstant,Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)

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
        
        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('Anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiTAtIP = &
               XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiPAtIP   = &
               XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
        CASE DEFAULT ! Hartikainen model
          deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,GasConstant,dw1,dw2,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP)
          B1AtIP = B1(deltaInElement,deltaGAtIP,ew,Mw,GasConstant,TemperatureAtIP)
          B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
          Xi0Tilde = GetXi0Tilde(Xi0,mu0,PorosityAtIP)
          XiAtIP = GetXi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0Tilde)
          XiTAtIP= XiT(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,p0,Mw,ew,&
               deltaInElement,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,TemperatureAtIP,PressureAtIP)
          XiPAtIP= XiP(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,ew,&
               deltaInElement,rhow0,rhoi0,GasConstant,TemperatureAtIP)
        END SELECT

        fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
        fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
        !fTildewpATIP = 0.0_dp !!!!!!!!!!!!!!!!!
        KgwAtIP = 0.0_dp
        KgwpTAtIP = 0.0_dp
        KgwppAtIP = 0.0_dp        
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
        !PRINT *,"KgwAtIP",KgwAtIP
        KgwpTAtIP = GetKgwpT(rhow0,fTildewTATIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(rhow0,fTildewpATIP,KgwAtIP)
        Weight = IntegStuff % s(t) * detJ

        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
          END DO
        END DO

        FluxpAtIp = 0.0_dp
        fluxTAtIP = 0.0_dp
        fluxgAtIp = 0.0_dp
        DO i=1,dim
          !fluxpAtIP(i) = 1.0_dp * i
          fluxpAtIP(i) = -1.0_dp * SUM(KgwppAtIP(i,1:DIM)*gradpAtIP(1:DIM))
          !fluxTAtIP(i) =  -1.0_dp * SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
          !fluxgAtIp(i) = ( (1.0_dp - SalinityAtIP) * rhow0  + SalinityAtIP * rhoc0)&
          !SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM)
          fluxgAtIp(i) =    rhow0 * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
        END DO
        !PRINT *,"Flux",fluxpAtIP,"gradp",gradpAtIP(1:DIM),"K",KgwppAtIP(1:DIM,1:DIM)
        DO i=1,dim
          !Coeff = 1.0_dp * Weight * (FluxpAtIP(i) + fluxgAtIp(i) + fluxTAtIP(i))
          Coeff = 1.0_dp * Weight * (FluxpAtIP(i)+ fluxgAtIp(i))
          !PRINT *,i,Coeff,FluxpAtIP(i),fluxgAtIp(i),fluxTAtIP(i)
          FORCE(i,1:nd) = FORCE(i,1:nd) + Coeff * Basis(1:nd)
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
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux_Init"
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  dim = CoordinateSystemDimension()

  IF( dim < 2 .OR. dim > 3 ) THEN
    CALL Fatal('PermafrostGroundwaterFlux_init','Flux computation makes sense only in 2D and 3D')
  END IF


  VarName = TRIM('GroundWater')

  IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    EqName = ListGetString( SolverParams,'Equation')
    CALL ListAddString( SolverParams, 'Variable','-nooutput '//TRIM(EqName)//'_temp' )
  END IF

  FluxName = TRIM(VarName)//' Flux'
  CALL Info('PermafrostGroundwaterFlux_init','Saving flux to: '//TRIM(FluxName), Level=3) 
  IF(dim == 2) THEN
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':2]')
  ELSE IF(dim == 3) THEN
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':3]')
  END IF
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
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,&
       NoGWflux=.TRUE.,ComputeGWFlux=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,NumberOfRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux
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
        NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )        
        IF (NumberOfRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
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

      CALL LocalMatrixHTEQ(  Element, n, nd+nb, NodalTemperature, NodalPressure, &
           NodalPorosity, NodalSalinity, NodalGWflux, ComputeGWFlux, &
           NoGWflux, CurrentRockMaterial, PhaseChangeModel)
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

  SUBROUTINE ReadVarsHTEQ(N)
    INTEGER :: N
    NodalPressure(1:N) = 0.0_dp
    NodalSalinity(1:N) = 0.0_dp
    NodalGWflux(1:3,1:N) = 0.0_dp
    NodalPorosity(1:N) = 0.0_dp
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
    IF (.NOT.NoPressure) THEN
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
    END IF
    IF (.NOT.NoSalinity) THEN
      NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    END IF
    IF (.NOT.NoGWflux) THEN
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
      !PRINT *,"Nodal-up1",NodalGWflux(1,1:N)
      !PRINT *,"Nodal-up2",NodalGWflux(2,1:N)
      !PRINT *,"Nodal-up3",NodalGWflux(3,1:N)
    END IF
   ! PRINT *, "ph",NodalPorosity(1:N),"P", NodalPressure(1:N),"S",&
   !      NodalSalinity(1:N),"J", NodalGWflux(1,1:N)
  END SUBROUTINE ReadVarsHTEQ

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
         'Groundwater Flux Variable', Found )
    IF (.NOT.Found) THEN
      NoGWflux = .TRUE.
    ELSE
      WRITE(Message,'(A,A)') "'Groundwater flux Variable' found and set to: ", GWfluxName
      CALL INFO(SolverName,Message,Level=9)
      !      NoGWflux = .FALSE.
      !    END IF
      NoGWflux = .FALSE.
      GWfluxVar1 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 1")
      IF (.NOT.ASSOCIATED(GWfluxVar1)) THEN
        PRINT *, TRIM(GWfluxName) // " 1", " not found"
        NoGWflux = .TRUE.
      END IF
      IF (DIM > 1) THEN
        GWfluxVar2 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 2")
        IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
          PRINT *, TRIM(GWfluxName) // " 2", " not found"
          NoGWflux = .TRUE.
        END IF
        IF (DIM > 2) THEN
          GWfluxVar3 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 3")
          IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
            PRINT *, TRIM(GWfluxName) // " 3", " not found"
            NoGWflux = .TRUE.
          END IF
        END IF
      END IF
    END IF
    IF (NoGWFlux) THEN  
      IF (NoPressure) THEN
        CALL WARN(SolverName,'Neither Pressure nor Groundwater Flux variable found. No convection will be computed')
        ComputeGWFlux = .FALSE.
      ELSE
        CALL INFO(SolverName,'Groundwater flux Variable not found. Using Pressure and Temperature to compute flux',Level=9)
        ComputeGWFlux = .TRUE.
      END IF
    ELSE
      GWflux1 => GWfluxVar1 % Values
      GWfluxPerm1 => GWfluxVar1 % Perm
      IF (DIM > 1) THEN
        GWflux2 => GWfluxVar2 % Values
        GWfluxPerm2 => GWfluxVar2 % Perm
        IF (DIM > 3) THEN
          GWflux3 => GWfluxVar3 % Values
          GWfluxPerm3 => GWfluxVar3 % Perm
        END IF
      END IF
      NoGWflux = .FALSE.
      ComputeGWFlux = .FALSE.
      CALL INFO(SolverName,'Groundwater flux Variable found. Using this as prescribed groundwater flux',Level=9)
    END IF
  END SUBROUTINE AssignVarsHTEQ

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHTEQ( Element, n, nd, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, NodalGWflux, ComputeGWFlux, NoGWFlux, &
       CurrentRockMaterial,PhaseChangeModel )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL :: ComputeGWFlux,NoGWFlux
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,As0  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, Mc, DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc    ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         GWfluxAtIP(3),StiffPQ, meanfactor
    REAL(KIND=DP) :: gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),fluxPAtIP(3),fluxgAtIP(3),Gravity(3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
    !------------------------------------------------------------------------------   
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      dim = CoordinateSystemDimension()
      ConstantsRead = &
           ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
           l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity)

    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Heat Source', Found )

    ! read variable material parameters from CurrentRockMateria
    CALL ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
         ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,&
         GasConstant,Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)

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
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))
            
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,GasConstant,dw1,dw2,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP)
        B1AtIP = B1(deltaInElement,deltaGAtIP,ew,Mw,GasConstant,TemperatureAtIP)
        B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
        Xi0Tilde = GetXi0Tilde(Xi0,mu0,PorosityAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0Tilde)
        XiTAtIP= XiT(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,p0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,GasConstant,TemperatureAtIP)
      END SELECT

      !Materialproperties needed at IP
      rhosAtIP = rhos(rhos0,TemperatureAtIP,PressureAtIP)  !!! NEW
      rhowAtIP = rhow(rhow0,TemperatureAtIP,PressureAtIP)  !!! NEW
      rhoiAtIP = rhow(rhoi0,TemperatureAtIP,PressureAtIP)  !!! NEW
      rhocAtIP = rhoc(rhoc0,TemperatureAtIP,PressureAtIP)  !!! NEW
      csAtIP   = cs(cs0,TemperatureAtIP,PressureAtIP)  !!! NEW
      cwAtIP   = cw(cw0,TemperatureAtIP,PressureAtIP)  !!! NEW
      ciAtIP   = ci(ci0,TemperatureAtIP,PressureAtIP)  !!! NEW
      ccAtIP   = cc(cc0,TemperatureAtIP,PressureAtIP)  !!! NEW
      
      ! heat conductivity at IP
      ksthAtIP = GetKalphath(ks0th,bs,T0,TemperatureAtIP)
      kwthAtIP = GetKalphath(kw0th,bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(ki0th,bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(kc0th,bc,T0,TemperatureAtIP)      
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP,&
           SalinityATIP,PorosityAtIP,meanfactor)
      
      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,cwAtIP,ciAtIP,csAtIP,ccAtIP,l0,&
           PorosityAtIP,SalinityAtIP)
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,SalinityAtIP)

      ! compute groundwater flux for advection term
      IF (.NOT.ComputeGWFlux .AND. .NOT.NoGWFlux) THEN
        JgwDAtIP = 0.0_dp
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N))
        END DO
      ELSE IF(ComputeGWFlux .AND. NoGWFlux) THEN        
        fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
        fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
        KgwpTAtIP = GetKgwpT(rhow0,fTildewTATIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(rhow0,fTildewpATIP,KgwAtIP)
        gradTAtIP = 0.0_dp
        gradPAtIP = 0.0_dp
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO

        DO i=1,DIM
          fluxTAtIP(i) =  -1.0_dp * SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
          fluxgAtIP(i) = ( (1.0_dp - SalinityAtIP) * rhow0  + SalinityAtIP * rhoc0) * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
          fluxPAtIP(i) =  -1.0_dp * SUM(KgwppAtIP(i,1:DIM)*gradPAtIP(1:DIM))
          !JgwDAtIP(i) = fluxgAtIP(i) - fluxTAtIP(i) - fluxPAtIP(i)
          JgwDAtIP(i) = fluxgAtIP(i) + fluxPAtIP(i)          
        END DO
      ELSE ! nothing at all is computed or read in
        JgwDAtIP(1:DIM) = 0.0_dp
      END IF

      Weight = IP % s(t) * DetJ

      DO p=1,nd
        DO q=1,nd
          ! diffusion term (KGTTAtIP.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          ! advection term (CgwTT * (Jgw.grad(u)),v)
          ! -----------------------------------
          IF (.NOT.NoGWFlux .OR. ComputeGWFlux) &
               STIFF (p,q) = STIFF(p,q) + Weight * &
               CgwTTAtIP * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * CGTTAtIP * Basis(q) * Basis(p)
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
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),WaterContentPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),WaterContent(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.   
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostUnfrozenWaterContent'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName, PhaseChangeModel

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,NumberOfRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature
  !------------------------------------------------------------------------------
  Params => GetSolverParams()

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
      NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )
      IF (NumberOfRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
        FirstTime = .FALSE.
      END IF
      dim = CoordinateSystemDimension()
    END IF
    CALL ReadVarsXi(N)
    CALL LocalMatrixXi(  Element, n, NodalTemperature, NodalPressure, &
         NodalPorosity, NodalSalinity, CurrentRockMaterial,PhaseChangeModel)
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
CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixXi( Element, n, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, PhaseChangeModel )
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,As0  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, Mc, DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc    ! constants read only once
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         GWfluxAtIP(3),StiffPQ, meanfactor
    REAL(KIND=DP) :: gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),fluxPAtIP(3),fluxgAtIP(3),Gravity(3)
    REAL(KIND=dp) :: MASS(n,n), STIFF(n,n), FORCE(n), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixXi)'
    !------------------------------------------------------------------------------   
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,&
         Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      dim = CoordinateSystemDimension()
      ConstantsRead = &
           ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
           l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! read variable material parameters from CurrentRockMateria
    CALL ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
         ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,&
         GasConstant,Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)

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

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('Anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,rhow0,rhos0,T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,GasConstant,dw1,dw2,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP)
        B1AtIP = B1(deltaInElement,deltaGAtIP,ew,Mw,GasConstant,TemperatureAtIP)
        B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
        Xi0Tilde = GetXi0Tilde(Xi0,mu0,PorosityAtIP)
        XiAtIP = GetXi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0Tilde)
        XiTAtIP= XiT(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,p0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,TemperatureAtIP,PressureAtIP)
        XiPAtIP= XiP(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,ew,&
             deltaInElement,rhow0,rhoi0,GasConstant,TemperatureAtIP)
      END SELECT

      Weight = IP % s(t) * DetJ

      DO p=1,n
        DO q=1,n
          Stiff(p,q) = Stiff(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:n) = FORCE(1:n) + Weight * XiAtIP * Basis(1:n)
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
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE PermafrostSalinity( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,TemperatureVar,TemperatureTimeDerVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20
  INTEGER,POINTER :: TemperaturePerm(:), TemperatureTimeDerPerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:),TemperatureTimeDer(:), Pressure(:), Porosity(:),&
       Salinity(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalTemperatureTimeDer(:),NodalGWflux(:,:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoPressure=.TRUE.,&
       NoGWflux=.TRUE.,ComputeGWFlux=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostSalinity'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, TemperatureName,TemperatureTimeDerName, GWfluxName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,NumberOfRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,&
       NodalTemperatureTimeDer,NodalGWflux
  !------------------------------------------------------------------------------

  CALL DefaultStart()


  Params => GetSolverParams()


  CALL AssignVarsSalinity()

  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    WRITE (Message,'(I5,A,I5,A)') iter,' out of ', maxiter,' max non-linear iterations'
    CALL INFO(SolverName,Message,Level=6)

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      ! Read Material information
      Material => GetMaterial(Element)
      IF (FirstTime) THEN
        NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )
        IF (NumberOfRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost rock material read',Level=3)
          FirstTime = .FALSE.
        END IF
      END IF
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      CALL ReadVarsSalinity(N)

      CALL LocalMatrixSalinity(  Element, n, nd+nb, NodalTemperature, NodalTemperatureTimeDer,NodalPressure, &
           NodalPorosity, NodalSalinity, NodalGWflux, ComputeGWFlux, NoGWflux, CurrentRockMaterial)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixSalinityBC(  Element, n, nd+nb )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    !CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()

CONTAINS

  SUBROUTINE ReadVarsSalinity(N)
    INTEGER :: N
    ! Nodal variable dependencies
    NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
    NodalTemperatureTimeDer(1:N) = TemperatureTimeDer(TemperatureTimeDerPerm(Element % NodeIndexes(1:N)))
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
    NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
    IF (NoGWflux) THEN
      NodalGWflux(1:3,1:N) = 0.0_dp
    ELSE
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
    END IF
  END SUBROUTINE ReadVarsSalinity

  SUBROUTINE AssignVarsSalinity()

    IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (AllocationsDone) &
           DEALLOCATE(NodalTemperature,NodalTemperatureTimeDer,NodalPorosity,NodalPressure,&
           NodalSalinity,NodalGWflux)
      ALLOCATE(NodalTemperature(N),NodalTemperatureTimeDer(N),NodalPorosity(N),NodalPressure(N),&
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
      CALL FATAL(SolverName,'Temperature Variable not associtated.')
    ELSE
      Temperature => TemperatureVar % Values
      TemperaturePerm => TemperatureVar % Perm
    END IF


    TemperatureTimeDerName = ListGetString(Params, &
         'TemperatureTimeDer Variable', Found )
    IF (.NOT.Found) THEN
      WRITE(TemperatureTimeDerName,'(A)') TRIM(TemperatureName),' Velocity'
      WRITE(Message,'(A,A)') " 'Temperature Time Derivative Variable' not found. Using default ",&
           TRIM(TemperatureTimeDerName)
      CALL WARN(SolverName,Message)
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Time Derivative Variable' found and set to: ",&
           TemperatureTimeDerName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    TemperatureTimeDerVar => VariableGet(Solver % Mesh % Variables,TemperatureTimeDerName)
    IF (.NOT.ASSOCIATED(TemperatureTimeDerVar)) THEN
      CALL FATAL(SolverName,'Temperature Time Derivative Variable not associtated.')
    ELSE
      TemperatureTimeDer => TemperatureTimeDerVar % Values
      TemperatureTimeDerPerm => TemperatureTimeDerVar % Perm
    END IF

    GWfluxName = ListGetString(Params, &
         'Groundwater Flux Variable', Found )
    IF (.NOT.Found) THEN
      NoGWflux = .TRUE.
    ELSE
      WRITE(Message,'(A,A)') "'Groundwater flux Variable' found and set to: ", GWfluxName
      CALL INFO(SolverName,Message,Level=9)
      !      NoGWflux = .FALSE.
      !    END IF
      NoGWflux = .FALSE.
      GWfluxVar1 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 1")
      IF (.NOT.ASSOCIATED(GWfluxVar1)) THEN
        PRINT *, TRIM(GWfluxName) // " 1", " not found"
        NoGWflux = .TRUE.
      END IF
      IF (DIM > 1) THEN
        GWfluxVar2 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 2")
        IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
          PRINT *, TRIM(GWfluxName) // " 2", " not found"
          NoGWflux = .TRUE.
        END IF
        IF (DIM > 2) THEN
          GWfluxVar3 => VariableGet(Solver % Mesh % Variables,TRIM(GWfluxName) // " 3")
          IF (.NOT.ASSOCIATED(GWfluxVar2)) THEN
            PRINT *, TRIM(GWfluxName) // " 3", " not found"
            NoGWflux = .TRUE.
          END IF
        END IF
      END IF
    END IF
    IF (NoGWFlux) THEN  
      IF (NoPressure) THEN
        CALL WARN(SolverName,'Neither Pressure nor Groundwater Flux variable found. No convection will be computed')
        ComputeGWFlux = .FALSE.
      ELSE
        CALL INFO(SolverName,'Groundwater flux Variable not found. Using Pressure and Temperature to compute flux',Level=9)
        ComputeGWFlux = .TRUE.
      END IF
    ELSE
      GWflux1 => GWfluxVar1 % Values
      GWfluxPerm1 => GWfluxVar1 % Perm
      IF (DIM > 1) THEN
        GWflux2 => GWfluxVar2 % Values
        GWfluxPerm2 => GWfluxVar2 % Perm
        IF (DIM > 3) THEN
          GWflux3 => GWfluxVar3 % Values
          GWfluxPerm3 => GWfluxVar3 % Perm
        END IF
      END IF
      NoGWflux = .FALSE.
      ComputeGWFlux = .FALSE.
      CALL INFO(SolverName,'Groundwater flux Variable found. Using this as prescribed groundwater flux',Level=9)
    END IF
  END SUBROUTINE AssignVarsSalinity

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSalinity( Element, n, nd, NodalTemperature, NodalTempTimeDer, NodalPressure, &
       NodalPorosity, NodalSalinity, NodalGWflux, ComputeGWFlux, NoGWFlux, CurrentRockMaterial )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalTempTimeDer(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL :: ComputeGWFlux,NoGWFlux
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,As0  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, Mc, DeltaT, T0,p0,rhow0,rhoi0,rhoc0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc    ! constants read only once
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP, TemperatureTimeDerAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         GWfluxAtIP(3),StiffPQ, meanfactor,AbsJgwDAtIP, gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),&
         fluxPAtIP(3),fluxgAtIP(3),Gravity(3),el(3),MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n),&
         KcAtIP(3,3), KcXcXcAtIP(3,3),fluxXcG(3),XixcAtIP,XietaAtIP,TimeWeight,ReactionWeight
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixSalinity)'
    !------------------------------------------------------------------------------   
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      dim = CoordinateSystemDimension()
      ConstantsRead = &
           ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,rhoc0,&
           l0,cw0,ci0,cc0,eps,kw0th,ki0th,kc0th,mu0,Dm0,dw1,dw2,dc0,dc1,bw,bi,bc,Gravity)
    END IF

    CALL GetElementNodes( Nodes )

    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Salinity Source', Found )

    ! read variable material parameters from CurrentRockMateria
    CALL ReadPermafrostRockMaterialVariables(Element,CurrentRockMaterial,meanfactor,MinKgw,ks0th,&
         ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0,qexp,alphaL,alphaT,As0,deltaInElement,D1InElement,D2InElement,&
         GasConstant,Mw, Mc, DeltaT, T0, p0, rhow0,rhoi0,l0,cw0,ci0,eps,kw0th,ki0th,mu0,Gravity,DIM)

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
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      TemperatureTimeDerAtIP = SUM( Basis(1:N) * NodalTemperatureTimeDer(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      ! functions at IP
      deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,Mc,l0,cw0,ci0,rhow0,rhoi0,GasConstant,dw1,dw2,&
           TemperatureAtIP,PressureAtIP,SalinityAtIP)
      B1AtIP = B1(deltaInElement,deltaGAtIP,ew,Mw,GasConstant,TemperatureAtIP)
      B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
      CGTTAtIP = GetCGTT(XiAtIP,XiTAtIP,rhos0,rhow0,rhoi0,rhoc0,cw0,ci0,cs0,cc0,l0,PorosityAtIP,SalinityAtIP)
      CgwTTAtIP = GetCgwTT(rhow0,rhoc0,cw0,cc0,SalinityAtIP)

      ! compute advection term
      IF (.NOT.ComputeGWFlux .AND. .NOT.NoGWFlux) THEN
        DO I=1,3
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N))
        END DO
      ELSE IF(ComputeGWFlux .AND. NoGWFlux) THEN        
        fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
        fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
        KgwpTAtIP = GetKgwpT(rhow0,fTildewTATIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(rhow0,fTildewpATIP,KgwAtIP)
        gradTAtIP = 0.0_dp
        gradPAtIP = 0.0_dp
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO

        DO i=1,DIM
          fluxTAtIP(i) =  -1.0_dp * SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
          fluxgAtIP(i) = ( (1.0_dp - SalinityAtIP) * rhow0  + SalinityAtIP * rhoc0) * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
          fluxPAtIP(i) =  -1.0_dp * SUM(KgwppAtIP(i,1:DIM)*gradPAtIP(1:DIM))
          !JgwDAtIP(i) = fluxgAtIP(i) - fluxTAtIP(i) - fluxPAtIP(i)
          JgwDAtIP(i) = fluxgAtIP(i) + fluxPAtIP(i)          
        END DO
      ELSE ! nothing at all is computed or read in
        JgwDAtIP(1:DIM) = 0.0_dp
      END IF

      AbsJgwDAtIP = SQRT(SUM(JgwDAtIP(1:DIM)*JgwDAtIP(1:DIM)))
      eL = 0.0_dp
      IF (AbsJgwDAtIP > 0.0_dp) &
           eL(1:DIM) = JgwDAtIP(1:DIM)/AbsJgwDAtIP
      
      Xi0Tilde = GetXi0Tilde(Xi0,mu0,PorosityAtIP)
      XiAtIP = GetXi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0Tilde)
      XixcAtIP = XiXc(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,Mc,ew,dw1,dw2,deltaInElement,GasConstant,SalinityAtIP)
      XietaAtIP = XiEta(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Xi0Tilde,eta0,PorosityAtIP)
      KcAtIP = 0.0_dp
      KcAtIP = GetKc(alphaL,alphaT,Dm0,XiAtIP,AbsJgwDAtIP,eL,PorosityAtIP)
      KcXcXcAtIP = GetKcXcXc(T0,rhoc0,dw1,dw2,dc0,dc1,KcAtIP,TemperatureAtIP,SalinityAtIP,PressureAtIP)

      DO i=1,DIM
        fluxXcG(i) = (1.0_dp - SalinityAtIP) * Mc /(GasConstant * T0 * dc0) * SUM(KcAtIP(1:DIM,i)*Gravity(1:DIM))
      END DO

      TimeWeight = PorosityAtIP * (XiAtIP  + (rhoi0/rhow0) * XixcAtIP * SalinityAtIP)
      ReactionWeight = (rhoi0/rhow0) * PorosityAtIP * XixcAtIP * TemperatureTimeDerAtIP
      Weight = IP % s(t) * DetJ

      ! remove that part
      KcXcXcAtIP = 0.0_dp
      KcXcXcAtIP(1,1) = 0.0001
      KcXcXcAtIP(2,2) = 0.0001
      XiAtIP = 1.0
      PorosityAtIP = 1.0
      TimeWeight = 1.0
      ! ---------------
      DO p=1,nd
        DO q=1,nd
          ! diffusion term (Xi * eta * KcXcXcAtIP.grad(u) - fluxXcG * u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + Weight * XiAtIP * PorosityAtIP * KcXcXcAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
              !PRINT *, "Diff: Stiff = XiAtIP * PorosityAtIP * KcXcXcAtIP(i,j)", Stiff(p,q),XiAtIP,PorosityAtIP,KcXcXcAtIP(i,j)
              !Stiff(p,q) = Stiff(p,q) - Weight * XiAtIP * PorosityAtIP * fluxXcG(i) * Basis(p) * dBasisdx(q,i)
            END DO
          END DO
          ! advection term ( JgwD.grad(u),v)
          ! -----------------------------------
          !IF (.NOT.NoGWFlux .OR. ComputeGWFlux) &
          !     STIFF (p,q) = STIFF(p,q) + Weight * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! reaction term (Rc*u,v)
          ! -----------------------------------
          !STIFF(p,q) = STIFF(p,q) + ReactionWeight * Weight *Basis(q) * Basis(p)

          ! time derivative (TimeWeight du/dt,v):
          ! ------------------------------
          !MASS(p,q) = MASS(p,q) +  TimeWeight * Weight * Basis(q) * Basis(p)
          MASS(p,q) = MASS(p,q) + 0.0001 * Weight * Basis(q) * Basis(p)
          !PRINT *,"Mass=TimeWeight * CGTTAtIP", MASS(p,q), TimeWeight
        END DO
      END DO

      !FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSalinity
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSalinityBC( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), ImposedSalinity, Salinity(n), F,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
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

    Flux(1:n)  = GetReal( BoundaryCondition,'Salinity Flux', FluxCondition )
    Salinity(1:n) = GetReal( BoundaryCondition,'Imposed Salinity', RobinCondition )

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

        ! Robin condition (C*(u-u_0)):
        ! ---------------------------     
        ImposedSalinity = SUM(Basis(1:n)*Salinity(1:n))

        IF (Robincondition) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * ImposedSalinity * Basis(1:nd)
        ELSE IF (Fluxcondition) THEN
          ! Given flux:
          ! -----------
          F = SUM(Basis(1:n)*flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) + Weight * F  * Basis(1:nd)
        END IF
      END DO
    END IF
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSalinityBC
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
END SUBROUTINE PermafrostSalinity
!------------------------------------------------------------------------------

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
  INTEGER :: DIM, i, j, k, NumberOfRecords,RockMaterialID,CurrentNode
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PorosityInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName
  LOGICAL :: Visited = .False., Found, GotIt

  SAVE Visited
  !,DIM,CurrentRockMaterial,NumberOfRecords

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
  PorosityValues = -9999.0_dp
  
  ! Loop over elements
  DO i = 1, Solver % NumberOFActiveElements
    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial()

    ! get RockMaterial pointer
    IF (.NOT.Visited) THEN
      NumberOfRecords =  ReadPermafrostRockMaterial( Material,CurrentRockMaterial )
      IF (NumberOfRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
      END IF
      dim = CoordinateSystemDimension()
      Visited=.True.
    END IF

    RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
    IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    
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
