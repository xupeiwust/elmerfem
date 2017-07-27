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
          cs0(:),Xi0(:),eta0(:),hs0(:),Kgwh0(:,:,:),qexp(:)
     REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         hw0,hi0,cw0,ci0,eps,kw0th,ki0th,mu0
     CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  END TYPE RockMaterial_t

CONTAINS

  FUNCTION ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
       NumberOfRecords,GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
       l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity) RESULT(Constantsread)
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    INTEGER :: DIM, NumberOfRecords
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity(3)
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
    Mw= CurrentRockMaterial % MW
    DeltaT= CurrentRockMaterial % DeltaT
    T0= CurrentRockMaterial % T0
    p0= CurrentRockMaterial % p0
    rhow0= CurrentRockMaterial % rhow0
    rhoi0= CurrentRockMaterial % rhoi0
    cw0= CurrentRockMaterial % cw0
    ci0= CurrentRockMaterial % ci0
    eps=CurrentRockMaterial % eps
    kw0th= CurrentRockMaterial % kw0th
    ki0th= CurrentRockMaterial % ki0th
    mu0 = CurrentRockMaterial % mu0
    l0= (CurrentRockMaterial % hw0) - (CurrentRockMaterial % hi0)
    CgwTT = rhow0*cw0
    ConstantsRead=.TRUE.

    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(FunctionName,"General Constants:", Level=9)
    WRITE(Message,'(A)') "GasConstant,Mw,DeltaT,T0,p0,rhow0,rhoi0,hw0,hi0,cw0,ci0,eps,kw0th,ki0th,mu0:"
    CALL INFO(FunctionName,Message,Level=9)
    WRITE(Message,'(15E12.5)') CurrentRockMaterial % GasConstant, &
         CurrentRockMaterial % Mw, CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0,&
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
      WRITE(Message,'(E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5,A,E10.5)') CurrentRockMaterial % Xi0(I),&
           ",",CurrentRockMaterial % eta0(I), ",", CurrentRockMaterial % Ks0th(I), "," ,&
           CurrentRockMaterial % ew(I), ",", CurrentRockMaterial % bs(I),",", CurrentRockMaterial % rhos0(I),&
           ",",CurrentRockMaterial % cs0(I)
      CALL INFO(FunctionName,Message,Level=9)
    END DO
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9) 

  END FUNCTION ReadPermafrostRockMaterialConstants

  FUNCTION ReadPermafrostRockMaterial(Params,CurrentRockMaterial ) RESULT(NumberOfRecords)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(RockMaterial_t), TARGET :: LocalRockMaterial
    Integer :: NumberOfRecords

    INTEGER :: i,j,k,l, t, active, DIM, ok,InitialNumberOfRecords
    INTEGER,parameter :: io=20,NumberOfEntries=10
    LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE., DataRead=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  MaterialFileName, NewMaterialFileName, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='ReadPermafrostRockMaterial'

    SAVE AllocationsDone,DataRead,InitialNumberOfRecords,LocalRockMaterial,MaterialFileName

    IF (DataRead) THEN
      NewMaterialFileName = GetString( Params, 'Rock Material File', Found )
      IF (NewMaterialFileName /= MaterialFileName) THEN
        WRITE (Message, '(A,A,A,A)') NewMaterialFileName,' does not match existing datafile ', MaterialFileName,'. Exiting!'
        CALL FATAL(FunctionName,Message)
      END IF
      NumberOfRecords = InitialNumberOfRecords
      CurrentRockMaterial => LocalRockMaterial
      RETURN
    ELSE
      DIM = CoordinateSystemDimension()
      !------------------------------------------------------------------------------
      ! Inquire and open file
      !------------------------------------------------------------------------------
      MaterialFileName = GetString( Params, 'Rock Material File', Found )
      IF (.NOT.Found) &
           CALL FATAL(FunctionName," 'Rock Material File' keyword not found - you have to specify an input file!")
      NumberOfRecords = GetInteger( Params, 'Records in File', Found )
      IF (.NOT.Found) &
           CALL FATAL(FunctionName," 'Records in File' keyword not found - you have to specify an integer!")
      OPEN(unit = io, file = TRIM(MaterialFileName), status = 'old',iostat = ok)
      IF (ok /= 0) THEN
        WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
        CALL FATAL(Trim(FunctionName),Trim(message))
      END IF
      WRITE (Message,'(A,I2,A,A)') "Attempting read ",NumberOfRecords," from data file ",TRIM(MaterialFileName)
      CALL INFO(FunctionName,Message,level=9)
      InitialNumberOfRecords = NumberOfRecords

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
           LocalRockMaterial % VariableBaseName(NumberOfRecords),&
           STAT=OK)
      AllocationsDone = .TRUE.
      DataRead = .TRUE.
      IF (OK /= 0) &
           CALL FATAL(FunctionName, 'Allocation Error of input data array')
      !------------------------------------------------------------------------------
      ! Read in information from material file
      ! General constants
      !  GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
      !     hw0,hi0,cw0,ci0,eps,kw0th,ki0th,mu0
      !------------------------------------------------------------------------------
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % GasConstant, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % Mw, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % DeltaT,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % T0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % p0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhow0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % rhoi0,Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % hw0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % hi0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % cw0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ci0, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % eps, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % kw0th, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % ki0th, Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % mu0, Comment

      ! for each material
      !       ks0th(:),ew(:),b(:),rhos0(:),cs0(:)
      DO I=1,NumberOfRecords
        READ (io, *, END=10, IOSTAT=OK, ERR=30) LocalRockMaterial % VariableBaseName(I)
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
        WRITE(Message,'(A)') "GasConstant,Mw,DeltaT,T0,p0,rhow0,rhoi0,hw0,hi0,cw0,ci0,eps,kw0th,ki0th,mu0:"
        CALL INFO(FunctionName,Message,Level=9)
        WRITE(Message,'(15E12.5)') CurrentRockMaterial % GasConstant, &
             CurrentRockMaterial % Mw, CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0,&
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

  RECURSIVE REAL FUNCTION delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant
    delta = 0.5_dp*l0*DeltaT/T0 &
         + (cw0 - ci0)*((T0 + 0.5_dp*DeltaT)*LOG(1.0_dp + 0.5_dp*DeltaT/T0) - 0.5_dp*DeltaT)
    delta = delta*(eps*(1.0_dp - eps)/(2.0_dp*eps -1.0_dp))* Mw/(GasConstant*(T0 + 0.5_dp*DeltaT))
  END FUNCTION delta

  RECURSIVE REAL FUNCTION deltaG(ew,eps,DeltaT,T0,p0,Mw,l0,cw0,ci0,rhow0,rhoi0,&
       GasConstant,Temperature,Pressure)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ew,eps,DeltaT,T0,p0,Mw,l0,cw0,ci0,rhow0,rhoi0,&
         GasConstant,Temperature,Pressure
    deltaG = -l0*(Temperature - T0)/T0 &  ! the first one is L Zero
         + ((1.0_dp/rhow0) + (1.0_dp/rhoi0))*(Pressure - p0) &
         - (cw0 - ci0)*(Temperature*LOG(Temperature/T0) - (Temperature - T0))
  END FUNCTION deltaG

  RECURSIVE REAL FUNCTION B1(delta,ew,Mw,GasConstant,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,ew,Mw,GasConstant,Temperature
    B1 = 1.0_dp/(ew + delta)
    B1 = B1*Mw/(GasConstant*Temperature)
  END FUNCTION B1

  RECURSIVE REAL FUNCTION B2(delta,deltaG,GasConstant,Mw,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,deltaG,Mw,GasConstant,Temperature
    B2 = Mw*deltaG/(GasConstant*Temperature*delta)
  END FUNCTION B2

  RECURSIVE REAL FUNCTION D1(delta,ew)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: delta,ew
    ! local
    D1 = delta/(ew + delta)
  END FUNCTION D1

  RECURSIVE REAL FUNCTION Xi(B1,B2,D1,D2,Xi0)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0
    Xi= Xi0/(1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1)) &
         + (1.0_dp - Xi0)/(1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))
  END FUNCTION Xi

  RECURSIVE REAL FUNCTION XiT(B1,B2,D1,D2,Xi0,p0,Mw,ew,delta,rhow0,rhoi0,cw0,ci0,&
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

  RECURSIVE REAL FUNCTION XiP(B1,B2,D1,D2,Xi0,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2
    IF (Temperature .LE. 0.0_dp) CALL FATAL("Permafrost (XiP)","(sub-)Zero Temperature detected")
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))**2.0_dp)
    XiP = (0.5_dp*Xi0*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0)*aux2/delta) &
         *((1.0_dp/rhoi0) - (1.0_dp/rhow0))* Mw/(GasConstant*Temperature)
  END FUNCTION XiP

  RECURSIVE REAL FUNCTION ksth(ks0th,bs,T0,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ks0th,bs,T0,Temperature
    ksth = ks0th/( 1.0_dp + bs*(Temperature - T0)/T0)
  END FUNCTION ksth

  RECURSIVE FUNCTION GetKGTT(ks0th,kw0th,ki0th,eta0,Xi,Temperature,Pressure,Porosity,meanfactor)RESULT(KGTT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ks0th,kw0th,ki0th,eta0,Xi,Temperature,Pressure,Porosity,meanfactor
    ! local
    REAL(KIND=dp) :: KGTT(3,3), factor,unittensor(3,3)
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    factor = (1.0_dp - meanfactor)*(1.0_dp/((1.0_dp - eta0)/ks0th + Xi*eta0/kw0th + (1.0_dp - Xi)*eta0/ki0th))
    factor = factor + meanfactor*((1.0_dp - eta0)*ks0th + Xi*eta0*kw0th + (1.0_dp - Xi)*eta0*ki0th)
    KGTT = unittensor*factor
  END FUNCTION GetKGTT

  RECURSIVE REAL FUNCTION CGTT(Xi,XiT,rhos0,rhow0,rhoi0,cw0,ci0,cs0,l0,eta0)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,XiT,rhos0,rhow0,rhoi0,cw0,ci0,cs0,l0,eta0
    CGTT = (1.0_dp - eta0)*rhos0*cs0 + Xi*eta0*rhow0*cw0 &
         + (1.0_dp - Xi)*eta0*rhoi0*ci0 &
         + rhoi0*l0*eta0*XiT
  END FUNCTION CGTT

  RECURSIVE REAL FUNCTION fTildewT(B1,Temperature,D1,delta,ew,l0,cw0,ci0,T0,Xi,Xi0)
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

  RECURSIVE REAL FUNCTION fTildewp(B1,D1,delta,ew,rhow0,rhoi0,Xi,Xi0)
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
  
  RECURSIVE FUNCTION GetKgw(mu0,mu,Xi,rhow0,qexp,Kgwh0,MinKgw)RESULT(Kgw)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: mu0,mu,Xi,rhow0,qexp,Kgwh0(3,3),MinKgw
    REAL(KIND=dp) :: Kgw(3,3), factor
    REAL(KIND=dp), PARAMETER :: gval=9.81 !hard coded, so match Kgwh0 with this value
    factor = (mu0/mu)*(Xi**qexp)/(rhow0*gval)
    Kgw = MAX(Kgwh0*factor,MinKgw)
    !PRINT *,"factor",factor,"gval",gval,"rhow0",rhow0,"Xi",Xi,"qexp",qexp
  END FUNCTION GetKgw
  
  RECURSIVE FUNCTION GetKgwpT(rhow0,fTildewT,Kgw)RESULT(KgwpT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,fTildewT,Kgw(3,3)
    REAL(KIND=dp) :: KgwpT(3,3), factor
    factor = rhow0*fTildewT
    KgwpT = Kgw*factor
  END FUNCTION GetKgwpT

  RECURSIVE FUNCTION GetKgwpp(rhow0,fTildewp,Kgw)RESULT(Kgwpp)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow0,fTildewp,Kgw(3,3)
    REAL(KIND=dp) :: Kgwpp(3,3), factor
    factor = (1.0_dp + rhow0*fTildewp)
    Kgwpp = Kgw*factor
  END FUNCTION GetKgwpp
END MODULE PermafrostMaterials



!-----------------------------------------------------------------------------
!> heat transfer equation for enhanced permafrost model
!-----------------------------------------------------------------------------
SUBROUTINE PermafrostHeatEquation( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,GWfluxVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, Active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20,NumberOfEntries=10
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),PorosityPerm(:),SalinityPerm(:),GWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),GWflux(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
        NodalTemperature(:),NodalGWflux(:,:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,NoGWflux=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,NumberOfRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  IF (FirstTime) THEN
    NumberOfRecords =  ReadPermafrostRockMaterial( Params,CurrentRockMaterial )
    IF (NumberOfRecords < 1) CALL FATAL(SolverName,'No Rock Material specified')
    FirstTime = .FALSE.
  END IF
  Temperature => Solver % Variable % Values
  TemperaturePerm => Solver % Variable % Perm
  IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
    N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    IF (AllocationsDone) &
         DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,NodalSalinity,NodalGWflux)
    ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),&
         NodalSalinity(N),NodalGWflux(3,N),STAT=istat )
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
  PressureName = ListGetString(Params, &
       'Pressure Variable', Found )
  IF (.NOT.Found) THEN
    CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
    WRITE(PressureName,'(A)') 'Pressure'
  ELSE
    WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
    CALL INFO(SolverName,Message,Level=3)
  END IF
  PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
  IF (.NOT.ASSOCIATED(PressureVar)) THEN
    NULLIFY(Pressure)
    NoPressure = .TRUE.
    WRITE(Message,'(A,A,A)') "'Pressure Variable ", PressureName, " not associated"
    CALL WARN(SolverName,Message)
  ELSE
    Pressure => PressureVar % Values
    PressurePerm => PressureVar % Perm
    NoPressure = .FALSE.
    WRITE(Message,'(A,A,A)') "'Pressure Variable ", PressureName, " associated"
    CALL INFO(SolverName,Message,Level=3)
  END IF

  PorosityName = ListGetString(Params, &
       'Porosity Variable', Found )
  IF (.NOT.Found) THEN
    CALL WARN(SolverName," 'Porosity Variable' not found. Using default 'Porosity' ")
    WRITE(PorosityName,'(A)') 'Porosity'
  ELSE
    WRITE(Message,'(A,A)') "'Porosity Variable' found and set to: ", PorosityName
    CALL INFO(SolverName,Message,Level=3)
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
    CALL INFO(SolverName,Message,Level=3)
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
    CALL WARN(SolverName," 'Groundwater flux Variable' not found. Using default 'Groundwater flux' ")
    WRITE(GWfluxName,'(A)') 'Groundwater flux'
  ELSE
    WRITE(Message,'(A,A)') "'Groundwater flux Variable' found and set to: ", GWfluxName
    CALL INFO(SolverName,Message,Level=3)
  END IF
  GWfluxVar => VariableGet(Solver % Mesh % Variables,GWfluxName)
  IF (.NOT.ASSOCIATED(GWfluxVar)) THEN
    CALL WARN(SolverName,'Groundwater flux Variable not found. Using Pressure and Temperature to compute flux')
    NoGWflux = .TRUE.
  ELSE
    GWflux => GWfluxVar % Values
    GWfluxPerm => GWfluxVar % Perm
    NoGWflux = .FALSE.
    CALL INFO(SolverName,'Groundwater flux Variable not found. Using this as prescribed groundwater flux',Level=4)
  END IF
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    !------------------------------------------------------------------------------
    Active = Solver % NumberOfActiveElements
    DO t=1,Active
      Element => GetActiveElement(t,Solver)
      IF (.NOT.ASSOCIATED(Element)) CYCLE
      ! cycle halo elements
      !-------------------
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      Material => GetMaterial()
      IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
        CALL FATAL(SolverName,Message)
      END IF
      N  = GetElementNOFNodes()
      ND = GetElementNOFDOFs()
      NB = GetElementNOFBDOFs()

      ! Nodal variable dependencies
      NodalTemperature(1:N) = Temperature(TemperaturePerm(Element % NodeIndexes(1:N)))
      !PRINT *, NodalTemperature(1:N), TemperaturePerm(Element % NodeIndexes(1:N)), Element % NodeIndexes(1:N)
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
      IF (NoGWflux) THEN
        NodalGWflux(1:3,1:N) = 0.0_dp
      ELSE
        DO I=1,GWfluxVar % DOFs
          NodalGWflux(I,1:N) = &
               GWflux(GWfluxVar % DOFs * GWfluxPerm(Element % NodeIndexes(1:N)) - I)
        END DO
      END IF

      CALL LocalMatrixHTEQ(  Element, N, ND+NB, NodalTemperature, NodalPressure, &
           NodalPorosity, NodalSalinity, NodalGWflux, NoGWflux, CurrentRockMaterial)
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
        CALL LocalMatrixBCHTEQ(  Element, n, nd+nb )
      END IF
    END DO
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT

  END DO

CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements for
  !    Permafrost Heat Transfer Equation

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHTEQ( Element, n, nd, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, NodalGWflux, ComputeGWFlux, CurrentRockMaterial )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL :: ComputeGWFlux
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT    ! constants read only once
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,GWfluxAtIP(3),StiffPQ
    REAL(KIND=DP) :: gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),pFluxAtIP(3),gFlux(3),Gravity(3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
    
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
           l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Heat source', Found )

    ! read variable material parameters from CurrentRockMaterial
    Material => GetMaterial()
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw .LE. 0.0_dp))  &
      MinKgw = 1.0D-14
      

    ! read element rock material specific parameters    
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    ew = CurrentRockMaterial % ew(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)

!    PRINT *, "ks0th", ks0th,"ew", ew, "bs",bs, "rhos0", rhos0, "cs0", cs0,&
!         "Xi0",Xi0, "eta0",eta0,"Kgwh0", Kgwh0(1:3,1:3)

    ! derive element rock material specific parameters
    deltaInElement = delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    D1InElement = D1(deltaInElement,ew)
    D2InElement = 1.0_dp

    !PRINT *, "deltaInElement",deltaInElement,"D1InElement",D1InElement,"D2InElement",D2InElement

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


      ! functions at IP
      deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,l0,cw0,ci0,rhow0,rhoi0,GasConstant,&
           TemperatureAtIP,PressureAtIP)
      B1AtIP = B1(deltaInElement,ew,Mw,GasConstant,TemperatureAtIP)
      B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
      XiAtIP = Xi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0)
      XiTAtIP= XiT(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,p0,Mw,ew,&
           deltaInElement,rhow0,rhoi0,cw0,ci0,l0,T0,GasConstant,TemperatureAtIP,PressureAtIP)
      XiPAtIP= XiP(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0,Mw,ew,&
           deltaInElement,rhow0,rhoi0,GasConstant,TemperatureAtIP)
      ksthAtIP = ksth(ks0th,bs,T0,TemperatureAtIP)
      KGTTAtIP = GetKGTT(ksthAtIP,kw0th,ki0th,eta0,XiAtIP,&
           TemperatureAtIP,PressureAtIP,PorosityAtIP,meanfactor)
      CGTTAtIP = CGTT(XiAtIP,XiTAtIP,rhos0,rhow0,rhoi0,cw0,ci0,cs0,l0,eta0)
      IF (.NOT.ComputeGWFlux) THEN
        DO I=1,3
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N))
        END DO
      ELSE
        fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
        fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
        KgwpTAtIP = GetKgwpT(rhow0,fTildewTATIP,KgwAtIP)
        !PRINT *,"KgwAtIP",KgwAtIP
        !PRINT *, fTildewTATIP, rhow0
        KgwppAtIP = GetKgwpp(rhow0,fTildewpATIP,KgwAtIP)
        gradTAtIP = 0.0_dp
        gradPAtIP = 0.0_dp
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO
        DO i=1,3
          fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:3)*gradTAtIP(1:3))
          gFlux(i) = rhow0 * SUM(Kgwh0(i,1:3)*Gravity(1:3))
          pFluxAtIP(i) =  SUM(KgwppAtIP(i,1:3)*gradPAtIP(1:3))
          !JgwDAtIP(i) = gFlux(i) - fluxTAtIP(i) - pFluxAtIP(i)
          JgwDAtIP(i) = gFlux(i) - pFluxAtIP(i)
          !JgwDAtIP(i) = 0.0_dp
        END DO
        !PRINT *,"pFluxAtIP(1:DIM) =", pFluxAtIP(1:3)
        !PRINT *,"KgwppAtIP=", KgwppAtIP(1:3,1:3)
!        PRINT *,"fluxTAtIP(1:DIM) =", fluxTAtIP(1:DIM)
!        PRINT *,"gFlux(1:DIM) =", gFlux(1:DIM) 
      END IF
      !fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiPAtIP,Xi0)
      !fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiPAtIP,Xi0)
      !KgwpTAtIP = KgwpT(rhow0,fTildewTATIP,Kgwh0)


      !PRINT *,"KGTTAtIP",KGTTAtIP
      !PRINT *,"CGTTAtIP",CGTTAtIP

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      DO p=1,nd
        DO q=1,nd
          StiffPQ = 0.0
          ! groundwater advection term (C*grad(u),v)
          ! C_GW^TT dT/dx_i J_gw^D_i
          ! -----------------------------------
          !StiffPQ = StiffPQ + &
          !     CGWTT * SUM(JgwDAtIP(1:DIM)*dBasisdx(q,1:DIM)) * Basis(p)

          ! diffusion term ( grad(u),grad(v))
          ! div(JGH) = d(KGTT_i,j dT/dx_j)/dx_i
          ! -----------------------------------
          DO i=1,DIM
            DO J=1,DIM
              StiffPQ = StiffPQ + KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ

          ! 

          ! time derivative (c*du/dt,v): !! THIS IS OK AS IS !!!
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * (CGTTAtIP) * Basis(q) * Basis(p)
        END DO
      END DO
      ! body force
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixHTEQ
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !  for Permafrost Heat Transfer equation
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCHTEQ( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,ConstantsRead=.FALSE.
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    
    SAVE Nodes
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'Heat Flux', Found )
    Coeff(1:n) = GetReal( BC,'Heat Transfer Coefficient', Found )
    IF (.NOT.Found) Coeff(1:n) = 0.0_dp
    Ext_t(1:n) = GetReal( BC,'External Temperature', Found )

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

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO
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
END SUBROUTINE PermafrostHeatEquation

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
  INTEGER,PARAMETER :: io=20,NumberOfEntries=10
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:)!,GWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:)!,GWflux(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
        NodalPressure(:) !NodalGWflux(:,:),
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.!, NoGWflux=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, VarName 

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,&
       NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure!,NodalGWflux,
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  IF (FirstTime) THEN
    NumberOfRecords =  ReadPermafrostRockMaterial( Params,CurrentRockMaterial )
    IF (NumberOfRecords < 1) CALL FATAL(SolverName,'No Rock Material specified')
    FirstTime = .FALSE.
  END IF
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
    CALL INFO(SolverName,Message,Level=3)
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
    CALL INFO(SolverName,Message,Level=3)
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
    CALL INFO(SolverName,Message,Level=3)
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
    !------------------------------------------------------------------------------
    Active = Solver % NumberOfActiveElements
    DO t=1,Active
      Element => GetActiveElement(t,Solver)
      IF (.NOT.ASSOCIATED(Element)) CYCLE
      ! cycle halo elements
      !-------------------
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      Material => GetMaterial()
      IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
        CALL FATAL(SolverName,Message)
      END IF
      N  = GetElementNOFNodes()
      ND = GetElementNOFDOFs()
      NB = GetElementNOFBDOFs()

      ! Nodal variable dependencies
      NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
      !PRINT *, NodalPressure(1:N), PressurePerm(Element % NodeIndexes(1:N)), Element % NodeIndexes(1:N)
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


      CALL LocalMatrixDarcy(  Element, N, ND+NB, NodalPressure, NodalTemperature, &
           NodalPorosity, NodalSalinity, CurrentRockMaterial)!, NodalGWflux
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

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT

    
  END DO

CONTAINS
  ! PermafrostGroundWaterFlow : Assembly of the matrix entries arising from the bulk elements 

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixDarcy( Element, n, nd,  NodalPressure, NodalTemperature, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
          NodalPorosity(:), NodalPressure(:)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),gFlux(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: fTildewTAtIP, fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT, Gravity(3)    ! constants read only once
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
         
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(Kind=dp) , POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow', &
         FunctionName='Permafrost (LocalMatrixDarcy)'

    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT, Gravity
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = ReadPermafrostRockMaterialConstants(Model, FunctionName, CurrentRockMaterial, DIM, &
           NumberOfRecords,GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
           l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) THEN
      LOAD(1:n) = GetReal( BodyForce,'Groundwater source', Found )   
    END IF

    ! read variable material parameters from CurrentRockMaterial
    Material => GetMaterial()
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw .LE. 0.0_dp)) MinKgw = 1.0D-15
    
    ! read element rock material specific parameters    
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    ew = CurrentRockMaterial % ew(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
    qexp = CurrentRockMaterial % qexp(RockMaterialID)

    !PRINT *, "ks0th", ks0th,"ew", ew, "bs",bs, "rhos0", rhos0, "cs0", cs0,&
    !     "Xi0",Xi0, "eta0",eta0,"Kgw", Kgwh0(1:3,1:3)

    ! derive element rock material specific parameters
    deltaInElement = delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    D1InElement = D1(deltaInElement,ew)
    D2InElement = 1.0_dp

    !PRINT *, "deltaInElement",deltaInElement,"D1InElement",D1InElement,"D2InElement",D2InElement
    
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
      !PRINT *, NodalTemperature(1:N),N
     
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))


      ! functions at IP
      deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,l0,cw0,ci0,rhow0,rhoi0,GasConstant,&
           TemperatureAtIP,PressureAtIP)
      B1AtIP = B1(deltaInElement,ew,Mw,GasConstant,TemperatureAtIP)
      !PRINT *,"B1(",deltaInElement,ew,Mw,GasConstant,TemperatureAtIP,")"
      B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
      XiAtIP = Xi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0)
      fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
      fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
      !DO i=1,DIM
      !  KgwAtIP(i,i) = 9.81d-04  ! REMOVE *************************************************************
      !END DO
      KgwpTAtIP = GetKgwpT(rhow0,fTildewTATIP,KgwAtIP)
      KgwppAtIP = GetKgwpp(rhow0,fTildewpATIP,KgwAtIP)
      !PRINT *,"KgwAtIP",KgwAtIP(1,1),'KgwpTAtIP',KgwpTAtIP(1,1),"KgwppAtIP",KgwppAtIP(1,1),&
      !     "fTildewpATIP",fTildewpATIP,"XiAtIP",XiAtIP,"qexp",qexp,"rhow0",rhow0,"Gravity",Gravity(2
      !PRINT *,"KgwAtIP",KgwAtIP,"Gravity",Gravity
      
      !PRINT *,"ftildewp(",B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0,")"

      !PRINT *,"Kgw",Kgw
      !PRINT *,"KgwppAtIP",KgwppAtIP

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
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
        END DO
      END DO
       
      ! body forcexs
      DO i=1,DIM
        gradTAtIP(i) =  SUM(NodalTemperature(1:n)*dBasisdx(1:n,i))
      END DO
      DO i=1,DIM
        fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
        gFlux(i) = rhow0 * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
      END DO
      DO p=1,nd     
        !FORCE(p) = FORCE(p) - Weight * SUM(fluxTAtIP(1:DIM)*dBasisdx(p,1:DIM))
        FORCE(p) = FORCE(p) + Weight * SUM(gFlux(1:DIM)*dBasisdx(p,1:DIM))
      END DO
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixDarcy
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCDarcy( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Pressure(n),Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,ConstantsRead=.FALSE.,FluxCondition,DirichletCondition
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    
    TYPE(Element_t),POINTER :: ParentElement
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3), MinKgw, KgwppAtIP(3,3), fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT    ! constants read only once
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         ngFlux
    REAL(KIND=DP) :: gradTAtIP(3),gradPAtIP(3),fluxTAtIP(3),pFluxAtIP(3),gFlux(3),Gravity(3),Normal(3)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: RockMaterialID,other_body_id

    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT, Gravity
    !------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'Groundwater Flux', FluxCondition )

    ! Need to access Parent Element to get Material properties
    other_body_id = Element % BoundaryInfo % outbody
    IF (other_body_id < 1) THEN ! only one body in calculation
      ParentElement => Element % BoundaryInfo % Right
      IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => Element % BoundaryInfo % Left
    ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
      ParentElement => Element %  BoundaryInfo % Right
      IF (ParentElement % BodyId == other_body_id) ParentElement =>  Element % BoundaryInfo % Left
    END IF
    Material => GetMaterial(ParentElement)
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw .LE. 0.0_dp))  &
         MinKgw = 1.0D-14
    Pressure(1:n) = GetReal( BC,'Imposed '// TRIM(VarName), DirichletCondition)
    
    ! read element rock material specific parameters    
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    ew = CurrentRockMaterial % ew(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 = CurrentRockMaterial % eta0(RockMaterialID)
    Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      IF (Fluxcondition) THEN
        ! Given flux:
        ! -----------
        F = SUM(Basis(1:n)*flux(1:n))
        FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
      END IF
      ! Additional terms to be set if there is an explicit Dirichlet condition
      !----------------------------------------------------------------------
      IF (DirichletCondition) THEN
        !CALL GetElementNodes( EdgeNodes, Element )
        Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .FALSE. )
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
        DO i=1,DIM
          gFlux(i) = rhow0 * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
        END DO
        ngFlux =  SUM(gFlux(1:DIM)*Normal(1:DIM))
        PressureAtIP = SUM(Pressure(1:n)*Basis(1:n))
        C = 10000.0_dp ! super high transfer coefficient to ensure pressure value in weak formulation
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
          FORCE(1:nd) = Force(1:nd) + Weight * C * PressureAtIP * Basis(1:nd)
          FORCE(1:nd) = Force(1:nd) + Weight * C * ngflux * Basis(1:nd)
        END DO
      END IF
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE)
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
  !LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry, GotIt
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
  END TYPE FieldTable_t
  TYPE(FieldTable_t) :: Fields(3)

  TYPE(Variable_t), POINTER :: PressureVar,TemperatureVar,PorosityVar,SalinityVar
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:)
  INTEGER :: NumberOfRecords
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:)
  LOGICAL :: ConstantPorosity, NoSalinity, FirstTime=.TRUE.,UnfoundFatal=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, PressureName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux"

  
  SAVE FirstTime,CurrentRockMaterial,NumberOfRecords
  
  CALL Info( SolverName, '-------------------------------------',Level=4 )
  CALL Info( SolverName, 'Computing the groundwater flux       ',Level=4 )
  CALL Info( SolverName, '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!  Check what needs to be computed
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()
  Dofs = Dim
  
  ! Read Variables
  CALL ReadVars()

  IF (FirstTime) THEN
    NumberOfRecords =  ReadPermafrostRockMaterial( SolverParams,CurrentRockMaterial )
    IF (NumberOfRecords < 1) CALL FATAL(SolverName,'No Rock Material specified')
    FirstTime = .FALSE.
  END IF
!-------------------------------------------------------------------------------
! If only one component is used use the scalar equation, otherwise use an
! auxiliary variable to store all the dimensions
!-------------------------------------------------------------------------------
  Varname = TRIM('Groundwater')

  i = 0
  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 1',UnFoundFatal=UnFoundFatal )
  Fields(1) % Values => FluxSol % Values

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 2',UnFoundFatal=UnFoundFatal )
  Fields(2) % Values => FluxSol % Values
    
  IF( dim == 3 ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 3',UnFoundFatal=UnFoundFatal )
    Fields(3) % Values => FluxSol % Values
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
    WRITE(Message,'(A,I1,A,I1)') "Working on Dof ",i," out of ",Dofs
    CALL INFO(SolverName,Message,Level=3)
    Solver % Matrix % RHS => ForceVector(:,i)
    UNorm = DefaultSolve()

    TotNorm = TotNorm + Unorm ** 2
    Fields(i) % Values = Solver % Variable % Values
  END DO
  
  DEALLOCATE( ForceVector )  
  Solver % Matrix % RHS => SaveRHS
  TotNorm = SQRT(TotNorm)
  Solver % Variable % Norm = Totnorm


!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( SolverName, Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( SolverName, Message, Level=4 )

  CALL Info( SolverName, 'All done',Level=6 )
  CALL Info( SolverName, '-------------------------------------',Level=6 )


  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
    !------------------------------------------------------------------------------

    INTEGER :: elem,t,i,j,k,p,q,n,nd, DIM,Rank, RockMaterialID
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,C(3,3),coeff,detJ,GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: LocalPotential(:)
    LOGICAL :: Found, ConstantsRead=.FALSE.
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), POINTER :: Conductivity(:,:,:)=>NULL()
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT, Gravity(3)    ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         fluxTAtIP(3),fluxpAtIP(3),gFlux(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: fTildewTAtIP, fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp  ! stuff comming from RockMaterial
    REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:), NodalSalinity(:), NodalPressure(:),NodalPorosity(:)
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP

    SAVE Conductivity, Nodes, ConstantsRead,DIM,&
         GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity

    n = 2*MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dofs,n) )
    ALLOCATE( LocalPotential(n), Basis(n), dBasisdx(n,3) )
    ALLOCATE( NodalPressure(N),NodalPorosity(N),NodalTemperature(N),NodalSalinity(N) )

    DO elem = 1,Solver % NumberOFActiveElements

      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      Material => GetMaterial()

      IF(.NOT.ConstantsRead) THEN
        ConstantsRead = ReadPermafrostRockMaterialConstants(Model, SolverName, CurrentRockMaterial, DIM, &
             NumberOfRecords,GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
             l0,cw0,ci0,eps,kw0th,ki0th,mu0,CgwTT,Gravity)
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
      RockMaterialID =&
           ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
      MinKgw = GetConstReal( Material, &
           'Hydraulic Conductivity Limit', GotIt)      
      IF (.NOT.GotIt .OR. (MinKgw .LE. 0.0_dp)) MinKgw = 1.0D-15

      ! read element rock material specific parameters    
      ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
      ew = CurrentRockMaterial % ew(RockMaterialID)
      bs = CurrentRockMaterial % bs(RockMaterialID)
      rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
      cs0 = CurrentRockMaterial % cs0(RockMaterialID)
      Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
      eta0 = CurrentRockMaterial % eta0(RockMaterialID)
      Kgwh0(1:3,1:3) = CurrentRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
      qexp = CurrentRockMaterial % qexp(RockMaterialID)

      ! derive element rock material specific parameters
      deltaInElement = delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
      D1InElement = D1(deltaInElement,ew)
      D2InElement = 1.0_dp

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

        ! functions at IP
        deltaGAtIP = deltaG(ew,eps,DeltaT,T0,p0,Mw,l0,cw0,ci0,rhow0,rhoi0,GasConstant,&
             TemperatureAtIP,PressureAtIP)
        B1AtIP = B1(deltaInElement,ew,Mw,GasConstant,TemperatureAtIP)
        B2AtIP = B2(deltaInElement,deltaGAtIP,GasConstant,Mw,TemperatureAtIP)
        XiAtIP = Xi(B1AtIP,B2AtIP,D1InElement,D2InElement,Xi0)
        fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0,XiAtIP,Xi0)
        fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0,XiAtIP,Xi0)

        KgwAtIP = 0.0_dp
        KgwpTAtIP = 0.0_dp
        KgwppAtIP = 0.0_dp        
        KgwAtIP = GetKgw(mu0,mu0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw)
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
        gFlux = 0.0_dp
        DO i=1,dim
          fluxpAtIP(i) = SUM(KgwppAtIP(i,1:DIM)*gradpAtIP(1:DIM))
          fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
          gFlux(i) = rhow0 * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
        END DO
 
        DO i=1,dim
          Coeff = -1.0_dp * Weight * (FluxpAtIP(i) - gFlux(i) + fluxTAtIP(i))
          FORCE(i,1:nd) = FORCE(i,1:nd) + Coeff * Basis(1:nd)
        END DO
      END DO

      !------------------------------------------------------------------------------
      !      Update global matrices from local matrices 
      !------------------------------------------------------------------------------
      !      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
      Solver % Matrix % Rhs => SaveRhs
      CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      !      END IF

      DO i=1,Dofs
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO

    END DO

    DEALLOCATE( LocalPotential, STIFF, FORCE, Basis, dBasisdx )

    !------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


  SUBROUTINE ReadVars()
    ! find variables for dependencies
    !--------------------------------
    TemperatureName = ListGetString(SolverParams, &
         'Temperature Variable', GotIt )
    IF (.NOT.GotIt) THEN
      CALL WARN(SolverName," 'Temperature Variable' not found. Using default 'Temperature' ")
      WRITE(TemperatureName,'(A)') 'Temperature'
    ELSE
      WRITE(Message,'(A,A)') "'Temperature Variable' found and set to: ", TemperatureName
      CALL INFO(SolverName,Message,Level=3)
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
      CALL INFO(SolverName,Message,Level=3)
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
      CALL INFO(SolverName,Message,Level=3)
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
      CALL INFO(SolverName,Message,Level=3)
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
  END SUBROUTINE ReadVars
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
  

