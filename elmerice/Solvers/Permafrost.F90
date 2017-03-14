MODULE PermafrostMaterials
  USE Types
  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE

  TYPE RockMaterial_t
     INTEGER :: NumberOfRecords
     REAL(KIND=dp), ALLOCATABLE :: ks0th(:),ew(:),bs(:),rhos0(:),cs0(:),Xi0(:),eta0(:),hs0(:),Kgw(:,:,:)
     REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         hw0,hi0,cw0,ci0,eps,kw0th,ki0th
     CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  END TYPE RockMaterial_t

CONTAINS

  FUNCTION ReadPermafrostRockMaterial( Element,Params,CurrentRockMaterial ) RESULT(NumberOfRecords)
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: Params
    TYPE(RockMaterial_t) :: CurrentRockMaterial
    Integer :: NumberOfRecords
    
    INTEGER :: i,j,k,l, t, active, DIM, ok
    INTEGER,parameter :: io=20,NumberOfEntries=10
    LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  MaterialFileName, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='ReadPermafrostRockMaterial'

    SAVE AllocationsDone

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
      WRITE(message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
      CALL FATAL(Trim(FunctionName),Trim(message))
    END IF
    WRITE (Message,'(A,I2,A,A)') "Attempting read ",NumberOfRecords," from data file ",TRIM(MaterialFileName)
    CALL INFO(FunctionName,Message,level=9)
    !------------------------------------------------------------------------------
    ! Allocate and read stuff
    !------------------------------------------------------------------------------
    !M = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) THEN
      DEALLOCATE(&
           CurrentRockMaterial % ks0th,&
           CurrentRockMaterial % ew,&
           CurrentRockMaterial % bs,&
           CurrentRockMaterial % rhos0,&
           CurrentRockMaterial % cs0,&
           CurrentRockMaterial % Xi0,&
           CurrentRockMaterial % eta0,&
           CurrentRockMaterial % hs0,&
           CurrentRockMaterial % VariableBaseName)
    END IF
    ALLOCATE(&
         CurrentRockMaterial % ks0th(NumberOfRecords),&
         CurrentRockMaterial % ew(NumberOfRecords),&
         CurrentRockMaterial % bs(NumberOfRecords),&
         CurrentRockMaterial % rhos0(NumberOfRecords),&
         CurrentRockMaterial % cs0(NumberOfRecords),&
         CurrentRockMaterial % Xi0(NumberOfRecords),&
         CurrentRockMaterial % eta0(NumberOfRecords),&
         CurrentRockMaterial % hs0(NumberOfRecords),&
         CurrentRockMaterial % Kgw(3,3,NumberOfRecords),&
         CurrentRockMaterial % VariableBaseName(NumberOfRecords),&
         STAT=OK)
    AllocationsDone = .TRUE.
    IF (OK /= 0) &
         CALL FATAL(FunctionName, 'Allocation Error of input data array')
    !------------------------------------------------------------------------------
    ! Read in information from material file
    ! General constants
    !  GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
    !     hw0,hi0,cw0,ci0,eps,kw0th,ki0th
    !------------------------------------------------------------------------------
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % GasConstant, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % Mw, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % DeltaT,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % T0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % p0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhow0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhoi0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % hw0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % hi0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % cw0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ci0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % eps, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % kw0th, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ki0th, Comment

    ! for each material
    !       ks0th(:),ew(:),b(:),rhos0(:),cs0(:)
    DO I=1,NumberOfRecords
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % VariableBaseName(I)
      WRITE(Message,'(A,I2,A,A)') "Input for Variable No.",I,": ", CurrentRockMaterial % VariableBaseName(I)
      CALL INFO(FunctionName,Message,Level=9)
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % Xi0(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % eta0(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ks0th(I), Comment          
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ew(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % bs(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhos0(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % cs0(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % hs0(I), Comment
      DO J=1,3
        DO K=1,3
          READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % Kgw(J,K,I), Comment
        END DO
      END DO
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % hs0(I), Comment      
    END DO
    WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRecords," rock material records from file ", TRIM(MaterialFileName)
    CALL INFO(FunctionName,Message,Level=1)
10  CLOSE(io)
    IF (I < NumberOfRecords) THEN
      WRITE(Message,'(I3,A,I3)') I,"records read, which is smaller than given number ", NumberOfRecords
      CALL FATAL(FunctionName,Message)
    ELSE
      WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRecords," rock material records from file ", TRIM(MaterialFileName)
      CALL INFO(FunctionName,Message,Level=1)
      CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
      CALL INFO(FunctionName,"General Constants:", Level=9)
      WRITE(Message,'(A)') "GasConstant,Mw,DeltaT,T0,p0,rhow0,rhoi0,hw0,hi0,cw0,ci0,eps,kw0th,ki0th:"
      CALL INFO(FunctionName,Message,Level=9)
      WRITE(Message,'(14E12.5)') CurrentRockMaterial % GasConstant, &
           CurrentRockMaterial % Mw, CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0,&
           CurrentRockMaterial % p0, CurrentRockMaterial % rhow0, CurrentRockMaterial % rhoi0,&           
           CurrentRockMaterial % hw0, CurrentRockMaterial % hi0, CurrentRockMaterial % cw0,&
           CurrentRockMaterial % ci0, CurrentRockMaterial % eps, CurrentRockMaterial % kw0th,&
           CurrentRockMaterial % ki0th
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
    deltaG = -l0*(Temperature - T0)/T0 &
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
    aux3 = (l0 + (cw0 - ci0)*(Temperature - T0) &
         + (-(1.0_dp/rhoi0) + (1.0_dp/rhow0))*(Pressure - p0))/Temperature
    XiT = (0.5_dp*Xi0*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0)*aux2/delta) *Mw*aux3/(T0*GasConstant*Temperature)
  END FUNCTION XiT

  RECURSIVE REAL FUNCTION XiP(B1,B2,D1,D2,Xi0,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0,Mw,ew,delta,rhow0,rhoi0,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/((1.0_dp + 0.5_dp*B1 + SQRT(0.25*B1*B1 + D1))**2.0_dp)
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/((1.0_dp + 0.5_dp*B2 + SQRT(0.25*B2*B2 + D2))**2.0_dp)
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
    !PRINT *,rhoi0*l0*eta0*XiT, rhoi0,l0,eta0,XiT   
  END FUNCTION CGTT
  
  RECURSIVE REAL FUNCTION fTildewT(B1,Temperature,D1,delta,ew,l0,cw0,ci0,T0)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,Temperature,D1,delta,ew,l0,cw0,ci0,T0
    fTildewT = 0.0_dp ! TBD
  END FUNCTION fTildewT
  
  RECURSIVE REAL FUNCTION  fTildewp(B1,D1,delta,ew,rhow0,rhoi0)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: B1,D1,delta,ew,rhow0,rhoi0 
    fTildewp = 0.0_dp ! TBD
  END FUNCTION fTildewp

  RECURSIVE REAL FUNCTION KgwpT(rhow0,fTildewTATIP,Kgw)
   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: rhow0,fTildewTATIP,Kgw(3,3)
       KgwpT = 0.0_dp ! TBD
  END FUNCTION KgwpT
  
  END MODULE PermafrostMaterials



!-----------------------------------------------------------------------------
!> 
!------------------------------------------------------------------------------
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar
  TYPE(RockMaterial_t) :: CurrentRockMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRecords, Active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=20,NumberOfEntries=10
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),PorosityPerm(:),SalinityPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:)
  REAL(KIND=dp),ALLOCATABLE :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:), NodalTemperature(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       NoDarcy=.FALSE.,ConstantPorosity=.FALSE., NoSalinity=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  IF (FirstTime) THEN
    NumberOfRecords =  ReadPermafrostRockMaterial( Element,Params,CurrentRockMaterial )
    IF (NumberOfRecords < 1) CALL FATAL(SolverName,'No Rock Material specified')
    FirstTime = .FALSE.
  END IF
  Temperature => Solver % Variable % Values
  TemperaturePerm => Solver % Variable % Perm
  IF ((.NOT.AllocationsDone) .OR. (Model % Mesh % Changed)) THEN
    N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    IF (AllocationsDone) DEALLOCATE(NodalTemperature,NodalPorosity,NodalPressure,NodalSalinity)
    ALLOCATE(NodalTemperature(N),NodalPorosity(N),NodalPressure(N),NodalSalinity(N),STAT=istat )
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
    CALL WARN(SolverName,'Pressure Variable not found. Switching Darcy flow terms off.')
    NoDarcy = .TRUE.
    NULLIFY(Pressure)
  ELSE
    Pressure => PressureVar % Values
    PressurePerm => PressureVar % Perm
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
      IF (NoDarcy) THEN
        NodalPressure(1:N) = 0.0_dp
      ELSE
        NodalPressure(1:N) = Pressure(PressurePerm(Element % NodeIndexes(1:N)))
      END IF
      IF (NoSalinity) THEN
        NodalSalinity(1:N) = 0.0_dp
      ELSE
        NodalSalinity(1:N) = Salinity(SalinityPerm(Element % NodeIndexes(1:N)))
      END IF
   
      CALL LocalMatrix(  Element, N, ND+NB, NodalTemperature, NodalPressure, &
           NodalPorosity, NodalSalinity, CurrentRockMaterial)
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
        CALL LocalMatrixBC(  Element, n, nd+nb )
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
! Assembly of the matrix entries arising from the bulk elements

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t) :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP,XiTAtIP,XiPAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: JgwDAtIP(3),KgwpTAtIP, KgwppAtIP, fTildewTAtIP,fTildewpAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,ew,bs,rhos0,cs0,Xi0,eta0,Kgw(3,3)  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0,p0,rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,CgwTT    ! constants read only once
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,&
         Weight,LoadAtIP,TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,StiffPQ
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Remove this Output:'
    
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kw0th,ki0th,CgwTT
    !------------------------------------------------------------------------------
    IF(.NOT.ConstantsRead) THEN
      DIM = CoordinateSystemDimension()
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
      l0= (CurrentRockMaterial % hw0) - (CurrentRockMaterial % hi0)
      CgwTT = rhow0*cw0
      ConstantsRead=.TRUE.

      CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
      CALL INFO(FunctionName,"General Constants:", Level=9)
      WRITE(Message,'(A)') "GasConstant,Mw,DeltaT,T0,p0,rhow0,rhoi0,hw0,hi0,cw0,ci0,eps,kw0th,ki0th:"
      CALL INFO(FunctionName,Message,Level=9)
      WRITE(Message,'(14E12.5)') CurrentRockMaterial % GasConstant, &
           CurrentRockMaterial % Mw, CurrentRockMaterial % DeltaT, CurrentRockMaterial % T0,&
           CurrentRockMaterial % p0, CurrentRockMaterial % rhow0, CurrentRockMaterial % rhoi0,&           
           CurrentRockMaterial % hw0, CurrentRockMaterial % hi0, CurrentRockMaterial % cw0,&
           CurrentRockMaterial % ci0, CurrentRockMaterial % eps, CurrentRockMaterial % kw0th,&
           CurrentRockMaterial % ki0th
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
    
    ! read element rock material specific parameters    
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    ew = CurrentRockMaterial % ew(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    Xi0 = CurrentRockMaterial % Xi0(RockMaterialID)
    eta0 =CurrentRockMaterial % eta0(RockMaterialID)
    Kgw(1:3,1:3) =CurrentRockMaterial % Kgw(1:3,1:3,RockMaterialID)

    !PRINT *, "ks0th", ks0th,"ew", ew, "bs",bs, "rhos0", rhos0, "cs0", cs0,&
    !     "Xi0",Xi0, "eta0",eta0,"Kgw", Kgw(1:3,1:3)
    
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
      fTildewTAtIP = fTildewT(B1AtIP,TemperatureAtIP,D1InElement,deltaInElement,ew,l0,cw0,ci0,T0)
      fTildewpAtIP = fTildewp(B1AtIP,D1InElement,deltaInElement,ew,rhow0,rhoi0)
      KgwpTAtIP = KgwpT(rhow0,fTildewTATIP,Kgw)
      JgwDAtIP = 0.0_dp ! TBD

      !PRINT *,"KGTTAtIP",KGTTAtIP
      !PRINT *,"CGTTAtIP",CGTTAtIP
      
      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      DO p=1,nd
        DO q=1,nd
          StiffPQ = 0.0
          ! advection term (C*grad(u),v)
          ! C_GW^TT dT/dx_i J_gw^D_i
          ! -----------------------------------
          StiffPQ = StiffPQ + &
               CGWTT * SUM(JgwDAtIP(1:DIM)*dBasisdx(q,1:DIM)) * Basis(p)

          ! diffusion term ( grad(u),grad(v))
          ! div(JGH) = d(KGTT_i,j dT/dx_j)/dx_i
          ! -----------------------------------
          DO i=1,DIM
            DO J=1,DIM
              StiffPQ = StiffPQ + KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
          

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
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'Heat Flux', Found )
    Coeff(1:n) = GetReal( BC,'Robin coefficient', Found )
    IF (.NOT.Found) Coeff(1:n) = 0.0_dp
    Ext_t(1:n) = GetReal( BC,'External field', Found )

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
  END SUBROUTINE LocalMatrixBC
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
      
    
