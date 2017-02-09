MODULE PermafrostMaterials
  USE Types
  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE

  TYPE RockMaterial_t
     INTEGER :: NumberOfRecords
     REAL(KIND=dp), ALLOCATABLE :: ks0th(:),ew(:),bs(:),rhos0(:),cs0(:)
     REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0, p0,T20, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kwth,kith,kcth, Xi0,eta0
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
           CurrentRockMaterial % VariableBaseName)
    END IF
    ALLOCATE(&
         CurrentRockMaterial % ks0th(NumberOfRecords),&
         CurrentRockMaterial % ew(NumberOfRecords),&
         CurrentRockMaterial % bs(NumberOfRecords),&
         CurrentRockMaterial % rhos0(NumberOfRecords),&
         CurrentRockMaterial % cs0(NumberOfRecords),&
         CurrentRockMaterial % VariableBaseName(NumberOfRecords),&
         STAT=OK)
    AllocationsDone = .TRUE.
    IF (OK /= 0) &
         CALL FATAL(FunctionName, 'Allocation Error of input data array')
    !------------------------------------------------------------------------------
    ! Read in information from material file
    ! General constants
    !  GasConstant, Mw, DeltaT, T0, p0,T20, rhow0,rhoi0,&
    !     l0,cw0,ci0,eps,kwth,kith,kcth, Xi0,eta0
    !------------------------------------------------------------------------------
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % GasConstant, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % Mw, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % DeltaT,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % T0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % p0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % T20, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhow0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhoi0,Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % l0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % cw0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ci0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % eps, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % kwth, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % kith, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % kwth, Comment   
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % Xi0, Comment
    READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % eta0, Comment
    ! for each material
    !       ks0th(:),ew(:),b(:),rhos0(:),cs0(:)
    DO I=1,NumberOfRecords
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % VariableBaseName(I)
      WRITE(Message,'(A,I2,A,A)') "Input for Variable No.",I,": ", CurrentRockMaterial % VariableBaseName(I)
      CALL INFO(FunctionName,Message,Level=9)
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ks0th(I), Comment          
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % ew(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % bs(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % rhos0(I), Comment
      READ (io, *, END=10, IOSTAT=OK, ERR=30) CurrentRockMaterial % cs0(I), Comment
      WRITE(Message,'(A)') "Ks0th,Xi0,ew,b,rhos0,cs0:"
      CALL INFO(FunctionName,Message,Level=9)
      WRITE(Message,'(E20.5,E20.5,E20.5,E20.5,E20.5)') CurrentRockMaterial % Ks0th(I), &
           CurrentRockMaterial % ew(I),CurrentRockMaterial % bs(I),CurrentRockMaterial % rhos0(I),&
           CurrentRockMaterial % cs0(I)
      CALL INFO(FunctionName,Message,Level=9)
      CALL INFO(FunctionName,"------------------------------",Level=9)        
    END DO
    WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRecords," rock material records from file ", TRIM(MaterialFileName)
    CALL INFO(FunctionName,Message,Level=1)
10  CLOSE(io)
    IF (I < NumberOfRecords) THEN
      WRITE(Message,'(I3,A,I3)') I,"records read, which is smaller than given number ", NumberOfRecords
      CALL WARN(FunctionName,Message)
    END IF
    RETURN
30  CALL FATAL(FunctionName,"I/O error")
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

  RECURSIVE REAL FUNCTION XiT(B1,B2,D1,D2,Xi0,ew,Temperature)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0
    !local
    REAL(KIND=dp) :: aux1, aux2
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/(1.0_dp + 0.5_dp*B1 + SQRT(0.25*B1*B1 + D1))
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/(1.0_dp + 0.5_dp*B2 + SQRT(0.25*B2*B2 + D2))
    XiT= (0.5_dp*Xi0*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0)*aux2/delta) &
         *Mw*(l0 + (cw0 - ci0)*(Temperature -T0))/(T0*R*Temperature)
  END FUNCTION XiT

  RECURSIVE REAL FUNCTION XiP(B1,B2,D1,D2,Xi0,ewrhow0,rhoi0,Temperature)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0
    !local
    REAL(KIND=dp) :: aux1, aux2
    aux1 = (1.0_dp + B1/SQRT(B1*B1 + 4.0_dp*D1))/(1.0_dp + 0.5_dp*B1 + SQRT(0.25*B1*B1 + D1))
    aux2 = (1.0_dp + B2/SQRT(B2*B2 + 4.0_dp*D2))/(1.0_dp + 0.5_dp*B2 + SQRT(0.25*B2*B2 + D2))
    XiT= (0.5_dp*Xi0*aux1/(ew + delta) + 0.5_dp*(1.0_dp - Xi0)*aux2/delta) &
         *Mw*((1.0_dp/rhoi0) - (1.0_dp/rhow0))/(R*Temperature)
  END FUNCTION XiT
  

  RECURSIVE REAL FUNCTION Xi(B1,B2,D1,D2,Xi0)
    REAL(KIND=dp), INTENT(IN) :: B1,B2,D1,D2,Xi0
    Xi= Xi0/(1.0_dp + 0.5_dp*B1 + SQRT(0.25_dp*B1*B1 + D1)) &
         + (1.0_dp - Xi0)/(1.0_dp + 0.5_dp*B2 + SQRT(0.25_dp*B2*B2 + D2))
  END FUNCTION Xi

  RECURSIVE REAL FUNCTION ksth(ks0th,bs,T0,Temperature)
    REAL(KIND=dp), INTENT(IN) :: ks0th,bs,T0,Temperature
    ksth = ks0th/( 1.0_dp + bs*(Temperature - T0)/T0)
  END FUNCTION ksth
  
  RECURSIVE FUNCTION GetKGTT(ksth,kwth,kith,eta0,Xi,Temperature,Pressure,Porosity)RESULT(KGTT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ksth,kwth,kith,eta0,Xi,Temperature,Pressure,Porosity
    ! local
    REAL(KIND=dp) :: KGTT(3,3), factor,unittensor(3,3)
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    factor = (1.0_dp - eta0)/ksth + Xi*eta0/kwth + (1.0_dp - Xi)*eta0/kith
    KGTT = unittensor/factor
  END FUNCTION GetKGTT

  RECURSIVE REAL FUNCTION CGTT(Xi,XiT,rhos0,rhow0,rhoi0,cw0,ci0,cs0,l0,eta0,T0,Temperature)
    REAL(KIND=dp), INTENT(IN) :: Xi,XiT,rhow0,rhoi0,cw0,ci0,cs0,l0,eta0,T0,Temperature
    CGTT = (1.0_dp - eta0)*rhos0*cs0 + Xi*eta0*rhow0*cw0 &
         + (1.0_dp - Xi)*eta0*rhoi0*ci0 &
         + rhoi0*l0*eta0*XiT
  END FUNCTION CGTT
  
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
  REAL(KIND=dp) :: Norm
  REAL(KIND=dp),POINTER :: Ks0TT(:,:,:),Xi0(:),ew(:),bs(:),rhos0(:),cs0(:), &
       Temperature(:), Pressure(:), Porosity(:), Salinity(:)
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
 
        
      CALL LocalMatrix(  Element, N, ND+NB, NodalTemperature,&
           NodalPorosity,NodalPressure, NodalSalinity,CurrentRockMaterial )
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
        !CALL LocalMatrixBC(  Element, n, nd+nb, CurrentRockMaterial )
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
  SUBROUTINE LocalMatrix( Element, n, nd, NodalTemperature, NodalSalinity,&
       NodalPorosity, NodalPressure, CurrentRockMaterial )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t) :: CurrentRockMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTT, CgwTT, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP !needed by XI
    REAL(KIND=dp) :: deltaInElement,D1InElement,D2InElement
    REAL(KIND=dp) :: ks0th,rhos0,bs,cs0,ew  ! stuff comming from RockMaterial
    REAL(KIND=dp) :: GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kwth,kith,kcth,Xi0,eta0     ! coonstants read only once
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,&
         Weight,LoadAtIP,TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,StiffPQ
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes

    
    SAVE Nodes, ConstantsRead, DIM, GasConstant, Mw, DeltaT, T0, p0, rhow0,rhoi0,&
         l0,cw0,ci0,eps,kwth,kith,kcth,Xi0,eta0
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
      l0= CurrentRockMaterial % l0
      cw0= CurrentRockMaterial % cw0
      ci0= CurrentRockMaterial % ci0
      eps=CurrentRockMaterial % eps
      kwth= CurrentRockMaterial % kwth
      kith= CurrentRockMaterial % kith
      kcth= CurrentRockMaterial % kcth
      Xi0= CurrentRockMaterial % Xi0
      eta0= CurrentRockMaterial % eta0
      ConstantsRead=.TRUE.
    END IF
    
    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         Load(1:n) = GetReal( BodyForce,'Heat source', Found )

    ! read variable material parameters from CurrentRockMaterial
    Material => GetMaterial()
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    
    ! read element rock material specific parameters
    ew = CurrentRockMaterial % ew(RockMaterialID)
    ks0th = CurrentRockMaterial % ks0th(RockMaterialID)
    rhos0 = CurrentRockMaterial % rhos0(RockMaterialID)
    cs0 = CurrentRockMaterial % cs0(RockMaterialID)
    bs = CurrentRockMaterial % bs(RockMaterialID)

    ! derive element rock material specific parameters
    deltaInElement = delta(ew,eps,DeltaT,T0,Mw,l0,cw0,ci0,GasConstant)
    D1InElement = D1(deltaInElement,ew)
    D2InElement = 1.0_dp
    
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
      ksthAtIP = ksth(ks0th,bs,T0,TemperatureAtIP)
      KGTTAtIP = GetKGTT(ksthAtIP,kwth,kith,eta0,XiAtIP,TemperatureAtIP,PressureAtIP,PorosityAtIP)
      
      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      StiffPQ = 0.0
      DO p=1,nd
        DO q=1,nd
          ! advection term (C*grad(u),v)
          ! -----------------------------------
          !IF (VeloFound) StiffPQ = StiffPQ + &
          !     CGWTT * SUM(Velo(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! diffusion term ( grad(u),grad(v))
          ! -----------------------------------
          !STIFF(p,q) = STIFF(p,q) + Weight * R*Basis(q) * Basis(p)
          DO i=1,DIM
            DO J=1,DIM
              StiffPQ = StiffPQ + KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ

          ! time derivative (rho*du/dt,v): !! THIS IS OK AS IS !!!
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * CGTT * Basis(q) * Basis(p)
        END DO
      END DO

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

    Flux(1:n)  = GetReal( BC,'field flux', Found )
    Coeff(1:n) = GetReal( BC,'robin coefficient', Found )
    Ext_t(1:n) = GetReal( BC,'external field', Found )

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
      
    
