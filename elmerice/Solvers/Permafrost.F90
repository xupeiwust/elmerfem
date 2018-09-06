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
! *       Original Date:  January 2017  -               
! * 
! *****************************************************************************
!>  Solvers for enhanced permafrost problem 
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> Solver for groundwater flow of the enhanced permafrost model
!    (i.e. Darcy Flow representing saturated aquifer)
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
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, DummyDtVar,SalinityDtVar,&
       DummyGWfluxVar,StressInvVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, Active,iter, maxiter, istat,StressInvDOFs
  INTEGER,PARAMETER :: io=22
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), DummyDtPerm(:),SalinityDtPerm(:),&
       StressInvPerm(:),DummyGWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), DummyDt(:),SalinityDt(:),&
       DummyGWflux(:),StressInv(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
       NodalPressure(:), DummyNodalGWflux(:,:), NodalStressInv(:),&
        NodalTemperatureDt(:), NodalSalinityDt(:),&
       NodalDummyDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.,GivenGWFlux,ElementWiseRockMaterial, DummyLog=.FALSE.,&
       InitializeSteadyState, ActiveMassMatrix, ComputeDeformation=.FALSE.,&
       StressInvAllocationsDone=.FALSE.,HydroGeo=.FALSE.,ComputeDt=.FALSE.,&
       TemperatureTimeDerExists=.FALSE.,SalinityTimeDerExists=.FALSE., FluxOutput=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, StressInvName, &
       VarName,PhaseChangeModel,ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure,NodalStressInv, &
       NodalTemperatureDt,NodalDummyDt,NodalSalinityDt,DummyNodalGWflux, &
       ElementWiseRockMaterial, ComputeDeformation, StressInvAllocationsDone
  !------------------------------------------------------------------------------
  CALL DefaultStart()

  Params => GetSolverParams()

  ! check, whether we assume steady state (despite transient run)
  ! this can come handy to produce a balance-pressure field at the
  ! start of the simulation
  !---------------------------------------------------------------
  InitializeSteadyState = GetLogical(Params,'Initialize Steady State',Found)
  IF (Found .AND. InitializeSteadyState .AND. (GetTimeStep() == 1)) THEN
    CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=1)
    ActiveMassMatrix = .FALSE.
  ELSE
    ActiveMassMatrix = .TRUE.
  END IF
  
  maxiter = ListGetInteger( Params, &
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  FluxOutput = GetLogical(Params,'Groundwater Flux Output',Found)
  
  ! find variables for dependencies
  !--------------------------------

  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux, &
       NodalTemperatureDt,NodalDummyDt,NodalSalinityDt,&
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, DummyDtVar, SalinityDtVar, &
       DummyGWfluxVar,DummyGWfluxVar,DummyGWfluxVar, &       
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, DummyDtPerm, SalinityDtPerm, &
       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, DummyDt, SalinityDt,&
       DummyGWflux,DummyGWflux,DummyGWflux, &       
       DummyLog, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt,SolverName)
    
  Pressure => Solver % Variable % Values
  PressurePerm => Solver % Variable % Perm
  VarName = Solver % Variable % Name

  StressInvName =  ListGetString(params,'Rock Stress Invariant Variable Name',ComputeDeformation)

  IF (ComputeDeformation) THEN
    CALL AssignSingleVar(Solver,Model,NodalStressInv,StressInvVar,StressInvName,StressInvDOFs,ComputeDeformation)
    IF (.NOT.ComputeDeformation) THEN
      WRITE (Message,*) '"Rock Stress Invariant Variable Name" assigned as ',TRIM(StressInvName),&
           ', but variable not found - no stress used'
      CALL WARN(SolverName,Message)
    END IF
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

      ! inquire what components have to computed/omitted
      HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
      IF (.NOT.Found) HydroGeo = .FALSE.
      IF(HydroGeo) ComputeDt = .FALSE.
      IF (ComputeDt) THEN
        CALL AssignSingleVarTimeDer(Solver,Model,Element,NodalTemperatureDt,&
             TemperatureDtVar,TemperatureTimeDerExists,dt)
        CALL AssignSingleVarTimeDer(Solver,Model,Element,NodalSalinityDt,&
             SalinityDtVar,SalinityTimeDerExists,dt)        
      END IF
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
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          PRINT *, "NumberOfRockRecords=", NumberOfRockRecords
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

      CALL ReadVars(N,Element,Model,Material,&
           NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux,&
           Temperature, Pressure, Porosity,Salinity,DummyGWFlux,DummyGWFlux,DummyGWFlux,&
           TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
           DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
           NoSalinity,.FALSE.,ConstantPorosity,GivenGWFlux,&
           PorosityName,SolverName,DIM)
      
      IF (ComputeDeformation) CALL ReadSingleVar(N,Element,StressInvPerm,NodalStressInv,StressInv,1)

      IF(ComputeDt) CALL ReadVarsDt(N,Element,Model,Material,&
       NodalTemperatureDt,NodalDummyDt,NodalSalinityDt,&
       TemperatureDtPerm, DummyDtPerm, SalinityDtPerm,&
       TemperatureDt, DummyDt, SalinityDt,&
       NoSalinity,.FALSE.,SolverName,DIM)
      
      ! compose element-wise contributions to matrix and R.H.S
      CALL LocalMatrixDarcy( Model, Element, t, N, ND+NB, Active, NodalPressure, NodalTemperature, &
           NodalPorosity, NodalSalinity, NodalTemperatureDt,NodalDummyDt,NodalSalinityDt, &
           CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
           PhaseChangeModel,ElementWiseRockMaterial, ActiveMassMatrix, &
           NodalStressInv,NoSalinity,ComputeDt,ComputeDeformation)
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
       NodalPorosity, NodalSalinity, NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
       CurrentRockMaterial, CurrentSoluteMaterial,CurrentSolventMaterial,&
       PhaseChangeModel, ElementWiseRockMaterial, ActiveMassMatrix,&
       NodalStress,NoSalinity,ComputeDt,ComputeDeformation )
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:),&
         NodalTemperatureDt(:),NodalPressureDt(:),NodalSalinityDt(:),&
         NodalStress(:)
    LOGICAL :: ElementWiseRockMaterial, ActiveMassMatrix,ComputeDt,ComputeDeformation,NoSalinity
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CgwppAtIP,CgwpTAtIP,CgwpYcAtIP,CgwpI1AtIP,KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),&
         meanfactor,MinKgw,gradTAtIP(3),gradPAtIP(3),gradYcAtIP(3),fluxTAtIP(3),fluxgAtIP(3) ! needed in equation
    REAL(KIND=dp) :: JgwDAtIP(3),JcFAtIP(3), DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3), fcAtIP(3), DispersionCoefficient ! from salinity transport
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP,bijAtIP(2,2),bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: ks0th,e1,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,RadGen,acs(0:5),&
         as0,aas(0:5),ks0,cks(0:5)  ! stuff comming from RockMaterial
    INTEGER :: acsl,aasl,cksl       ! stuff comming from RockMaterial
    REAL(KIND=dp) :: EGAtIP,nuGAtIP,kappaGAtIP ! bedrock deformation
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps,Gravity(3)! real constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,&
         rhowTAtIP,rhowPAtIP,rhowYcAtIP,&
         rhoiPAtIP,rhoiTAtIP,&
         rhocPAtIP,rhocTAtIP,rhocYcAtIP,&
         rhogwPAtIP,rhogwTAtIP,rhogwYcAtIP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,KPorosityAtIP,SalinityAtIP,PressureAtIP
    REAL(KIND=dp) :: TemperatureDtAtIP,SalinityDtAtIP,PressureDtAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp) , POINTER :: gWork(:,:)
    !REAL(KIND=dp) , ALLOCATABLE :: CgwpI1AtNodes(:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE., ConstantDispersion=.FALSE.,&
         CryogenicSuction=.FALSE.
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
    DispersionCoefficient = GetConstReal(Material,"Dispersion Coefficient", ConstantDispersion)
    CryogenicSuction = GetLogical(Material,"Compute Cryogenic Suction", Found)
    IF (.NOT.Found) CryogenicSuction = .FALSE.
   

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

 
      ! Variable gradients at IP
      DO i=1,DIM        
        gradTAtIP(i) =  SUM(NodalTemperature(1:n)*dBasisdx(1:n,i))
        gradYcAtIP(i) = SUM(NodalSalinity(1:n)*dBasisdx(1:n,i))
        gradpAtIP(i) = SUM(NodalPressure(1:n)*dBasisdx(1:n,i))
      END DO

      ! Time derivatives of other variables
      IF (ComputeDt) THEN
        TemperatureDtAtIP = SUM( Basis(1:N) * NodalTemperatureDt(1:N) )
        SalinityDtAtIP = SUM( Basis(1:N) * NodalSalinityDt(1:N) )
      END IF

      
      !Materialproperties needed at IP for Xi computation (anything else thereafter)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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

        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.FALSE.,.TRUE.,.FALSE.)
        IF (XiAtIP .NE. XiAtIP) THEN
          PRINT *, "Darcy: XiAtIP", B1AtIP,D1AtIP,Xi0Tilde
          PRINT *, "Darcy:  XiAtIP", deltaInElement, CurrentRockMaterial % e1(RockMaterialID), bijAtIP
          PRINT *, "Darcy:  XiAtIP", Xi0tilde,SalinityAtIP
          PRINT *, "Darcy: XiAtIP", B1AtIP*B1AtIP + D1AtIP !1.0/(1.0 + 0.5*B1AtIP + SQRT(B1AtIP*B1AtIP + D1AtIP)
          CALL FATAL(SolverName,"XiAtIP is NaN")
        END IF 
      END SELECT

      ! on Xi (directly or indirectly) dependent material parameters (incl. updates) at IP
      rhowAtIP  = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal) ! update
      rhowPAtIP = rhowP(CurrentSolventMaterial,rhowAtIP,p0,PressureAtIP) ! update with new rhowAtIP
      rhowTAtIP = rhowT(CurrentSolventMaterial,rhowAtIP,T0,TemperatureAtIP)
      rhoiPAtIP = rhoiP(CurrentSolventMaterial,rhoiAtIP,p0,PressureAtIP)
      rhoiTAtIP = rhoiT(CurrentSolventMaterial,rhoiAtIP,T0,TemperatureAtIP)
      IF (.NOT.NoSalinity) THEN
        rhocAtIP    = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        rhocPAtIP   = rhocP(CurrentSoluteMaterial,rhocAtIP,ConstVal)
        rhocTAtIP   = rhocT(CurrentSoluteMaterial,rhocAtIP,T0,TemperatureAtIP,ConstVal)
        rhocYcAtIP  = rhocYc(CurrentSoluteMaterial,rhocAtIP,XiAtIP,SalinityAtIP,ConstVal)
        rhowYcAtIP  = rhowYc(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP)
        rhogwYcAtIP = rhogwYc(rhowAtIP, rhocAtIP, rhowYcAtIP,rhocYcAtIP,XiAtIP,SalinityAtIP)
      ELSE
        rhocAtIP    = 0.0_dp
        rhocPAtIP   = 0.0_dp
        rhocTAtIP   = 0.0_dp
        rhocYcAtIP  = 0.0_dp
        rhowYcAtIP  = 0.0_dp
        rhogwYcAtIP = 0.0_dp      
      END IF
      rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
      rhogwpAtIP = rhogwP(rhowPAtIP,rhocPAtIP,XiAtIP,SalinityAtIP)
      rhogwTAtIP = rhogwT(rhowTAtIP,rhocTAtIP,XiAtIP,SalinityAtIP)
      !IF ((rhogwAtIP < 980.0_dp) .OR. (rhogwpAtIP > 1250.0_dp)) THEN
      !  PRINT *,"rhogwAtIP:",rhogwAtIP, rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP
      !END IF

      ! conductivities at IP
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP,MinKgw)
      KgwpTAtIP = 0.0_dp
      KgwppAtIP = 0.0_dp
      IF (CryogenicSuction) THEN
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
      ELSE
        KgwppAtIP = KgwAtIP
      END IF
      KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)

      ! capacities at IP
      EGAtIP = EG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiTAtIP,PorosityAtIP)
      nuGAtIP = nuG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiTAtIP,PorosityAtIP)
      kappaGAtIP = kappaG(EGAtIP,nuGAtIP)
      CgwppAtIP = GetCgwpp(rhogwAtIP,rhoiAtIP,rhogwPAtIP,rhoiPAtIP,kappaGAtIP,&
           XiAtIP,XiPAtIP,CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      CgwpTAtIP = GetCgwpT(rhogwAtIP,rhoiAtIP,rhogwTAtIP,rhoiTAtIP,XiAtIP,XiTAtIP,PorosityAtIP)
      IF (.NOT.NoSalinity) THEN
        CgwpYcAtIP = GetCgwpYc(rhogwAtIP,rhoiAtIP,rhogwYcAtIP,XiAtIP,XiYcAtIP,PorosityAtIP)        
      END IF
      IF (ComputeDeformation) THEN
        CgwpI1AtIP = GetCgwpI1(rhogwAtIP,rhoiAtIP,kappaGAtIP,XiAtIP,CurrentRockMaterial,RockMaterialID)
      END IF
      
     
      ! parameters for diffusion-dispersion flow
      r12AtIP = GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,rhowAtIP,rhocAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)
      DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)
      IF (ConstantDispersion) THEN
        KcAtIP = GetConstKC(DispersionCoefficient)
      ELSE
        JgwDAtIP = 0.0_dp
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,&
             Gravity,rhogwAtIP,DIM,CryogenicSuction)
        KcAtIP = GetKc(CurrentRockMaterial,RockMaterialID,DmAtIP,XiAtIP,JgwDAtIP,PorosityAtIP)
      END IF      
      IF (.NOT.NoSalinity) THEN        
        fcAtIP = GetFc(rhocAtIP,rhowAtIP,Gravity,r12AtIP,XiTAtIP,XiPAtIP,XiAtIP,gradPAtIP,gradTAtIP)
        KcYcYcAtIP = GetKcYcYc(KcAtIP,r12AtIP)
        JcFAtIP = GetJcF(KcYcYcAtIP,KcAtIP,fcAtIP,gradYcAtIP,SalinityAtIP)        
      END IF


      ! fluxes other than pressure induced at IP
      DO i=1,DIM
        fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
        fluxgAtIP(i) = rhogwAtIP * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))   !!
        ! insert missing JcF here
        IF ((fluxgAtIP(i) .NE. fluxgAtIP(i)) .OR. (fluxTAtIP(i) .NE. fluxTAtIP(i))) THEN
          PRINT *, "NaN in r.h.s. of Darcy fluxes"
          PRINT *, "flux(",i,")= Jgwg",fluxgAtIP(i),"+ JgwpT", fluxTAtIP(i)
          PRINT *, "KgwAtIP=",KgwAtIP(i,1:DIM)
          PRINT *, "KgwpTAtIP=",KgwpTAtIP(i,1:DIM)
          PRINT *, "rhowAtIP=",rhowAtIP," rhocAtIP=",rhocAtIP
          STOP
        END IF
      END DO

      ! composition of the matrix:
      ! -----------------------------------
      DO p=1,nd
        DO q=1,nd
          ! next term can be switched off
          ! time derivative (Cgwpp*dp/dt,v):
          ! ------------------------------
          IF (ActiveMassMatrix) &
               MASS(p,q) = MASS(p,q) + Weight * CgwppAtIP * Basis(q) * Basis(p) 

          ! advection term (still needs vstar, hence commented)
          !STIFF (p,q) = STIFF(p,q) + Weight * &
          !   CgwppAtIP * SUM(vstar(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! diffusion term ( Kgwpp * grad(p),grad(v))
          ! div(J_gwp) = d(Kgwpp_i,j dp/dx_j)/dx_i
          ! -----------------------------------
          StiffPQ = 0.0
          DO i=1,DIM
            DO j=1,DIM
              StiffPQ = StiffPQ +  rhogwAtIP * KgwppAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)              
            END DO
          END DO
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
        END DO
      END DO
      ! body forces
      DO p=1,nd     
        FORCE(p) = FORCE(p) + Weight * rhogwAtIP * SUM(fluxgAtIP(1:DIM)*dBasisdx(p,1:DIM))
        IF (CryogenicSuction) &
             FORCE(p) = FORCE(p) + Weight * rhogwAtIP * SUM(fluxTAtIP(1:DIM)*dBasisdx(p,1:DIM))
        IF (ComputeDt) THEN
          FORCE(p) = FORCE(p) + Weight * CgwpTAtIP * TemperatureDtAtIP !dT/dt + v* grad T
          FORCE(p) = FORCE(p) + Weight * CgwpYcAtIP * SalinityDtAtIP ! dyc/dt + v* grad yc
        END IF
        !FORCE(p) = FORCE(p) + &
        !     Weight * PorosityAtIP * (rhocAtIP - rhowAtIP)* SUM(JcFAtIP(1:DIM)*dBasisdx(p,1:DIM))
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
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
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
  TYPE(Variable_t), POINTER :: PressureVar,TemperatureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       DummyGWfluxVar
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       DummyGWfluxPerm(:)
  INTEGER :: NumberOfRockRecords
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       DummyGWflux(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
       NodalPressure(:), NodalGWflux(:,:),NodalTemperatureDt(:),NodalPressureDt(:),&
       NodalSalinityDt(:) ! all dummies
  LOGICAL :: AllocationsDone,ConstantPorosity, NoSalinity,GivenGWFlux=.FALSE.,UnfoundFatal=.TRUE.,ComputeDt=.FALSE.,DummyLog
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, PressureName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux"

  SAVE SaveRHS,AllocationsDone,NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt

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
  
  !ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  ! Read Variables
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       DummyGWfluxVar,DummyGWfluxVar,DummyGWfluxVar, &       
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       DummyGWflux,DummyGWflux,DummyGWflux, &       
       DummyLog, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt,SolverName)
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
    LOGICAL :: Found, ConstVal=.FALSE.,ConstantsRead=.FALSE., FirstTime=.TRUE., ElementWiseRockMaterial,&
         CryogenicSuction=.FALSE.
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: GasConstant, N0, meanfactor,DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         JgwDAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT
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
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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
        rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        !        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        rhowAtIP =  rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        
        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
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
          CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .TRUE.,.FALSE.,.TRUE.,.FALSE.)
        END SELECT
        rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
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
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        !PRINT *,"Flux:", rhogwAtIP,rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP
        Weight = IntegStuff % s(t) * detJ

        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
          END DO
        END DO

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
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
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3,DepthVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:),NodalDepth(:),&
        NodalTemperatureDt(:), NodalSalinityDt(:),&
       NodalPressureDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName
  TYPE(ValueHandle_t) :: Load_h
  
  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,DepthName,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,NodalDepth,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
       ElementWiseRockMaterial,ComputeDt!,NodalGWflux,NoGWflux
  !------------------------------------------------------------------------------
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, 'Computing heat transfer              ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )

  ! Handle to Heat Source (possible description of heat source at elements/IP's) 
  CALL ListInitElementKeyword( Load_h,'Body Force','Heat Source')

  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &       
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt, SolverName)
  
  IF (FirstTime) THEN
    DepthName = ListGetString(Params,'Depth Variable Name', Found)
    IF (.NOT.Found) THEN
      WRITE(DepthName,'(A)') 'Depth'
      CALL WARN(SolverName,' "Depth Variable Name" not found. Assuming default "Depth"')
    END IF
  END IF
  
  CALL AssignSingleVar(Solver, Model,NodalDepth,DepthVar,DepthName,DepthDOFs,DepthExists)
  
  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    WRITE(Message,*) "Nonlinear iteration ", iter, " out of ", maxiter
    CALL INFO( SolverName, Message, Level=3)
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
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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

      CALL ReadVars(N,Element,Model,Material,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
       Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
       PorosityName,SolverName,DIM)

      IF (DepthExists) CALL ReadSingleVar(N,Element,DepthVar % Perm,NodalDepth,DepthVar % Values,DepthVar % DOFs)
     ! IF (DepthExists) &
     !      NodalDepth(1:N) = DepthVar % Values(DepthVar % Perm(Element % NodeIndexes(1:N)))
      
      CALL LocalMatrixHTEQ(  Element, t, Active, n, nd+nb,&
           NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&           
           NodalGWflux, NodalDepth, GivenGWflux, DepthExists, &
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

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHTEQ(  Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, NodalDepth,GivenGWflux,DepthExists, &
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
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:),NodalDepth(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial,DepthExists
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: RefDepth,CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,&
         KgwppAtIP(3,3),fwAtIP,mugwAtIP,DtdAtIP(3,3)!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,CryogenicSuction=.FALSE.,HydroGeo=.FALSE.
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


    
    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementNo  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    IF (DepthExists) THEN
      RefDepth = GetConstReal(Material,'Radiogenic Reference Depth',Found)
      IF (Found) THEN
        DO I=1,N
          LOAD(I) = LOAD(I) + &
               RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth)
          !PRINT *,"HTEQ: RGEN",RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth), NodalDepth(I)
        END DO
 
      END IF
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
    IF (.NOT.Found) HydroGeo = .FALSE.
    
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
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      ! The heat soruce term
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
      !IF (LoadAtIP > 0.0_dp) PRINT *,"HTEQ:LoadAtIP", LoadAtIP
      ! Contribution from other heat source
      LoadAtIP = LoadAtIP + SUM( Basis(1:n) * LOAD(1:n) )
      
      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))
      
      !Materialproperties needed for computing Xi at IP

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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
        CALL  GetXiHartikainen (CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"HTEQ: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP
      
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
      !ci(ci0,aci,T0,TemperatureAtIP,PressureAtIP)  !!
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
      !IF (TemperatureAtIP > 419.00_dp) &
      !     PRINT *, "HTEQ: KGTTAtIP",KGTTAtIP(1,1),KGTTAtIP(1,2),KGTTAtIP(2,2),KGTTAtIP(2,1),"ksthAtIP",ksthAtIP

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      !IF ((ElementNo == 23739) .AND. (t == 1)) &
      !     PRINT *,"HTEQ:", CGTTAtIP, KGTTAtIP, TemperatureAtIP
      !IF (TemperatureAtIP > 419.0_dp) PRINT *,"HTEQ: CGTTAtIP",CGTTAtIP,csAtIP,rhosAtIP,csAtIP*rhosAtIP
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
      !IF (TemperatureAtIP > 419.0_dp) PRINT *,"HTEQ: CgwTTAtIP",CgwTTAtIP,rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,SalinityAtIP
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
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
         !PRINT *,"HTEQ: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO
        
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        !PRINT *,"HTEQ: JgwD=(",JgwDAtIP(1:DIM)*365.5*24.0*3600.0,")"
      END IF
      
      ! add thermal dispersion in Hydro-Geological Mode
      IF (HydroGeo) THEN
        DtdAtIP = GetDtd(CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP,JgwDAtIP)
        DO I=1,DIM
          DO J=1,DIM
            KGTTAtIP(I,J) = KGTTAtIP(I,J) + CGWTTAtIP * DtdAtIP(I,J)
          END DO
        END DO
      END IF

      Weight = IP % s(t) * DetJ
      !PRINT *, "Weight=", Weight
      DO p=1,nd
        DO q=1,nd
          ! diffusion term (KGTTAtIP.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
              !PRINT *,"cond", Weight," *", KGTTAtIP(i,j)," *", dBasisdx(p,j),"*", dBasisdx(q,i)
            END DO
          END DO
          ! advection term (CgwTT * (Jgw.grad(u)),v)
          ! -----------------------------------
          STIFF (p,q) = STIFF(p,q) +&
               Weight * CgwTTAtIP * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p) !
          !PRINT *,"adv", Weight," *", CgwTTAtIP," *", SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)
          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * (CGTTAtIP + 0.0_dp*CGTycAtIP) * Basis(q) * Basis(p) ! check CGTycAtIP
          !PRINT *,"storage", CGTTAtIP, "*",Basis(q) * Basis(p) 
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
          !FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
          FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
          !PRINT *,"LocalMatrixBCHTEQ:",F
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
!> solute (salt ions) transport equation for enhanced permafrost model
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
  TYPE(Variable_t), POINTER :: PressureVar,SalinityVar,PorosityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:), NodalGWflux(:,:),  NodalTemperatureDt(:),&
       NodalSalinityDt(:), NodalPressureDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE., ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostSoluteTransport'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, VarName, TemperatureName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       ElementWiseRockMaterial!,NodalGWflux,NoGWflux
  !------------------------------------------------------------------------------
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, 'Computing solute transport           ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()

  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt, SolverName) 
  
  Salinity => Solver % Variable % Values
  IF (.NOT.ASSOCIATED(Salinity)) THEN
    WRITE(Message,*) "Variable for solute fraction not associated"
    CALL FATAL(Solvername,Message)
  END IF
  SalinityPerm => Solver % Variable % Perm

  
  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    WRITE(Message,*) "Nonlinear iteration ", iter, " out of ", maxiter
    CALL INFO( SolverName, Message, Level=3)
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
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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

      CALL ReadVars(N,Element,Model,Material,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
       Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
       PorosityName,SolverName,DIM)
!      CALL LocalMatrixSolute2(  Element, t, n, nd+nb )
      
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
        CALL ReadVars(N,Element,Model,Material,&
             NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
             Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
             TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
             GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
             NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
             PorosityName,SolverName,DIM)
        CALL LocalMatrixBCSolute(  Element, t, Active, n, nd+nb,&
           NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
           NodalGWflux,GivenGWflux,&
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
        !PRINT *,t,"of",Active,":",n, nb
      END IF
    END DO
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    
    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    DO I=1,Solver % Mesh % NumberOfNodes
      Salinity(SalinityPerm(I)) = MAX(0.0,Salinity(SalinityPerm(I)))
      Salinity(SalinityPerm(I)) = MIN(0.3,Salinity(SalinityPerm(I)))
    END DO
    
    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
    
  END DO
  
  CALL DefaultFinish()
CONTAINS
    
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
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),TemperatureTimeDer,PressureTimeDer,JgwDAtIP(3),&
         KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,KgwppAtIP(3,3),fwAtIP,mugwAtIP!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP, & ! material properties at IP
         CcYcTAtIP, CcYcPAtIP, CcYcYcAtIP, rhocPAtIP, rhocYcAtIP, rhocTAtIP,& ! material properties at IP
         DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3), fcAtIP(3), extforceFlux(3), DispersionCoefficient ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor, DummyTensor(3,3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,j,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,ConstantDispersion,CryogenicSuction=.FALSE.
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

    DispersionCoefficient = GetConstReal(Material,"Dispersion Coefficient", ConstantDispersion)
    
    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
!      PRINT *,"Here0"
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
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT
      
      ! solute and rock densities and derivatives
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      rhocPAtIP = rhocP(CurrentSoluteMaterial,rhocAtIP,ConstVal)
      rhocYcAtIP = rhocYc(CurrentSoluteMaterial,rhocAtIP,XiAtIP,SalinityAtIP,ConstVal)
      rhocTAtIP = rhocT(CurrentSoluteMaterial,rhocAtIP,T0,TemperatureAtIP,ConstVal)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!

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
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        !PRINT *, "Solute: Kgw", KgwAtIP(1,1)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF       
        !PRINT *,"Solute: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        !IF (SalinityAtIP > 0.2_dp) THEN
        !  PRINT *,"Solute: rhogw", rhogwAtIP, rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP          
        !END IF
        ! gradT and gradP have been moved upwards, as needed elsewhere
        
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,&
             Gravity,rhogwAtIP,DIM,CryogenicSuction)

      END IF
      
      ! parameters for diffusion-dispersion flow
      !r12AtIP = GetR(CurrentSoluteMaterial,GasConstant,rhocAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)
      r12AtIP = GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,rhowAtIP,rhocAtIP,XiAtIP,TemperatureAtIP,SalinityAtIP)
      !IF (r12AtIP(2) > 1.2_dp) PRINT *,"Salinity: R2", r12AtIp(2)
      
      !DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)
      DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,CurrentSoluteMaterial % rhoc0,CurrentSolventMaterial % muw0,TemperatureAtIP)
      !PRINT *, "Solute: SalinityAtIP", SalinityAtIP
      !PRINT *, "Solute: Dm", DmAtIP,CurrentSoluteMaterial % rhoc0,CurrentSolventMaterial % muw0,TemperatureAtIP
      IF (ConstantDispersion) THEN
        KcAtIP = GetConstKC(DispersionCoefficient)
        !PRINT *,"DispersionCoefficient",KcAtIP
      ELSE
        DummyTensor = GetConstKC(3.565d-06)
        KcAtIP = GetKc(CurrentRockMaterial,RockMaterialID,DmAtIP,XiAtIP,JgwDAtIP,PorosityAtIP)
        !PRINT *,"Solute: Kc", KcAtIP(1,1), DummyTensor(1,1), DmAtIP,XiAtIP,JgwDAtIP(1:2),PorosityAtIP
      END IF
      KcYcYcAtIP = GetKcYcYc(KcAtIP,r12AtIP)
      !PRINT *,"Solute: KcYcYc", KcYcYcAtIP(1,1)
      fcAtIP = GetFc(rhocAtIP,rhowAtIP,Gravity,r12AtIP,XiTAtIP,XiPAtIP,XiAtIP,gradPAtIP,gradTAtIP) 
      
      Weight = IP % s(t) * DetJ
      !PRINT *,"Solute:",DIM,Weight
      !PRINT *, "Solute: rhoc, Porosity", rhocAtIP,PorosityAtIP
      
      !PRINT *, rhocAtIP
      DO p=1,nd
        DO q=1,nd
          ! diffusion term (Porosity rhoc KcYcYc.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + &
                   Weight * PorosityAtIP * rhocAtIP * KcYcYcAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          !PRINT *, "Solute:  KcYcYcAtIP", KcYcYcAtIP(1,1)
          ! advection term (CcYcYc * (v* .grad(u)),v) ! V* not implemented yet and set to zero
          ! -----------------------------------
          !STIFF (p,q) = STIFF(p,q) &
          !     + Weight * CcYcYcAtIP * SUM(vstarAtIP(1:DIM)*dBasisdx(q,1:dim)) * Basis(p)

          ! left overs from partial integration of fluxes
          ! (rhoc/Xi) (u, Jgw.grad(v))
          !JgwDAtIP(1) = 0.0_dp
          !JgwDAtIP(2) = -1.0d-07
          STIFF (p,q) = STIFF(p,q) &
               - Weight * rhocAtIP * Basis(q) * SUM(JgwDAtIP(1:DIM) * dBasisdx(p,1:DIM))/XiAtIP
!                   PRINT *, "LocalMatrixSolute: ", JgwDAtIP(1:DIM), XiAtIP, rhocAtIP
          
          ! porosity rhoc  (u,(Kc.fc).grad(v))         
          DO i=1,DIM
            extforceFlux =  SUM(KcAtIP(i,1:DIM)*fcAtIP(1:DIM))
          END DO          
          STIFF (p,q) = STIFF(p,q) &
               - Weight * PorosityAtIP * rhocAtIP * Basis(q) * SUM(extforceFlux(1:DIM) * dBasisdx(p,1:DIM))
          
          ! time derivative (CcYcYc*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * CcYcYcAtIP  * Basis(q) * Basis(p)
          !PRINT *, MASS(p,q)
        END DO
      END DO

      
      LoadAtIP = LoadAtIP !+ TemperatureTimeDer * CcYcTAtIP + PressureTimeDer * CcYcPAtIP 
      

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
  SUBROUTINE LocalMatrixBCSolute(Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, GivenGWflux,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    USE DefUtils
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
    REAL(KIND=dp) :: Flux(n), Coeff(n), ImposedSalinity(n), JgwDN(n),F,JgwDNAtIP,Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP, gradTAtIP(3),gradPAtIP(3)       
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    LOGICAL :: Stat,Found,ConstVal,FluxCondition,GWFluxCondition,WeakDirichletCond,ConstantsRead=.FALSE.
    INTEGER :: i,t,p,q,DIM,body_id, other_body_id, material_id, RockMaterialID
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition, ParentMaterial
    TYPE(Element_t), POINTER ::  ParentElement
    TYPE(Nodes_t) :: Nodes    
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostSoluteTransport (LocalMatrixBCSolute)'
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde, XiTAtIP, XiPAtIP, XiYcAtIP,XiEtaAtIP
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP,&
         rhowAtIP, rhoiAtIP, rhocAtIP !needed by XI
    SAVE Nodes,DIM, ConstantsRead, N0, DeltaT, T0, p0, eps, GasConstant, Gravity
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      !PRINT *, "BCSolute: (Constantsread) ", GasConstant, N0, DeltaT, T0, p0, eps, Gravity, ConstantsRead
    END IF

   
    ! inquire parent element and material
    other_body_id = Element % BoundaryInfo % outbody
    IF (other_body_id < 1) THEN ! only one body in calculation
      ParentElement => Element % BoundaryInfo % Right
      IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => Element % BoundaryInfo % Left
    ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
      ParentElement => Element % BoundaryInfo % Right
      IF (ParentElement % BodyId == other_body_id) ParentElement => Element % BoundaryInfo % Left
    END IF
    ! all the above was just so we can get the material properties of the parent element...
    body_id = ParentElement % BodyId
    material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', Found)
    IF (.NOT.Found) CALL FATAL(FunctionName,'Parent Material ID not found')
    
    ParentMaterial => Model % Materials(material_id) % Values
    IF (.NOT. ASSOCIATED(ParentMaterial)) THEN
      WRITE(Message,*)&
           'No material values found for body no ', body_id,&
           ' under material id ', material_id
      CALL FATAL(FunctionName,Message)
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(ParentElement)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementNo  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(ParentMaterial,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(ParentMaterial,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF
    
    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    !Check, whether we have a prescribed solute flow
    Flux(1:n)  = GetReal( BoundaryCondition,'Solute Flow', FluxCondition )
    ! check, whether we have a prescribed groundwater flux
    JgwDN(1:n)  = GetReal( BoundaryCondition,'Groundwater Flux', GWFluxCondition )
    ! Check, whether we have a weakly imposed Dirichlet condition
    ImposedSalinity(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)

    ! if none of the above, we can call it a day
    IF (.NOT.(FluxCondition .OR. GWFluxCondition .OR. WeakDirichletCond)) RETURN

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    !PRINT *,"BCSolute:",deltaInElement,eps,DeltaT,T0,GasConstant

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! we need XiAtIP and rhocAtIP only if we have a non-zero groundwater flux or solute flow
      IF (FluxCondition .OR. GWFluxCondition) THEN ! else spare us the computation

        ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
        TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
        PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
        PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
        SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))


        !Materialproperties needed at IP

        ! water/ice densitities

        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal) !!
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          ! NB: XiTAtIP, XiPAtIP not needed
        CASE DEFAULT ! Hartikainen model
          CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .FALSE., .FALSE., .FALSE.,.FALSE.)
          ! NB: XiTAtIP, XiPAtIP, XiYcAtIP not needed
        END SELECT
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        
        ! Check, whether we have a weakly imposed Dirichlet condition
        IF (GivenGWFlux) THEN
          JgwDNAtIP = SUM(Basis(1:n)*JgwDN(1:n))
          ! contribution from partial integration of groundwater flux term (always on)
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) &
                 + Weight * rhocAtIP * Basis(q) * Basis(p) * JgwDNAtIP/XiAtIP
            END DO
          END DO
        END IF

        !PRINT *,"BCSolute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

        ! Given flux:
        ! -----------
        IF (Fluxcondition) THEN
          F = SUM(Basis(1:n)*Flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) - Weight * PorosityAtIP * rhocAtIP * F * Basis(1:nd)
          !PRINT *,"Salinity BC: Flux:", F, PorosityAtIP , rhocAtIP , Weight, Weight * PorosityAtIP * rhocAtIP * F * Basis(1:nd)
        END IF
      END IF
      
      ! Given salinity, weakly imposed
      !----------------------------------------------------------------------
      IF (WeakDirichletCond) THEN
        SalinityAtIP = SUM(Salinity(1:n)*Basis(1:n))
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * C * SalinityAtIP * Basis(1:nd)
      END IF
    END DO
    !PRINT *, "Salinity BC: Flux:",SUM(FORCE(1:nd))
    CALL DefaultUpdateEquations(STIFF,FORCE)
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
!> output of unfrozen water content as variable (post-processing)
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE PermafrostUnfrozenWaterContentOld( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       DummyGWfluxVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       WaterContentPerm(:),DummyGWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       DummyGWflux(:),WaterContent(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),DummyNodalGWflux(:,:),&
       NodalTemperatureDt(:),NodalPressureDt(:),NodalSalinityDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeDt=.FALSE.,ComputeXiT=.FALSE., DummyLog,GivenGWFlux,ElementWiseRockMaterial
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostUnfrozenWaterContentOld'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName,&
       PhaseChangeModel, ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NumberOfRockRecords,ElementWiseRockMaterial,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,DummyNodalGWflux,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  ComputeXiT = GetLogical(Params,"Compute XiT",Found)
  IF (.NOT.Found) ComputeXiT=.FALSE.
  
  CALL DefaultInitialize()
  
  ! Assign output variables
  WaterContent => Solver % Variable % Values
  WaterContentPerm => Solver % Variable % Perm
  
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  ! Read Variables
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       DummyGWfluxVar,DummyGWfluxVar,DummyGWfluxVar, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       DummyGWflux,DummyGWflux,DummyGWflux, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt,SolverName)

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
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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
      dim = CoordinateSystemDimension()
    END IF
    
    CALL ReadVars(N,Element,Model,Material,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux,&
       Temperature, Pressure, Porosity,Salinity,DummyGWflux,DummyGWflux,DummyGWflux,&
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm,&
       NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
       PorosityName,SolverName,DIM)
    CALL LocalMatrixXi(  Element, n, t, NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
          CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
          PhaseChangeModel,ComputeXiT, ElementWiseRockMaterial)
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
  SUBROUTINE LocalMatrixXi( Element, n, elem, NodalTemperature, NodalPressure, &
       NodalPorosity, NodalSalinity, CurrentRockMaterial, CurrentSoluteMaterial,&
       CurrentSolventMaterial, PhaseChangeModel, ComputeXit, ElementWiseRockMaterial)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, elem
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    LOGICAL :: ComputeXiT, ElementWiseRockMaterial
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT
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
      RockMaterialID = elem  ! each element has it's own set of parameters
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
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      ! Material properties at IP
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)


      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.FALSE.,.TRUE.,.FALSE.)
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


END SUBROUTINE PermafrostUnfrozenWaterContentOld

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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       DummyGWfluxVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       WaterContentPerm(:),DummyGWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       DummyGWflux(:),WaterContent(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),DummyNodalGWflux(:,:),&
       NodalTemperatureDt(:),NodalPressureDt(:),NodalSalinityDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeDt=.FALSE.,ComputeXiT=.FALSE., DummyLog,GivenGWFlux,ElementWiseRockMaterial
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostUnfrozenWaterContent'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName,&
       PhaseChangeModel, ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NumberOfRockRecords,ElementWiseRockMaterial,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,DummyNodalGWflux,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  ComputeXiT = GetLogical(Params,"Compute XiT",Found)
  IF (.NOT.Found) ComputeXiT=.FALSE.
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  CALL DefaultInitialize()
  
  ! Assign output variables
  WaterContent => Solver % Variable % Values
  WaterContentPerm => Solver % Variable % Perm
  
  
  ! Read Variables
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       DummyGWfluxVar,DummyGWfluxVar,DummyGWfluxVar, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       DummyGWflux,DummyGWflux,DummyGWflux, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt,SolverName)

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
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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
      dim = CoordinateSystemDimension()
    END IF
    
    CALL ReadVars(N,Element,Model,Material,&
         NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,DummyNodalGWflux,&
         Temperature, Pressure, Porosity,Salinity,DummyGWflux,DummyGWflux,DummyGWflux,&
         TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
         DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm,&
         NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
         PorosityName,SolverName,DIM)
    CALL LocalSetValue(  Element, n, t, NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
         WaterContent, WaterContentPerm, &
         CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
         PhaseChangeModel,ComputeXiT, ElementWiseRockMaterial)
  END DO
  

  ! Trim values into [0,1]
  IF (.NOT.ComputeXiT) THEN
    DO I=1,Model % NumberOfNodes
      WaterContent(WaterContentPerm(I)) =  MAX(WaterContent(WaterContentPerm(I)),0.0_dp)
      WaterContent(WaterContentPerm(I)) =  MIN(WaterContent(WaterContentPerm(I)),1.0_dp)
    END DO
  END IF

  CALL INFO("SolverName","Computation of unfrozen water content (Xi) for post-processing done",Level=1)

CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalSetValue( Element, n, elem, NodalTemperature, NodalPressure, &       
       NodalPorosity, NodalSalinity, WaterContent, WaterContentPerm, &
       CurrentRockMaterial, CurrentSoluteMaterial,&
       CurrentSolventMaterial, PhaseChangeModel, ComputeXit, ElementWiseRockMaterial)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, elem
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp),  POINTER ::WaterContent(:),&
         NodalTemperature(:), NodalSalinity(:),&
         NodalPorosity(:), NodalPressure(:)
    INTEGER, POINTER :: WaterContentPerm(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    LOGICAL :: ComputeXiT, ElementWiseRockMaterial
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtNode, CgwTTAtNode, KGTTAtNode(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtNode, Xi0Tilde,XiTAtNode,XiPAtNode,XiYcAtNode,XiEtaAtNode,ksthAtNode  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtNode,B2AtNode,DeltaGAtNode, bijAtNode(2,2), bijYcAtNode(2,2),&
         gwaAtNode, giaAtNode, gwaTAtNode,giaTAtNode,gwapAtNode,giapAtNode !needed by XI
    REAL(KIND=dp) :: JgwDAtNode(3),KgwAtNode(3,3),KgwpTAtNode(3,3), MinKgw, KgwppAtNode(3,3), fwAtIp, mugwAtNode !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtNode,D2AtNode
    REAL(KIND=dp) :: GasConstant, N0,DeltaT, T0, p0,eps,Gravity(3)
    REAL(KIND=dp) :: rhowAtNode, rhoiAtNode, rhosAtNode, rhocAtNode
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),U, V, W,DetJ,Weight,LoadAtNode,&
         StiffPQ, meanfactor
    REAL(KIND=DP) :: gradTAtNode(3),gradPAtNode(3),fluxTAtNode(3),fluxPAtNode(3),fluxgAtNode(3)
    REAL(KIND=dp) :: MASS(n,n), STIFF(n,n), FORCE(n), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE.
    !TYPE(GaussIntegrationPoints_t) :: IP
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
      RockMaterialID = elem  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF
    
    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    ! Numerical integration:
    !-----------------------
    DO t=1,N
 
       ! get local coordinates of the point t inside the element
       U = Element % Type % NodeU(t)
       V = Element % Type % NodeV(t)
       W = Element % Type % NodeW(t)

       ! get local information on test-functions and derivatives of the point t
       stat = ElementInfo( Element,Nodes,U,V,W,detJ, &
            Basis,dBasisdx)  
  

      ! Material properties at IP
      rhowAtNode = rhow(CurrentSolventMaterial,T0,p0,NodalTemperature(t),NodalPressure(t),ConstVal)
      rhoiAtNode = rhoi(CurrentSolventMaterial,T0,p0,NodalTemperature(t),NodalPressure(t),ConstVal)
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,NodalPorosity(t))


      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtNode = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,NodalTemperature(t),NodalPressure(t),NodalPorosity(t))
        XiTAtNode = &
             XiAndersonT(XiAtNode,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,NodalTemperature(t),NodalPressure(t),NodalPorosity(t))
        XiPAtNode   = &
             XiAndersonP(XiAtNode,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,NodalTemperature(t),NodalPressure(t),NodalPorosity(t))        
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen (CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             NodalTemperature(t),NodalPressure(t),NodalSalinity(t),NodalPorosity(t),&
             Xi0tilde,deltaInElement,rhowAtNode,rhoiAtNode,&
             GasConstant,p0,T0,&
             XiAtNode,XiTAtNode,XiYcAtNode,XiPAtNode,XiEtaAtNode,&
             .TRUE.,.FALSE.,.TRUE.,.FALSE.)
      END SELECT
      IF (ComputeXiT) THEN
        WaterContent(WaterContentPerm(Element % NodeIndexes(t))) = XiTAtNode
      ELSE
         WaterContent(WaterContentPerm(Element % NodeIndexes(t))) = XiAtNode
      END IF

    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalSetValue


END SUBROUTINE PermafrostUnfrozenWaterContent



!==============================================================================
!>  initialization of Porosity to given reference value in material
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
    Material => GetMaterial(CurrentElement)
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (.NOT.Visited) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      !PRINT *,"PorosityInit:",TRIM(ElementRockMaterialName),ElementWiseRockMaterial
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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

!==============================================================================
!>  initialization of arbitrary scalar given by nodal file
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
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER, POINTER :: NodalVariablePerm(:)
  INTEGER,PARAMETER :: io=26
  INTEGER, ALLOCATABLE :: GlobalToLocalPerm(:)
  REAL(KIND=dp), POINTER :: NodalVariableValues(:)
  REAL(KIND=dp) :: InputField, ValueOffset
  INTEGER :: DIM, i, j, CurrentNode, NumberOfNodes, MaxNumberOfGNodes, MinNumberOfGNodes,&
       OK,  counter, localGlobalRange
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="NodalVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: NodalVariableName,NodalVariableFileName
  LOGICAL :: Visited = .False., Found, Parallel, GotIt

  !SAVE Visited
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------
  
  ! Execute solver only once at beginning
  !if (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing variable to reference ', Level=1)
  CALL Info(SolverName, 'levels in material file            ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

 
  DIM = CoordinateSystemDimension()
  Parallel = (ParEnv % PEs > 1)
  Mesh => GetMesh()
  
  ! Get variable to fill in
  SolverParams => GetSolverParams()

  NodalVariableName = ListGetString(SolverParams, &
       'Nodal Variable', GotIt )
  IF (.NOT.GotIt) THEN
    CALL FATAL(SolverName, ' "Nodal Variable" not found')
  END IF
  NodalVariable => VariableGet( Mesh % Variables, NodalVariableName,GotIt )
  IF (.NOT.GotIt) CALL FATAL(SolverName,"Variable not found")
  
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

  ValueOffset = GetConstReal(SolverParams,'Variable Offset',GotIt)
  IF (.NOT.GotIt) THEN
    ValueOffset = 0.0_dp
  ELSE
    WRITE (Message,*) ' "Variable Offset" found and set to: ', ValueOffset
    CALL INFO(SolverName,Message,Level=3)
  END IF

  NumberOfNodes = Mesh % NumberOfNodes
  
  IF (Parallel) THEN
    MaxNumberOfGNodes = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
    MinNumberOfGNodes = MINVAL(Mesh % ParallelInfo % GlobalDOFs)
    !localGlobalRange = MaxNumberOfGNodes - MinNumberOfGNodes
    IF (MaxNumberOfGNodes <= MinNumberOfGNodes) CALL FATAL(SolverName,"No nodes in parallel domain")
    ALLOCATE(GlobalToLocalPerm(MinNumberOfGNodes:MaxNumberOfGNodes), STAT=OK)
    IF (OK /= 0) CALL FATAL(SolverName,"Allocation error of GlobalToLocalPerm")
    GlobalToLocalPerm = 0
    DO I=1,NumberOfNodes
      GlobalToLocalPerm(Mesh % ParallelInfo % GlobalDOFs(I)) = I
    END DO
    PRINT *, TRIM(SolverName),": ParENV:",ParEnv % MyPE,".  Global Nodal Numbers from",&
         MinNumberOfGNodes,"to",MaxNumberOfGNodes
  ELSE
    MinNumberOfGNodes = 1
    MaxNumberOfGNodes = NumberOfNodes   
  END IF
  
  OPEN(unit = io, file = TRIM(NodalVariableFileName), status = 'old',action='read',iostat = ok)
  IF (ok /= 0) THEN
    WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(NodalVariableFileName)
    CALL FATAL(TRIM(SolverName),TRIM(message))
  ELSE
    !------------------------------------------------------------------------------
    ! Read in the number of records ordered in global node-numbering
    ! in file (first line integer)
    !------------------------------------------------------------------------------
    DO J=1,MaxNumberOfGNodes ! all or in parallel up to max global index
      READ (io, *, END=70, IOSTAT=OK, ERR=80) counter, InputField
      IF (counter .NE. J) CALL FATAL(SolverName,'No concecutive numbering in file')
      IF (J < MinNumberOfGNodes) CYCLE
      
      IF (Parallel) THEN
        I = GlobalToLocalPerm(J)
        IF (I == 0) CYCLE ! point in range, but not in partition        
      ELSE
        I=J
      END IF
      !IF ((NodalVariablePerm(I)<1) .OR. (NodalVariablePerm(I)>NumberOfNodes)) THEN
      !  PRINT *, "NodalVariableInit:", ParEnv % myPE, "NodalVariablePerm(",I,")=",&
      !       NodalVariablePerm(I),">",NumberOfNodes
      !  CALL FATAL(SolverName,'No corresponding entry of target variable')
      !END IF
      NodalVariableValues(NodalVariablePerm(I)) = InputField + ValueOffset
     ! PRINT *,i,counter
    END DO
    !PRINT *, "END", i,counter
70  IF (J-1 .NE. MaxNumberOfGNodes) THEN
      WRITE (Message,*) 'Number of records ',i,' in file ',&
           TRIM(NodalVariableFileName),' does not match number of nodes ',&
           NumberOfNodes, ' in mesh'
      CALL FATAL(SolverName,Message)
    END IF
    CLOSE (io)
    IF (Parallel) &
         DEALLOCATE(GlobalToLocalPerm)
    RETURN
80  CALL FATAL(SolverName,"I/O error")
  END IF
END SUBROUTINE NodalVariableInit

!==============================================================================
!> output of material parameter
!==============================================================================
SUBROUTINE PermafrostMaterialOutput( Model,Solver,dt,TransientSimulation )
!==============================================================================
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       GWFluxVar1,GWFluxVar2,GWFluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,TensorComponent(2)
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       OutputPropertyPerm(:),GWfluxPerm1(:),GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       OutputProperty(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:), NodalGwFlux(:,:),&
       NodalTemperatureDt(:),NodalPressureDt(:),NodalSalinityDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeDt=.FALSE.,ComputeXiT=.FALSE.,ElementWiseRockMaterial,GivenGWflux
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostMaterialOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName,&
       PhaseChangeModel, VariableName,ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NumberOfRockRecords,NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGwFlux,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       ElementWiseRockMaterial
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  VariableName = ListGetString(Params,"Output Property",Found)
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  
  IF (.NOT.Found) &
       CALL INFO(SolverName, ' "Output Property" not found')

  TensorComponent(1:2) = ListGetInteger(Params,"Output Property Component",Found)
  IF (.NOT.Found) TensorComponent = 1

  CALL DefaultInitialize()

  ! Assign output variables
  OutputProperty => Solver % Variable % Values
  OutputPropertyPerm => Solver % Variable % Perm
  
  ! Read Variables
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux,DIM,ComputeDt,SolverName)
  
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
      dim = CoordinateSystemDimension()
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
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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

    CALL ReadVars(N,Element,Model,Material,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
       Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
       GWfluxPerm1,GWfluxPerm2,GWfluxPerm3,&
       NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
       PorosityName,SolverName,DIM)

    CALL LocalMatrixMaterialOutput( Element, t, Active, n, nd,&
         NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
         NodalGWflux, GivenGWflux,&
         CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
         NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, &
         VariableName,TensorComponent)
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
  SUBROUTINE LocalMatrixMaterialOutput( Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, GivenGWflux,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, &
       VariableName,TensorComponent)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, ElementNo, NoElements, NumberOfRockRecords, TensorComponent(2)
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel

    !    INTEGER :: n,TensorComponent(2)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
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
    REAL(KIND=DP) :: PropertyAtIP
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,CryogenicSuction=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName, VariableName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixMaterialOutput)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, GasConstant, N0,DeltaT, T0, p0,eps,Gravity
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

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal) !!
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.FALSE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)

      ! heat capacities
      csAtIP   = cs(CurrentRockMaterial,RockMaterialID,&
           T0,TemperatureAtIP,ConstVal)
      cwAtIP   = cw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
      ciAtIP   = ci(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      ccAtIP   = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)

      ! latent heat
      hiAtIP = hi(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      hwAtIP = hw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)

      ! heat conductivity at IP
      ksthAtIP = GetKalphath(CurrentRockMaterial % ks0th(RockMaterialID),&
           CurrentRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
      kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)      
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP,&
           SalinityATIP,PorosityAtIP,meanfactor)

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
      CGTpAtIP = GetCGTp(rhoiAtIP,hiAtIP,hwAtIP,XiPAtIP,PorosityAtIP) !NEW
      CGTycAtIP = GetCGTyc(rhoiAtIP,hiAtIP,hwAtIP,XiYcAtIP,PorosityAtIP) !NEW
      
      ! groundwater flux at IP
      JgwDAtIP = 0.0_dp
      IF (GivenGWFlux) THEN
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N)) 
        END DO
      ELSE        
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
      END IF

      ! select parameter name for output
      SELECT CASE(VariableName)
      CASE('kgw')
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP =  GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        IF (KgwAtIP(TensorComponent(1),TensorComponent(2)) <= 0.0_dp) STOP
        PropertyAtIP = KgwAtIP(TensorComponent(1),TensorComponent(2))
        !IF (PropertyAtIP < 0.0_dp) THEN
        !  PRINT *, "KgwAtIP < 0 =", KgwAtIP, muw0,muw0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw!
        ! STOP
        ! ELSE
        !   PRINT *, "KgwAtIP=", PropertyAtIP
        ! END IF
      CASE('rhogw')       
        PropertyAtIP = rhogwAtIP
      CASE('gwa')
        PropertyAtIP = gwaAtIP
      CASE('mugw')
        PropertyAtIP = mugwAtIP
      CASE('kgwpt')
        PropertyAtIP = KgwpTAtIP(TensorComponent(1),TensorComponent(2))
      CASE('kgwpp')
        PropertyAtIP = KgwppAtIP(TensorComponent(1),TensorComponent(2))
      CASE('CgwTT')
        PropertyAtIP = CGTTAtIP
      CASE('KGTT')
        PropertyAtIP = KGTTAtIP(TensorComponent(1),TensorComponent(2))
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
END SUBROUTINE PermafrostMaterialOutput


!==============================================================================
!> output of material parameter at IP
!==============================================================================
SUBROUTINE PermafrostIPOutput( Model,Solver,dt,TransientSimulation )
  !==============================================================================
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
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3,DepthVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:),NodalDepth(:),&
       NodalTemperatureDt(:), NodalSalinityDt(:),&
       NodalPressureDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostIPOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName
  TYPE(ValueHandle_t) :: Load_h

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,DepthName,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,NodalDepth,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
       ElementWiseRockMaterial,ComputeDt

  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, ' Assignment of IP variables          ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )

  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &       
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt, SolverName)

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
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
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

    CALL ReadVars(N,Element,Model,Material,&
         NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
         Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
         TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
         GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
         NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
         PorosityName,SolverName,DIM)
    PRINT *,"SetIPValues"
    CALL SetIPValues(  Element, t, Active, n, nd+nb,&
           NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&           
           NodalGWflux, NodalDepth, GivenGWflux, DepthExists, &
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
    
  END DO

CONTAINS
  SUBROUTINE SetIPValues(   Element, ElementNo, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, NodalDepth,GivenGWflux,DepthExists, &
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
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:),NodalDepth(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial,DepthExists
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: RefDepth,CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,&
         KgwppAtIP(3,3),fwAtIP,mugwAtIP,DtdAtIP(3,3)!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,CryogenicSuction=.FALSE.,HydroGeo=.FALSE.
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

    PRINT *,"inside"

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementNo  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    IF (DepthExists) THEN
      RefDepth = GetConstReal(Material,'Radiogenic Reference Depth',Found)
      IF (Found) THEN
        DO I=1,N

!               RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth)
          !PRINT *,"HTEQ: RGEN",RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth), NodalDepth(I)
        END DO

      END IF
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
    IF (.NOT.Found) HydroGeo = .FALSE.

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

    ! Loop all Gauss-points
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      ! The heat soruce term
      !LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
      !IF (LoadAtIP > 0.0_dp) PRINT *,"HTEQ:LoadAtIP", LoadAtIP
      ! Contribution from other heat source
      !LoadAtIP = LoadAtIP + SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      !Materialproperties needed for computing Xi at IP

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
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
        CALL  GetXiHartikainen (CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"HTEQ: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

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
      !ci(ci0,aci,T0,TemperatureAtIP,PressureAtIP)  !!
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
      !IF (TemperatureAtIP > 419.00_dp) &
      !     PRINT *, "HTEQ: KGTTAtIP",KGTTAtIP(1,1),KGTTAtIP(1,2),KGTTAtIP(2,2),KGTTAtIP(2,1),"ksthAtIP",ksthAtIP

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      !IF ((ElementNo == 23739) .AND. (t == 1)) &
      !     PRINT *,"HTEQ:", CGTTAtIP, KGTTAtIP, TemperatureAtIP
      !IF (TemperatureAtIP > 419.0_dp) PRINT *,"HTEQ: CGTTAtIP",CGTTAtIP,csAtIP,rhosAtIP,csAtIP*rhosAtIP
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
      !IF (TemperatureAtIP > 419.0_dp) PRINT *,"HTEQ: CgwTTAtIP",CgwTTAtIP,rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,SalinityAtIP
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
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        !PRINT *,"HTEQ: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)

        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        PRINT *,"IPOutput: JgwD=(",JgwDAtIP(1:DIM)*365.5*24.0*3600.0,")"
      END IF

      ! add thermal dispersion in Hydro-Geological Mode
      IF (HydroGeo) THEN
        DtdAtIP = GetDtd(CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP,JgwDAtIP)
        DO I=1,DIM
          DO J=1,DIM
            KGTTAtIP(I,J) = KGTTAtIP(I,J) + CGWTTAtIP * DtdAtIP(I,J)
          END DO
        END DO
      END IF

      
    END DO
END SUBROUTINE SetIPValues
END SUBROUTINE PermafrostIPOutput











