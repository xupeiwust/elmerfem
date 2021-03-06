!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
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
! *
! * Solve the mixed formulation of the Poisson equation by using div-conforming
! * (face) finite elements of degree k = 1
! *
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Feb 13, 2019
! *
!******************************************************************************

!------------------------------------------------------------------------------
SUBROUTINE MixedPoisson(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found, InitHandles, SecondFamily

  INTEGER :: dim, n, nb, nd, t, istat, active

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), Load(:,:)
  REAL(KIND=dp) :: Norm


  SAVE STIFF, FORCE, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()

  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs  ! just big enough
    ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'MixedPoisson', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !-----------------------
  ! System assembly:
  !----------------------
  SecondFamily = GetLogical(GetSolverParams(), 'Second Kind Basis', Found)

  active = GetNOFActive()
  CALL DefaultInitialize()

  DO t=1,active
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() ! Nodes count corresponding to the background mesh
    nd = GetElementNOFDOFs()  ! The total number of degrees of freedom
    nb = SIZE(Element % BubbleIndexes(:)) ! The number of elementwise degrees 
                                          ! of freedom. NOTE: GetElementNOFBDOFs()
                                          ! doesn't return the right value here 

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(STIFF, FORCE, Element, n, nd, nb, dim, SecondFamily)
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(STIFF, FORCE)

  END DO

  CALL DefaultFinishBulkAssembly()

  InitHandles = .TRUE.
  active = GetNOFBoundaryElements()
  DO t=1,active
    Element => GetBoundaryElement(t)
    IF (ActiveBoundaryElement()) THEN
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrixBC(Element, Mesh, n, nd, SecondFamily, InitHandles)
    END IF
  END DO

  CALL DefaultFinishBoundaryAssembly()

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()  

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(STIFF, FORCE, Element, n, nd, nb, dim, SecondFamily)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: nd  ! The total count of DOFs (nodal, facial and elementwise)
    INTEGER :: nb  ! The number of elementwise DOFs (for the scalar unknown)
    INTEGER :: dim
    LOGICAL :: SecondFamily
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat, Found, EvaluateMatPar, EvaluateSource

    INTEGER :: t, i, j, p, q, np

    REAL(KIND=dp) :: Load(n), a_parameter(n), MatPar, f
    REAL(KIND=dp) :: FaceBasis(nd-nb,3), DivFaceBasis(nd-nb)
    REAL(KIND=dp) :: Basis(n), DetJ, s
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0

    !------------------------------------------------------------------------
    ! The reference element is chosen to be that used for p-approximation,
    ! so we need to switch to using a quadrature which would not be used 
    ! otherwise
    !------------------------------------------------------------------------
    SELECT CASE( GetElementFamily(Element) )
    CASE(3)
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    CASE(5)
      IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
    CASE DEFAULT
      CALL Fatal('MixedPoisson', 'A simplicial mesh assumed currently')
    END SELECT

    !----------------------------------------------------------------
    ! A material parameter and source:
    !----------------------------------------------------------------
    Material => GetMaterial()
    a_parameter(1:n) = GetReal(Material, 'Material Parameter', EvaluateMatPar)
    IF (.NOT. EvaluateMatPar) MatPar = 1.0_dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
        Load(1:n) = GetReal(BodyForce, 'Source Field', EvaluateSource)

    ! Set np = n, if nodal dofs are employed; otherwise set np = 0:
    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)    

    DO t=1,IP % n
      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)

      IF (EvaluateMatPar) MatPar = 1.0_dp/SUM(Basis(1:n) * a_parameter(1:n))
      IF (EvaluateSource) f = SUM(Basis(1:n) * Load(1:n))

      s = detJ * IP % s(t)

      !----------------------------------------------------------------
      ! The following branch could be used to produce the 
      ! Galerkin projection of the pressure for visualization.
      !------------------------------------------------------------------
      IF (np > 0) THEN
        DO p = 1,n
          DO q = 1,n       
            STIFF(p,q) = STIFF(p,q) + Basis(p) * Basis(q) * s    
          END DO

          DO q = nd-nb+1,nd
            STIFF(p,q) = STIFF(p,q) - Basis(p) * 1.0d0 * s            
          END DO
        END DO
      END IF

      !--------------------------------------------------------------
      ! The contribution from the variation with the flux variable q
      !---------------------------------------------------------------
      DO p = 1,nd-np-nb
        i = np + p
        DO q = 1,nd-np-nb
          j = np + q
          STIFF(i,j) = STIFF(i,j) + MatPar * &
              SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * s
        END DO

        DO q = nd-nb+1,nd
          STIFF(i,q) = STIFF(i,q) + 1.0d0 * DivFaceBasis(p) * s
        END DO
      END DO

      !--------------------------------------------------
      ! The contribution from the constraint div q = -f
      !--------------------------------------------------
      DO p = nd-nb+1,nd
        DO q = 1,nd-np-nb
          j = np + q
          STIFF(p,j) = STIFF(p,j) + 1.0d0 * DivFaceBasis(q) * s
        END DO
      END DO
      
      IF (EvaluateSource) THEN
        DO p = nd-nb+1,nd
          FORCE(p) = FORCE(p) - f * 1.0d0 * s 
        END DO
      END IF

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, Mesh, n, nd, SecondFamily, InitHandles)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: nd  ! The total count of DOFs (nodal and facial)
    LOGICAL :: SecondFamily
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Parent, Edge
    TYPE(ValueHandle_t) :: ScalarField
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: RevertSign(6), stat, AssembleForce
    INTEGER :: j, k, l, t, np, p, ActiveFaceId, Family, matches
    REAL(KIND=dp) :: FORCE(nd), Basis(n), TraceBasis(nd), detJ, s, w, g

    SAVE ScalarField, Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT. ASSOCIATED(BC)) RETURN

    IF (InitHandles) THEN
      CALL ListInitElementKeyword(ScalarField, 'Boundary Condition', &
          'Scalar Field')
      InitHandles = .FALSE.
    END IF
    IF (ScalarField % NotPresentAnywhere) RETURN

    ! 
    ! The sign reversion of basis will be checked via the parent element:
    ! 
    Parent => Element % BoundaryInfo % Left
    IF (.NOT. ASSOCIATED(Parent)) THEN
      Parent => Element % BoundaryInfo % Right
    END IF
    IF (.NOT. ASSOCIATED(Parent)) RETURN

    Family = GetElementFamily(Element)
    SELECT CASE(Family)
    CASE(2)
      IF (.NOT. ASSOCIATED(Parent % EdgeIndexes)) RETURN
      !
      ! Pick the edge (face) of the parent element corresponding to the boundary element:
      ! (TO CONSIDER: write a subroutine for this routine task)
      !
      DO ActiveFaceId=1,Parent % TYPE % NumberOfEdges
        Edge => Mesh % Edges(Parent % EdgeIndexes(ActiveFaceId))
        matches = 0
        DO k=1,Element % TYPE % NumberOfNodes
          DO l=1,Edge % TYPE % NumberOfNodes
            IF (Element % NodeIndexes(k)==Edge % NodeIndexes(l)) matches=matches+1
          END DO
        END DO
        IF (matches==Element % TYPE % NumberOfNodes) EXIT
      END DO
      IF (matches /= Element % TYPE % NumberOfNodes) RETURN

      !
      ! Use the parent element to check whether sign reversions are needed:
      !
      CALL FaceElementOrientation(Parent, RevertSign, ActiveFaceId)

    CASE DEFAULT
      CALL Warn('ModelMixedPoisson', 'Neumann BCs in 3-D have not been implemented yet')
      RETURN
    END SELECT

    np = n * Solver % Def_Dofs(GetElementFamily(Parent), Parent % BodyId, 1)

    IF (RevertSign(ActiveFaceId)) THEN
      s = -1.0d0
    ELSE
      s = 1.0d0
    END IF

    CALL GetElementNodes(Nodes)
    IP = GaussPoints(Element)

    Force = 0.0d0
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), DetJ, Basis)
      
      !
      ! The normal traces of face basis functions - now this is available only for RT_0.
      ! NOTE: Here the effect of the Piola transformation is taken into account
      !       such that the multiplication with DetJ is not needed
      ! TO CONSIDER: Get the traces of vector-values basis functions 
      !              by calling a subroutine
      !
      IF (SecondFamily) THEN
        CALL Fatal('ModelMixedPoisson', 'Cannot yet set natural BCS for 2nd family')
      ELSE
        TraceBasis(1) = s * 0.5d0
      END IF

      w = IP % s(t) ! NOTE: No need to multiply with DetJ
      g = ListGetElementReal(ScalarField, Basis, Element, AssembleForce)

      IF (AssembleForce) THEN
        DO p = 1,nd-np
          j = np + p
          FORCE(j) = FORCE(j) + g * TraceBasis(p) * w 
        END DO
      END IF
    END DO

    IF (AssembleForce) CALL DefaultUpdateForce(Force)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE MixedPoisson
!------------------------------------------------------------------------------
