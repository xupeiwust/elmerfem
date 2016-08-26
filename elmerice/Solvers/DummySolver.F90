!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *   Authors: Thomas Zwinger, Martina Schäfer
! *
! ******************************************************************************
!>  As the name DummySolver suggests, this is a Sovler that can be used 
!>  for setting values of variables without actually doing anything
 

RECURSIVE SUBROUTINE DummySolver( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE


  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), Pointer :: BC, BodyForce
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: VarPerm(:)
  INTEGER :: VarDOFs, i, j, k, N, t
  REAL(KIND=dp), POINTER :: VarValues(:)
  LOGICAL :: GotIt

  CALL INFO("DummySolver", "DummySolver", Level=1)

  Var => Solver % Variable
  IF (ASSOCIATED(Var)) THEN    
    VarPerm => Var % Perm
    VarDOFs =  Var % DOFs
    VarValues => Var % Values
  ELSE
    CALL FATAL('DummySolver','No Variable associated')
  END IF
  VarValues = 0.0_dp
  DO t = 1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    N = GetElementNOFNodes(Element)
    BodyForce => GetBodyForce(Element)
    DO j=1,VarDOFs      
      VarValues(VarDOFs*(VarPerm(Element % Nodeindexes(1:N)) - 1)+j) =&
           ListGetReal(BodyForce,TRIM(Solver % Variable % Name),N,Element % NodeIndexes(1:N),GotIt)
      IF (.NOT.GotIt) VarValues(VarDOFs*(VarPerm(Element % Nodeindexes(1:N)) - 1)+j) =0.0_dp             
    END DO
  END DO
  ! VarValues = 0.0_dp
  DO t=1, Solver % Mesh % NumberOfBoundaryElements
    ! get element information
    Element => GetBoundaryElement(t)
    IF ( .NOT.ActiveBoundaryElement() ) CYCLE
    BC => GetBC()
    N = GetElementNOFNodes(Element)
    DO j=1,VarDOFs
      VarValues(VarDOFs*(VarPerm(Element % NodeIndexes(1:N)) - 1) + j) =&
           ListGetReal(BC,TRIM(Solver % Variable % Name),N,Element % NodeIndexes,GotIt)
    END DO
  END DO

END SUBROUTINE DummySolver
