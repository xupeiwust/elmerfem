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
! *
! *  Authors: Denis Cohen, Thomas Zwinger
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! Contains functions for material properties of ice:
! - IceConductivity
! - IceConductivity_m_Mpa_a
! - IceCapacity
! - IceCapacity_m_MPa_a
! - IcePressureMeltingPoint
! 
!==============================================================================
FUNCTION IceConductivity(Model, Node, temp) RESULT(cond)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2
  REAL(KIND=dp) :: Height

  !-----------------------------------------------
  ! Get the height of ice layer (either upper or plain layer)
  !-----------------------------------------------
  Height = 100 ! Assume there is ice
  HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL('IceConductivity', 'Cound not find depth or upper depth')
    END IF
  END IF

  If (Height > 10.100) THEN  ! This is ice
    cond = 9.828*exp(-5.7E-03*temp)
  ELSE                       ! A very conductive layer
    cond = 20
  ENDIF

END FUNCTION IceConductivity

!==============================================================================
FUNCTION IceConductivity_m_Mpa_a(Model, Node, temp) RESULT(cond)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2
  REAL(KIND=dp) :: Height

  !-----------------------------------------------
  ! Get the height of ice layer (either upper or plain layer)
  !-----------------------------------------------
  Height = 100 ! Assume there is ice
  HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL('IceConductivity', 'Cound not find depth or upper depth')
    END IF
  END IF

  If (Height > 10.100) THEN  ! This is ice
    cond = 9.828*exp(-5.7E-03*temp)
  ELSE                       ! A very conductive layer
    cond = 20
  ENDIF
  cond = cond * 31557600 * 1.0e-06 ! From SI to m-MPa-a

END FUNCTION IceConductivity_m_MPa_a

!==============================================================================
FUNCTION IceCapacity(Model, Node, temp) RESULT(capac)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2
  REAL(KIND=dp) :: Height

  !-----------------------------------------------
  ! Get the height of ice layer (either upper or plain layer)
  !-----------------------------------------------
  Height = 100 ! Assume there is ice
  HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL('IceCapacity', 'Cound not find depth or max depth')
    END IF
  END IF

  If (Height > 10.100) THEN ! This is ice
    capac = 146.3+(7.253*temp)
  ELSE                      ! A low heat capacity layer
    capac = 200
  ENDIF

END FUNCTION IceCapacity

!==============================================================================
FUNCTION IceCapacity_m_MPa_a(Model, Node, temp) RESULT(capac)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2
  REAL(KIND=dp) :: Height

  !-----------------------------------------------
  ! Get the height of ice layer (either upper or plain layer)
  !-----------------------------------------------
  Height = 100 ! Assume there is ice
  HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL('IceCapacity', 'Cound not find depth or max depth')
    END IF
  END IF

  If (Height > 10.100) THEN ! This is ice
    capac = 146.3+(7.253*temp)
  ELSE                      ! A low heat capacity layer
    capac = 200
  ENDIF
  capac = capac * 31557600.0**2

END FUNCTION IceCapacity_m_Mpa_a

!==============================================================================
FUNCTION IcePressureMeltingPoint(Model, Node, press) RESULT(Tpmp)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Tpmp, press

  INTEGER :: N
  REAL(KIND=dp) :: ClausiusClapeyron
  TYPE(ValueList_t), POINTER :: Constants
  LOGICAL :: FirstTime = .TRUE., GotIt

  SAVE FirstTime, ClausiusClapeyron

  IF (FirstTime) THEN
    FirstTime = .FALSE.
    Constants => GetConstants()
    IF (.NOT.ASSOCIATED(Constants)) CALL FATAL("IcePressureMeltingPoint","No Constants associated.")
    ClausiusClapeyron = GetConstReal( Constants, 'Clausius Clapeyron Constant', GotIt)
    IF (.NOT.GotIt) THEN
      ClausiusClapeyron = 9.8d-08
      CALL INFO("IcePressureMeltingPoint","No entry found for >Clausius Clapeyron Constant<. Setting to 9.8d-08 (SI units)")
    END IF
  END IF

  Tpmp = 273.15 - ClausiusClapeyron*MAX(press, 0.0_dp)

END FUNCTION IcePressureMeltingPoint

!==============================================================================
