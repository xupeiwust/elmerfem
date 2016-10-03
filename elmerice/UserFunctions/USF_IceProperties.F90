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
! - IceConductivity_SI
! - IceConductivity_m_Mpa_a
! - IceCapacity_SI
! - IceCapacity_m_MPa_a
! - IcePressureMeltingPoint
!
MODULE IceProperties
  USE DefUtils
  USE Types
  
  IMPLICIT None

  CONTAINS

  !==============================================================================
  FUNCTION IceConductivity(Model, temp) RESULT(cond)
  !==============================================================================

    USE DefUtils

    IMPLICIT None

    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: temp, cond

    cond = 9.828*exp(-5.7E-03*temp)

  END FUNCTION IceConductivity
  
  !==============================================================================
  FUNCTION IceCapacity(Model, temp) RESULT(capac)
  !==============================================================================

    USE DefUtils

    IMPLICIT None

    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: temp, capac



    capac = 146.3_dp + (7.253_dp * temp)
  END FUNCTION IceCapacity
  
  !==============================================================================
  FUNCTION IcePressureMeltingPoint(Model, ClausiusClapeyron, press) RESULT(Tpmp)
  !==============================================================================

    USE DefUtils
    
    IMPLICIT None

    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: Tpmp, ClausiusClapeyron, press

    Tpmp = 273.15 - ClausiusClapeyron*MAX(press, 0.0_dp)

  END FUNCTION IcePressureMeltingPoint

END MODULE IceProperties

!==============================================================================
FUNCTION IceConductivity_SI(Model, Node, temp) RESULT(cond)
!==============================================================================
  USE IceProperties
  
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond
  
  cond = IceConductivity(Model, temp)
  
END FUNCTION IceConductivity_SI
!==============================================================================
FUNCTION IceConductivity_m_Mpa_a(Model, Node, temp) RESULT(cond)
!==============================================================================

  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond

  cond = IceConductivity(Model,temp) * 31557600.0_dp * 1.0d-06 ! From SI to m-MPa-a

END FUNCTION IceConductivity_m_MPa_a

!==============================================================================
FUNCTION IceCapacity_SI(Model, Node, temp) RESULT(capac)
!==============================================================================
  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  capac = IceCapacity(Model,temp)

END FUNCTION IceCapacity_SI

!==============================================================================
FUNCTION IceCapacity_m_MPa_a(Model, Node, temp) RESULT(capac)
!==============================================================================

  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  capac = IceCapacity(Model,temp) * (31557600.0_dp**2.0_dp)

END FUNCTION IceCapacity_m_Mpa_a

!==============================================================================
FUNCTION IcePressureMeltingPoint_K_SI(Model, Node, press) RESULT(Tpmp)
!==============================================================================

  USE IceProperties

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
      CALL INFO("IcePressureMeltingPoint","No entry found for >Clausius Clapeyron Constant<.",Level=9)
      CALL INFO("IcePressureMeltingPoint","Setting to 9.8d-08 (SI units)",Level=9)
    END IF
  END IF

  Tpmp = IcePressureMeltingPoint(Model,ClausiusClapeyron,press)

END FUNCTION IcePressureMeltingPoint_K_SI
!==============================================================================
FUNCTION IcePressureMeltingPoint_C_SI(Model, Node, press) RESULT(Tpmp)
!==============================================================================

  USE IceProperties

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
      CALL INFO("IcePressureMeltingPoint","No entry found for >Clausius Clapeyron Constant<.",Level=9)
      CALL INFO("IcePressureMeltingPoint","Setting to 9.8d-08 (SI units)",Level=9)
    END IF
  END IF

  Tpmp = IcePressureMeltingPoint(Model,ClausiusClapeyron,press) - 273.15_dp

END FUNCTION IcePressureMeltingPoint_C_SI
!==============================================================================
