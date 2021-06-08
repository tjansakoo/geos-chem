!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dust_mod.F90
!
! !DESCRIPTION: Module DUST\_MOD contains routines for computing dust aerosol
!  emissions, chemistry, and optical depths.
!\\
!\\
! !INTERFACE:
!
MODULE DUST_MOD
!
! !USES:
!
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEMDUST
#ifdef TOMAS
  PUBLIC  :: SETTLEDUST
#endif
  PUBLIC  :: RDUST_ONLINE
!  PUBLIC  :: GET_DUST_ALK
  PUBLIC  :: INIT_DUST
  PUBLIC  :: CLEANUP_DUST
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DRY_SETTLING
#ifdef APM
  PRIVATE :: DRY_SETTLINGBIN
#endif
!
!  !REVISION HISTORY:
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER               :: id_DST1,  id_DST2,  id_DST3,   id_DST4
  INTEGER               :: id_DAL1,  id_DAL2,  id_DAL3,   id_DAL4
  INTEGER               :: id_DUST1, id_NK1

  ! Arrays
  REAL(fp), ALLOCATABLE :: FRAC_S(:)
  REAL(fp), ALLOCATABLE :: SRCE_FUNC(:,:,:)

#ifdef TOMAS
  ! To replicate the obsolete Input_Opt%IDDEP field
  ! Set to a large placeholder value
  INTEGER               :: TOMAS_IDDEP(1000)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemdust
!
! !DESCRIPTION: Subroutine CHEMDUST is the interface between the GEOS-Chem
!  main program and the dust chemistry routines that mostly calculates dust
!  dry deposition.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMDUST( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Diag_Mod,     ONLY : DgnState
#ifdef APM
    USE HCO_State_GC_Mod,   ONLY : HcoState
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
    USE APM_INIT_MOD,       ONLY : APMIDS
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Non-SAVEd scalars
    LOGICAL            :: LDRYD
    LOGICAL            :: LDUST
    LOGICAL            :: LDSTUP
    LOGICAL            :: prtDebug

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc
#ifdef APM
    INTEGER       :: I,J,N,IDDST
    REAL(fp)      :: A_M2, E_DST, DTSRCE
    REAL(fp), POINTER :: Spc(:,:,:,:)
#endif

    !=================================================================
    ! CHEMDUST begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at CHEMDUST (in module GeosCore/dust_mod.F)'

    ! Copy fields from INPUT_OPT to local variables for use below
    LDRYD  = Input_Opt%LDRYD
    LDUST  = Input_Opt%LDUST
    LDSTUP = Input_Opt%LDSTUP

    ! Set a flag for debug output
    prtDebug = ( Input_Opt%LPrt .and. Input_Opt%amIRoot )

    ! Execute on first call only
    IF ( FIRST ) THEN

       ! Stop w/ error if dust species flags are undefined
       ! NOTE: Ind_() returns -1 for species not found, so for the
       ! algorithm below to work, we need to make sure all id's are
       ! not less than zero (bmy, 7/5/16)
       IF ( MAX( id_DST1, 0 )  + &
            MAX( id_DST2, 0 )  + &
            MAX( id_DST3, 0 )  + &
            MAX( id_DST4, 0 ) == 0 ) THEN
          IF ( LDUST ) THEN
             ErrMsg = 'LDUST=T but dust species are undefined!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

#ifdef APM
    ! Initialize GEOS-Chem tracer array [kg] from Chemistry State
    ! object
    ! (mpayer, 12/6/12)
    Spc => State_Chm%Species

    ! Emission timestep
    DTSRCE = HcoState%TS_EMIS

    IDDST   = HCO_GetHcoID( 'DST1',   HcoState )-1

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( N, J, I, A_M2, E_DST )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

       DO N=1,HcoState%nDust
          ! Get emissions [kg/m2/s] and convert to [kg/box]
          E_DST = HcoState%Spc(IDDST+N)%Emis%Val(I,J,1) * A_M2 * DTSRCE

          Spc(I,J,1,APMIDS%id_DSTBIN1)    = Spc(I,J,1,APMIDS%id_DSTBIN1)    &
                                            +E_DST*1.5083735055071434d-006
          Spc(I,J,1,APMIDS%id_DSTBIN1+1)  = Spc(I,J,1,APMIDS%id_DSTBIN1+1)  &
                                            +E_DST*4.7465319694587000d-005
          Spc(I,J,1,APMIDS%id_DSTBIN1+2)  = Spc(I,J,1,APMIDS%id_DSTBIN1+2)  &
                                            +E_DST*6.5459286253814980d-004
          Spc(I,J,1,APMIDS%id_DSTBIN1+3)  = Spc(I,J,1,APMIDS%id_DSTBIN1+3)  &
                                            +E_DST*3.6946414378027244d-003
          Spc(I,J,1,APMIDS%id_DSTBIN1+4)  = Spc(I,J,1,APMIDS%id_DSTBIN1+4)  &
                                            +E_DST*8.0323357006383631d-003
          Spc(I,J,1,APMIDS%id_DSTBIN1+5)  = Spc(I,J,1,APMIDS%id_DSTBIN1+5)  &
                                            +E_DST*1.5922384531190468d-002
          Spc(I,J,1,APMIDS%id_DSTBIN1+6)  = Spc(I,J,1,APMIDS%id_DSTBIN1+6)  &
                                            +E_DST*3.7324808601695736d-002
          Spc(I,J,1,APMIDS%id_DSTBIN1+7)  = Spc(I,J,1,APMIDS%id_DSTBIN1+7)  &
                                            +E_DST*0.1144517534356116
          Spc(I,J,1,APMIDS%id_DSTBIN1+8)  = Spc(I,J,1,APMIDS%id_DSTBIN1+8)  &
                                            +E_DST*0.1880182758312848
          Spc(I,J,1,APMIDS%id_DSTBIN1+9)  = Spc(I,J,1,APMIDS%id_DSTBIN1+9)  &
                                            +E_DST*0.2448371224641443
          Spc(I,J,1,APMIDS%id_DSTBIN1+10) = Spc(I,J,1,APMIDS%id_DSTBIN1+10) &
                                            +E_DST*0.1452025570524453
          Spc(I,J,1,APMIDS%id_DSTBIN1+11) = Spc(I,J,1,APMIDS%id_DSTBIN1+11) &
                                            +E_DST*0.1954553759504486
          Spc(I,J,1,APMIDS%id_DSTBIN1+12) = Spc(I,J,1,APMIDS%id_DSTBIN1+12) &
                                            +E_DST*4.0358202390254526d-002
          Spc(I,J,1,APMIDS%id_DSTBIN1+13) = Spc(I,J,1,APMIDS%id_DSTBIN1+13) &
                                            +E_DST*5.2032641469382107d-003
          Spc(I,J,1,APMIDS%id_DSTBIN1+14) = Spc(I,J,1,APMIDS%id_DSTBIN1+14) &
                                            +E_DST*7.2351157743278424d-004
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()
#endif

    !=================================================================
    ! Do dust settling & deposition
    !=================================================================

    !-------------------------------------------------------------
    ! Dust settling for regular dust species
    !
    ! NOTE: Assumes that id_DST1, id_DST2, id_DST3, id_DST4 are
    ! contiguous in the species list.  This is usually correct,
    ! but we may have to update this later. (bmy, 7/5/16)
    !-------------------------------------------------------------
    CALL DRY_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                       State_Met, id_DST1,   id_DST4,    RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in call to "DRY_SETTLING"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------------------------------------
    ! Dust settling for acid uptake on dust species
    !
    ! NOTE: Assumes that id_DAL1, id_DAL2, id_DAL3, id_DAL4 are
    ! contiguous in the species list.  This is usually correct,
    ! but we may have to update this later. (bmy, 7/5/16)
    !
    ! ALSO NOTE: Dry settling of DSTNIT, DSTSO4 occurs in
    ! GRAV_SETTLING, called from CHEMSULFATE in SULFATE_MOD.
    !-------------------------------------------------------------
    IF ( LDSTUP ) THEN
       CALL DRY_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                          State_Met, id_DAL1,   id_DAL4,    RC )

#ifdef APM
       ! Call APM size-resolved settling algorithm
       CALL DRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, State_Met, RC )
#endif

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in call to "DRY_SETTLING"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Debug print
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### CHEMDUST: a DRY_SETTLING' )
    ENDIF

  END SUBROUTINE CHEMDUST
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: settledust
!
! !DESCRIPTION: Subroutine SETTLEDUST is the interface between the
!  size-resolved dry deposition subroutine AERO\_DRYDEP and the dust module.
!  This is to call only gravitational settling and deals with removal of
!  aerosol number with the dust mass.  (win, 7/17/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SETTLEDUST( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD             ! ND44
    USE DIAG_MOD,           ONLY : AD44
#endif
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AVO
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TOMAS_MOD,          ONLY : IBINS, Xk, SRTDUST
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  17 Jul 2009 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE          :: FIRST = .TRUE.

    ! Non-SAVEd scalars
    INTEGER                :: N, BIN, I, J, L, NN, IDISP
    REAL(fp)               :: DU0(State_Grid%NX,State_Grid%NY,State_Grid%NZ,IBINS)
    REAL(fp)               :: DIF, FLUXN, FLUXD
    REAL(fp)               :: DT_SETTL, AREA_CM2

    !debug
    integer                :: ii, jj , ix, jx, bb
    data                      ii,jj, ix, jx, bb /37, 24, 58, 34, 30 /
    CHARACTER(LEN=255)     :: MSG, LOC

    ! Pointers
    REAL(fp),      POINTER :: Spc(:,:,:,:)
    TYPE(Species), POINTER :: ThisSpc

    !=================================================================
    ! SETTLEDUST begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine SETTLEDUST in dust_mod.F'
       CALL GC_Error( MSG, RC, LOC )
    ENDIF

    ! Set pointers
    Spc     => State_Chm%Species  ! [kg], for now
    ThisSpc => NULL()

    !=================================================================
    ! Do dust settling & deposition
    !=================================================================

    ! Dust settling timestep [s]
    DT_SETTL = GET_TS_CHEM()

    ! Save initial dust mass
    DO BIN = 1, IBINS

       ! Species ID
       N = id_DUST1 - 1 + BIN

       DO L   = 1, State_Grid%NZ
       DO J   = 1, State_Grid%NY
       DO I   = 1, State_Grid%NX
          DU0(I,J,L,BIN) = Spc(I,J,L,N)
       ENDDO
       ENDDO
       ENDDO
    ENDDO

    ! Dust settling
    CALL DRY_SETTLING( Input_Opt, State_Chm,  State_Diag, State_Grid, &
                       State_Met, id_DUST1,   id_DUST1-1+IBINS, RC )

    ! Calculate change in number to correspond with dust redistribution
    ! by gravitational settling
    DO BIN = 1, IBINS

       ! Species ID
       N       =  id_DUST1 - 1 + BIN

       ! Look up this species in the species database
       ThisSpc => State_Chm%SpcData(N)%Info

       ! Drydep index (??? -- replace w/ species DB info?)
       NN      = State_Chm%nDryDep + (SRTDUST-1)*IBINS + BIN

       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Surface area [cm2]
          AREA_CM2 = State_Grid%Area_M2(I,J) * 1e+4_fp

          DIF   = 0.0_fp
          FLUXD = 0e+0_fp
          FLUXN = 0e+0_fp

          DO L = 1, State_Grid%NZ
             DIF = DU0(I,J,L,BIN) - Spc(I,J,L,id_DUST1-1+BIN)

             Spc(I,J,L,id_NK1-1+BIN) = Spc(I,J,L,id_NK1-1+BIN) - &
                                       DIF/(SQRT( Xk(BIN)*Xk(BIN+1)))

             ! Convert flux from [kg/s] to [molec/cm2/s]
             FLUXD = FLUXD + DIF / DT_SETTL * AVO / &
                     ( 1.e-3_fp * ThisSpc%MW_g ) / AREA_CM2


             FLUXN = FLUXN + DIF/(SQRT( Xk(BIN)*Xk(BIN+1))) / &
                     DT_SETTL * AVO / &
                     ( 1.e-3_fp * ThisSpc%MW_g ) / AREA_CM2
          ENDDO

#ifdef BPCH_DIAG
          !===========================================================
          !  Dry deposition diagnostic [#/cm2/s] (bpch)
          !===========================================================
          IF ( ND44 > 0 ) THEN

             !--------------------------------------------------------
             ! ND44 DIAGNOSTIC (bpch)
             ! Dry deposition flux loss [#/cm2/s]
             !
             ! NOTE: Bpch diagnostics are being phased out.
             !--------------------------------------------------------
             !%%% NOTE: Now use TOMAS_IDDEP, which replicates the
             !%%% since-removed Input_Opt%IDDEP field (bmy, 3/17/17)
             AD44(I,J,TOMAS_IDDEP(BIN),1) = AD44(I,J,TOMAS_IDDEP(BIN),1) +  &
                                            FLUXN
             AD44(I,J,NN,1) = AD44(I,J,NN,1) + FLUXD
          ENDIF
#endif

       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE SETTLEDUST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissdust
!
! !DESCRIPTION: Subroutine EMISSDUST is the driver routine for the dust
!  emission module.  You may call either the GINOUX or the DEAD dust source
!  function.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSDUST( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD             ! ND59
    USE DIAG_MOD,           ONLY : AD59_DUST, AD59_NUMB  !(win, 7/17/09)
#endif
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TOMAS_MOD,          ONLY : IBINS, XK             !(win, 7/17/09)
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL,  SAVE    :: FIRST = .TRUE.

    ! Non-SAVEd scalars
    LOGICAL           :: LDEAD
    LOGICAL           :: LDUST
    LOGICAL           :: prtDebug
    LOGICAL           :: LINTERP
    INTEGER           :: I, J, K
    REAL(fp)          :: MEMIS
    REAL(fp)          :: MINIT(State_Grid%NX,State_Grid%NY,1,IBINS)

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! EMISSDUST begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at EMISSDUST (in module GeosCore/dust_mod.F)'

    ! Initialize
    LDEAD     =  Input_Opt%LDEAD
    LDUST     =  Input_Opt%LDUST
    prtDebug  =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Execute on first-call only
    IF ( FIRST ) THEN

       ! Return if dust ID flags are not defined
       IF ( id_DUST1 == 0 ) THEN
          IF ( LDUST ) THEN
             ErrMsg = 'LDUST=T but dust species are undefined!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    !=================================================================
    ! For TOMAS microphysics
    !=================================================================
    IF ( id_NK1 > 0 .and. id_DUST1 > 0 ) THEN

       MINIT(:,:,1,1:IBINS) = Spc(:,:,1,id_DUST1:id_DUST1-1+IBINS)

       IF ( LDEAD ) THEN
          ! still didn't figure out why run would crash w/ this option
          ! (win, 7/17/09)
          PRINT*,'Currently the TOMAS code does not work with ', &
                 'dust DEAD emission yet!  Switch to GINOUX for now'
          STOP

       ELSE
          IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')
       ENDIF

#ifdef BPCH_DIAG
       IF ( ND59 > 0 ) THEN
          DO K = 1, IBINS
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             MEMIS = Spc(I,J,1,id_DUST1-1+K) - MINIT(I,J,1,K)
             IF ( MEMIS == 0.e+0_fp ) GOTO 10

             AD59_DUST(I,J,1,K) = AD59_DUST(I,J,1,K) + MEMIS ! kg ????
             Spc(I,J,1,id_NK1-1+K) = Spc(I,J,1,id_NK1-1+K) + &
                                     MEMIS / (sqrt(Xk(K)*Xk(K+1)))

             AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + &
                                  MEMIS / (sqrt(Xk(K)*Xk(K+1)))
10           CONTINUE
          ENDDO
          ENDDO
          ENDDO
       ENDIF
#endif
    ENDIF

    ! Free pointers
    Spc => NULL()

  END SUBROUTINE EMISSDUST
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settling
!
! !DESCRIPTION: Subroutine DRY\_SETTLING computes the dry settling of
!  dust species.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DRY_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                           State_Met, Ind0,      Ind1,       RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: Ind0        ! Starting dust species #
    INTEGER,        INTENT(IN)    :: Ind1        ! Ending   dust species #
!
! !INPUT/OUTPUT PARAMETERS
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I,        J,     L
    INTEGER                :: N,        NA,    ND,   S
    REAL(fp)               :: DT_SETTL, DELZ,  DELZ1
    REAL(fp)               :: REFF,     DEN,   CONST
    REAL(fp)               :: NUM,      LAMDA, FLUX
    REAL(fp)               :: AREA_CM2, TC0(State_Grid%NZ)
    REAL(fp)               :: TOT1,     TOT2

    ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
    REAL(fp)               :: P

    ! Diameter of aerosol [um]
    REAL(fp)               :: Dp

    ! Pressure * DP
    REAL(fp)               :: PDp

    ! Temperature (K)
    REAL(fp)               :: TEMP

    ! Slip correction factor
    REAL(fp)               :: Slip

    ! Viscosity of air (Pa s)
    REAL(fp)               :: Visc

    ! Settling velocity of particle (m/s)
    REAL(fp)               :: VTS(State_Grid%NZ)

    ! Pointers
    TYPE(Species), POINTER :: ThisSpc
    REAL(fp),      POINTER :: TC(:,:,:,:)

!
! !DEFINED PARAMETERS:
!
    REAL(fp),  PARAMETER   :: C1 =  0.7674e+0_fp
    REAL(fp),  PARAMETER   :: C2 =  3.079e+0_fp
    REAL(fp),  PARAMETER   :: C3 =  2.573e-11_fp
    REAL(fp),  PARAMETER   :: C4 = -1.424e+0_fp

    !=================================================================
    ! DRY_SETTLING begins here!
    !=================================================================

    ! Initialize
    RC        =  GC_SUCCESS
    ThisSpc   => NULL()
    TC        => State_Chm%Species

    ! Dust settling timestep [s]
    DT_SETTL  = GET_TS_CHEM()

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    ) &
    !$OMP PRIVATE( CONST, AREA_CM2, VTS,  TEMP,  P,    PDP,  SLIP  ) &
    !$OMP PRIVATE( VISC,  TC0,      DELZ, DELZ1, TOT1, TOT2, FLUX  ) &
    !$OMP PRIVATE( NA,    ThisSpc,  ND,   S                        )

    ! Loop over only the advected dust species
    DO NA = Ind0, Ind1

       ! Look up this species in the species database
       ! NOTE: The advected species are listed first in the main
       ! species list, so it's OK to use the advected species ID
       ! to query the Species Database (bmy, 3/16/17)
       ThisSpc => State_Chm%SpcData(NA)%Info

       ! Get the drydep ID corresponding to this species
       ND      =  ThisSpc%DryDepId

       ! Density [kg/m3] and radius [m]
       DEN     =  ThisSpc%Density
       REFF    =  ThisSpc%Radius

       ! Dp [um] = particle diameter
       DP      =  2e+0_fp * REFF * 1.e+6_fp
       CONST   =  2e+0_fp * DEN * REFF**2 * g0 / 9e+0_fp

       ! Loop over grid latitude and longitude
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Surface area [cm2]
          AREA_CM2 = State_Grid%AREA_M2(I,J) * 1e+4_fp

          ! Initialize settling velocity
          DO L = 1, State_Grid%NZ
             VTS(L) = 0e+0_fp
          ENDDO

          ! Loop over levels
          DO L = 1, State_Grid%NZ

             ! Get P [kPa], T [K], and P*DP
             ! Use moist air pressure for mean free path (ewl, 3/2/15)
             P    = State_Met%PMID(I,J,L) * 0.1e+0_fp
             TEMP = State_Met%T(I,J,L)
             PDP  = P * DP

             !=====================================================
             ! # air molecule number density
             ! num = P * 1d3 * 6.023d23 / (8.314 * Temp)
             !
             ! # gas mean free path
             ! lamda = 1.d6 /
             !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
             !
             ! # Slip correction
             ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * &
             !        exp( -1.1 * Dp / (2. * lamda))) / Dp
             !=====================================================
             ! NOTE, Slip correction factor calculations following
             !       Seinfeld, pp464 which is thought to be more
             !       accurate but more computation required.
             !=====================================================

             ! Slip correction factor as function of (P*dp)
             SLIP = 1e+0_fp + ( 15.60e+0_fp + 7.0e+0_fp * &
                    EXP(-0.059e+0_fp*PDP) ) / PDP

             !=====================================================
             ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
             ! which produce slip correction factor with small
             ! error compared to the above with less computation.
             !=====================================================

             ! Viscosity [Pa s] of air as a function of temp (K)
             VISC = 1.458e-6_fp * (TEMP)**(1.5e+0_fp) / &
                    ( TEMP + 110.4e+0_fp )

             ! Settling velocity [m/s]
             VTS(L) = CONST * SLIP / VISC

          ENDDO

          ! Method is to solve bidiagonal matrix
          ! which is implicit and first order accurate in Z
          DO L = 1, State_Grid%NZ
             TC0(L) = TC(I,J,L,NA)
          ENDDO

          ! We know the boundary condition at the model top
          L           = State_Grid%MaxChemLev
          DELZ        = State_Met%BXHEIGHT(I,J,L)
          TC(I,J,L,NA) = TC(I,J,L,NA) / &
                         ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

          DO L = State_Grid%MaxChemLev-1, 1, -1
             DELZ         = State_Met%BXHEIGHT(I,J,L)
             DELZ1        = State_Met%BXHEIGHT(I,J,L+1)
             TC(I,J,L,NA) = 1.e+0_fp / &
                            ( 1.e+0_fp + DT_SETTL * VTS(L)   / DELZ ) &
                            * (TC(I,J,L,NA)  + DT_SETTL * VTS(L+1) / DELZ1 &
                            *  TC(I,J,L+1,NA) )
          ENDDO

          !========================================================
          ! Dry deposition diagnostic [molec/cm2/s]
          !========================================================
          IF ( State_Diag%Archive_DryDepChm .OR. &
               State_Diag%Archive_DryDep        ) THEN
             !-----------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             ! Dry deposition flux loss [molec/cm2/s]
             !
             ! NOTE: Eventually think about converting this
             ! diagnostic to more standard units [kg/m2/s]
             !-----------------------------------------------------
             ! Initialize
             TOT1 = 0e+0_fp
             TOT2 = 0e+0_fp

             ! Compute column totals of TCO(:) and TC(I,J,:,N)
             DO L = 1, State_Grid%NZ
                TOT1 = TOT1 + TC0(L)
                TOT2 = TOT2 + TC(I,J,L,NA)
             ENDDO

             ! Convert dust flux from [kg/s] to [molec/cm2/s]
             FLUX = ( TOT1 - TOT2 ) / DT_SETTL
             FLUX = FLUX * AVO * ( AIRMW / ThisSpc%Mw_g    ) / &
                    ( AIRMW * 1.e-3_fp ) / AREA_CM2

             ! Drydep flux in chemistry only
             S = State_Diag%Map_DryDepChm%id2slot(ND)
             IF ( S > 0 ) THEN
                State_Diag%DryDepChm(I,J,S) = FLUX
             ENDIF
          ENDIF
       ENDDO
       ENDDO

       ! Nullify pointer
       ThisSpc => NULL()
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    ThisSpc => NULL()
    TC      => NULL()

  END SUBROUTINE DRY_SETTLING
#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settlingbin
!
! !DESCRIPTION: Subroutine DRY\_SETTLINGBIN computes the dry settling of
!  aerosol tracers. Modified for APM simulation. (G. Luo)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                              State_Grid, State_Met, RC)
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    USE PRESSURE_MOD,   ONLY : GET_PCENTER
    USE TIME_MOD,       ONLY : GET_TS_CHEM
    USE PhysConstants
    USE APM_INIT_MOD,   ONLY : APMIDS
    USE APM_INIT_MOD,   ONLY : NCTDST,NDSTB
    USE APM_INIT_MOD,   ONLY : RDST, DENDST
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)   :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)   :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)   :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Aug 2011 - G. Luo - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Local variables
    INTEGER               :: I, J, L, N, K
    INTEGER               :: IDTEMP
    REAL*8                :: DT_SETTL, DELZ,  DELZ1
    REAL*8                :: REFF,     DEN,   CONST
    REAL*8                :: NUM,      LAMDA, FLUX
    REAL*8                :: AREA_CM2, TC0(State_Grid%NZ)
    REAL*8                :: TOT1,     TOT2

    ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
    REAL*8                :: P

    ! Diameter of aerosol [um]
    REAL*8                :: Dp

    ! Pressure * DP
    REAL*8                :: PDp

    ! Temperature (K)
    REAL*8                :: TEMP

    ! Slip correction factor
    REAL*8                :: Slip

    ! Viscosity of air (Pa s)
    REAL*8                :: Visc

    ! Settling velocity of particle (m/s)
    REAL*8                :: VTS(State_Grid%NZ)
    REAL*8                :: MASS(State_Grid%NZ)
    REAL*8                :: OLD(State_Grid%NZ,NCTDST)

    ! Make a pointer to the tracer array
    REAL*8, POINTER       :: Spc(:,:,:,:)

    !=================================================================
    ! DRY_SETTLINGBIN begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Point to Spc
    Spc => State_Chm%species

    ! Aerosol settling timestep [s]
    DT_SETTL = GET_TS_CHEM()

    IDTEMP = APMIDS%id_DSTBIN1+NDSTB-1

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, K, DEN, REFF, DP )       &
    !$OMP PRIVATE( CONST, VTS, TEMP, P, PDP, SLIP )     &
    !$OMP PRIVATE( MASS, OLD, VISC, TC0, DELZ, DELZ1  ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       DO L = 1, State_Grid%NZ
          MASS(L) = SUM(Spc(I,J,L,APMIDS%id_DSTBIN1:IDTEMP))
          DO K = 1, NCTDST
             OLD(L,K) = Spc(I,J,L,(APMIDS%id_CTDST+K-1))
             Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = 0.D0
          ENDDO
       ENDDO

       ! Loop over aerosol bins
       DO N = 1, NDSTB
          DO L = 1, State_Grid%NZ
             TC0(L) = Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1))

             IF(TC0(L)>1.D-30)THEN
                ! Initialize
                DEN   = DENDST(N)*1.d3
                REFF  = RDST(N)
                DP    = 2D0 * REFF * 1.D6 ! Dp [um] = particle diameter
                CONST = 2D0 * DEN * REFF**2 * g0 / 9D0

                ! Get P [kPa], T [K], and P*DP
                P    = GET_PCENTER(I,J,L) * 0.1d0
                TEMP = State_Met%T(I,J,L)
                PDP  = P * DP

                ! Slip correction factor as function of (P*dp)
                SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP

                ! Viscosity [Pa s] of air as a function of temp (K)
                VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                ! Settling velocity [m/s]
                VTS(L) = CONST * SLIP / VISC
             ELSE
                VTS(L) = 0.D0
             ENDIF

          ENDDO

          ! Method is to solve bidiagonal matrix
          ! which is implicit and first order accurate in Z
          L = State_Grid%NZ
          IF(MASS(L)>1.D-30)THEN
             DELZ = State_Met%BXHEIGHT(I,J,L)

             Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) = &
                  Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

             DO K = 1, NCTDST
                Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTDST+K-1))+ &
                     OLD(L,K)*TC0(L)/MASS(L) / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
             ENDDO
          ENDIF

          DO L = State_Grid%NZ-1, 1, -1
             IF((MASS(L)*MASS(L+1))>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) = 1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * (Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) &
                     + DT_SETTL * VTS(L+1) / DELZ1 * TC0(L+1) )

                DO K = 1, NCTDST
                   Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = &
                        Spc(I,J,L,(APMIDS%id_CTDST+K-1))+ &
                        1.e+0_fp / ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                        * (OLD(L,K)*TC0(L)/MASS(L) &
                        + DT_SETTL * VTS(L+1) / DELZ1 &
                        * OLD(L+1,K)*TC0(L+1)/MASS(L+1) )
                ENDDO
             ELSE IF(MASS(L)>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) = 1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * (Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) &
                     + DT_SETTL * VTS(L+1) / DELZ1 * TC0(L+1) )

                DO K = 1, NCTDST
                   Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = &
                        Spc(I,J,L,(APMIDS%id_CTDST+K-1))+ &
                        1.e+0_fp / ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                        * OLD(L,K)*TC0(L)/MASS(L)
                ENDDO
             ELSE IF(MASS(L+1)>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) = 1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * (Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) &
                     + DT_SETTL * VTS(L+1) / DELZ1 * TC0(L+1) )

                DO K = 1, NCTDST
                   Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = &
                        Spc(I,J,L,(APMIDS%id_CTDST+K-1))+ &
                        1.e+0_fp / ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                        * DT_SETTL * VTS(L+1) / DELZ1 &
                        * OLD(L+1,K)*TC0(L+1)/MASS(L+1)
                ENDDO
             ENDIF

          ENDDO

       ENDDO

       DO L = 1, State_Grid%NZ
          DO K = 1, NCTDST
             Spc(I,J,L,(APMIDS%id_CTDST+K-1)) = &
                  MAX(1.d-30,Spc(I,J,L,(APMIDS%id_CTDST+K-1)))
          ENDDO
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Clear the pointer
    Spc => NULL()

  END SUBROUTINE DRY_SETTLINGBIN
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdust_online
!
! !DESCRIPTION: Subroutine RDUST\_ONLINE reads global mineral dust
!  concentrations as determined by P. Ginoux.  Calculates dust optical
!  depth at each level for the FAST-J routine "set\_prof.f".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RDUST_ONLINE(  Input_Opt, State_Chm, State_Diag, State_Grid, &
                            State_Met, DUST,      ODSWITCH,   RC )
!
! !USES:
!
    USE CMN_FJX_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid    ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met     ! Meteorology State object
    REAL(fp),       INTENT(IN)  :: DUST(State_Grid%NX,State_Grid%NY, &
                                        State_Grid%NZ,NDUST) !Dust [kg/m3]
    INTEGER,        INTENT(IN)  :: ODSWITCH
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Martin, R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL           :: LINTERP
    INTEGER           :: I,       J,       L,       N
    INTEGER           :: NOUT,    S,       W,       IWV
    INTEGER           :: IIWV,    NWVS,    IDST,    IWV1
    INTEGER           :: IWV2
    REAL(fp)          :: XTAU,    XRH,     CRITRH,  ACOEF
    REAL(fp)          :: BCOEF,   LOGTERM

    ! Arrays
    LOGICAL           :: LINTERPARR(Input_Opt%NWVSELECT)
    REAL(fp)          :: MSDENS(NDUST)
    REAL(fp)          :: tempOD(State_Grid%NX,State_Grid%NY,               &
                                State_Grid%NZ,NDUST,         3)

    ! Pointers
    REAL(fp), POINTER :: ERADIUS(:,:,:,:)
    REAL(fp), POINTER :: TAREA(:,:,:,:)
    REAL(fp), POINTER :: WERADIUS(:,:,:,:)
    REAL(fp), POINTER :: WTAREA(:,:,:,:)

    !=================================================================
    ! RDUST_ONLINE begins here!
    !=================================================================

    ! Assume success
    RC        =  GC_SUCCESS

    ! Initialize pointers
    ERADIUS   => State_Chm%AeroRadi     ! Aerosol Radius     [cm]
    TAREA     => State_Chm%AeroArea     ! Aerosol Area       [cm2/cm3]
    WERADIUS  => State_Chm%WetAeroRadi  ! Wet Aerosol Radius [cm]
    WTAREA    => State_Chm%WetAeroArea  ! Wet Aerosol Area   [cm2/cm3]

    ! Index for dust in ODAER and LUT arrays
    ! Dust properties are saved to different indices in RD_AOD for
    ! UCX vs tropchem simulations
    IF ( Input_Opt%LUCX ) THEN
       IDST   = 8
    ELSE
       IDST   = 6
    ENDIF

    ! Dust density
    MSDENS(1) = 2500.0_fp
    MSDENS(2) = 2500.0_fp
    MSDENS(3) = 2500.0_fp
    MSDENS(4) = 2500.0_fp
    MSDENS(5) = 2650.0_fp
    MSDENS(6) = 2650.0_fp
    MSDENS(7) = 2650.0_fp

    ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
    CRITRH = 35.0e+0_fp   ! [%]

    !=================================================================
    ! Convert concentration [kg/m3] to optical depth [unitless].
    !
    ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) /
    !           ( MSDENS * RAA * 1e-6 )
    ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
    !
    !  Units ==> DUST     [ kg/m3    ]
    !            MSDENS   [ kg/m3    ]
    !            RAA      [ um       ]
    !            BXHEIGHT [ m        ]
    !            QAA      [ unitless ]
    !            ODMDUST  [ unitless ]
    !
    ! NOTES:
    ! (1) Do the calculation at QAA(IND999,:) (i.e. 999 nm).
    !=================================================================
    ! DAR Oct 2012
    ! if the call is from chemistry (ODSWITCH=0) only need one wavelength
    ! if radiation on (LRAD=TRUE) then cycle over all LUT wavelengths
    ! if radiation off but call is for AOD then cycle over the
    ! wavelength required or the wavelengths between which to
    ! interpolate

    IF (ODSWITCH .EQ. 0) THEN
       NWVS   = 1
    ELSE
       IF ( Input_Opt%LRAD ) THEN !Loop over all RT wavelengths (30)
          ! plus any required for calculating the AOD
          NWVS   = NWVAART+NWVREQUIRED
       ELSE                       !Loop over wavelengths needed (from RD_AOD)
          NWVS   = NWVREQUIRED
       ENDIF
    ENDIF

    DO IIWV = 1, NWVS
       ! get current wavelength index
       ! (specified by user or all wavelengths for RRTMG)
       ! (1000nm always used for FAST-J and stored in first index)
       IF (ODSWITCH .EQ. 0) THEN
          ! only doing for 1000nm i.e. IWV=10 in LUT
          ! N.B. NWVS is fixed to 1 above - only one wavelength
          IWV=IWV1000
       ELSE
          IF ( Input_Opt%LRAD ) THEN
             ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
             ! but add on any others required
             IF (IIWV.LE.NWVAART) THEN
                !index of RRTMG wavelengths starts after the standard NWVAA0
                !(currently NWVAA0=11, set in CMN_FJX_mod based on the .dat
                !LUT)
                IWV = IIWV+NWVAA0
             ELSE
                !now we calculate at wvs for the requested AOD
                !offset index by NWVAART i.e. start from 1
                IWV = IWVREQUIRED(IIWV-NWVAART)
             ENDIF
          ELSE
             ! IWVREQUIRED lists the index of requires standard wavelengths
             IWV = IWVREQUIRED(IIWV)
          ENDIF
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, N )
       DO N = 1, NDUST
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! dust stored in the IDST species bin of LUT variables
          ODMDUST(I,J,L,IWV,N) = 0.75e+0_fp * &
                                 State_Met%BXHEIGHT(I,J,L) * &
                                 DUST(I,J,L,N) * QQAA(IWV,N,IDST)  / &
                                 ( MSDENS(N) * RDAA(N,IDST) * 1.0e-6_fp)

#ifdef RRTMG
          !add dust optics to the RT code arrays
          !SSA and ASYM copying seems a little redundant...
          !will keep this way for uniformity for now but
          !possibly could deal with SSA and ASYM in RT module
          RTODAER(I,J,L,IWV,NAER+2+N) = ODMDUST(I,J,L,IWV,N)
          RTSSAER(I,J,L,IWV,NAER+2+N) = SSAA(IWV,N,IDST)
          RTASYMAER(I,J,L,IWV,NAER+2+N) = ASYMAA(IWV,N,IDST)
#endif

       ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO !wavelength loop

    !==============================================================
    ! Calculate Dust Surface Area
    !
    ! Units ==> DUST     [ kg dust/m^3 air    ]
    !           MSDENS   [ kg dust/m^3 dust   ]
    !           RAA      [ um                 ]
    !           TAREA    [ cm^2 dust/cm^3 air ]
    !           ERADIUS  [ cm                 ]
    !
    ! NOTE: first find volume of dust (cm3 dust/cm3 air), then
    !       multiply by 3/radius to convert to surface area in cm2
    !
    ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for
    ! the NDUST FAST-J dust wavelength bins (read into DUST)
    !==============================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, XRH )
    DO N = 1, NDUST
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ERADIUS(I,J,L,N) = RDAA(N,IDST) * 1.0e-4_fp

       TAREA(I,J,L,N)   = 3.e+0_fp / ERADIUS(I,J,L,N) * &
                          DUST(I,J,L,N) / MSDENS(N)

       ! Archive WTAREA and WERADIUS when RH > 35%  (tmf, 6/13/07)
       ! Get RH
       XRH                = State_Met%RH( I, J, L )  ! [%]
       WTAREA(I,J,L, N)   = 0.e+0_fp
       WERADIUS(I,J,L, N) = 0.e+0_fp

       IF ( XRH >= CRITRH ) THEN
          WTAREA(I,J,L, N)   = TAREA(I,J,L, N)
          WERADIUS(I,J,L, N) = ERADIUS(I,J,L, N)
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( State_Diag%Archive_AOD .AND. ODSWITCH .EQ. 1 ) THEN

       ! Initialize local arrays
       tempOD = 0.0_fp
       LINTERPARR(:) = .FALSE.

       ! Interpolate?
       DO W = 1, Input_Opt%NWVSELECT
          IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
             LINTERPARR(W) =.FALSE.
          ELSE
             LINTERPARR(W) =.TRUE.
          ENDIF
       ENDDO

       ! Loop over dust bins, # of wavelengths, and all grid cells
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, N, W, NOUT, LINTERP, S )
       DO W = 1, Input_Opt%NWVSELECT
       DO N = 1, NDUST
       DO L = 1, State_Grid%MaxChemLev
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          IF ( .not. LINTERPARR(W) ) THEN

             tempOD(I,J,L,N,W) = ODMDUST(I,J,L,IWVSELECT(1,W),N)

          ELSEIF ((ODMDUST(I,J,L,IWVSELECT(1,W),N).GT.0).AND. &
                  (ODMDUST(I,J,L,IWVSELECT(2,W),N).GT.0)) THEN

             ! Interpolate using angstrom exponent between
             ! closest available wavelengths (coefs pre-calculated
             ! in CALC_AOD)
             ! AOD sometimes zero (if Q zero), must catch this!
             tempOD(I,J,L,N,W) = &
                  ODMDUST(I,J,L,IWVSELECT(2,W),N)* &
                  ACOEF_WV(W)**(BCOEF_WV(W)*LOG( &
                  ODMDUST(I,J,L,IWVSELECT(1,W),N)/ &
                  ODMDUST(I,J,L,IWVSELECT(2,W),N)))

          ENDIF

          !---------------------------------------------------
          ! Set size-resolved dust optical depth diagnostic
          !---------------------------------------------------
          IF ( State_Diag%Archive_AODDustWL1 .AND. W == 1 ) THEN
             S = State_Diag%Map_AODDustWL1%id2slot(N)
             IF ( S > 0 ) THEN
                State_Diag%AODDustWL1(I,J,L,S) = tempOD(I,J,L,N,1)
             ENDIF
          ENDIF

          IF ( State_Diag%Archive_AODDustWL2 .AND. W == 2 ) THEN
             S = State_Diag%Map_AODDustWL2%id2slot(N)
             IF ( S > 0 ) THEN
                State_Diag%AODDustWL2(I,J,L,S) = tempOD(I,J,L,N,2)
             ENDIF
          ENDIF

          IF ( State_Diag%Archive_AODDustWL3 .AND. W == 3 ) THEN
             S = State_Diag%Map_AODDustWL3%id2slot(N)
             IF ( S > 0 ) THEN
                State_Diag%AODDustWL3(I,J,L,S) = tempOD(I,J,L,N,3)
             ENDIF
          ENDIF

       ENDDO ! longitude
       ENDDO ! latitude
       ENDDO ! level
       ENDDO ! bin
       ENDDO ! wavelength
       !$OMP END PARALLEL DO

       !---------------------------------------------------
       ! Set dust optical depth diagnostic
       !---------------------------------------------------
       IF ( State_Diag%Archive_AODDust ) THEN
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L )
          DO L = 1, State_Grid%MaxChemLev
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             State_Diag%AODDust(I,J,L) = SUM( tempOD(I,J,L,:,:) )
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF

    ENDIF

    ! Archive total dust surface area (sum across all bins)
    IF ( State_Diag%Archive_AerSurfAreaDust ) THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%MaxChemLev
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Diag%AerSurfAreaDust(I,J,L) = SUM( TAREA(I,J,L,1:NDUST) )
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Free pointers
    ERADIUS  => NULL()
    TAREA    => NULL()
    WERADIUS => NULL()
    WTAREA   => NULL()

  END SUBROUTINE RDUST_ONLINE
!EOC!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_dust
!
! !DESCRIPTION: Subroutine INIT\_DUST allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_DUST( Input_Opt,  State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : PI, AIRMW
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
#ifdef TOMAS
    USE TOMAS_MOD,          ONLY : IBINS, Xk
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  30 Mar 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, N
#ifdef TOMAS
    INTEGER                :: I
    INTEGER                :: BIN
#endif

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: ThisLoc

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! INIT_DUST begins here!
    !=================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  '-> at INIT_DUST (in module GeosCore/dust_mod.F)'
    SpcInfo => NULL()

    ! Exit immediately if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    !=================================================================
    ! INITIALIZATION SECTION FOR STANDARD GEOS-CHEM SIMULATION
    !=================================================================

    !-----------------------------------------------------------------
    ! DST1 - DST4 species (most fullchem + aerosol simulations)
    !
    ! NOTE: We can consider removing DUSTDEN and DUSTREFF from
    ! the Input_Opt object at a later point.  These can be directly
    ! obtained from the species database object now. (bmy, 6/20/16)
    !-----------------------------------------------------------------
    id_DST1 =  Ind_('DST1'    )
    id_DST2 =  Ind_('DST2'    )
    id_DST3 =  Ind_('DST3'    )
    id_DST4 =  Ind_('DST4'    )

    !-----------------------------------------------------------------
    ! DAL1 - DAL4 and SO4d1 - SO4d4 species (acid uptake sims only)
    !-----------------------------------------------------------------
    id_DAL1 =  Ind_('DSTAL1'  )
    id_DAL2 =  Ind_('DSTAL2'  )
    id_DAL3 =  Ind_('DSTAL3'  )
    id_DAL4 =  Ind_('DSTAL4'  )

    ! Error check the dust uptake species
    IF ( Input_Opt%LDSTUP ) THEN
       IF ( id_DAL1 < 0 .or. id_DAL2 < 0   .or. &
            id_DAL3 < 0 .or. id_DAL4 < 0 ) THEN
          ErrMsg = 'Dust uptake species are undefined!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#ifdef TOMAS
    !=================================================================
    ! INITIALIZATION SECTION FOR TOMAS MICROPHYSICS
    !=================================================================

    ! TOMAS species ID flags
    id_DUST1 = Ind_('DUST1')
    id_NK1   = Ind_('NK1'  )

    !----------------------------------
    ! Set up FRAC_S (only for Ginoux)
    !----------------------------------
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% NOTE: LDEAD has not been set by HCOI_GC_INIT yet. This code needs %
    !% to be moved or modified accordingly (mps, 4/9/15)                 %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% NOTE: This code replicates the obsolete Input_Opt%IDDEP field,    %
    !% which has since been removed.  TOMAS TEAM please take note.       %
    !% (bmy, 3/17/17)                                                    %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF ( .not. Input_Opt%LDEAD ) THEN

       ! Initialize
       TOMAS_IDDEP = 0

       DO BIN = 1, IBINS
          DO N   = 1, State_Chm%nDryDep
             IF ( State_Chm%Map_DryDep(N) == ( id_NK1-1+BIN ) )THEN
                TOMAS_IDDEP(BIN) = N
                GOTO 100
             ENDIF
          ENDDO
100       CONTINUE
       ENDDO

    ENDIF
#endif

#if defined( ESMF_ ) || defined( TOMAS )
    ! EXPERIMENTAL: For archiving the dust source for GCHP
    !
    ! Changed to use the ESMF_ flag and not EXTERNAL_GRID/EXTERNAL_FORCING, as
    ! WRF-GC which uses these flags does not require SRCE_FUNC (hplin, 1/22/19)
    ALLOCATE( SRCE_FUNC(State_Grid%NX,State_Grid%NY,3), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRCE_FUNC' )
    SRCE_FUNC = 0.e+0_fp
#endif

  END SUBROUTINE INIT_DUST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_dust
!
! !DESCRIPTION: Subroutine CLEANUP\_DUST deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_DUST
!
! !REVISION HISTORY:
!  30 Mar 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_DUST begins here!
    !=================================================================
    IF ( ALLOCATED( FRAC_S    ) ) DEALLOCATE( FRAC_S    )
    IF ( ALLOCATED( SRCE_FUNC ) ) DEALLOCATE( SRCE_FUNC )

  END SUBROUTINE CLEANUP_DUST
!EOC
END MODULE DUST_MOD
