! THIS MODULE INCLUDES THREE GROUP OF GLOBAL VARIABLES THAT ARE USED FOR
! GEOMETRY PROJECTED METHOD
      module user_global_variable_FE
        REAL, ALLOCATABLE, SAVE :: CENTROIDS(:),ELEM_R, ELEM_VOL
        INTEGER, SAVE :: NO_EL, NONODE
        INTEGER, SAVE :: DIM, N_BAR_DOFS
        REAL, SAVE :: RHO_MIN = 0.01
        REAL, ALLOCATABLE :: PENALIZED_ELEM_DENS, ELEM_DENS        
      end module 
C
      module user_global_variable_GEOM
      ! BAR, DESIGN PARAMETER
        INTEGER, SAVE :: N_POINT, N_BAR
        REAL, ALLOCATABLE, SAVE :: POINT_MATRIX(:,:), BAR_MATRIX(:,:)
        REAL, SAVE :: MIN_BAR_RADIUS, MAX_BAR_RADIUS
        REAL, ALLOCATABLE :: CURRENT_POINT_MATRIX(:,:),
     +    CURRENT_BAR_MATRIX(:,:) 
      end module
C
      module user_global_variable_OPT
        CHARACTER, SAVE :: SMOOTH_MAX_SCHEME*18, PENALIZATION_SCHEME*6
        INTEGER, SAVE :: SMOOTH_MAX_PARAM, PENALIZATION_PARAM
      ! SCALING VARIABLE
        REAL, ALLOCATABLE :: POINT_SCALE(:)
        REAL :: POINT_MIN, RADIUS_SCALE, RADIUS_MIN
      ! DV VARIABLE
        INTEGER, SAVE :: N_DV
        INTEGER, ALLOCATABLE :: POINT_DV(:), SIZE_DV(:), RADIUS_DV(:) 
        LOGICAL, SAVE :: DV_SCALING
      ! GEOMETRY PROJECTED DENSITY
C        REAL, ALLOCATABLE :: ELEM_DENS(:)
C
      ! SENSITIVITIES OF DENSITY GEOMETRY PROJECTED 
        REAL, ALLOCATABLE :: DELEM_DENS_DDV(:),
     +                       DPENALIZED_ELEM_DENS_DDV(:)              
        ! PARAMETER OF ELEMENT
       ! REAL, ALLOCATABLE, SAVE :: ELEM_R(:) 
      end module

C-----------------------------------------------------------------------
C*********************  INPUT_FILE *************************************
C-----------------------------------------------------------------------
C
C---------------------------------------------------------------------------------------- 
C USING SUBROUTINE UEXTERNALDB TO IMPORT AND EXPORT DATA IN ABAQUS
C----------------------------------------------------------------------------------------

      SUBROUTINE UEXTERNALDB(LOP,ANALYSISPOS, LRESTART, TIME, DTIME, 
     +       KSTEP, KINC)
C----------------------------------------------------------------------------------------
C THIS SUBROUTINE IS USED TO READ INPUT DATA AND OUPUT RESULTS 
C FIRST MISSION: INPORT THE COORDINATE AND DATA OF BARS
C----------------------------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!     ###############################           
!     Declare SUBROUTINE variables
!     ###############################
C
      DIMENSION TIME(2)
      INTEGER :: ANALYSISPOS, KSTEP
      CHARACTER*6O
      INTEGER ERROR      
!     ###############################           
!     Declare LOCAL variables
!     ###############################
C
      INTEGER LENJOBNAME,LENOUTDIR,NDIM,NINT,NINTS
      CHARACTER*256 JOBNAME,OUTDIR,FILENAME 

C----------------------------------------------------------------------------------------
C   EXECUTION AREA
C----------------------------------------------------------------------------------------
      
      IF (LOP.EQ.0) THEN

      CALL GETJOBNAME(JOBNAME,LENJOBNAME)
      CALL GETOUTDIR(OUTDIR, LENOUTDIR)
      ! AS FILEPATH THE ACTUAL ABAQUS OUTPUT DIRECTORY WILL BE USED.
      ! GET ABAQUS OUTPUT DIRECTORY FOR ANALYSIS
      OPEN(UNIT = 3,
     1 FILE=
     2 OUTDIR(:LENOUTDIR) // '/' // 'ZDATA.txt'
     3 ,STATUS = 'UNKNOWN', FORM = 'FORMATTED', IOSTAT = ERROR)

      OPEN(UNIT = 4,
     1 FILE=
     2 OUTDIR(:LENOUTDIR) // '/' // 'ZRESULTS.txt'
     3  ,STATUS = 'REPLACE')

      WRITE(4,*) ERROR

c      CALL get_elem_con_2D
     

      IF (ERROR.EQ.0) THEN
      CALL INDATA
      END IF
      
      
      END IF
      
            
     
      RETURN
      
      END SUBROUTINE UEXTERNALDB 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine UVARM
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     + NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     + JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      CHARACTER*80 CMNAME,ORNAME
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION JMAC(*),JMATYP(*),COORD(*)
!-----Data from ABAQUS
      CHARACTER*3 FLGRAY(15)
      DIMENSION ARRAY(15),JARRAY(15)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      REAL*8 EVOLC,EVOLF
      COMMON/EVOL/ EVOLC,EVOLF
!-----------------------------------------------------------------------
!     Access integration point volume
!-----------------------------------------------------------------------
C
      CALL GETVRM('IVOL',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
C
      EVOLC = ARRAY(1)
      EVOLF = EVOLC      
!-----------------------------------------------------------------------
!     Update user-defined variables
!-----------------------------------------------------------------------
      UVAR(1) = EVOLF
      

!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE UVARM     
      
C---------------------------------------------------------------------------------------
C Start of Base code, DO NOT change
C---------------------------------------------------------------------------------------
C  
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C      
      use user_global_variable_FE
      use user_global_variable_OPT      
      use user_global_variable_GEOM
C
!     COPYRIGHT (2012,2020,2021) EVER J. BARBERO, ALL RIGHTS RESERVED
!     Ex. 3.13 Finite Element Analysis of Composite Materials with Abaqus
C
      INCLUDE 'ABA_PARAM.INC'
C
!
      CHARACTER*80 CMNAME   !CMNAME<=MATERL in R10
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C  
C
C     ! LOCAL VARIABLE DECLARATION
      REAL X,Y,Z
      REAL DSTRESS(NTENS), DSTRAIN(NTENS)      
      REAL DC_DDV(N_DV)
C
      REAL V_E, V, DVOLFRAC_DDV(N_DV)
      REAL myE, myNu, myLambda, myMu
C
C     GET ELEMENT INFORMATION
      ! COMMON VARIABLE OF ELEMENT INFORMATION
C      
      INTERFACE
         SUBROUTINE COMPUTE_ELEM_PARAMETER(NOEL)
           INTEGER, INTENT(IN) :: NOEL
         END SUBROUTINE COMPUTE_ELEM_PARAMETER
      END INTERFACE       

      INTEGER Q

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0) 
C
C      CALL COMPUTE_ELEM_PARAMETER(NOEL)            
      X = CENTROIDS(1)
      Y = CENTROIDS(2)
      Z = CENTROIDS(3)

C------------------------------------------------------------------------------
C     APPLYING THE GEOMETRY PROJECTED DENSITY METHOD    
C------------------------------------------------------------------------------

C      CALL PROJECT_ELEMENT_DENSITY

C---------------------------------------------------------------------------------------

C Start of USER code

C---------------------------------------------------------------------------------------
C           PROPS(NPROPS)
C           User-specified array of material constants associated with this user material.
C
C           NPROPS
C           User-defined number of material constants associated with this user material.
C
C           PROPS(1) = Elasticity Modulus
C           PROPS(2) = Poisson's Ratio
C---------------------------------------------------------------------------------------

       E=PROPS(1)*PENALIZED_ELEM_DENS
       v=PROPS(2)
       myE=PROPS(1)
       myNu=PROPS(2)       

      IF (DIM.EQ.2) THEN  
      DO I = 1, NTENS
       DO J = 1, NTENS
       DDSDDE(I,J)= ZERO
       END DO
      END DO       


C           DDSDDE(NTENS,NTENS)
C           Jacobian matrix of the constitutive model

C           NTENS
C           Size of the stress or strain component array (NDI + NSHR).

C           NDI
C           Number of direct stress components at this point.

C           NSHR
C           Number of engineering shear stress components at this point.
                           
C          Stiffness Matrix for Plane Stress
C		
      DDSDDE(1,1) = E/(ONE-v*v)
      DDSDDE(2,2) = E/(ONE-v*v)
      DDSDDE(3,3) = E*(ONE-v)/(TWO*(ONE-v*v))
         
      DDSDDE(1,2) = E*v/(ONE-v*v)
      DDSDDE(2,1) = E*v/(ONE-v*v)

C     Calculation of STRESSES 
C
C     STRESS(NTENS)
C     This array is passed in as the stress tensor at the beginning of the increment and must be updated
C     in this routine to be the stress tensor at the end of the increment.
C      
C     DSTRAN(NTENS)
C     Array of strain increments.
C
       DO I=1, NTENS
       DO J=1, NTENS
C            STRESS(J) = 0.0
           STRESS(J)=STRESS(J)+DDSDDE(J, I)*DSTRAN(I)
       END DO
       END DO
       ELSE
         myMu= myE / (TWO * (ONE + myNu)) 
         myLambda= (myE * myNu) / ((ONE + myNu) * (ONE - TWO * myNu))

C     Stiffness Matrix for Isotropic Isothermal Elasticity with Lamé parameters         
         DO i=1, NDI
             DO j=1, NDI
                 DDSDDE(j, i)= myLambda
             END DO
         DDSDDE(i, i)= TWO * myMu + myLambda
         END DO
C        
         DO i=NDI+1, NTENS
             DDSDDE(i ,i)= myMu
         END DO
C           Calculation of STRESSES 
C
C           STRESS(NTENS)
C           This array is passed in as the stress tensor at the beginning of the increment and must be updated
C           in this routine to be the stress tensor at the end of the increment.
C      
C           DSTRAN(NTENS)
C           Array of strain increments.
         DO i=1, NTENS
             DO j=1, NTENS
                 STRESS(j)=STRESS(j) + DDSDDE(j, i) * DSTRAN(i)
             END DO
         END DO                              
       ENDIF
C       DO I = 1, NTENS
C            DSTRESS(I) = STRESS(I)
C            DSTRAIN(I) = DSTRAN(I)
C      END DO      
C
C     ! ELASTIC STRAIN ENGERGY
C     
      DEE = 0.0
      DO I = 1, NTENS
            DEE = DEE + 0.5*STRESS(I)*DSTRAN(I)
      END DO

      SSE = DEE
      Q = 1

      DC_PENALIZED_ELEM_DENS = -2.0*DEE*(ELEM_VOL)
     1  *(Q*PENALIZED_ELEM_DENS**(Q-1)/PENALIZED_ELEM_DENS)

      DO N = 1, N_DV
        DC_DDV(N) = DC_PENALIZED_ELEM_DENS*DPENALIZED_ELEM_DENS_DDV(N)
      END DO

C     STORE COMPLIANCE AND THE SENSITIVITIES OF COMPLIANCE
C
      DO N = 1, N_DV
        STATEV(N) = DC_DDV(N)
      END DO
      STATEV(N_DV+1) = 2.0*DEE*(ELEM_VOL)

C      STORE COMBINED PROJECTED DENSITY

      STATEV(N_DV+2) = PENALIZED_ELEM_DENS

C     COMPUTE THE VOLUME FRACTION
      V_E = ELEM_VOL      
      V = (ELEM_VOL)*ELEM_DENS
      STATEV(500+N_DV+1) = V
      STATEV(500+N_DV+2) = V_E
C     COMPUTE THE DESIGN SENSITIVITY

      DO N = 1, N_DV
        DVOLFRAC_DDV(N) = V_E*DELEM_DENS_DDV(N)      
      END DO
C     STORE THE DESIGN SENSITIVITY
      DO N = 1, N_DV
        STATEV(500+N) = DVOLFRAC_DDV(N)
      END DO


   50 FORMAT( I10)
C   60 FORMAT(10F10.4)
      RETURN
      END SUBROUTINE UMAT
C-----------------------------------------------------------------------------------      
C-----------------------------------------------------------------------           |
C*********************  INPUT_FILE *************************************           |
C-----------------------------------------------------------------------           |      
      SUBROUTINE INDATA                                                            
C                                                                                  |
C**********************************************************************            |
C*                                                                    *            |
C** SUBPROGRAM TO INPUT DATA FOR BARS **                      *                    |
C*                                                                    *            |
C**********************************************************************            |
C     USE MODULE
      use user_global_variable_FE
      use user_global_variable_OPT      
      use user_global_variable_GEOM
C     LOCAL VARIABLE DECLARATION SECTION                 
      CHARACTER TITLE(7)*10, CODE*20, CODE2*16,
     1 CODE3*19, CODE4*18
C     EXECUTION SECTION      
      READ(3,10) TITLE
C   
      READ(3,11) N_POINT, N_BAR
C 
      WRITE(4,10) TITLE
      WRITE(4,22) NN_POINT, N_BAR
C     ALLOCATE POINT_MATRIX AND BAR_MATRIX
      IF(.NOT. ALLOCATED(POINT_MATRIX)) 
     +     ALLOCATE(POINT_MATRIX(N_POINT,3))
          IF(.NOT. ALLOCATED(BAR_MATRIX)) 
     +     ALLOCATE(BAR_MATRIX(N_BAR,5))            
      CALL READ_DATA
C     ALLOCATE MEMORY GLOBAL VARIABLE      
      IF(.NOT. ALLOCATED(DELEM_DENS_DDV)) 
     +      ALLOCATE(DELEM_DENS_DDV(N_DV))
      IF(.NOT. ALLOCATED(DPENALIZED_ELEM_DENS_DDV))
     +      ALLOCATE(DPENALIZED_ELEM_DENS_DDV(N_DV))      
C      
   10 FORMAT(7A10)
   11 FORMAT(2I5)
   22 FORMAT(/27H NUMB OF POINTS - N_POINT =, I5,
     1        24H, NUMB OF BARS - N_BAR =, I5/)
C
      RETURN           
      END SUBROUTINE INDATA
      
C------------------------------------------------------------------------------
      SUBROUTINE READ_DATA
c------------------------------------------------------------------------------      
      use user_global_variable_FE
      use user_global_variable_OPT      
      use user_global_variable_GEOM
c                   
      CHARACTER TITLE(7)*10, CODE*20,CODE2*16,
     1 CODE3*19, CODE4*18, DIMENSION*19

      REAL NU_VOID

C     ! EXECUTION REGION      
      READ(3,50) (POINT_MATRIX(N,1), POINT_MATRIX(N,2), POINT_MATRIX(N,3)
     1 , N = 1,N_POINT)
C 
      READ(3,60) (BAR_MATRIX(N,1), BAR_MATRIX(N,2), BAR_MATRIX(N,3),
     1 BAR_MATRIX(N,4), BAR_MATRIX(N,5), N = 1, N_BAR)
C
      READ(3,2) DIM
      READ(3,1) MIN_BAR_RADIUS
      READ(3,1) MAX_BAR_RADIUS
      READ(3,*) SMOOTH_MAX_SCHEME
      READ(3,2) SMOOTH_MAX_PARAM       
      READ(3,*) PENALIZATION_SCHEME
      READ(3,2) PENALIZATION_PARAM          
C
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,*) 'COORDINATE OF POINTS'
      WRITE(4,*) '-----------------------------------------------------'   
      WRITE(4,50)(POINT_MATRIX(N,1), POINT_MATRIX(N,2), POINT_MATRIX(N,3)
     1 , N = 1,N_POINT)
C
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,*) 'DESIGN VARIABLE OF BARS' 
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,60)(BAR_MATRIX(N,1), BAR_MATRIX(N,2), BAR_MATRIX(N,3),
     1 BAR_MATRIX(N,4), BAR_MATRIX(N,5), N = 1,N_BAR) 
C      

C
C MATERIAL INFORMATION
C
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,*) 'INITIAL DESGIN (GEOMETRY)'
      WRITE(4,*) '-----------------------------------------------------'     
      ! SPECIFY DIMENSION                
C      DIM = 2
      N_BAR_DOFS = 2*DIM
      WRITE(4,20)'DIMENSION', DIM
C     INITIAL DESIGN (GEOMETRY)
C     YOU MUST SPECIFY THE LOWER AND UPPER BOUNDS ON THE BAR RADIUS. IF 
C     YOU WANT A DESIGN WITH FIXED BAR RADII, SIMPLY SET BOTH FIELDS TO
C     THE SAME VALUE.      
C      MIN_BAR_RADIUS = 0.5
C      MAX_BAR_RADIUS = 0.501
      WRITE(4,10)'MIN_BAR_RADIUS', MIN_BAR_RADIUS
      WRITE(4,10)'MAX_BAR_RADIUS', MAX_BAR_RADIUS
C
      ! SPECIFY THE YOUNG'S MODULUS AND POISSON RATIO 
      RHO_MIN = 0.01 ! MINIMUM DENSITY IN VOID REGION
      WRITE(4,40) 'MINIMUM_DENSITY', RHO_MIN
      NU_VOID = 0.3  ! POISSON RATIO OF THE DESIGN MATERIAL
C-----------------------------------------------------------------------
C     GEOMETRY PROJECTION PARAMETERS
C     SAMPLE WINDOW RADIUS. BY DEFAULT, THE ROUTINE INIT_OPTIMIZATION
C     COMPUTE THIS RADIUS AS SQRT(DIM)/2 * ELEM_VOL^(1/DIM), WHICH
C     CORRESPONDS TO THE RADIUS OF THE CIRCLE (IN 2-D) OR SPHERE (IN 3D)
C     THAT CIRCUMSCRIBES A SQUARE OF CUBE ELEMENT. IF YOU UNCOMMENT THAT
C     FOLLOWING LINE, YOU CAN OVERRIDE THIS VALUE (IF YOU DO, THEN VALUE
C     CORRESPOINGS TO THE ACTUAL RADIUS DIMENSION).
C
C     TYPE THE SMOOTH MAXIMUM FUNCTION FOR AGGREGATION (BOOLEAN UNION)
C     OPTION ARE: 'MOD_P-NORM', 'KS', AND 'KS_UNDER'
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,*) 'GEOMETRY PROJECTION PARAMETERS'
      WRITE(4,*) '-----------------------------------------------------'
C      SMOOTH_MAX_SCHEME = 'MOD_P-NORM'
      WRITE(4,33)'SMOOTH_MAX_SCHEME', SMOOTH_MAX_SCHEME
C     PARAMETER TO BE USED FOR THE SMOOTH_MAX FUNCTION
C      SMOOTH_MAX_PARAM = 8
      WRITE(4,54) 'SMOOTH_MAX_PARAM', SMOOTH_MAX_PARAM
C     PENALIZATION SCHEME
C      PENALIZATION_SCHEME = 'SIMP' 
      WRITE(4,66) 'PENALIZATION_SCHEME', PENALIZATION_SCHEME    
C     PARAMETER TO BE USED FOR THE PENALIZATION
C      PENALIZATION_PARAM = 3 
      WRITE(4,99)'PENALIZATION_PARAM', PENALIZATION_PARAM
C-----------------------------------------------------------------------
C     OPTIMIZATION PARAMETERS      
C     WHETHER TO SCALE THE DESIGN VARIABLES TO RANGE [0,1]
      WRITE(4,*) '-----------------------------------------------------'
      WRITE(4,*) 'OPTIMIZATION PARAMETER'
      WRITE(4,*) '-----------------------------------------------------'
      DV_SCALING = .TRUE.
      WRITE(4,100) 'DV_SCALING', DV_SCALING
C     OPT GLOBAL VARIABLE FROM INT_OPTIMIZATION
      CALL INT_OPTIMIZATION
C-----------------------------------------------------------------------
C     FINITE DIFFERENCE CHECK OF SENSITIVITIES
C     ! IF REQUESTED      
C      IF (MAKE_FD_CHECK) THEN
C            CALL RUN_FINITE_DIFFERENCE_CHECK
C      ENDIF      

C
    1 FORMAT(F10.3)
    2 FORMAT(I10) 
C
   10 FORMAT(A19,5X,F10.3)
   20 FORMAT(A19,5X,I10)
   33 FORMAT(A19,5X,A15)
   54 FORMAT(A19,5X,I10)
   66 FORMAT(A19,5X,A12)
   99 FORMAT(A19,5X,I10) 
  100 FORMAT(A19,5X,L10)
  110 FORMAT(A19,5X,E10.3)
  120 FORMAT(A19,5X,L10) 

   40 FORMAT(A19,5X,F10.3)   

   50 FORMAT(3F20.15)     
   60 FORMAT(5F20.15) 

  200 FORMAT(A15,I10)            
      
      RETURN
      END SUBROUTINE READ_DATA
C********************************************************************************
C******************************** END INPUTE_FILE********************************
C********************************************************************************      
C-----------------------------------------------------------------------------------      
C-----------------------------------------------------------------------           |
C*********************  GEOMETRY PROJECTION FILE************************           |
C-----------------------------------------------------------------------           | 
      SUBROUTINE COMPUTE_BAR_ELEMENT_DISTANCE(DIST,DDIST_DBAR_ENDS)
C------------------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES AN ARRAY DIST OF DIMENSION N BAR AND N ELEMENT
C     WITH THE DISTANCE FROM THE CENTROID OF EACH ELEMENT TO EACH BAR'S MEDIAL
C     AXIS.
C
C     DDIST_DBAR_ENDS IS A 3-DIMENSIONAL ARRAY OF DIMENSIONS N_BAR_DOFS
C     X N_BAR X N_ELEM THAT CONTAINS THE SENSITIVITIES OF THE SIGNED
C     DISTANCES IN DIST WITH RESPECT TO EACH OF THE N_BAR_DOFS COORDS OF
C     THE MEDIAL AXIS END POINTS     
C------------------------------------------------------------------------------
C     USE MODULE
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM 
C     DATA PARAMETER DECLARATION SECTION     
      REAL, INTENT(OUT):: DDIST_DBAR_ENDS(N_BAR_DOFS,N_BAR), DIST(N_BAR)          
C     LOCAL VARIABLE DECLARATION SECTION
!      INTEGER N_BAR_DOFS
      REAL POINTS(N_POINT,DIM)
      REAL X_1B(DIM,N_BAR), X_2B(DIM,N_BAR), A_B(DIM,N_BAR),
     +  L_BE(N_BAR), R_BE(N_BAR) 
      REAL X_E_1B(DIM,N_BAR), X_E_2B(DIM,N_BAR), VEC_R_BE(DIM,N_BAR)
      REAL X_E(DIM), L_B(N_BAR) 
      REAL NORM_X_E_1B(N_BAR), NORM_X_E_2B(N_BAR)
C     CONVERTING NUMERICAL-LOGICAL VARIABLE
      INTEGER BRANCH1(N_BAR), BRANCH2(N_BAR),
     +  BRANCH3(N_BAR) 
      ! THE SENSITITY OF SIGNED DISTANCE
      REAL DD_BE_DX_1B, DD_BE_DX_2B, D_INV, L_BE_OVER_L_B
      DIMENSION DD_BE_DX_1B(DIM,N_BAR),DD_BE_DX_2B(DIM,N_BAR),
     +  D_INV(N_BAR), L_BE_OVER_L_B(N_BAR)
C-----------------------------------------------------------------------
C     EXECUTION SECTION
C-----------------------------------------------------------------------      
C-----------------------------------------------------------------------
C     CREATE INPUT GEOMETRY OF BAR      
C-----------------------------------------------------------------------
      ! SET PARAMETERS
      TOL = 1E-12          !TOLERANCE ON THE LENGTH OF A BAR
      ! COMPUTE VARIABLE LIMITS FOR EQ. (32)
      DIM = DIM 
      N_BAR  = N_BAR
C
C     THE FOLLOWING CODE IS VECTORIZED OVER THE BARS AND THE ELEMENTS.
C     WE HAVE TO BE CONSISTENT WITH THE ORDER OF INDICES TO PERFORM 
C     ELEMENT-WISE ARRAY OPERATIONS. THE ORDER OF INDICES IS (DIM, BAR,
C     ELEMENT)
C      
      POINTS(:,:) = POINT_MATRIX(:,2:(DIM+1))
C      
      DO N = 1, N_BAR
        X_1B(1,N) = POINTS((2*N-1),1)
        X_1B(2,N) = POINTS((2*N-1),2)      ! (I,B)
        X_2B(1,N) = POINTS((2*N),1)
        X_2B(2,N) = POINTS((2*N),2)        ! (I,B)      
      END DO
C
      A_B = X_2B - X_1B  ! NUMERATOR OF EQ. (11)
      L_B = SQRT(SUM(A_B**2, DIM = 1))  ! LENGTH OF THE BARS, EQ. (10)
      ! TO AVOID DIVISION BY ZERO
      DO N = 1,N_BAR
        IF (L_B(N).LE.TOL) THEN 
          L_B = 1.0         
        END IF
      END DO
      A_B(1,:) = A_B(1,:)/L_B(:)
      A_B(2,:) = A_B(2,:)/L_B(:)        ! NORMALIZE THE BAR DIRECTION TO UNIT VECTOR, EQ. (11)
C 
C                        
      DO J = 1, N_BAR
        X_E_1B(1,J) = CENTROIDS(1) - X_1B(1,J)
        X_E_1B(2,J) = CENTROIDS(2) - X_1B(2,J)
        X_E_2B(1,J) = CENTROIDS(1) - X_2B(1,J)
        X_E_2B(2,J) = CENTROIDS(2) - X_2B(2,J)
      END DO
C      
      NORM_X_E_1B(:) = SQRT(SUM(X_E_1B(:,:)**2, DIM = 1))     ! (B,E)
      NORM_X_E_2B(:) = SQRT(SUM(X_E_2B(:,:)**2, DIM = 1))     ! (B,E)     
      L_BE(:) = SUM(X_E_1B(:,:)*A_B(:,:), DIM = 1)            ! (B,E), EQ. (12)
      VEC_R_BE(1,:) = X_E_1B(1,:) - A_B(1,:)*L_BE(:)      
      VEC_R_BE(2,:) = X_E_1B(2,:) - A_B(2,:)*L_BE(:)        ! (I,B,E) 
      R_BE (:) = SQRT(SUM(VEC_R_BE(:,:)**2, DIM = 1))         ! (B,E), EQ. (13)
C     CONVERTING NUMERICAL-LOGICAL VARIABLE
      BRANCH1(:) = (L_BE(:) .LE. 0.0)           ! (B,E)
      BRANCH2(:) = (L_BE(:) .GT. L_B(:))
      BRANCH3(:) = ((L_BE(:).GT.0.0).AND.(L_BE(:).LE.L_B(:)))
C
      BRANCH1(:)  = -BRANCH1(:)
      BRANCH2(:)  = -BRANCH2(:)
      BRANCH3(:)  = -BRANCH3(:)
C
      DIST = BRANCH1*NORM_X_E_1B + BRANCH2*NORM_X_E_2B + BRANCH3*R_BE             
C
C------------------------------------------------------------------------------
C     COMPUTE SENSITIVITIES
C------------------------------------------------------------------------------
C      
      DD_BE_DX_1B = 0
      DD_BE_DX_2B = 0
C
      D_INV = DIST**(-1) ! THIS CAN RENDER A DIVISION BY ZERO(IF POINT
      ! LIES ON MEDIAL AXIS), AND SO WE NOW FIX IT
      L_BE_OVER_L_B = L_BE/L_B
C     THE SENSITIVITIES BELOW ARE OBATAINED FROM EQ. (30)
!     SENSITIVITIES TO X_1B AND X_2B
C
      DO N = 1, N_BAR              
        IF (L_BE(N).LE.0.0) THEN
          !X_1B    
          DD_BE_DX_1B(:,N) = -X_E_1B(:,N)*D_INV(N)
          !X_2B
          DD_BE_DX_2B(:,N) = 0
        ELSE IF (L_BE(N).GT.L_B(N)) THEN
          !X_1B
          DD_BE_DX_1B(:,N) = 0  
          !X_2B
          DD_BE_DX_2B(:,N) = -X_E_2B(:,N)*D_INV(N)
        ELSE 
          !X_1B
          DD_BE_DX_1B(:,N) = -VEC_R_BE(:,N)*D_INV(N)*
     +                        (1.0 - L_BE_OVER_L_B(N)) 
          !X_2B
          DD_BE_DX_2B(:,N) = -VEC_R_BE(:,N)*D_INV(N)*L_BE_OVER_L_B(N) 
        END IF
      END DO       
!     ASSEMBLE THE SENSITIVITIES TO THE BAR DESIGN PARAMETERS (SCALED)
      DO I = 1, DIM
        DDIST_DBAR_ENDS(I,:) = DD_BE_DX_1B(I,:)*POINT_SCALE(I) 
        DDIST_DBAR_ENDS(I+DIM,:) = DD_BE_DX_2B(I,:)*POINT_SCALE(I)
      END DO

      ! CHECK VALUE OF VARIABLE
C      WRITE(4,1001) ID, CENTROIDS(ID,1), CENTROIDS(ID,2)
C      WRITE(4,1060) (DDIST_DBAR_ENDS(1,N), N = 1, N_BAR)
C      WRITE(4,1060) (DDIST_DBAR_ENDS(2,N), N = 1, N_BAR)
C      WRITE(4,*) '----------------------------------------------------'
C      WRITE(4,1060) (DDIST_DBAR_ENDS(3,N), N = 1, N_BAR)
C      WRITE(4,1060) (DDIST_DBAR_ENDS(4,N), N = 1, N_BAR)
C      WRITE(4,*) '----------------------------------------------------'

C      DO I = 1, N_ELEM
!        WRITE(4,1001) ID, CENTROIDS(ID,1), CENTROIDS(ID,2)             
!      WRITE(4,'(5F15.5)') (DIST(N,ID), N = 1, N_BAR)    
C      END DO
 1001 FORMAT(I10,2F15.5)
 1060 FORMAT(5F15.5)
C
      RETURN
      END SUBROUTINE COMPUTE_BAR_ELEMENT_DISTANCE
      
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
      SUBROUTINE PROJECT_ELEMENT_DENSITY
C
C     THIS SUBROUTINE COMPUTES THE COMBINED UNPENALIZED DENSITIES (USED TO 
C     COMPUTE THE VOLUME) AND PENALIZED DENSITIES (USED TO COMPUTE THE ERSATZ
C     MATERIAL FOR THE ANALYSIS)
C     
C     IT ALSO COMPUTES THE DERIVATIVES OF THE UNPENALIZED AND PENALIZED
C     DENSITIES WITH RESPECT TO THE DESIGN PARAMETERS, AND SAVES THEM IN THE
C     COMMON VARIABE DELE_DENS_DDV AND DPENALIZED_ELE_DENS_DDV
C------------------------------------------------------------------------------ 
C
C     USING MODULE
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM
C-----------------------------------------------------------------------
C     LOCAL VARIABLE DECLARATION SECTION
C-----------------------------------------------------------------------      
      REAL D_BE(N_BAR), X_BE(N_BAR)
      REAL R_E, R_B(N_BAR), ALPHA_B(N_BAR)
      REAL DD_BE_DBAR_ENDS(N_BAR_DOFS,N_BAR)
C     CONVERTING NUMERICAL-LOGICAL VARIABLE
      LOGICAL BRANCH1(N_BAR), BRANCH2(N_BAR)
C     PROJECT DENSITY VARIABLE     
      REAL RHO_BE(N_BAR)
C     WITHOUT PENALIZED VARIABLE
      REAL RHO_E, DRHO_E_DX_BE(N_BAR)
C     WITH PENALIZED VARIALBE
      REAL :: PENAL_X_BE(N_BAR), DPENAL_X_BE_DX_BE(N_BAR),
     + DPENAL_RHO_E_DPENAL_X_BE(N_BAR), PENAL_RHO_E
      ! SENSITIVITY OF PROJECTED DENSITY VARIABLE
      REAL DRHO_BE_DBAR_RADII(N_BAR),DRHO_BE_DBAR_ENDS(N_BAR_DOFS,N_BAR)
      REAL DRHO_BE_DX_BE(N_BAR)
      REAL DX_BE_DBAR_RADII(N_BAR), DX_BE_DBAR_ENDS(N_BAR_DOFS,N_BAR),
     +     DX_BE_DBAR_SIZE(N_BAR)
      REAL DRHO_E_DBAR_SIZE(N_BAR), DRHO_E_DBAR_RADII(N_BAR),
     +     DRHO_E_DBAR_ENDS(N_BAR_DOFS, N_BAR),
     +     DPENAL_X_BE_DBAR_SIZE(N_BAR), DPENAL_X_BE_DBAR_RADII(N_BAR),  
     +     DPENAL_X_BE_DBAR_ENDS(N_BAR_DOFS, N_BAR),
     +     DPENAL_RHO_E_DBAR_SIZE(N_BAR),DPENAL_RHO_E_DBAR_RADII(N_BAR),
     +     DPENAL_RHO_E_DBAR_ENDS(N_BAR_DOFS, N_BAR)
      REAL DRHO_E_DDV(N_DV), 
     +     DPENAL_RHO_E_DDV(N_DV)
      INTEGER DDV     
C-----------------------------------------------------------------------
C     EXECUTION SECTION
C----------------------------------------------------------------------- 
      PI = 4.D0*ATAN(1.D0)
C DISTANCES FROM THE ELEMENT CENTROIDS TO THE MEDIAL SEGMENT OF EACH BAR
      CALL COMPUTE_BAR_ELEMENT_DISTANCE(D_BE,DD_BE_DBAR_ENDS)             
      R_B(:) = BAR_MATRIX(:,5)     ! BAR RADII 
      R_E = ELEM_R
      ! X_BE IS \PHI_B/R IN EQ.(2). NOTE THAT THE NUMERATOR CORRESPONDS TO
      ! THE SIGNED DISTANCE OF EQ. (8).
      X_BE = (R_B - D_BE)/R_E 
C      
      ! PROJECTED DENSITY
C     CONVERTING NUMERICAL-LOGICAL VARIABLE
        BRANCH1 = (ABS(X_BE).LT. 1.0)           ! (B,E)
        BRANCH2 = (X_BE.GE.1.0)
C          
      IF (DIM .EQ. 2) THEN       ! 2D
        DO N = 1, N_BAR
          IF (BRANCH1(N)) THEN
            RHO_BE(N) = ((PI - ACOS(X_BE(N)) + X_BE(N)
     +         *SQRT(1.0 - X_BE(N)**2))/PI)
            DRHO_BE_DX_BE(N) = ((SQRT(1.0-X_BE(N)**2.0)*2.0)/PI)
          ELSE IF (BRANCH2(N)) THEN
            RHO_BE(N) = 1.0
            DRHO_BE_DX_BE(N) = 0
          ELSE
            RHO_BE(N) = 0.0
            DRHO_BE_DX_BE(N) = 0
c         ((PI - ACOS(X_BE) + X_BE*SQRT(1.0 - X_BE**2))/PI)
c     +            *BRANCH1 + 1.0*BRANCH2        ! EQS (2) AND (3)
          END IF
        END DO
      END IF
      ! SENSITIVITIES OF RAW PROJECTED DENSITIES, EQS. (27) AND (29)
      DRHO_BE_DBAR_RADII = DRHO_BE_DX_BE*(1/ELEM_R)*RADIUS_SCALE
      DO I = 1, N_BAR_DOFS
        DRHO_BE_DBAR_ENDS(I,:) = DRHO_BE_DX_BE(:)*(-1/ELEM_R)*
     +                           DD_BE_DBAR_ENDS(I,:)                                    
      END DO

      ! GET SIZE VARIABLES
      ALPHA_B(:) = BAR_MATRIX(:,4)               ! BAR SIZE
C 
C     WITHOUT PENALIZATION:
C=========================
      ! X_BE HERE IS \HAT{\RHO}_B IN EQ.(4) WITH THE VALUE OF Q SUCH THAT
      ! THERE IS NO PENALIZATION (E.Q., Q = 1 IN SIMP)
        X_BE = RHO_BE*ALPHA_B  
      ! SENSITIVITIES OF UNPENALIZED EFFECTIVE DENSITIES, EQ. (26) WITH  
      ! \PARTIAL \MU / \PARTIAL (\ALPHA_B \RHO_{BE}) = 1
      DX_BE_DBAR_SIZE = RHO_BE
      DX_BE_DBAR_RADII = DRHO_BE_DBAR_RADII * ALPHA_B
      DO I = 1, N_BAR_DOFS      
        DX_BE_DBAR_ENDS(I,:) = DRHO_BE_DBAR_ENDS(I,:)*ALPHA_B(:)
      END DO

      ! COMBINED DENSITY OF EQ.(5).     
      CALL SMOOTH_MAX(X_BE,SMOOTH_MAX_PARAM,SMOOTH_MAX_SCHEME,
     +       RHO_MIN, RHO_E, DRHO_E_DX_BE)
      ! SENSITIVITIES OF COMBINED DENSITIES, EQ. (25)
      DRHO_E_DBAR_SIZE = DRHO_E_DX_BE*DX_BE_DBAR_SIZE
      DRHO_E_DBAR_RADII = DRHO_E_DX_BE*DX_BE_DBAR_RADII
      DO I = 1, N_BAR_DOFS
        DRHO_E_DBAR_ENDS(I,:) = DRHO_E_DX_BE(:)*DX_BE_DBAR_ENDS(I,:)
      END DO       
C
      ! STACK TOGETHER SENSITIVITIES WITH RESPECT TO DIFFERENT DESIGN
      ! VARIABLES INTO A SINGLE VECTOR PER ELEMENT
C
      DO M = 1, N_BAR
      DO N = 1, N_BAR_DOFS
            DRHO_E_DDV(N+(M-1)*N_BAR_DOFS) = DRHO_E_DBAR_ENDS(N,M)
      END DO
            DRHO_E_DDV(4*N_BAR+M) = DRHO_E_DBAR_SIZE(M)
            DRHO_E_DDV(5*N_BAR+M) = DRHO_E_DBAR_RADII(M)
      END DO 
C     WITH PENALIZATION:
C=======================================================================
      ! IN THIS CASE X_BE *IS* PENALIZED (EQ. (4)).
      CALL PENALIZED(X_BE,PENALIZATION_PARAM,PENALIZATION_SCHEME,
     +       PENAL_X_BE, DPENAL_X_BE_DX_BE)
      ! SENSITIVITIES OF EFFECTIVE (PENALIZED) DENSITIES, EQ. (26)
      DPENAL_X_BE_DBAR_SIZE = DPENAL_X_BE_DX_BE*DX_BE_DBAR_SIZE
      DPENAL_X_BE_DBAR_RADII = DPENAL_X_BE_DX_BE*DX_BE_DBAR_RADII
      DO I = 1, N_BAR_DOFS
        DPENAL_X_BE_DBAR_ENDS(I,:) = DPENAL_X_BE_DX_BE(:)
     +                               *DX_BE_DBAR_ENDS(I,:)
      END DO                  
C
      ! COMBINED DENSITY OF EQ. (5).      
      CALL SMOOTH_MAX(PENAL_X_BE, SMOOTH_MAX_PARAM,
     +     SMOOTH_MAX_SCHEME, RHO_MIN, PENAL_RHO_E,
     +     DPENAL_RHO_E_DPENAL_X_BE)
      ! SENSITIVITIES OF COMBINED DENSITIES, EQ. (25)
      DPENAL_RHO_E_DBAR_SIZE = DPENAL_RHO_E_DPENAL_X_BE
     +                        *DPENAL_X_BE_DBAR_SIZE
      DPENAL_RHO_E_DBAR_RADII = DPENAL_RHO_E_DPENAL_X_BE
     +                         *DPENAL_X_BE_DBAR_RADII
      DO I = 1, N_BAR_DOFS
        DPENAL_RHO_E_DBAR_ENDS(I,:) = DPENAL_RHO_E_DPENAL_X_BE(:)
     +                             *DPENAL_X_BE_DBAR_ENDS(I,:)
      END DO          
C
      ! STACK TOGETHER SENSITIVITIES WITH RESPECT TO DIFFERENT DESIGN
      ! VARIABLES INTO A SINGLE VECTOR PER ELEMENT
C 
      DO M = 1, N_BAR
      DO N = 1, N_BAR_DOFS
            DPENAL_RHO_E_DDV(N+(M-1)*N_BAR_DOFS) 
     1       = DPENAL_RHO_E_DBAR_ENDS(N,M)
      END DO
            DPENAL_RHO_E_DDV(4*N_BAR+M) = DPENAL_RHO_E_DBAR_SIZE(M)
            DPENAL_RHO_E_DDV(5*N_BAR+M) = DPENAL_RHO_E_DBAR_RADII(M)
      END DO
C     STORING THE ELEMENT DENSITIES WITH & WITHOU PENALTY
      ELEM_DENS =  RHO_E
      PENALIZED_ELEM_DENS = PENAL_RHO_E      
      ! STORING THE SENSITIVITIES OF DENSITIES TO GLOBAL VARIABLE OPT
      DELEM_DENS_DDV = DRHO_E_DDV
      DPENALIZED_ELEM_DENS_DDV = DPENAL_RHO_E_DDV        
C      WRITE(4,'(I10,2F15.5)') ID, CENTROIDS(1), CENTROIDS(2)
C      DO I = 1,6
C        WRITE(4,1060) (DPENAL_RHO_E_DDV((I-1)*N_BAR + N), N = 1, N_BAR)
C      END DO   
 1001 FORMAT(I10,2F15.5)
 1060 FORMAT(8F15.5)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SMOOTH_MAX(X, P, FORM_DEF, X_MIN, S, DSDX)
C     THIS SUBROUTINE COMPUTES A SMOOTH APPROXIMATION OF THE MAXIMUM OF X. THE
C     TYPE OF SMOOTH APPROXIMATION IS GIVEN BY THE ARGUMENT FORM DEF, AND
C     THE CORRESPONDING APPROXIMATION PARAMETER IS GIVEN BY P. X_MIN IS A LOWER
C     BOUND TO THE SMOOTH APPROXIMATION FOR THE MODIFIED P-NORM AND MODIFIED 
C     P-MEAN APPROXIMATION.
C 
C     THE OPTIONAL THIRD ARGUMENT IS A STRING THAT INDICATES THE WAY THE
C     APPROXIMATION IS DEFINED, POSSIBLE VALUES ARE:
C     'MOD_P-NORM'     : OVERESTIMATE USING MODIFIED P-NORM (SUPPORTS X=0) 
C     'MOD_P-MEAN'     : UNDERESTIMATE USING MODIFIED P-NORM (SUPPORTS X=0)
C     'KS'             : KREISSELMEIER-STEINHAUSER, OVERESTIMATE
C     'KS_UNDER'       : KREISSELMEIER-STEINHAUSER, UNDERRESTIMATE        
C------------------------------------------------------------------------------ 
C     USE MODULE
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM
C     INPUT ARGUMENT VARIABLES
      REAL, INTENT(IN) :: X(N_BAR)
      INTEGER, INTENT(IN) :: P
      CHARACTER, INTENT(IN) :: FORM_DEF*18
      REAL, INTENT(IN) :: X_MIN                          
C     OUTPUT ARGUMENT VARIABLES      
      REAL, INTENT(OUT) :: DSDX(N_BAR)
      REAL, INTENT(OUT) :: S
C
C     DECLARATION SECTION
C
      REAL SUM_PENAL_RHO, SUM_EPX_RHO, XP(N_BAR) 
C            
C     EXECUTABLE SECTION
      ! COMPUTING SUMMATION OF DENSITY
      SUM_PENAL_RHO = SUM(X**P, DIM = 1)
      ! COMPUTING SUMMATION OF EXP DENSITY        
      SUM_EPX_RHO = SUM(EXP(P*X), DIM = 1)
C
      SELECT CASE(FORM_DEF)
      CASE ("MOD_P-NORM") 
      ! EQ. 6
      ! IN THIS CASE, WE ASSUME X >= 0
      S = (X_MIN**P + (1.0 - X_MIN**P)*SUM_PENAL_RHO)**(1.0/P)
      DSDX = (1 - X_MIN**P)*(X/S)**(P-1)
      CASE ("MOD_P-MEAN")
      ! IN THIS CASE, WE ASSUME X >= 0      
      S = (X_MIN**P + (1.0 - X_MIN**P)*SUM_PENAL_RHO/N_BAR)**(1.0/P)
      DSDX = (1 - X_MIN**P)*(1/N)*(X/S)**(P-1)
      CASE ("KS")
      S = X_MIN + (1.0 - X_MIN)*LOG(SUM_EPX_RHO)/P
      DSDX = (1 - X_MIN)*EXP(P*X)/SUM_EPX_RHO
      CASE ("KS_UNDER")
      ! NOTE: CONVERGENCE MIGHT BE FIXED WITH EULER-GAMMA
      S = X_MIN + (1.0 - X_MIN)*LOG(SUM_EPX_RHO/N_BAR)/P 
      DSDX = (1-X_MIN)*EXP(P*X)/SUM_EPX_RHO                
      CASE DEFAULT
      WRITE(4,1120) FORM_DEF
      END SELECT      
 1120 FORMAT(//37H SMOOTH_MAX RECEIVED INVALID FORM_DEF , A18)
C
      RETURN
      END SUBROUTINE SMOOTH_MAX
C------------------------------------------------------------------------------
      SUBROUTINE PENALIZED(X, Q, FORM_DEF, P, DPDX) 
C
C     THIS SUBROUTINE COMPUTES THE COMBINED PENALIZED DENSITIES 
C     [P, DPDX] = PENALIZE(X, Q, PENAL_SCHEME)
C     PENALIZE(X) ASSUME X IN [0,1] AND DECREASES THE INTERMEDIATE VALUES
C
C     FOR A SINGLE INPUT, THE INTERPOLATION IS SIMP WITH P =3
C
C     THE OPTIONAL SECOND ARGUMENT IS THE PARAMETER VALUE P.
C
C     THE OPTIONAL THIRD ARGUMENT IS A STRING THAT INDICATES THE WAY THE
C     INTERPOLATION IS DEFINED, POSSIBLE VALUES ARE:
C      'SIMP'     : DEFAULT
C      'RAMP'     :
C
C------------------------------------------------------------------------------    
C     USE MODULE
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM
C     DECLARATION SECTION
      CHARACTER, INTENT(IN):: FORM_DEF*6
      INTEGER, INTENT(IN):: Q
      REAL, INTENT(IN):: X(N_BAR)
      REAL, INTENT(OUT):: P(N_BAR), DPDX(N_BAR)
C     EXECUTION SECTION
C     CONSIDER INPUT
      SELECT CASE(FORM_DEF)
      CASE ("SIMP")
      P = X**Q
      DPDX = Q*X**(Q-1)       
      CASE ("RAMP")      
      P = X/(1.0 + Q*(1.0 - X))
      DPDX = (1+Q)/(1+Q*(1-X))**2
      CASE DEFAULT
      WRITE(4,1120), FORM_DEF
      END SELECT
 1120 FORMAT(//36H PENALIZE RECEIVED INVALID FORM_DEF , A18)
C      DO I = 1, NO_EL
C        WRITE(4,'(I10,2F15.5)') I, CENTROIDS(I,1), CENTROIDS(I,2)    
C        WRITE(4,'(F15.5,E15.5)') S(I), SUM_PENAL_RHO(I)
C        WRITE(4,'(4F15.5)') (P(N,I), N = 1, N_BAR)             
C      END DO 
      RETURN
      END SUBROUTINE PENALIZED
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C                  
C********************************************************************************
C************************* END GEOMETRY_PROJECTION FILE**************************
C********************************************************************************        
C-----------------------------------------------------------------------------

C-----------------------------------------------------------------------           |
C****************************  OPTIMIZATION FILE************************           |
C-----------------------------------------------------------------------           |       
C-----------------------------------------------------------------------------
      SUBROUTINE INT_OPTIMIZATION
C     INITIALIZE FUNCTIONS TO COMPUTE
C     CONCATENATE LIST OF FUNCTIONS TO BE COMPUTED
C-----------------------------------------------------------------------      
C     USING MODULE
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM
C     LOCAL VARIABLE DECLARATION SECTION     
      REAL DELTA_RADIUS
      REAL X_1B_ID(N_BAR), X_2B_ID(N_BAR)
C     ALLOCATE MEMORY FOR VARIABLES
C     ALLOCATE MEMORY FOR NODE COORDINATES AND CONNECTION TABLE       
      IF(.NOT. ALLOCATED(POINT_SCALE)) ALLOCATE(POINT_SCALE(DIM))
      IF(.NOT. ALLOCATED(CURRENT_POINT_MATRIX))
     +                   ALLOCATE(CURRENT_POINT_MATRIX(N_POINT,(3))) 
      IF(.NOT. ALLOCATED(CURRENT_POINT_MATRIX))
     +                        ALLOCATE(CURRENT_BAR_MATRIX(N_BAR,5))
      IF(.NOT. ALLOCATED(POINT_DV)) ALLOCATE(POINT_DV(DIM*N_POINT))
      IF(.NOT. ALLOCATED(SIZE_DV)) ALLOCATE(SIZE_DV(N_BAR))
      IF(.NOT. ALLOCATED(RADIUS_DV)) ALLOCATE(RADIUS_DV(N_BAR)) 
C     EXECUTION SECTION
C     ! INITILIZE THE DESIGN VARIABLE AND ITS INDEXING SCHEMES
C     WE ARE DESIGNING THE POINTS, THE SIZE VARIABLES, AND THE RADII OF 
C     THE BARS:
      N_DV = DIM*N_POINT + 2*N_BAR
      ! POINT_DV
      DO I = 1, (DIM*N_POINT)
        POINT_DV(I) = I 
      END DO
      ! SIZE_DV
      DO I = 1, N_BAR
        SIZE_DV(I) =  POINT_DV(DIM*N_POINT)+I
      END DO  
      ! RADIUS_DV
      DO I = 1, N_BAR
        RADIUS_DV(I) = SIZE_DV(N_BAR) + I
      END DO 

      IF (DV_SCALING) THEN
      ! COMPUTE VARIABLE LIMITS FOR EQ. (32)
      POINT_SCALE(1) = 20.0
      POINT_SCALE(2) = 10.0 
      POINT_MIN = 0
      ! CONSIDER POSSIBILITY THAT MAX_BAR_RADIUS AND MIN_BAR_RADIUS ARE 
      ! THE SAME (WHEN BARS ARE OF FIXED RADIUS)
      DELTA_RADIUS = MAX_BAR_RADIUS - MIN_BAR_RADIUS
            IF (DELTA_RADIUS.LT.1E-12) THEN
                  RADIUS_SCALE = 1.0
            ELSE 
                  RADIUS_SCALE = DELTA_RADIUS
            ENDIF
      RADIUS_MIN = MIN_BAR_RADIUS 
      ELSE
      POINT_SCALE = 1.0
      POINT_MIN = 0.0
      RADIUS_SCALE = 1.0
      RADIUS_MIN = 0.0
      ENDIF
C     FILL IN DESIGN VARIABLE VECTOR BASED ON THE INITIAL DESIGN
C     SET THE CURRENT DESIGN TO THE INITIAL DESIGN:
      CURRENT_POINT_MATRIX = POINT_MATRIX
      CURRENT_BAR_MATRIX = BAR_MATRIX
C     CONSIDER THE BAR DESIGN VARIABLES
C     X_1B, X_2B, ALPHA_B, R_B
      X_1B_ID = CURRENT_BAR_MATRIX(:,2)
      X_2B_ID = CURRENT_BAR_MATRIX(:,3) 

      RETURN
      END SUBROUTINE INT_OPTIMIZATION
C********************************************************************************
C***************************** END UNTILITY FILE*********************************
C******************************************************************************** 

      subroutine get_elem_con_2D
C    
      include 'aba_param.inc' 
c       
      character*256 jobname,outdir,filename,input
      integer ele(20000,5)
      real node(40000,3)
      common /elem_inform/ ele, node
C     LOCAL VARIABLE
      integer non
      INTEGER NELCNT, NBNODE
c           
c           Open INP file for reading.
c
      call getjobname(jobname,lenjobname)
      call getoutdir(outdir,lenoutdir)
      filename=outdir(1:lenoutdir)//'/'//
     +        jobname(1:lenjobname)//'.inp' 
      open(unit=17,file=filename(1:lenoutdir+                       
     +        lenjobname+5),status='unknown')
c           
      ncnt=0
      nelcnt=0
      nbnode=0
c           
c           Skip down to *Node in INP File
c      
      do while (index(input,'*Node')==0)
            read(17,'(a)')input
      end do
c
c           Read number of node and node coordinate
c
      do while (.true.)
            read(17,*) input
            if(index(input,'*')==0)then             
                  backspace(17)
                  nbnode = nbnode +1
                  read(17,*) non, node(nbnode,1), node(nbnode,2)
c             write(4,*) input
C                  write(4,*) non, node(nbnode,1), node(nbnode,2)                          
c              read(17,*) input, node(nbnode,1), node(nbnode,2)
c              write(4,*) nbnode, node(nbnode,1), node(nbnode,2)       
            else
                  exit
            endif                  
      end do

c                                         
c           Read in Element Nodal Connectivity  
c     
      do while(.true.)
            read(17,*)input                                       
            if(index(input,'*')==0)then
                  backspace(17)
                  nelcnt=nelcnt+1   
                  read(17,*)ele(nelcnt,1),ele(nelcnt,2),
     +                   ele(nelcnt,3),ele(nelcnt,4),ele(nelcnt,5)
            else
                  exit
            endif
C            WRITE(4,*)ele(1,1),ele(1,2),
C     +                   ele(1,3),ele(1,4),ele(1,5)
      end do
c
      return
      end subroutine get_elem_con_2D
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      SUBROUTINE COMPUTE_ELEM_PARAMETER(NOEL)
      use user_global_variable_FE
      use user_global_variable_OPT
      use user_global_variable_GEOM      
      INTEGER, INTENT(IN) :: NOEL
C     THIS SUBROUTINE COMPUTES ELEMENT VOLUMES, ELEMENT RADIUS
C     ELEMENT CENTROID LOCATIONS AND MAX/MIN NODAL COORDINATES VALUE FOR
C     MESH. IT ASSUMES THE FE STRUCTURE HAS ALREADY BEEN POPULATED.     
      ! SUBROUTINE VARIALBE
      integer ele(20000,5)
      real node(40000,3)
      common /elem_inform/ ele, node
      ! SUBROUTINE VARIABLE
      REAL N1(DIM), N2(DIM), N3(DIM), N4(DIM), DIAG1(DIM), DIAG2(DIM),
     +  V1(DIM)
C
      INTERFACE
         FUNCTION CROSS(X,Y)
         REAL, DIMENSION(3) :: CROSS, X, Y
         END FUNCTION CROSS
      END INTERFACE
C
C     ALLOCATE MEMORY FOR NODE COORDINATES AND CONNECTION TABLE
      IF(.NOT. ALLOCATED(CENTROIDS)) ALLOCATE(CENTROIDS(DIM))      
C
C     EXECUTION AREA
C
      DIAG1 = 0.0
      DIAG2 = 0.0
      V1 = 0.0   
      ELEM_VOL = 0.0
      ! COORDINATE OF NODE OF ELEMENT
      N1(1) = NODE(ELE(NOEL,2),1)
      N1(2) = NODE(ELE(NOEL,2),2)
C
      N2(1) = NODE(ELE(NOEL,3),1)
      N2(2) = NODE(ELE(NOEL,3),2)
C
      N3(1) = NODE(ELE(NOEL,4),1)
      N3(2) = NODE(ELE(NOEL,4),2)
C
      N4(1) = NODE(ELE(NOEL,5),1)
      N4(2) = NODE(ELE(NOEL,5),2) 
      ! CENTROID COORDINATE
      CENTROIDS = (N1 + N2 + N3 + N4)/4.0
C
      ! HERE WE CAN TAKE ADVANTAGE OF THE PLANAR QUADRILATERALS AND USE
      ! BRETSCHNEIDER'S FORMULA (CROSS PRODUCT OF DIAGONALS):
      DIAG1 = N3 - N1
      DIAG2 = N4 - N2
      !NOTE: SINCE ONLY THE 3RD COMPONENT IS NONZERO, WE USE IT INSTEAD
      ! OF, 0.5*SQRT(SUM(CROSS(DIAG1,DIAG2)**2))
      V1 = CROSS(DIAG1, DIAG2)
      ELEM_VOL = 0.5*(SQRT(V1(3)**2))
      ! COMPUTE SAMPLING RADIUS
      ! THE RADIUS CORRESPONDS TO THE CIRCLE (OR SPHERE) THAT
      ! CIRCUMSCRIBES A SQUARE (OR CUBE) THAT HAS THE EDGE LENGTH OF
      ! ELEM_SIZE.
      DIM_2D = REAL(DIM)
      ELEM_R = SQRT(DIM_2D)/2 * ELEM_VOL**(1/DIM_2D)

C      WRITE(4,*) NOEL,ele(NOEL,2), ele(NOEL,3), ele(NOEL,4), ele(NOEL,5)
C      WRITE(4,*) NOEL, CENTROID(1), CENTROID(2)
C      WRITE(4,*) NOEL, ELEM_VOL, RSAMP
      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------      
C     UTILITY FUNCTION
C-----------------------------------------------------------------------
C----------------------------------------------------------------------- 
      FUNCTION CROSS(X,Y)
        REAL, DIMENSION(3) :: CROSS, X, Y 
        CROSS(1) = X(2)*Y(3) - X(3)*Y(2)
        CROSS(2) =-X(1)*Y(3) + X(3)*Y(1)
        CROSS(3) = X(1)*Y(2) - X(2)*Y(1)
      END FUNCTION CROSS

