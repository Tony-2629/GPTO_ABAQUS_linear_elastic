       E=PROPS(1)*PENALIZED_ELEM_DENS
       v=PROPS(2)
       myE=PROPS(1)!*PENALIZED_ELEM_DENS
       myNu=PROPS(2)       
C
      myMu= myE / (TWO * (ONE + myNu)) 
      myLambda= (myE * myNu) / ((ONE + myNu) * (ONE - TWO * myNu))
C
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J)= ZERO
        END DO
      END DO
C            
      IF(DIM.EQ.2) THEN
        DDSDDE(1,1) = E/(ONE-v*v)
        DDSDDE(2,2) = E/(ONE-v*v)
        DDSDDE(3,3) = E*(ONE-v)/(TWO*(ONE-v*v))
C         
        DDSDDE(1,2) = E*v/(ONE-v*v)
        DDSDDE(2,1) = E*v/(ONE-v*v)
C
      ELSEIF (DIM.EQ.3) THEN        
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
      ELSE
        WRITE(4,*) 'ERROR, DIMENSION SHOULD BE 2D OR 3D'
      ENDIF
C               
      DO i=1, NTENS
       DO j=1, NTENS
           STRESS(j)=STRESS(j) + DDSDDE(j, i) * DSTRAN(i)
       END DO
      END DO 
C
C     ! ELASTIC STRAIN ENGERGY
C     
      DEE = 0.0
      DO I = 1, NTENS
            DEE = DEE + 0.5*STRESS(I)*DSTRAN(I)
      END DO
C
      SSE = DEE
C
!      Q = 1

!      DC_PENALIZED_ELEM_DENS = -2.0*DEE*(ELEM_VOL)
!     1  *(Q*PENALIZED_ELEM_DENS**(Q-1)/PENALIZED_ELEM_DENS)

!      DO N = 1, N_DV
!        DC_DDV(N) = DC_PENALIZED_ELEM_DENS*DPENALIZED_ELEM_DENS_DDV(N)
!      END DO

C     STORE COMPLIANCE AND THE SENSITIVITIES OF COMPLIANCE
C
!      DO N = 1, N_DV
!        STATEV(N) = DC_DDV(N)
!      END DO
!      STATEV(N_DV+1) = 2.0*DEE*(ELEM_VOL)

C      STORE COMBINED PROJECTED DENSITY

!      STATEV(N_DV+2) = PENALIZED_ELEM_DENS

C     COMPUTE THE VOLUME FRACTION
!      V_E = ELEM_VOL      
!      V = (ELEM_VOL)*ELEM_DENS
!      STATEV(500+N_DV+1) = V
!      STATEV(500+N_DV+2) = V_E
C     COMPUTE THE DESIGN SENSITIVITY

!      DO N = 1, N_DV
!        DVOLFRAC_DDV(N) = V_E*DELEM_DENS_DDV(N)      
!      END DO
C     STORE THE DESIGN SENSITIVITY
!      DO N = 1, N_DV
!        STATEV(500+N) = DVOLFRAC_DDV(N)
!      END DO
!   50 FORMAT( I10)
C   60 FORMAT(10F10.4)