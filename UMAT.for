        SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
       INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

        DIMENSION DSTRESS(NTENS), TEMP(NTENS,NTENS)
        CHARACTER*8 msg
        INTEGER*4 status,stat
 
C     Open output file         
c        OPEN (unit = 15, file = 'c:\repo\UEL\DSTRAIN.txt', status='replace')
c        write(15,'(I0)') JSTEP,KINC
c        write(15,'(4X,E24.16)') DSTRAN
c        CLOSE(15)


        CALL EXECUTE_COMMAND_LINE('python c:\repo\UEL\py.py',wait=.true.
     &, exitstat=status, CMDSTAT=stat, CMDMSG=msg)
        print *, status
        print *, stat
        print *, msg
        CALL XIT
        
        EMOD=2.0e8
        ENU=0.3

        EBULK3=EMOD/(1-2*ENU)
        EG2=EMOD/(1+ENU)
        EG=EG2/2 ! G
        EG2=EMOD/(1+ENU)
        EG3=3*EG
        ELAM=(EBULK3-EG2)/3 ! λ
c       initialize DDSDDE
        DO K1=1,NTENS
            DO K2=1,NTENS
                DDSDDE(K1,K2)=0.0
            END DO
        END DO

c          ⎡ λ + 2G    λ       λ                      ⎢
c          ⎢ λ         λ + 2G  λ                      ⎢
c     J=   ⎢ λ         λ       λ + 2G                 ⎢
c         ⎢                             G            ⎢
c         ⎢                                   G      ⎢
c         ⎣                                       G  ⎢
        DO K1=1,NDI
            DO K2=1,NDI
                DDSDDE(K2,K1)=ELAM
            END DO
            DDSDDE(K1,K1)=EG2+ELAM
        END DO
        DO K1=NDI+1,NTENS
            DDSDDE(K1,K1)=EG
        END DO
        
C
C     DETERMINE STRESS INCREMENT
C
        DO K1 = 1,NTENS
            TEMP(:,K1)=DDSDDE(:,K1)*DSTRAN(K1)
        END DO
        DSTRESS=SUM(TEMP,DIM=2)
C
C      UPDATE STRESS
C
        DO K1 = 1,NTENS
            STRESS(K1) = STRESS(K1) + DSTRESS(K1)
        END DO


       RETURN
       END