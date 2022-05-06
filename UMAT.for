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

        DIMENSION DSTRESS(NTENS), TEMP(NTENS,NTENS),DS(3,NTENS)
        CHARACTER*8 msg
        INTEGER*4 status,stat
          
c        print *,'JSETP     ', JSTEP
        print *,'KINC      ', KINC
        print *,'NOEL, NPT ', NOEL, NPT
        
        OPEN (unit = 15, file = 'c:\repo\UEL\DSTRAIN.txt', status='replace')
        write(15,'(I0)') JSTEP,KINC, NOEL, NPT
        write(15,'(4X,ES24.16)') DSTRAN
        CLOSE(15)
        
        CALL EXECUTE_COMMAND_LINE('python c:\repo\UEL\py.py',wait=.true.
     &, exitstat=status, CMDSTAT=stat)
        if (status==0) then 
c            print *, 'UMAT INFO: RVE script COMPLETED'
            OPEN (unit = 15, file = 'c:\repo\UEL\DDSDDE.OUT', status='OLD', ACTION='READ')
            OPEN (unit = 16, file = 'c:\repo\UEL\DSTRESS.OUT', status='OLD', ACTION='READ')
            DO  i=1,NTENS
                READ(15,*) DDSDDE(1:NTENS,i)
            END DO
            DO  i=1,NTENS
                READ(16,*) DS(1:3,i)
            END DO
            DSTRESS(1)=DS(1,1)
            DSTRESS(2)=DS(2,2)
            DSTRESS(3)=DS(3,3)
            DSTRESS(4)=DS(1,4)
            DSTRESS(5)=DS(1,5)
            DSTRESS(6)=DS(2,6)

            CLOSE(15)
            CLOSE(16)

            DO K1 = 1,NTENS
                STRESS(K1) = STRESS(K1) + DSTRESS(K1)
            END DO
        else
            print *, 'UMAT ERROR: RVE script ERROR exit'
        end if
c        CALL XIT

c        EMOD=2.0e8
c        ENU=0.3
c
c        EBULK3=EMOD/(1-2*ENU)
c        EG2=EMOD/(1+ENU)
c        EG=EG2/2 ! G
c        EG2=EMOD/(1+ENU)
c        EG3=3*EG
c        ELAM=(EBULK3-EG2)/3 ! λ
cc       initialize DDSDDE
c        DO K1=1,NTENS
c            DO K2=1,NTENS
c                DDSDDE(K1,K2)=0.0
c            END DO
c        END DO
cc
cc          ⎡ λ + 2G    λ       λ                      ⎢
cc          ⎢ λ         λ + 2G  λ                      ⎢
cc     J=   ⎢ λ         λ       λ + 2G                 ⎢
cc         ⎢                             G            ⎢
cc         ⎢                                   G      ⎢
cc         ⎣                                       G  ⎢
c        DO K1=1,NDI
c            DO K2=1,NDI
c                DDSDDE(K2,K1)=ELAM
c            END DO
c            DDSDDE(K1,K1)=EG2+ELAM
c        END DO
c        DO K1=NDI+1,NTENS
c            DDSDDE(K1,K1)=EG
c        END DO
cC
cC     DETERMINE STRESS INCREMENT
cC
c        DO K1 = 1,NTENS
c            TEMP(:,K1)=DDSDDE(:,K1)*DSTRAN(K1)
c        END DO
c        DSTRESS=SUM(TEMP,DIM=2)
cC
cC      UPDATE STRESS
cC
c        DO K1 = 1,NTENS
c            STRESS(K1) = STRESS(K1) + DSTRESS(K1)
c        END DO
c
c
       RETURN
       END