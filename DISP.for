        SUBROUTINE DISP (U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
        INCLUDE 'ABA_PARAM.INC'
C
        DIMENSION U(3),TIME(2),COORDS(3)
        INTEGER i=1

        DO WHILE (i <=121)
                IF(NODE.EQ.i.AND.JDOF.EQ.1) THEN
                        U(1)=1e-6
                ENDIF
                i=i+1
        END DO

        RETURN
        END