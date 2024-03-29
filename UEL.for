!DEC$ ATTRIBUTES ALIAS:"uel"::UEL
c
c     ** MATERIALS
c     **      cn=PROPS(1)
c     **      ct=PROPS(2)
c     **      Gn=PROPS(3)
c     **      Gt=PROPS(4)
c     **      Qn0=PROPS(5)
c     **      Qt0=PROPS(6)
c             tempK=PROPS(6) : adjust the cohesive constitution
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1 PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2 DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3 PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4 NJPROP, PERIOD)
      INCLUDE 'ABA_PARAM.INC'
c     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
cccccc The variations are defined in a seperate file as follow
      INCLUDE 'VAR_define.for'
cccccc TODO: to delete non-using VARs and move VARS to VAR FILE
      Real(Kind=8) :: Du_up(3,LAYERNODE),Du_low(3,LAYERNODE),!increasment of node displacement on surface
     1 Du_nodal(3,LAYERNODE*2),Du_seperation(3,LAYERNODE),!increasement of open displacement & mid-plane displacement
     2 U_up(3,LAYERNODE),U_low(3,LAYERNODE),!total node displacement
     3 U_nodal(3,LAYERNODE),U_mid(3,LAYERNODE),U_total_nodal(3,LAYERNODE*2)!total open & mid-plane displacement
     4 ,U_total_nodal_up(3,LAYERNODE),U_total_nodal_low(3,LAYERNODE)
     5 ,Seperation_total_nodal(3,LAYERNODE),Seperation_total_Gauss(3,LAYERNODE)
     9 ,Temp(2,LAYERNODE),Tensor_Product_T(3,3)
c
      Real(Kind=8) :: PD_I1_g1(3),PD_I2_g1(3),PD_I3_g1(3)
     5,PD_I1_g2(3),PD_I2_g2(3),PD_I3_g2(3)
     6,PD_sai0n_I1,PD_sai0n_I2,PD_sai0n_I3,PD_sai0t_I1,PD_sai0t_I2,PD_sai0t_I3
     7,PD_sai0n_u(3),PD_sai0t_u(3)
     8,PD_sai0n_g1(3),PD_sai0t_g1(3)
     9,PD_sai0n_g2(3),PD_sai0t_g2(3)
     3,PD_sain_u(3)
     4,R_g1_nodal(3,LAYERNODE),R_g2_nodal(3,LAYERNODE)
     4,R_g1(3,LAYERNODE),R_g2(3,LAYERNODE)
     4,R_u_low_nodal(3,LAYERNODE),R_u_up_nodal(3,LAYERNODE)
     4,R_u_low(3,LAYERNODE),R_u_up(3,LAYERNODE)
     4,R_up(3,LAYERNODE),R_low(3,LAYERNODE)
     4,J_up(3,3),J_low(3,3)
     4,detJ_up(1,LAYERNODE),detJ_low(1,LAYERNODE)

      CHARACTER*80 FNAME     
      DIMENSION  ARRAY(513),JRRAY(NPRECD,513),LRUNIT(2,1),S(6),E(6),DISP(6)     
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
      INTEGER*8 iVar1,iVar2
      CHARACTER*8 cVar1,cVar2
      Real*8 rVar1

c ================== VAR of calculating points =============
c all the data can be calculate at one 2D mid-surface
ccccccccc    NOT using ccccccccccccccc
cccc ALL_POINTS(2,9) structure:
ccccc               node 1~6, Gauss points 1~3
ccccc        r | elememt node : Gauss points |
cccccc       s | elememt node : Gauss points |
ccccccccccccccccccccccccccc
      interface
            function TensorProduct(g,gc)
                  real(Kind=8) :: g(3),gc(3),TensorProduct(3,3)
            end function
            function Dot_Matrix(A,B,NODE)
                  real(Kind=8) :: A(3,NODE),B(3,NODE),Dot_Matrix(1,NODE)
                  INTEGER :: NODE
            end function
            function Dot_3_1(A,B)
                  real(Kind=8) :: A(3,3),B(3),Dot_3_1(3)
            end function
            function Shape_poly(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly(3,6)
            end function
            function Shape_poly_PD_r(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_r(3,6)
            end function
            function Shape_poly_PD_s(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_s(3,6)
            end function
      end interface 
c ===========  K matrix output ================
c ==== every step the K matrix is being overwrited
c ==== expecting: every matrix should be stored in this file
      open(110, file = 'c:\repo\UEL\KK.dat', status = 'OLD')
c ================================================
      write (7,*) '======================================'
      write (7,*) '=====           UEL START        ====='
      write (7,*) '=====       Information Pass In  ====='
      write (7,*) '=====    ELEMENT    :', JELEM
      write (7,*) '=====    STEP       :', KSTEP 
      write (7,*) '=====    INCRESEMENT:', KINC
      write (7,*) '======================================'
C      write (7,*) '1. Original coordinates:'
C      write (7,*) COORDS
C      write (7,*) '2. Properties, total number is',NPROPS
C      write (7,*) PROPS(1:NPROPS)
c      write (7,*) '3. SVARS, total number is',NSVARS
c      write (7,*) SVARS(1:NSVARS)
C      write (7,*) '4. Total number of distributed loads and/or fluxes:',MDLOAD
C      write (7,*) '5. (JDLTYP)Distributed load types:'
C      write (7,*) JDLTYP(MDLOAD,1),MDLOAD
C      write (7,*) '6. (ADLMAG)Total load magnitude of the K1th distributed load at the end of the current increment' 
C      write (7,*) ADLMAG(MDLOAD,1)
C      write (7,*) '7. (DDLMAG)Increments in the magnitudes of the distributed loads that are currently active on this element'
C      write (7,*) DDLMAG(MDLOAD,1)
C      write (7,*) '8. (MLVARX)Dimensioning parameter of RHS:'
C      write (7,*)  MLVARX
c      write (7,*) '                 Information about calculating'
c      write (7,*) '1. (U)Total values of the variables'
c      write (7,*) U
c      write (7,*) '2. (DU)Incremental values of the variables for the current'
c      write (7,*) DU(:,1)
c ========== define the unit matrix ==========
      Matrix_1(:,1)=(/1.0,0.0,0.0/);Matrix_1(:,2)=(/0.0,1.0,0.0/);Matrix_1(:,3)=(/0.0,0.0,1.0/);
c ========== initialize parameters  ==========
      PD_I1_g1=0.0;PD_I1_g2=0.0
c=========== define element shape  =======================
c====== WHEN triangle LAYERNODE=6 ======================
c====== WHEN quadrilateral LAYERNODE=8 ======================
c      LAYERNODE=PROPS(7)
c      write (7,*) 'SVARS passin is'
c      write (7,*) SVARS(25),SVARS(26),SVARS(27),SVARS(28),SVARS(29),SVARS(30)
c========= define the local coordinates WHEN the element is triangle ==========
            Gauss_point(:,1)=(/0.16666667,0.16666667/)
            Gauss_point(:,2)=(/0.66666667,0.16666667/)
            Gauss_point(:,3)=(/0.16666667,0.66666667/)
            weight=0.33333333
c  local coordinates in calculating K, ranged in (r,s) sequence
            Nodal_Local_coord(:,1)=(/0.0,0.0/)
            Nodal_Local_coord(:,2)=(/1.0,0.0/)
            Nodal_Local_coord(:,3)=(/0.0,1.0/)
            Nodal_Local_coord(:,4)=(/0.5,0.0/)
            Nodal_Local_coord(:,5)=(/0.5,0.5/)
            Nodal_Local_coord(:,6)=(/0.0,0.5/)
c
      U_total_nodal=RESHAPE(U,(/3,12/))
      U_total_nodal_up=U_total_nodal(:,1:LAYERNODE)
      U_total_nodal_low=U_total_nodal(:,LAYERNODE+1:LAYERNODE*2)
      Coord_t_up=COORDS(:,1:LAYERNODE)+U_total_nodal_up
      Coord_t_low=COORDS(:,LAYERNODE+1:LAYERNODE*2)+U_total_nodal_low
c      write (7,*) '3. U OF UP'
c      write (7,*) U_total_nodal_up
c      write (7,*) '4. U OF low'
c      write (7,*) U_total_nodal_low
      Du_nodal=RESHAPE(DU(:,1),(/3,12/))
      Du_up=Du_nodal(:,1:LAYERNODE)
      Du_low=Du_nodal(:,LAYERNODE+1:LAYERNODE*2)
      Du_seperation=Du_up-Du_low
c      write (7,*) '3. Du_seperation'
c      write (7,*) Du_seperation
c      write (7,*) '3. DU OF UP'
c      write (7,*) Du_up
c      write (7,*) '4. DU OF low'
c      write (7,*) Du_up
c-----------------------------------------
c---------- caculate the czm model geometry -------------
      Coord_mid=0.5*(Coord_t_low+Coord_t_up)! x_m
      Seperation_total_nodal=Coord_t_up-Coord_t_low! [u]
c      write (7,*) '*** Seperation at node is'
c      write (7,*) Seperation_total_nodal
c =============  ================
      do i=1,LAYERNODE
c=========== calculate the geometric invariations =================
      Seperation_nodal=Seperation_total_nodal(:,i)
c   - Jacobian(local_coordinate,global_coord,output,detJ) ------
      call Jacobian(Nodal_Local_coord(:,i),Coord_mid,J_mid,detJ(:,i))
c            write (7,*) '*** Jacobian midsurface at point',i,'is'
c            write (7,*) J_mid
c   - g_beta at nodal point --------------------------
      g1=J_mid(:,1)
      g2=J_mid(:,2)
      call Base_g3(g1,g2,g3)
      call Base_gc(g1,g2,g3,gc1,gc2,gc3)
c   - I(u,g)
c   --- I1 (1,6) scaler
      I1=DOT_PRODUCT(Seperation_nodal,Seperation_nodal)
C   --- I2 (1,6) scaler
      call Product_I(g1,gc1,Seperation_nodal,I2)
C   --- I3 (1,6) scaler
      call Product_I(g2,gc2,Seperation_nodal,I3)
C   --- I4 (1,6) scaler
      call Product_I(g3,gc3,Seperation_nodal,I4)
      if ( (I1-I2-I3-I4)*Qn0 > 1e-5) then
            write (7,*) 'Seperation_nodal is: ',Seperation_nodal
            write (7,*) '!!!!!! I1I2I3I4 DISequal '
      end if
c===== damage parameter passin from last increasment =====
      kn_n0=SVARS(i)
      kt_n0=SVARS(6+i)
c      tempK=PROPS(7)
c      tempK2=PROPS(8)
c===== decide if the normal seperation is minus ===============
      if (DOT_PRODUCT(Du_seperation(:,i),g3)<0.) then ! if the crack is closing
c        --- jugde inject ---
            if (DOT_PRODUCT(Seperation_nodal,g3)<0.) then ! inject
                  write (7,*) 'WARN !! : OVERLAP! Seperation_nodal is MINUS: ', DOT_PRODUCT(Seperation_nodal,g3)
                  cn=PROPS(1)
                  Gn=PROPS(3)
                  Qn0=PROPS(5)
                  kn0=0.5*PROPS(5)*PROPS(5)/PROPS(1)
                  dnn=1-EXP((kn0-kn_n0)*PROPS(5)/PROPS(3))
                  dtn=dnn
                  Sai0n=0.5*cn*(I1-I2-I3) ! ψ0(I1,I2,I3)~(u,g) (3,6) scaler
            else ! not inject
                  write (7,*) 'WARN !: THE CRACK IS CLOSING: ', DOT_PRODUCT(Du_seperation(:,i),g3)
                  cn=PROPS(1)
                  Gn=PROPS(3)
                  Qn0=PROPS(5)
                  Sai0n=0.5*cn*(I1-I2-I3) ! ψ0(I1,I2,I3)~(u,g) (3,6) scaler
c              --- calculate the normal damage parameters (3,6) scaler
                  kn0=0.5*PROPS(5)*PROPS(5)/PROPS(1)
                  dnn=1-EXP((kn0-kn_n0)*PROPS(5)/PROPS(3))
                  if (dnn<0) then
                        dnn=0.
                  end if
                  dtn=dnn
            end if
c         --- keep the normal damage parameter still ---
      else ! if the crack is not closing
            cn=PROPS(1)
            Gn=PROPS(3)
            Qn0=PROPS(5)
            Sai0n=0.5*cn*(I1-I2-I3) ! ψ0(I1,I2,I3)~(u,g) (3,6) scaler
c        --- calculate the normal damage parameters (3,6) scaler
            kn0=0.5*Qn0*Qn0/cn
            dnn=1-EXP((kn0-kn_n0)*Qn0/Gn)
            if (dnn<0) then
                  dnn=0.
            end if
            dtn=dnn
c         --- update the normal damage parameter ---
            SVARS(i)=MAX(kn_n0,Sai0n,kn0)
      END if
c===== get the shear material parameters =====
      ct=PROPS(2)
      Gt=PROPS(4)
      Qt0=PROPS(6)
      Sai0t=0.5*ct*(I2+I3) ! ψ0(I1,I2,I3)~(u,g) (3,6) scaler
      kt0=0.5*Qt0*Qt0/ct
      dtt=1-EXP((kt0-kt_n0)*Qt0/Gt)
      if (dtt<0.) then
            dtt=0.
      end if
      dnt=dtt
      SVARS(i+6)=MAX(kt_n0,Sai0t,kt0)

c      write (7,*) '3. (k0) kn0 and kt0 of',i,' is'
c      write (7,*) kn0,kt0
c      write (7,*) '*** k of',i,'at last step is'
c      write (7,*) 'kn_n0', kn_n0
c      write (7,*) 'kt_n0', kt_n0
c      write (7,*) '*** Sai0n of',i, 'is',Sai0n

c ==== ∂ψ0/∂I
      PD_sai0n_I1=0.5*cn
      PD_sai0n_I2=-0.5*cn
      PD_sai0n_I3=-0.5*cn
      PD_sai0t_I1=0.0
      PD_sai0t_I2=0.5*ct
      PD_sai0t_I3=0.5*ct

c     ψ=(1-d)ψ0 scaler
      Sain=(1-dnn)*(1-dnt)*Sai0n
      Sait=(1-dtn)*(1-dtt)*Sai0t

c     [u].gc scaler
         Param_ugc1=DOT_PRODUCT(Seperation_nodal,gc1)
         Param_ugc2=DOT_PRODUCT(Seperation_nodal,gc2)
         Param_ugc3=DOT_PRODUCT(Seperation_nodal,gc3)
c     [u].g scaler
            Param_ug1=DOT_PRODUCT(Seperation_nodal,g1)
            Param_ug2=DOT_PRODUCT(Seperation_nodal,g2)
c     gc12 scaler
            Param_gc12=DOT_PRODUCT(gc1,gc2)
c     gc11 scaler
            Param_gc11=DOT_PRODUCT(gc1,gc1)
c     gc22 scaler
            Param_gc22=DOT_PRODUCT(gc2,gc2)
c     -([u].gc1)gc1+gc11([u].gc3)gc3 vector
            Param_u_gc11=-Param_ugc1*gc1+Param_gc11*Param_ugc3*gc3
c     -([u].gc2)gc1+gc12([u].gc3)gc3 vector
            Param_u_gc21=-Param_ugc2*gc1+Param_gc12*Param_ugc3*gc3
c     -([u].gc1)gc2+gc12([u].gc3)gc3 vector
            Param_u_gc12=-Param_ugc1*gc2+Param_gc12*Param_ugc3*gc3
c     -([u].gc2)gc2+gc22([u].gc3)gc3 vector
            Param_u_gc22=-Param_ugc2*gc2+Param_gc22*Param_ugc3*gc3
c     -gc1×gc1+gc11(gc3×gc3) (3,3) √
            Param_gc1_gc1=-TensorProduct(gc1,gc1)+Param_gc11*TensorProduct(gc3,gc3)
c     -gc1×gc2+gc12(gc3×gc3) (3,3) √
            Param_gc1_gc2=-TensorProduct(gc1,gc2)+Param_gc12*TensorProduct(gc3,gc3)
c     -gc2×gc1+gc12(gc3×gc3) (3,3) √
            Param_gc2_gc1=-TensorProduct(gc2,gc1)+Param_gc12*TensorProduct(gc3,gc3)
c     -gc2×gc2+gc22(gc3×gc3) (3,3) √
            Param_gc2_gc2=-TensorProduct(gc2,gc2)+Param_gc22*TensorProduct(gc3,gc3)
c     -(gc2.gc2)gc1-(gc1.gc2)gc2
            Param_dot_gc2=-Param_gc22*gc1-Param_gc12*gc2
c     -(gc1.gc2)gc1-(gc1.gc1)gc2
            Param_dot_gc1_2=-Param_gc12*gc1-Param_gc11*gc2
c     ∂I/∂u (3,1)
            PD_I1_u=2*Seperation_nodal   
            PD_I2_u=Param_ugc1*g1+Param_ug1*gc1
            PD_I3_u=Param_ugc2*g2+Param_ug2*gc2
c     ∂I/∂g_beta (3,1)
            PD_I2_g1=Param_ugc1*Seperation_nodal+Param_ug1*Param_u_gc11
            PD_I2_g2=-Param_ug1*Param_u_gc12
            PD_I3_g1=Param_ug2*Param_u_gc12
            PD_I3_g2=Param_ugc2*Seperation_nodal-Param_ug2*Param_u_gc22
c     ∂(∂I/∂u)/∂u (3,3)
            PD_I1_u_u=2*Matrix_1
            PD_I2_u_u=TensorProduct(gc1,g1)+TensorProduct(g1,gc1)
            PD_I3_u_u=TensorProduct(g2,gc2)+TensorProduct(gc2,g2)
c     ∂(∂I/∂u)∂g (3,3)
            PD_I2_u_g1=TensorProduct(Seperation_nodal,gc1)
     &+Param_ugc1*Matrix_1+TensorProduct(Param_u_gc11,g1)+Param_ug1*Param_gc1_gc1
            PD_I2_u_g2=TensorProduct(-Param_u_gc21,g1)-Param_ug1*Param_gc1_gc2
            PD_I3_u_g1=TensorProduct(Param_u_gc12,g2)+Param_ug2*Param_gc2_gc1
            PD_I3_u_g2=TensorProduct(Seperation_nodal,gc2)
     &+Param_ugc2*Matrix_1+TensorProduct(-Param_u_gc22,g2)-Param_ug2*Param_gc2_gc2
c     ∂(∂I/∂g)∂g (3,3)
            PD_I2_g1_g1=TensorProduct(Seperation_nodal,Param_u_gc11)
     &+TensorProduct(Param_u_gc11,Seperation_nodal)+2*Param_ug1
     &*(Param_ugc1*(-Param_gc1_gc1)-Param_gc11*Param_ugc3*
     &(TensorProduct(gc1,gc3)+TensorProduct(gc3,gc1)))
            PD_I2_g2_g2=Param_ug1*(TensorProduct(gc1,-Param_u_gc22)-Param_ugc2*Param_gc2_gc1
     &-Param_gc12*Param_ugc2*TensorProduct(gc3,gc3)+Param_ugc3*TensorProduct(gc3,Param_dot_gc2)
     &-Param_gc12*Param_ugc3*TensorProduct(gc2,gc3))
            PD_I2_g1_g2=TensorProduct(-Param_u_gc21,Seperation_nodal)+Param_ug1*
     &(TensorProduct(Param_u_gc21,gc1)+TensorProduct(Param_u_gc11,gc2)
     &+Param_gc12*Param_ugc1*TensorProduct(gc3,gc3)
     &+Param_gc11*Param_ugc2*TensorProduct(gc3,gc3)+2*Param_ugc3*Param_gc12*TensorProduct(gc1,gc3))
            PD_I3_g1_g1=Param_ug2*(TensorProduct(gc2,-Param_u_gc11)-Param_ugc1*Param_gc1_gc2
     &-Param_gc12*Param_ugc1*TensorProduct(gc3,gc3)+Param_ugc3*TensorProduct(gc3,Param_dot_gc1_2)
     &-Param_gc12*Param_ugc3*TensorProduct(gc1,gc3))
            PD_I3_g2_g2=TensorProduct(Seperation_nodal,-Param_u_gc22)
     &+TensorProduct(-Param_u_gc22,Seperation_nodal)+2*Param_ug2
     &*(Param_ugc2*(-Param_gc2_gc2)-Param_gc22*Param_ugc3*
     &(TensorProduct(gc2,gc3)+TensorProduct(gc3,gc2)))
            PD_I3_g1_g2=TensorProduct(Seperation_nodal,Param_u_gc12)+Param_ug2*
     &(TensorProduct(Param_u_gc21,gc2)+TensorProduct(Param_u_gc12,gc2)
     &+Param_gc22*Param_ugc1*TensorProduct(gc3,gc3)
     &+Param_gc12*Param_ugc2*TensorProduct(gc3,gc3)
     &+Param_ugc3*(Param_gc12*TensorProduct(gc2,gc3)+Param_gc22*TensorProduct(gc1,gc3)))
c     ∂ψ0/∂u
      PD_sai0n_u=PD_sai0n_I1*PD_I1_u+PD_sai0n_I2*PD_I2_u+PD_sai0n_I3*PD_I3_u
      PD_sai0t_u=PD_sai0t_I1*PD_I1_u+PD_sai0t_I2*PD_I2_u+PD_sai0t_I3*PD_I3_u
c     ∂ψ0/∂g
      PD_sai0n_g1=PD_sai0n_I1*PD_I1_g1+PD_sai0n_I2*PD_I2_g1+PD_sai0n_I3*PD_I3_g1
      PD_sai0n_g2=PD_sai0n_I1*PD_I1_g2+PD_sai0n_I2*PD_I2_g2+PD_sai0n_I3*PD_I3_g2
      PD_sai0t_g1=PD_sai0t_I1*PD_I1_g1+PD_sai0t_I2*PD_I2_g1+PD_sai0t_I3*PD_I3_g1
      PD_sai0t_g2=PD_sai0t_I1*PD_I1_g2+PD_sai0t_I2*PD_I2_g2+PD_sai0t_I3*PD_I3_g2
c      write (7,*) '*** Determinate of Jacobian at point',i,'is'
c      write (7,*) detJ
c      write (7,*) '*** g1 is'
c      write (7,*) g1
c      write (7,*) '*** g2 is'
c      write (7,*) g2
c      write (7,*) '*** g3 is'
c      write (7,*) g3     
c      write (7,*) '*** gc1 is'
c      write (7,*) gc1
c      write (7,*) '*** gc2 is'
c      write (7,*) gc2
c      write (7,*) '*** gc3 is'
c      write (7,*) gc3
c      write(7,*) '*** I'
c      write(7,*) '      I1'
c      write(7,*) I1
c      write(7,*) '      I2'
c      write(7,*) I2
c      write(7,*) '      I3'
c      write(7,*) I3
c      write(7,*) '*** ψ0'
c      write(7,*) 'ψn0',Sai0n
c      write(7,*) 'ψt0',Sai0t
c      write(7,*) '*** ψ'
c      write(7,*) 'ψn',Sain
c      write(7,*) 'ψt',Sait
c      write (7,*) '*** k of',i,'at end is'
c            write (7,*) 'kn', SVARS(i)
c            write (7,*) 'kt', SVARS(6+i)
      write(7,*) '*** damage parameters'
      write(7,*) '      dnn',dnn
      write(7,*) '      dtt',dtt
c      write(7,*) '      dtn',dtn
c      write(7,*) '      dtt',dtt
c      write(7,*)'***  ∂I1/∂u'
c      write(7,*) PD_I1_u
c      write(7,*)'***  ∂I2/∂u'
c      write(7,*) PD_I2_u
c      write(7,*) '*** ∂I2/∂g1'
c      write(7,*) PD_I2_g1
c      write(7,*) '*** ∂I3/∂g1'
c      write(7,*) PD_I3_g1
c      write(7,*)'***  ∂I3/∂u'
c      write(7,*) PD_I3_u
c      write(7,*) '*** ∂I3/∂g1'
c      write(7,*) PD_I3_g1
c      write(7,*) '*** ∂I2/∂u/∂u'
c      write(7,*) PD_I2_u_u
c      write(7,*) '*** ∂I2/∂u/∂g1'
c      write(7,*) PD_I2_u_g1
c      write(7,*) '*** ∂I2/∂u/∂g2'
c      write(7,*) PD_I2_u_g2
c      write(7,*) '*** ∂I2/∂g1/∂g1'
c      write(7,*) PD_I2_g1_g1
c      write(7,*) '*** ∂I2/∂g2/∂g2'
c      write(7,*) PD_I2_g2_g2
c      write(7,*) '*** ∂I2/∂g1/∂g2'
c      write(7,*) PD_I2_g1_g2
c      write(7,*) '*** ∂I3/∂u/∂u'
c      write(7,*) PD_I3_u_u
c      write(7,*) '*** ∂I3/∂u/∂g1'
c      write(7,*) PD_I3_u_g1
c      write(7,*) '*** ∂I3/∂u/∂g2'
c      write(7,*) PD_I3_u_g2
c      write(7,*) '*** ∂I3/∂g1/∂g1'
c      write(7,*) PD_I3_g1_g1
c      write(7,*) '*** ∂I3/∂g2/∂g2'
c      write(7,*) PD_I3_g2_g2
c      write(7,*) '*** ∂I3/∂g1/∂g2'
c      write(7,*) PD_I3_g1_g2      
c      write(7,*) '*** PD_sai0n_u'
c      write(7,*) PD_sai0n_u
c      write(7,*) '*** PD_sai0t_u'
c      write(7,*) PD_sai0t_u
c      write(7,*) '*** PD_sai0t_g1      '
c      write(7,*) PD_sai0t_g1
c      write(7,*) '*** PD_sai0t_g2      '
c      write(7,*) PD_sai0t_g2
c      write (7,*) 'k update is'
c      write (7,*) SVARS(25),SVARS(26),SVARS(27),SVARS(28),SVARS(29),SVARS(30)
c      write (7,*) '*** I1,I2,I3 at Gauss point is'
c      write (7,*) I1
c      write (7,*) I2
c      write (7,*) I3
c================= K calculating =======================================
c      write(7,*) '*** start calculate ∂(∂ψ0/∂u)∂u ∂(∂ψ0/∂u)∂g_β ∂(∂ψ0/∂g_β)∂g_β'
c      ∂(∂ψ0/∂u)∂u (3,3)
      PD_sai0n_u_u=PD_sai0n_I1*PD_I1_u_u+PD_sai0n_I2*PD_I2_u_u-PD_sai0n_I3*PD_I3_u_u
      PD_sai0t_u_u=PD_sai0t_I1*PD_I1_u_u+PD_sai0t_I2*PD_I2_u_u-PD_sai0t_I3*PD_I3_u_u
c      ∂(∂ψ0/∂u)∂g1 (3,3)
      PD_sai0n_g1_u=PD_sai0n_I2*PD_I2_u_g1+PD_sai0n_I3*PD_I3_u_g1
      PD_sai0t_g1_u=PD_sai0t_I2*PD_I2_u_g1+PD_sai0t_I3*PD_I3_u_g1
c      ∂(∂ψ0/∂u)∂g2 (3,3)
      PD_sai0n_g2_u=PD_sai0n_I2*PD_I2_u_g2+PD_sai0n_I3*PD_I3_u_g2
      PD_sai0t_g2_u=PD_sai0t_I2*PD_I2_u_g2+PD_sai0t_I3*PD_I3_u_g2
c      ∂(∂ψ0/∂g1)∂g1 (3,3)
      PD_sai0n_g1_g1=PD_sai0n_I2*PD_I2_g1_g1+PD_sai0n_I3*PD_I3_g1_g1
      PD_sai0t_g1_g1=PD_sai0t_I2*PD_I2_g1_g1+PD_sai0t_I3*PD_I3_g1_g1
c      ∂(∂ψ0/∂g1)∂g2 (3,3)
      PD_sai0n_g2_g1=PD_sai0n_I2*PD_I2_g1_g2+PD_sai0n_I3*PD_I3_g1_g2
      PD_sai0t_g2_g1=PD_sai0t_I2*PD_I2_g1_g2+PD_sai0t_I3*PD_I3_g1_g2
c      ∂(∂ψ0/∂g2)∂g2 (3,3)
      PD_sai0n_g2_g2=PD_sai0n_I2*PD_I2_g2_g2+PD_sai0n_I3*PD_I3_g2_g2
      PD_sai0t_g2_g2=PD_sai0t_I2*PD_I2_g2_g2+PD_sai0t_I3*PD_I3_g2_g2
c      ∂(∂ψ0/∂g1)∂u
      PD_sai0n_u_g1=TRANSPOSE(PD_sai0n_g1_u)
      PD_sai0t_u_g1=TRANSPOSE(PD_sai0t_g1_u)
c      ∂(∂ψ0/∂g2)∂u
      PD_sai0n_u_g2=TRANSPOSE(PD_sai0n_g2_u)
      PD_sai0t_u_g2=TRANSPOSE(PD_sai0t_g2_u)
c      ∂(∂ψ0/∂g2)∂g1
      PD_sai0n_g1_g2=TRANSPOSE(PD_sai0n_g2_g1)
      PD_sai0t_g1_g2=TRANSPOSE(PD_sai0n_g2_g1)
c      write(7,*) '-------------------------'
c      write(7,*) PD_I1_u_u
c      write(7,*) PD_I2_u_u
c      write(7,*) PD_I3_u_u
c      write(7,*) '-------------------------'
c      write(7,*) '∂(∂ψ0/∂u)∂u PD_sai0n_u_u'
c      write(7,*) PD_sai0n_u_u
c      write(7,*) '∂(∂ψ0/∂u)∂u PD_sai0t_u_u'
c      write(7,*) PD_sai0t_u_u
c      write(7,*) '∂(∂ψ0/∂u)∂g1 PD_sai0n_u_g1'
c      write(7,*) PD_sai0n_u_g1
c      write(7,*) '∂(∂ψ0/∂u)∂g2 PD_sai0n_u_g2'
c      write(7,*) PD_sai0n_u_g2
c      write(7,*) '∂(∂ψ0/∂g1)∂g2 PD_sai0n_g1_g2'
c      write(7,*) PD_sai0n_g1_g2
c====================================
c     ============ N(ξ) & ∂N/∂ξ  (3,6)===========================================
      N=Shape_poly(Nodal_Local_coord(:,i))
      PD_N_r=Shape_poly_PD_r(Nodal_Local_coord(:,i))
      PD_N_s=Shape_poly_PD_s(Nodal_Local_coord(:,i))
c      write(7,*) '*** N'
c      write(7,*) N
c      write(7,*) '*** PD_N_r'
c      write(7,*) PD_N_r
c      write(7,*) '*** PD_N_s'
c      write(7,*) PD_N_s
c     ----------------------------------
         do j=1,6
c ============  Δ[u] & Δg (3)================
            do m=1,3
      delta_U=0.0;delta_g1=0.0;delta_g2=0.0
      delta_U(m)=N(1,j)
      delta_g1(m)=PD_N_r(1,j)*0.5
      delta_g2(m)=PD_N_s(1,j)*0.5
c      write(7,*) '      delta_u'
c      write(7,*) delta_U
c      write(7,*) '      delta_g1'
c      write(7,*) delta_g1
c      write(7,*) '      delta_g2'
c      write(7,*) delta_g2
c      write(7,*) '∂(∂ψ0/∂u)∂u PD_sai0n_u_u'
c      write(7,*) PD_sai0n_u_u
c      write(7,*) '∂(∂ψ0/∂u)∂u PD_sai0t_u_u'
c      write(7,*) PD_sai0t_u_u
c===============   Δ∂ψ0/∂u (3,6) ====================================
c      write(7,*) '*** start calculate Δ∂ψ0/∂u '
      delta_PD_sai0n_u_up=Dot_3_1(PD_sai0n_u_u,delta_U)+Dot_3_1(PD_sai0n_u_g1,delta_g1)
     &+Dot_3_1(PD_sai0n_u_g2,delta_g2)
      delta_PD_sai0t_u_up=Dot_3_1(PD_sai0t_u_u,delta_U)+Dot_3_1(PD_sai0t_u_g1,delta_g1)
     &+Dot_3_1(PD_sai0t_u_g2,delta_g2)
      delta_PD_sai0n_u_low=-Dot_3_1(PD_sai0n_u_u,delta_U)+Dot_3_1(PD_sai0n_u_g1,delta_g1)
     &+Dot_3_1(PD_sai0n_u_g2,delta_g2)
      delta_PD_sai0t_u_low=-Dot_3_1(PD_sai0t_u_u,delta_U)+Dot_3_1(PD_sai0t_u_g1,delta_g1)
     &+Dot_3_1(PD_sai0t_u_g2,delta_g2)
      delta_PD_sai0n_g1_up=Dot_3_1(PD_sai0n_g1_g2,delta_g2)+Dot_3_1(PD_sai0n_g1_u,delta_U)
      delta_PD_sai0t_g1_up=Dot_3_1(PD_sai0t_g1_g2,delta_g2)+Dot_3_1(PD_sai0t_g1_u,delta_U)
      delta_PD_sai0n_g2_up=Dot_3_1(PD_sai0n_g2_g1,delta_g1)+Dot_3_1(PD_sai0n_g2_u,delta_U)
      delta_PD_sai0t_g2_up=Dot_3_1(PD_sai0t_g2_g1,delta_g1)+Dot_3_1(PD_sai0t_g2_u,delta_U)
      delta_PD_sai0n_g1_low=Dot_3_1(PD_sai0n_g1_g2,delta_g2)-Dot_3_1(PD_sai0n_g1_u,delta_U)
      delta_PD_sai0t_g1_low=Dot_3_1(PD_sai0t_g1_g2,delta_g2)-Dot_3_1(PD_sai0t_g1_u,delta_U)
      delta_PD_sai0n_g2_low=Dot_3_1(PD_sai0n_g2_g1,delta_g1)-Dot_3_1(PD_sai0n_g2_u,delta_U)
      delta_PD_sai0t_g2_low=Dot_3_1(PD_sai0t_g2_g1,delta_g1)-Dot_3_1(PD_sai0t_g2_u,delta_U)
c      write(7,*) '*** delta_PD_sai0n_u'
c      write(7,*) delta_PD_sai0n_u_low
c      write(7,*) '*** delta_PD_sai0n_g1'
c      write(7,*) delta_PD_sai0n_g1_low
c      write(7,*) '*** delta_PD_sai0n_g2'
c      write(7,*) delta_PD_sai0n_g2_low
c      write(7,*) '*** delta_PD_sai0y_u'
c      write(7,*) delta_PD_sai0t_u_low
c      write(7,*) '*** delta_PD_sai0t_g1'
c      write(7,*) delta_PD_sai0t_g1_up
c      write(7,*) '*** delta_PD_sai0t_g2'
c      write(7,*) delta_PD_sai0t_g2_up
c================= Δd    =====================================
c      write(7,*) '*** start calculate Δd '
c      ∂ψ0n scaler
      delta_sai0n_up=DOT_PRODUCT(PD_sai0n_u,delta_U)+DOT_PRODUCT(PD_sai0n_g1,delta_g1)
     &+DOT_PRODUCT(PD_sai0t_g1,delta_g1)
      delta_sai0n_low=-DOT_PRODUCT(PD_sai0n_u,delta_U)+DOT_PRODUCT(PD_sai0n_g1,delta_g1)
     &+DOT_PRODUCT(PD_sai0t_g1,delta_g1)
c      ∂ψ0t scaler
      delta_sai0t_up=DOT_PRODUCT(PD_sai0t_u,delta_U)+DOT_PRODUCT(PD_sai0n_g2,delta_g2)
     &+DOT_PRODUCT(PD_sai0t_g2,delta_g2)
      delta_sai0t_low=-DOT_PRODUCT(PD_sai0t_u,delta_U)+DOT_PRODUCT(PD_sai0n_g2,delta_g2)
     &+DOT_PRODUCT(PD_sai0t_g2,delta_g2)
c      write(7,*) '*** Δψ0'
c      write(7,*) '      Δψ0n'
c      write(7,*) delta_sai0n_up
c      write(7,*) '      Δψ0t'
c      write(7,*) delta_sai0t_up
      delta_dnn_up=(1-dnn)*Qn0/Gn*delta_sai0n_up !scaler
      delta_dnt_up=(1-dnt)*Qt0/Gt*delta_sai0t_up !scaler
      delta_dtn_up=(1-dtn)*Qn0/Gn*delta_sai0n_up !scaler
      delta_dtt_up=(1-dtt)*Qt0/Gt*delta_sai0t_up !scaler
      delta_dnn_low=(1-dnn)*Qn0/Gn*delta_sai0n_low !scaler
      delta_dnt_low=(1-dnt)*Qt0/Gt*delta_sai0t_low !scaler
      delta_dtn_low=(1-dtn)*Qn0/Gn*delta_sai0n_low !scaler
      delta_dtt_low=(1-dtt)*Qt0/Gt*delta_sai0t_low !scaler
c      write(7,*) '*** Δd'
c      write(7,*) '      Δdnn'
c      write(7,*) delta_dnn_low
c      write(7,*) '      Δdnt'
c      write(7,*) delta_dnt_low
c      write(7,*) '      Δdtn'
c      write(7,*) delta_dtn_low
c      write(7,*) '      Δdtt'
c      write(7,*) delta_dtt_low
c     K--
      K=(1-dnn)*(1-dnt)*(-N(1,i)*delta_PD_sai0n_u_low+0.5*PD_N_r(1,i)*delta_PD_sai0n_g1_low+0.5*PD_N_s(1,i)*delta_PD_sai0n_g2_low)
     &-delta_dnn_low*(1-dnt)*(-N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &-(1-dnn)*delta_dnt_low*(-N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &+(1-dtn)*(1-dtt)*(-N(1,i)*delta_PD_sai0t_u_low+0.5*PD_N_r(1,i)*delta_PD_sai0t_g1_low+0.5*PD_N_s(1,i)*delta_PD_sai0t_g2_low)
     &-delta_dtn_low*(1-dtt)*(-N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
     &-(1-dtn)*delta_dtt_low*(-N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
      AMATRX((j-1)*3+19:(j-1)*3+21,(i-1)*3+m+18)=K
c      write(7,*) '*** K--, i',i,'j',j
c      write(7,*) K
c     K+-
      K=(1-dnn)*(1-dnt)*(N(1,i)*delta_PD_sai0n_u_low+0.5*PD_N_r(1,i)*delta_PD_sai0n_g1_low+0.5*PD_N_s(1,i)*delta_PD_sai0n_g2_low)
     &-delta_dnn_low*(1-dnt)*(N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &-(1-dnn)*delta_dnt_low*(N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &+(1-dtn)*(1-dtt)*(N(1,i)*delta_PD_sai0t_u_low+0.5*PD_N_r(1,i)*delta_PD_sai0t_g1_low+0.5*PD_N_s(1,i)*delta_PD_sai0t_g2_low)
     &-delta_dtn_low*(1-dtt)*(N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
     &-(1-dtn)*delta_dtt_low*(N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
      AMATRX((j-1)*3+1:(j-1)*3+3,(i-1)*3+m+18)=K
c      write(7,*) '*** K+-, i',i,'j',j
c      write(7,*) K
c     K-+
      K=(1-dnn)*(1-dnt)*(-N(1,i)*delta_PD_sai0n_u_up+0.5*PD_N_r(1,i)*delta_PD_sai0n_g1_up+0.5*PD_N_s(1,i)*delta_PD_sai0n_g2_up)
     &-delta_dnn_up*(1-dnt)*(-N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &-(1-dnn)*delta_dnt_up*(-N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &+(1-dtn)*(1-dtt)*(-N(1,i)*delta_PD_sai0t_u_up+0.5*PD_N_r(1,i)*delta_PD_sai0t_g1_up+0.5*PD_N_s(1,i)*delta_PD_sai0t_g2_up)
     &-delta_dtn_up*(1-dtt)*(-N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
     &-(1-dtn)*delta_dtt_up*(-N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
      AMATRX((j-1)*3+19:(j-1)*3+21,(i-1)*3+m)=K
c      write(7,*) '*** K-+, i',i,'j',j
c      write(7,*) K
c     K++
      K=(1-dnn)*(1-dnt)*(N(1,i)*delta_PD_sai0n_u_up+0.5*PD_N_r(1,i)*delta_PD_sai0n_g1_up+0.5*PD_N_s(1,i)*delta_PD_sai0n_g2_up)
     &-delta_dnn_up*(1-dnt)*(N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &-(1-dnn)*delta_dnt_up*(N(1,i)*PD_sai0n_u+0.5*PD_N_r(1,i)*PD_sai0n_g1+0.5*PD_N_s(1,i)*PD_sai0n_g2)
     &+(1-dtn)*(1-dtt)*(N(1,i)*delta_PD_sai0t_u_up+0.5*PD_N_r(1,i)*delta_PD_sai0t_g1_up+0.5*PD_N_s(1,i)*delta_PD_sai0t_g2_up)
     &-delta_dtn_up*(1-dtt)*(N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
     &-(1-dtn)*delta_dtt_up*(N(1,i)*PD_sai0t_u+0.5*PD_N_r(1,i)*PD_sai0t_g1+0.5*PD_N_s(1,i)*PD_sai0t_g2)
      AMATRX((j-1)*3+1:(j-1)*3+3,(i-1)*3+m)=K
c      write(7,*) '*** K++, i',i,'j',j
      
c==============================================================
            end do !m
            end do !j
c =========== parameters using in numerical integration ==================
c     ∂ψ/∂u (3,6)
      R_u_up_nodal(:,i)=(1-dnn)*(1-dnt)*PD_sai0n_u+(1-dtn)*(1-dtt)*PD_sai0t_u
      R_u_low_nodal(:,i)=R_u_up_nodal(:,i)*(-1)
c     ∂ψ/∂g (3,6)
      R_g1_nodal(:,i)=(1-dnn)*(1-dnt)*PD_sai0n_g1+(1-dtn)*(1-dtt)*PD_sai0t_g1
      R_g2_nodal(:,i)=(1-dnn)*(1-dnt)*PD_sai0n_g2+(1-dtn)*(1-dtt)*PD_sai0t_g2
c      write(7,*) '*** R_u_low_nodal'
c      write(7,*) R_u_low_nodal
c ===========================================================================
      end do !i
c      write(110,*) '*** KK'
c      write(110,"(36f12.1)") AMATRX
c ===============================================
c      write(7,*) '************************'
c      write(7,*) '***** RHS ASSEMBLE *****'
c      write(7,*) '************************'
      call Integration(weight,R_u_low_nodal,R_u_low,Gauss_point,detJ)
      call Integration(weight,R_u_up_nodal,R_u_up,Gauss_point,detJ)
      call Integration_g1(weight,R_g1_nodal,R_g1,Gauss_point,detJ)
      call Integration_g2(weight,R_g2_nodal,R_g2,Gauss_point,detJ)
c      write(7,*) '*** R_u_low'
c      write(7,*) R_u_low
c      write(7,*) '*** R_u_up'
c      write(7,*) R_u_up
c      write(7,*) '*** R_g1'
c      write(7,*) R_g1
c      write(7,*) '*** R_g2'
c      write(7,*) R_g2
      R_up=R_u_up+R_g1+R_g2
      R_low=R_u_low+R_g1+R_g2
      do i=1,6
            RHS((i-1)*3+1:(i-1)*3+3,1)=-R_up(:,i)
            RHS((i-1)*3+19:(i-1)*3+21,1)=-R_low(:,i)
      end do
c      write(7,*) '*** RHS'
c      write(7,*) RHS(:,1)
      write(7,*) '*** T_bar'
      write(7,*) R_u_up
      write(7,*) '*** DIV P_bar'
      write(7,*) R_g1
      write(7,*) R_g2
c      write(7,*) '*** R_low'
c      write(7,*) R_low
c      write(110,*) '*** RHS'
c      write(110,"(36f12.1)") RHS(:,1)
      write(7,*) '*************************'
      write(7,*) '***** END     CYCLE *****'
      write(7,*) '*************************'
      close(110)
c ============================================
c ============== END of MAIN =================
c ============================================
      return
      END
c======================================================================
c=============================SUBROUTINES==============================
c======================================================================
c tensor product
      function TensorProduct(g,gc)
      real(Kind=8) :: g(3),gc(3),TensorProduct(3,3),k
c      write(7,*) '------ SUBROUTIE TensorProduct ------'
c      write(7,*) '  A'
c      write(7,*) g
c      write(7,*) '  B'
c      write(7,*) gc
      do i=1,3
            TensorProduct(:,i)=g
c            write(7,*) '  colomn',i,'is'
c            write(7,*) TensorProduct(:,i)
            TensorProduct(:,i)=TensorProduct(:,i)*gc(i)
c            write(7,*) '  product',i,'is'
c            write(7,*) TensorProduct(:,i)
      end do
c      write(7,*) '  the whole Product is'
c      write(7,*) TensorProduct
      return
      end function
c-------------------------------------
      function Dot_Matrix(A,B,NODE)
        real(Kind=8) :: A(3,NODE),B(3,NODE),Dot_Matrix(3,NODE),Temp(3,NODE),Temp2(1,NODE)
        INTEGER :: NODE
        Temp=A*B
        Temp2(1,:)=SUM(Temp,DIM=1)
        do i=1,NODE
            Dot_Matrix(:,i)=Temp2(1,i)
        end do
c        write(7,*) '------ SUBROUTIE Dot_Matrix ------'
c        write(7,*) '  A',A
c        write(7,*) '  B',B
c        write(7,*) '  Product',Dot_Matrix
      return
      end function
c----------- ??? -------------------------
      function Dot_3_1(A,B)
      real(Kind=8) :: A(3,3),B(3),Dot_3_1(3),TEMP(3,3)
      do i=1,3
            TEMP(:,i)=A(:,i)*B(i)
      end do
      Dot_3_1=SUM(TEMP,DIM=2)
c      write(7,*) '------ SUBROUTIE Dot_3_1 ------'
c      write(7,*) '  A'
c      write(7,*) A
c      write(7,*) '  B'
c      write(7,*) B
c      write(7,*) '  Product'
c      write(7,*) Dot_3_1
      return
      end function
c---------------------------------------------------------------------
      function Shape_poly(local_coordinate)
      real(Kind=8) :: local_coordinate(2,1), Shape_poly(3,6)
     &,r,s,h1,h2,h3,h4,h5,h6
      r=local_coordinate(1,1);s=local_coordinate(2,1)
      h4=4*r*(1-r-s)
      h5=4*r*s
      h6=4*s*(1-r-s)
      h1=(1-r-s)-0.5*h4-0.5*h6
      h2=r-0.5*h4-0.5*h5
      h3=s-0.5*h5-0.5*h6
c
      Shape_poly(:,1)=h1
      Shape_poly(:,2)=h2
      Shape_poly(:,3)=h3
      Shape_poly(:,4)=h4
      Shape_poly(:,5)=h5
      Shape_poly(:,6)=h6
      return
      end function
c
      function Shape_poly_PD_r(local_coordinate)
      real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_r(3,6)
     &,r,s,h1,h2,h3,h4,h5,h6
      r=local_coordinate(1,1);s=local_coordinate(2,1)
      h4=4-8*r-4*s
      h5=4*s
      h6=-4*s
      h1=-3+4*(r+s)
      h2=-1+4*r
      h3=0
c
      Shape_poly_PD_r(:,1)=h1
      Shape_poly_PD_r(:,2)=h2
      Shape_poly_PD_r(:,3)=h3
      Shape_poly_PD_r(:,4)=h4
      Shape_poly_PD_r(:,5)=h5
      Shape_poly_PD_r(:,6)=h6
      return
      end function
c
      function Shape_poly_PD_s(local_coordinate)
      real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_s(3,6)
     &,r,s,h1,h2,h3,h4,h5,h6
      r=local_coordinate(1,1);s=local_coordinate(2,1)
      h4=-4*r
      h5=4*r
      h6=4*(1-r-2*s)
      h1=4*(r+s)-3
      h2=0
      h3=4*s-1
c
      Shape_poly_PD_s(:,1)=h1
      Shape_poly_PD_s(:,2)=h2
      Shape_poly_PD_s(:,3)=h3
      Shape_poly_PD_s(:,4)=h4
      Shape_poly_PD_s(:,5)=h5
      Shape_poly_PD_s(:,6)=h6
      return
      end function
c--------------------------------------
c               (g1,gc1,Seperation_total_nodal,I2
      subroutine Product_I(g,gc,u,output)
      real(Kind=8) :: g(3),gc(3),u(3),output
      output=
     &u(1)*(g(1)*gc(1)*u(1) + g(2)*gc(1)*u(2) + g(3)*gc(1)*u(3)) 
     &+u(2)*(g(1)*gc(2)*u(1) + g(2)*gc(2)*u(2) + g(3)*gc(2)*u(3)) 
     &+u(3)*(g(1)*gc(3)*u(1) + g(2)*gc(3)*u(2) + g(3)*gc(3)*u(3))
c          write(7,*) '------ SUBROUTIE Product_I '
c          write(7,*) 'g'
c          write(7,*) g
c          write(7,*) 'gc'
c          write(7,*) gc
c          write(7,*) 'seperation'
c          write(7,*) u
c          write(7,*) 'I'
c          write(7,*) output
      return
      end
c================= Jacobian caculate at one point ===================================
      subroutine Jacobian(local_coordinate,global_coord,output,detJ)
c         ∂r ∂s      
c     x |       |
c     y |       |
c     z |       |
c         g1 g2
c      INCLUDE 'ABA_PARAM.INC'
      real(Kind=8) :: local_coordinate(2,1), output(3,2)
     &,NODE,shape_p(3,6),global_coord(3,6),Temp(3,6),detJ(3,1)
      interface
            function Shape_poly_PD_r(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_r(3,6)
            end function
            function Shape_poly_PD_s(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_s(3,6)
            end function
      end interface
c      write(7,*) '------- Subroutine Jacobian -------'
      shape_p=Shape_poly_PD_r(local_coordinate)!∂Ni/∂r
      output(:,1)=SUM(shape_p*global_coord,DIM=2)!  ∑∂Ni/∂r*x
c      write(7,*) '      output(:,1)'
c      write(7,*) output(:,1)
c
      shape_p=Shape_poly_PD_s(local_coordinate)
      output(:,2)=SUM(shape_p*global_coord,DIM=2)!
c      write(7,*) '      output(:,2)'
c      write(7,*) output(:,2)
c
      call det(output,detJ)
      return
      end
c--------------------------------
      subroutine det(matrix,output)
c      INCLUDE 'ABA_PARAM.INC'
      real(Kind=8) :: matrix(3,2),output(3,1),TEMP(3)
      TEMP(1)=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
      TEMP(2)=matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2)
      TEMP(3)=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
      output=SQRT(SUM(TEMP*TEMP,DIM=1))
      return
      end
c===================================
       subroutine Integration(weight,Var_nodal,output,Guass_coordinate,detJ)
      interface
      function Shape_poly(local_coordinate)
            real(Kind=8) :: local_coordinate(2,1), Shape_poly(3,6)
      end function
      end interface
c     ---------------------------------
      real(Kind=8) :: Shape_fun(3,6),weight,Var_nodal(3,6)
     & ,output(3,6),Guass_coordinate(2,3),detJ(3,6)
      Shape_fun=0.0
      do i=1,3
            Shape_fun=Shape_fun+Shape_poly(Guass_coordinate(:,i))
      end do
      output=Shape_fun*Var_nodal*weight*detJ*0.5
      return
      END
c------------------------
       subroutine Integration_g1(weight,Var_nodal,output,Guass_coordinate,detJ)
      interface
            function Shape_poly_PD_r(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_r(3,6)
            end function
            function Shape_poly_PD_s(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_s(3,6)
            end function
      end interface
c     ---------------------------------
      real(Kind=8) :: Shape_fun(3,6),weight,Var_nodal(3,6),shape_p(3,6)
     & ,output(3,6),Guass_coordinate(2,3),detJ(3,6)
      Shape_fun=0.0
      do i=1,3
            shape_p=Shape_poly_PD_r(Guass_coordinate(:,i))
            Shape_fun=Shape_fun+shape_poly_PD_r(Guass_coordinate(:,i))
      end do
      output=Shape_fun*Var_nodal*weight*detJ*0.5
      return
      END
c-------------------------------
       subroutine Integration_g2(weight,Var_nodal,output,Guass_coordinate,detJ)
      interface
            function Shape_poly_PD_r(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_r(3,6)
            end function
            function Shape_poly_PD_s(local_coordinate)
                  real(Kind=8) :: local_coordinate(2,1), Shape_poly_PD_s(3,6)
            end function
      end interface
c     ---------------------------------
      real(Kind=8) :: Shape_fun(3,6),weight,Var_nodal(3,6)
     & ,output(3,6),Guass_coordinate(2,3),detJ(3,6)
      Shape_fun=0.0
      do i=1,3
           Shape_fun=Shape_fun+Shape_poly_PD_s(Guass_coordinate(:,i))
      end do
      output=Shape_fun*Var_nodal*weight*detJ*0.5
      return
      END
c ============ caculate the g base==============
c-------------------------------------------------
      subroutine Base_cross(g1,g2,output)
      real(Kind=8) :: g1(3), g2(3),output(3)
      output(1)=-g1(3)*g2(2)+g1(2)*g2(3)
      output(2)=g1(3)*g2(1)-g1(1)*g2(3)
      output(3)=-g1(2)*g2(1)+g1(1)*g2(2)
      return
      END
c-------------------------------------------------
      subroutine Base_g3(g1,g2,g3)
      real(Kind=8) :: g1(3), g2(3),g3(3),Norm
c
      Norm=SQRT(ABS(-g1(2)*g2(1)+g1(1)*g2(2))**(2)
     &+ABS(g1(3)*g2(1)-g1(1)*g2(3))**(2)
     &+ABS(-g1(3)*g2(2)+g1(2)*g2(3))**(2))
      call Base_cross(g1,g2,g3)
      g3(1)=g3(1)/Norm
      g3(2)=g3(2)/Norm
      g3(3)=g3(3)/Norm
      return
      END
c-------------------------------------------------
      subroutine Base_gc(g1,g2,g3,gc1,gc2,gc3)
      real(Kind=8) :: g1(3), g2(3),g3(3)
     &,gc1(3),gc2(3),gc3(3),T(3),V
c
      gc3=g3
c
      call Base_cross(g2,g3,T)
      call Base_cross(g1,g3,gc2)
      V=DOT_PRODUCT(g1,T)
      gc2=1/V*gc2
c
      call Base_cross(g2,g3,gc1)
      gc1=1/V*gc1    
c      write(7,*) 'gc2 '
c      write(7,*) gc2
c
      return
      END
c=====================================================









      
