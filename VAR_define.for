        INTEGER,PARAMETER :: SURFACENODE=9,LAYERNODE=6,GaussNODE=3!node number on one side surface
c
        Real(Kind=8) :: Matrix_1(3,3),ALL_POINTS(2,SURFACENODE) ! local coordinates at the mid-surface
     &,Gauss_point(2,GaussNODE)
     &,Seperation_global_nodal(3,LAYERNODE) ! nodal seperation at global coordinate
     &,Seperation_ALL(3,9) ! seperation transfered by local coordinate 
     &,Seperation_nodal(3)
     &,J_mid(3,2),detJ(3,LAYERNODE)
     &,g1(3),g2(3),g3(3) !mid-plane base
     &,gc1(3),gc2(3),gc3(3) !mid-plane contra-base
     &,I1,I2,I3,I4
     &,Sai0n,Sai0t,Sain,Sait
     &,dnn,dnt,dtn,dtt
     &,kn_n0,kt_n0
     &,PD_I1_u(3),PD_I2_u(3),PD_I3_u(3)!partial derivate of I and u
     &,PD_I1_u_u(3,3),PD_I2_u_u(3,3),PD_I3_u_u(3,3)  !∂(∂I/∂u)/∂u
     &,PD_I1_u_g1(3,3),PD_I2_u_g1(3,3),PD_I3_u_g1(3,3)  !∂(∂I/∂u)∂g
     &,PD_I1_u_g2(3,3),PD_I2_u_g2(3,3),PD_I3_u_g2(3,3)  !∂(∂I/∂u)∂g
     &,PD_I1_g1_g1(3,3),PD_I2_g1_g1(3,3),PD_I3_g1_g1(3,3) !∂(∂I/∂g)∂g
     &,PD_I1_g1_g2(3,3),PD_I2_g1_g2(3,3),PD_I3_g1_g2(3,3) !∂(∂I/∂g)∂g
     &,PD_I1_g2_g2(3,3),PD_I2_g2_g2(3,3),PD_I3_g2_g2(3,3) !∂(∂I/∂g)∂g
     &,Param_u_gc11(3),Param_u_gc12(3),Param_u_gc21(3),Param_u_gc22(3)
     &,Param_gc1_gc1(3,3),Param_gc1_gc2(3,3),Param_gc2_gc1(3,3),Param_gc2_gc2(3,3)
     &,Param_dot_gc2(3),Param_dot_gc1_2(3)

c ================== VAR of material properties =============
        INTEGER :: i,j,m
        Real(Kind=8) :: Coord_t_up(3,LAYERNODE),Coord_t_low(3,LAYERNODE),Coord_mid(3,LAYERNODE)
     1,U_Gauss(3,LAYERNODE),du_Gauss(3,LAYERNODE)
     2,Param_ug1,Param_ugc1,Param_ug2,Param_ugc2,Param_ug3,Param_ugc3
     &,Param_gc11,Param_gc12,Param_gc22
c
        Real(Kind=8) :: cn,ct,Gn,Gt,Qn0,Qt0,kn0,kt0,T1,T2,weight(1,3),TEMP3(3,LAYERNODE)
c     &,TEMP4(3,LAYERNODE),TEMP5(3,LAYERNODE),TEMP6(3,LAYERNODE),TEMP7(3,LAYERNODE),TEMP8(3,LAYERNODE)
c
        Real(Kind=8) :: Nodal_Local_coord(2,6),delta_U(3),delta_g1(3),delta_g2(3)
     &,T_delta_g1(3,6,12),T_delta_g2(3,6,12)
     &,N(3,6),PD_N_r(3,6),PD_N_s(3,6),K(3),KK(36,36)
     &,delta_sai0n_up,delta_sai0t_up!Δψ0n,Δψ0t
     &,delta_sai0n_low,delta_sai0t_low!Δψ0n,Δψ0t
     &,delta_PD_sai0n_u_up(3),delta_PD_sai0t_u_up(3)!Δ∂ψ0/∂u
     &,delta_PD_sai0n_g1_up(3),delta_PD_sai0n_g2_up(3) !Δ∂ψ0n/∂g
     &,delta_PD_sai0t_g1_up(3),delta_PD_sai0t_g2_up(3) !Δ∂ψ0t/∂g
     &,delta_PD_sai0n_u_low(3),delta_PD_sai0t_u_low(3)!Δ∂ψ0/∂u
     &,delta_PD_sai0n_g1_low(3),delta_PD_sai0n_g2_low(3) !Δ∂ψ0n/∂g
     &,delta_PD_sai0t_g1_low(3),delta_PD_sai0t_g2_low(3) !Δ∂ψ0t/∂g
     &,PD_sai0n_u_u(3,3),PD_sai0t_u_u(3,3)
     &,PD_sai0n_u_g1(3,3),PD_sai0t_u_g1(3,3),PD_sai0n_u_g2(3,3),PD_sai0t_u_g2(3,3)
     &,PD_sai0n_g1_u(3,3),PD_sai0t_g1_u(3,3),PD_sai0n_g2_u(3,3),PD_sai0t_g2_u(3,3)
     &,PD_sai0n_g1_g1(3,3),PD_sai0t_g1_g1(3,3),PD_sai0n_g1_g2(3,3)
     &,PD_sai0t_g1_g2(3,3),PD_sai0n_g2_g2(3,3),PD_sai0t_g2_g2(3,3)
     &,PD_sai0n_g2_g1(3,3),PD_sai0t_g2_g1(3,3)
     &,delta_dnn_up,delta_dnt_up,delta_dtn_up,delta_dtt_up
     &,delta_dnn_low,delta_dnt_low,delta_dtn_low,delta_dtt_low