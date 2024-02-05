# -- bcp : boundary condition parameters 
bcp1=0; bcp2=0; bcp3=0;
# -- icp : initial  condition parameters 
icp1=0; icp2=0; icp3=0;
# -- vp  :      potential     parameters
vp1 =0; vp2 =0; vp3 =0; 
# -- time & space limits 
ti  = 0; tf  = 0;   # 
xti = 0; xtf = 0;   # xt := "x-tortoise"
# -- t & x* slices : Nt, Nxt 
Nt  = 0; Nxt = 0;   # Nxt := "Nx-tortoise"
Npl = 0;            # draw 2d plots of the solution if i%Npl = 0
# -- t \& xt steps 
Dt  = 0; Dxt = 0; 
#----------------------------------------
lam = 0;