function M = npend_alphas_DAE_M_1(I1,L1,d1,g,m1,omega1,theta1,in8)
%NPEND_ALPHAS_DAE_M_1
%    M = NPEND_ALPHAS_DAE_M_1(I1,L1,D1,G,M1,OMEGA1,THETA1,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Dec-2018 02:30:11

t2 = cos(theta1);
t3 = sin(theta1);
M = reshape([m1,0.0,0.0,1.0,0.0,0.0,m1,0.0,0.0,1.0,0.0,0.0,-I1,-d1.*t2,-d1.*t3,-1.0,0.0,-d1.*t2,0.0,0.0,0.0,-1.0,-d1.*t3,0.0,0.0],[5,5]);