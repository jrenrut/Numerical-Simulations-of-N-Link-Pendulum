function b = npend_alphas_Lagrange_b_1(I1,L1,d1,g,m1,omega1,theta1)
%NPEND_ALPHAS_LAGRANGE_B_1
%    B = NPEND_ALPHAS_LAGRANGE_B_1(I1,L1,D1,G,M1,OMEGA1,THETA1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:32:27

b = [-d1.*conj(g).*conj(m1).*sin(theta1)];
