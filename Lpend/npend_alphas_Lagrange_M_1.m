function M = npend_alphas_Lagrange_M_1(I1,L1,d1,g,m1,omega1,theta1)
%NPEND_ALPHAS_LAGRANGE_M_1
%    M = NPEND_ALPHAS_LAGRANGE_M_1(I1,L1,D1,G,M1,OMEGA1,THETA1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:32:27

t2 = cos(theta1);
t3 = conj(d1);
t4 = sin(theta1);
M = [I1+m1.*(d1.*t2.^2.*t3.*2.0+d1.*t3.*t4.^2.*2.0).*(1.0./2.0)];