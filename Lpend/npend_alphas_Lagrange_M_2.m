function M = npend_alphas_Lagrange_M_2(in1,in2,in3,g,in5,in6,in7)
%NPEND_ALPHAS_LAGRANGE_M_2
%    M = NPEND_ALPHAS_LAGRANGE_M_2(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:32:28

I1 = in1(1,:);
I2 = in1(2,:);
L1 = in2(1,:);
d1 = in3(1,:);
d2 = in3(2,:);
m1 = in5(1,:);
m2 = in5(2,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
t2 = cos(theta1);
t3 = conj(L1);
t4 = sin(theta1);
t5 = t2.^2;
t6 = conj(d1);
t7 = t4.^2;
t8 = cos(theta2);
t9 = conj(d2);
t10 = sin(theta2);
t11 = L1.*t2.*t8.*t9;
t12 = d2.*t2.*t3.*t8;
t13 = L1.*t4.*t9.*t10;
t14 = d2.*t3.*t4.*t10;
t15 = t11+t12+t13+t14;
t16 = m2.*t15.*(1.0./2.0);
M = reshape([I1+m2.*(L1.*t3.*t5.*2.0+L1.*t3.*t7.*2.0).*(1.0./2.0)+m1.*(d1.*t5.*t6.*2.0+d1.*t6.*t7.*2.0).*(1.0./2.0),t16,t16,I2+m2.*(d2.*t8.^2.*t9.*2.0+d2.*t9.*t10.^2.*2.0).*(1.0./2.0)],[2,2]);