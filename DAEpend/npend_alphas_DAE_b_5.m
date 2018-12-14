function b = npend_alphas_DAE_b_5(in1,in2,in3,g,in5,in6,in7,in8)
%NPEND_ALPHAS_DAE_B_5
%    B = NPEND_ALPHAS_DAE_B_5(IN1,IN2,IN3,G,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Dec-2018 02:30:17

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
d5 = in3(5,:);
d_a1 = in8(:,1);
d_a2 = in8(:,2);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
m4 = in5(4,:);
m5 = in5(5,:);
omega1 = in6(1,:);
omega2 = in6(2,:);
omega3 = in6(3,:);
omega4 = in6(4,:);
omega5 = in6(5,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
theta3 = in7(3,:);
theta4 = in7(4,:);
theta5 = in7(5,:);
t2 = omega1.^2;
t3 = sin(theta1);
t4 = cos(theta1);
t5 = d1.*t2.*t4;
t6 = omega2.^2;
t7 = sin(theta2);
t8 = cos(theta2);
t9 = d2.*t6.*t8;
t10 = omega3.^2;
t11 = sin(theta3);
t12 = cos(theta3);
t13 = d3.*t10.*t12;
t14 = omega4.^2;
t15 = sin(theta4);
t16 = cos(theta4);
t17 = d4.*t14.*t16;
t18 = omega5.^2;
b = [0.0;-g.*m1;0.0;-g.*m2;0.0;-g.*m3;0.0;-g.*m4;0.0;-g.*m5;0.0;0.0;0.0;0.0;0.0;d_a1-d1.*t2.*t3;d_a2+t5;-L1.*t2.*t3+d1.*t2.*t3-d2.*t6.*t7;-t5+t9+L1.*t2.*t4;-L2.*t6.*t7+d2.*t6.*t7-d3.*t10.*t11;-t9+t13+L2.*t6.*t8;-L3.*t10.*t11+d3.*t10.*t11-d4.*t14.*t15;-t13+t17+L3.*t10.*t12;-d5.*t18.*sin(theta5)-L4.*t14.*t15+d4.*t14.*t15;-t17+d5.*t18.*cos(theta5)+L4.*t14.*t16];
