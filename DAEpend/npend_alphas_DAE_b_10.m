function b = npend_alphas_DAE_b_10(in1,in2,in3,g,in5,in6,in7,in8)
%NPEND_ALPHAS_DAE_B_10
%    B = NPEND_ALPHAS_DAE_B_10(IN1,IN2,IN3,G,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Dec-2018 02:30:35

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
L6 = in2(6,:);
L7 = in2(7,:);
L8 = in2(8,:);
L9 = in2(9,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
d5 = in3(5,:);
d6 = in3(6,:);
d7 = in3(7,:);
d8 = in3(8,:);
d9 = in3(9,:);
d10 = in3(10,:);
d_a1 = in8(:,1);
d_a2 = in8(:,2);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
m4 = in5(4,:);
m5 = in5(5,:);
m6 = in5(6,:);
m7 = in5(7,:);
m8 = in5(8,:);
m9 = in5(9,:);
m10 = in5(10,:);
omega1 = in6(1,:);
omega2 = in6(2,:);
omega3 = in6(3,:);
omega4 = in6(4,:);
omega5 = in6(5,:);
omega6 = in6(6,:);
omega7 = in6(7,:);
omega8 = in6(8,:);
omega9 = in6(9,:);
omega10 = in6(10,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
theta3 = in7(3,:);
theta4 = in7(4,:);
theta5 = in7(5,:);
theta6 = in7(6,:);
theta7 = in7(7,:);
theta8 = in7(8,:);
theta9 = in7(9,:);
theta10 = in7(10,:);
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
t19 = sin(theta5);
t20 = cos(theta5);
t21 = d5.*t18.*t20;
t22 = omega6.^2;
t23 = sin(theta6);
t24 = cos(theta6);
t25 = d6.*t22.*t24;
t26 = omega7.^2;
t27 = sin(theta7);
t28 = cos(theta7);
t29 = d7.*t26.*t28;
t30 = omega8.^2;
t31 = sin(theta8);
t32 = cos(theta8);
t33 = d8.*t30.*t32;
t34 = omega9.^2;
t35 = sin(theta9);
t36 = cos(theta9);
t37 = d9.*t34.*t36;
t38 = omega10.^2;
b = [0.0;-g.*m1;0.0;-g.*m2;0.0;-g.*m3;0.0;-g.*m4;0.0;-g.*m5;0.0;-g.*m6;0.0;-g.*m7;0.0;-g.*m8;0.0;-g.*m9;0.0;-g.*m10;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;d_a1-d1.*t2.*t3;d_a2+t5;-L1.*t2.*t3+d1.*t2.*t3-d2.*t6.*t7;-t5+t9+L1.*t2.*t4;-L2.*t6.*t7+d2.*t6.*t7-d3.*t10.*t11;-t9+t13+L2.*t6.*t8;-L3.*t10.*t11+d3.*t10.*t11-d4.*t14.*t15;-t13+t17+L3.*t10.*t12;-L4.*t14.*t15+d4.*t14.*t15-d5.*t18.*t19;-t17+t21+L4.*t14.*t16;-L5.*t18.*t19+d5.*t18.*t19-d6.*t22.*t23;-t21+t25+L5.*t18.*t20;-L6.*t22.*t23+d6.*t22.*t23-d7.*t26.*t27;-t25+t29+L6.*t22.*t24;-L7.*t26.*t27+d7.*t26.*t27-d8.*t30.*t31;-t29+t33+L7.*t26.*t28;-L8.*t30.*t31+d8.*t30.*t31-d9.*t34.*t35;-t33+t37+L8.*t30.*t32;-d10.*t38.*sin(theta10)-L9.*t34.*t35+d9.*t34.*t35;-t37+d10.*t38.*cos(theta10)+L9.*t34.*t36];
