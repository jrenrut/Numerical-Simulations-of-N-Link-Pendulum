function M = nmodeM_4(in1,in2,in3,g,in5,in6,in7)
%NMODEM_4
%    M = NMODEM_4(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    11-Dec-2018 19:33:44

I1 = in1(1,:);
I2 = in1(2,:);
I3 = in1(3,:);
I4 = in1(4,:);
L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
m4 = in5(4,:);
t2 = conj(L1);
t3 = conj(L2);
t4 = L1.*t3;
t5 = L2.*t2;
t6 = t4+t5;
t7 = conj(d2);
t8 = L1.*t7;
t9 = d2.*t2;
t10 = t8+t9;
t11 = m2.*t10.*(1.0./2.0);
t12 = m3.*t6.*(1.0./2.0);
t13 = m4.*t6.*(1.0./2.0);
t14 = t11+t12+t13;
t15 = conj(d3);
t16 = conj(L3);
t17 = conj(d4);
t18 = L1.*t15;
t19 = d3.*t2;
t20 = t18+t19;
t21 = m3.*t20.*(1.0./2.0);
t22 = L1.*t16;
t23 = L3.*t2;
t24 = t22+t23;
t25 = m4.*t24.*(1.0./2.0);
t26 = t21+t25;
t27 = L2.*t15;
t28 = d3.*t3;
t29 = t27+t28;
t30 = m3.*t29.*(1.0./2.0);
t31 = L2.*t16;
t32 = L3.*t3;
t33 = t31+t32;
t34 = m4.*t33.*(1.0./2.0);
t35 = t30+t34;
t36 = L1.*t17;
t37 = d4.*t2;
t38 = t36+t37;
t39 = m4.*t38.*(1.0./2.0);
t40 = L2.*t17;
t41 = d4.*t3;
t42 = t40+t41;
t43 = m4.*t42.*(1.0./2.0);
t44 = L3.*t17;
t45 = d4.*t16;
t46 = t44+t45;
t47 = m4.*t46.*(1.0./2.0);
M = reshape([I1+L1.*m2.*t2+L1.*m3.*t2+L1.*m4.*t2+d1.*m1.*conj(d1),t14,t26,t39,t14,I2+L2.*m3.*t3+L2.*m4.*t3+d2.*m2.*t7,t35,t43,t26,t35,I3+L3.*m4.*t16+d3.*m3.*t15,t47,t39,t43,t47,I4+d4.*m4.*t17],[4,4]);
