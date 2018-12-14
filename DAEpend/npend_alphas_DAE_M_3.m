function M = npend_alphas_DAE_M_3(in1,in2,in3,g,in5,in6,in7,in8)
%NPEND_ALPHAS_DAE_M_3
%    M = NPEND_ALPHAS_DAE_M_3(IN1,IN2,IN3,G,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Dec-2018 02:30:13

I1 = in1(1,:);
I2 = in1(2,:);
I3 = in1(3,:);
L1 = in2(1,:);
L2 = in2(2,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
theta3 = in7(3,:);
t2 = cos(theta1);
t3 = sin(theta1);
t4 = L1-d1;
t5 = cos(theta2);
t6 = sin(theta2);
t7 = L2-d2;
t8 = cos(theta3);
t9 = sin(theta3);
M = reshape([m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,-I1,0.0,0.0,-d1.*t2,-d1.*t3,-L1.*t2+d1.*t2,-L1.*t3+d1.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I2,0.0,0.0,0.0,-d2.*t5,-d2.*t6,-L2.*t5+d2.*t5,-L2.*t6+d2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I3,0.0,0.0,0.0,0.0,-d3.*t8,-d3.*t9,-1.0,0.0,0.0,0.0,0.0,0.0,-d1.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,-t2.*t4,-d2.*t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,-t5.*t7,-d3.*t8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,-d1.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,-t3.*t4,-d2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,-t6.*t7,-d3.*t9,0.0,0.0,0.0,0.0,0.0,0.0],[15,15]);
