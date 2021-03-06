function b = npend_alphas_DAE_b_30(in1,in2,in3,g,in5,in6,in7,in8)
%NPEND_ALPHAS_DAE_B_30
%    B = NPEND_ALPHAS_DAE_B_30(IN1,IN2,IN3,G,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Dec-2018 02:33:26

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
L6 = in2(6,:);
L7 = in2(7,:);
L8 = in2(8,:);
L9 = in2(9,:);
L10 = in2(10,:);
L11 = in2(11,:);
L12 = in2(12,:);
L13 = in2(13,:);
L14 = in2(14,:);
L15 = in2(15,:);
L16 = in2(16,:);
L17 = in2(17,:);
L18 = in2(18,:);
L19 = in2(19,:);
L20 = in2(20,:);
L21 = in2(21,:);
L22 = in2(22,:);
L23 = in2(23,:);
L24 = in2(24,:);
L25 = in2(25,:);
L26 = in2(26,:);
L27 = in2(27,:);
L28 = in2(28,:);
L29 = in2(29,:);
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
d11 = in3(11,:);
d12 = in3(12,:);
d13 = in3(13,:);
d14 = in3(14,:);
d15 = in3(15,:);
d16 = in3(16,:);
d17 = in3(17,:);
d18 = in3(18,:);
d19 = in3(19,:);
d20 = in3(20,:);
d21 = in3(21,:);
d22 = in3(22,:);
d23 = in3(23,:);
d24 = in3(24,:);
d25 = in3(25,:);
d26 = in3(26,:);
d27 = in3(27,:);
d28 = in3(28,:);
d29 = in3(29,:);
d30 = in3(30,:);
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
m11 = in5(11,:);
m12 = in5(12,:);
m13 = in5(13,:);
m14 = in5(14,:);
m15 = in5(15,:);
m16 = in5(16,:);
m17 = in5(17,:);
m18 = in5(18,:);
m19 = in5(19,:);
m20 = in5(20,:);
m21 = in5(21,:);
m22 = in5(22,:);
m23 = in5(23,:);
m24 = in5(24,:);
m25 = in5(25,:);
m26 = in5(26,:);
m27 = in5(27,:);
m28 = in5(28,:);
m29 = in5(29,:);
m30 = in5(30,:);
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
omega11 = in6(11,:);
omega12 = in6(12,:);
omega13 = in6(13,:);
omega14 = in6(14,:);
omega15 = in6(15,:);
omega16 = in6(16,:);
omega17 = in6(17,:);
omega18 = in6(18,:);
omega19 = in6(19,:);
omega20 = in6(20,:);
omega21 = in6(21,:);
omega22 = in6(22,:);
omega23 = in6(23,:);
omega24 = in6(24,:);
omega25 = in6(25,:);
omega26 = in6(26,:);
omega27 = in6(27,:);
omega28 = in6(28,:);
omega29 = in6(29,:);
omega30 = in6(30,:);
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
theta11 = in7(11,:);
theta12 = in7(12,:);
theta13 = in7(13,:);
theta14 = in7(14,:);
theta15 = in7(15,:);
theta16 = in7(16,:);
theta17 = in7(17,:);
theta18 = in7(18,:);
theta19 = in7(19,:);
theta20 = in7(20,:);
theta21 = in7(21,:);
theta22 = in7(22,:);
theta23 = in7(23,:);
theta24 = in7(24,:);
theta25 = in7(25,:);
theta26 = in7(26,:);
theta27 = in7(27,:);
theta28 = in7(28,:);
theta29 = in7(29,:);
theta30 = in7(30,:);
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
t39 = sin(theta10);
t40 = cos(theta10);
t41 = d10.*t38.*t40;
t42 = omega11.^2;
t43 = sin(theta11);
t44 = cos(theta11);
t45 = d11.*t42.*t44;
t46 = omega12.^2;
t47 = sin(theta12);
t48 = cos(theta12);
t49 = d12.*t46.*t48;
t50 = omega13.^2;
t51 = sin(theta13);
t52 = cos(theta13);
t53 = d13.*t50.*t52;
t54 = omega14.^2;
t55 = sin(theta14);
t56 = cos(theta14);
t57 = d14.*t54.*t56;
t58 = omega15.^2;
t59 = sin(theta15);
t60 = cos(theta15);
t61 = d15.*t58.*t60;
t62 = omega16.^2;
t63 = sin(theta16);
t64 = cos(theta16);
t65 = d16.*t62.*t64;
t66 = omega17.^2;
t67 = sin(theta17);
t68 = cos(theta17);
t69 = d17.*t66.*t68;
t70 = omega18.^2;
t71 = sin(theta18);
t72 = cos(theta18);
t73 = d18.*t70.*t72;
t74 = omega19.^2;
t75 = sin(theta19);
t76 = cos(theta19);
t77 = d19.*t74.*t76;
t78 = omega20.^2;
t79 = sin(theta20);
t80 = cos(theta20);
t81 = d20.*t78.*t80;
t82 = omega21.^2;
t83 = sin(theta21);
t84 = cos(theta21);
t85 = d21.*t82.*t84;
t86 = omega22.^2;
t87 = sin(theta22);
t88 = cos(theta22);
t89 = d22.*t86.*t88;
t90 = omega23.^2;
t91 = sin(theta23);
t92 = cos(theta23);
t93 = d23.*t90.*t92;
t94 = omega24.^2;
t95 = sin(theta24);
t96 = cos(theta24);
t97 = d24.*t94.*t96;
t98 = omega25.^2;
t99 = sin(theta25);
t100 = cos(theta25);
t101 = d25.*t98.*t100;
t102 = omega26.^2;
t103 = sin(theta26);
t104 = cos(theta26);
t105 = d26.*t102.*t104;
t106 = omega27.^2;
t107 = sin(theta27);
t108 = cos(theta27);
t109 = d27.*t106.*t108;
t110 = omega28.^2;
t111 = sin(theta28);
t112 = cos(theta28);
t113 = d28.*t110.*t112;
t114 = omega29.^2;
t115 = sin(theta29);
t116 = cos(theta29);
t117 = d29.*t114.*t116;
t118 = omega30.^2;
b = [0.0;-g.*m1;0.0;-g.*m2;0.0;-g.*m3;0.0;-g.*m4;0.0;-g.*m5;0.0;-g.*m6;0.0;-g.*m7;0.0;-g.*m8;0.0;-g.*m9;0.0;-g.*m10;0.0;-g.*m11;0.0;-g.*m12;0.0;-g.*m13;0.0;-g.*m14;0.0;-g.*m15;0.0;-g.*m16;0.0;-g.*m17;0.0;-g.*m18;0.0;-g.*m19;0.0;-g.*m20;0.0;-g.*m21;0.0;-g.*m22;0.0;-g.*m23;0.0;-g.*m24;0.0;-g.*m25;0.0;-g.*m26;0.0;-g.*m27;0.0;-g.*m28;0.0;-g.*m29;0.0;-g.*m30;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;d_a1-d1.*t2.*t3;d_a2+t5;-L1.*t2.*t3+d1.*t2.*t3-d2.*t6.*t7;-t5+t9+L1.*t2.*t4;-L2.*t6.*t7+d2.*t6.*t7-d3.*t10.*t11;-t9+t13+L2.*t6.*t8;-L3.*t10.*t11+d3.*t10.*t11-d4.*t14.*t15;-t13+t17+L3.*t10.*t12;-L4.*t14.*t15+d4.*t14.*t15-d5.*t18.*t19;-t17+t21+L4.*t14.*t16;-L5.*t18.*t19+d5.*t18.*t19-d6.*t22.*t23;-t21+t25+L5.*t18.*t20;-L6.*t22.*t23+d6.*t22.*t23-d7.*t26.*t27;-t25+t29+L6.*t22.*t24;-L7.*t26.*t27+d7.*t26.*t27-d8.*t30.*t31;-t29+t33+L7.*t26.*t28;-L8.*t30.*t31+d8.*t30.*t31-d9.*t34.*t35;-t33+t37+L8.*t30.*t32;-L9.*t34.*t35+d9.*t34.*t35-d10.*t38.*t39;-t37+t41+L9.*t34.*t36;-L10.*t38.*t39+d10.*t38.*t39-d11.*t42.*t43;-t41+t45+L10.*t38.*t40;-L11.*t42.*t43+d11.*t42.*t43-d12.*t46.*t47;-t45+t49+L11.*t42.*t44;-L12.*t46.*t47+d12.*t46.*t47-d13.*t50.*t51;-t49+t53+L12.*t46.*t48;-L13.*t50.*t51+d13.*t50.*t51-d14.*t54.*t55;-t53+t57+L13.*t50.*t52;-L14.*t54.*t55+d14.*t54.*t55-d15.*t58.*t59;-t57+t61+L14.*t54.*t56;-L15.*t58.*t59+d15.*t58.*t59-d16.*t62.*t63;-t61+t65+L15.*t58.*t60;-L16.*t62.*t63+d16.*t62.*t63-d17.*t66.*t67;-t65+t69+L16.*t62.*t64;-L17.*t66.*t67+d17.*t66.*t67-d18.*t70.*t71;-t69+t73+L17.*t66.*t68;-L18.*t70.*t71+d18.*t70.*t71-d19.*t74.*t75;-t73+t77+L18.*t70.*t72;-L19.*t74.*t75+d19.*t74.*t75-d20.*t78.*t79;-t77+t81+L19.*t74.*t76;-L20.*t78.*t79+d20.*t78.*t79-d21.*t82.*t83;-t81+t85+L20.*t78.*t80;-L21.*t82.*t83+d21.*t82.*t83-d22.*t86.*t87;-t85+t89+L21.*t82.*t84;-L22.*t86.*t87+d22.*t86.*t87-d23.*t90.*t91;-t89+t93+L22.*t86.*t88;-L23.*t90.*t91+d23.*t90.*t91-d24.*t94.*t95;-t93+t97+L23.*t90.*t92;-L24.*t94.*t95+d24.*t94.*t95-d25.*t98.*t99;-t97+t101+L24.*t94.*t96;-L25.*t98.*t99+d25.*t98.*t99-d26.*t102.*t103;-t101+t105+L25.*t98.*t100;-L26.*t102.*t103+d26.*t102.*t103-d27.*t106.*t107;-t105+t109+L26.*t102.*t104;-L27.*t106.*t107+d27.*t106.*t107-d28.*t110.*t111;-t109+t113+L27.*t106.*t108;-L28.*t110.*t111+d28.*t110.*t111-d29.*t114.*t115;-t113+t117+L28.*t110.*t112;-d30.*t118.*sin(theta30)-L29.*t114.*t115+d29.*t114.*t115;-t117+d30.*t118.*cos(theta30)+L29.*t114.*t116];
