function b = npend_alphas_Lagrange_b_5(in1,in2,in3,g,in5,in6,in7)
%NPEND_ALPHAS_LAGRANGE_B_5
%    B = NPEND_ALPHAS_LAGRANGE_B_5(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:32:40

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
d5 = in3(5,:);
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
t2 = cos(theta1);
t3 = conj(L1);
t4 = sin(theta1);
t5 = conj(L2);
t6 = conj(d3);
t7 = cos(theta3);
t8 = cos(theta2);
t9 = sin(theta2);
t10 = sin(theta3);
t11 = conj(L4);
t12 = sin(theta4);
t13 = cos(theta4);
t14 = conj(d4);
t15 = omega1.*t2.*t3;
t16 = omega1.*t3.*t4;
t17 = conj(d2);
t18 = L1.*omega1.*t2;
t19 = L1.*omega1.*t4;
t20 = conj(L3);
t21 = L1.*omega3.*t4.*t7.*t20;
t22 = L3.*omega3.*t3.*t4.*t7;
t23 = t21+t22-L3.*omega3.*t2.*t3.*t10-L1.*omega3.*t2.*t10.*t20;
t24 = omega2.*t5.*t8;
t25 = omega2.*t5.*t9;
t26 = conj(d5);
t27 = cos(theta5);
t28 = L2.*omega2.*t8;
t29 = L2.*omega2.*t9;
t30 = sin(theta5);
t31 = omega3.*t6.*t7;
t32 = t15+t24+t31;
t33 = omega3.*t6.*t10;
t34 = t16+t25+t33;
t35 = d3.*omega3.*t7;
t36 = t18+t28+t35;
t37 = d3.*omega3.*t10;
t38 = t19+t29+t37;
t39 = d2.*omega2.*t8;
t40 = t18+t39;
t41 = d2.*omega2.*t9;
t42 = t19+t41;
t43 = omega2.*t8.*t17;
t44 = t15+t43;
t45 = omega2.*t9.*t17;
t46 = t16+t45;
t47 = d5.*omega5.*t27;
t48 = L3.*omega3.*t7;
t49 = L4.*omega4.*t13;
t50 = t18+t28+t47+t48+t49;
t51 = L3.*omega3.*t10;
t52 = L4.*omega4.*t12;
t53 = d5.*omega5.*t30;
t54 = t19+t29+t51+t52+t53;
t55 = omega3.*t7.*t20;
t56 = omega4.*t11.*t13;
t57 = omega5.*t26.*t27;
t58 = t15+t24+t55+t56+t57;
t59 = omega3.*t10.*t20;
t60 = omega4.*t11.*t12;
t61 = omega5.*t26.*t30;
t62 = t16+t25+t59+t60+t61;
t63 = L1.*omega2.*t2.*t5.*t9;
t64 = L2.*omega2.*t2.*t3.*t9;
t66 = L1.*omega2.*t4.*t5.*t8;
t67 = L2.*omega2.*t3.*t4.*t8;
t65 = t63+t64-t66-t67;
t68 = omega4.*t12.*t14;
t69 = t16+t25+t59+t68;
t70 = d4.*omega4.*t13;
t71 = t18+t28+t48+t70;
t72 = omega4.*t13.*t14;
t73 = t15+t24+t55+t72;
t74 = d4.*omega4.*t12;
t75 = t19+t29+t51+t74;
t76 = conj(m2);
t77 = conj(g);
t78 = conj(m3);
t79 = conj(m4);
t80 = conj(m5);
t81 = L2.*omega3.*t7.*t9.*t20;
t82 = L3.*omega3.*t5.*t7.*t9;
t83 = t81+t82-L3.*omega3.*t5.*t8.*t10-L2.*omega3.*t8.*t10.*t20;
t84 = L1.*omega1.*t2.*t5.*t9;
t85 = L2.*omega1.*t2.*t3.*t9;
t87 = L1.*omega1.*t4.*t5.*t8;
t88 = L2.*omega1.*t3.*t4.*t8;
t86 = t84+t85-t87-t88;
t89 = L1.*omega1.*t4.*t7.*t20;
t90 = L3.*omega1.*t3.*t4.*t7;
t91 = t89+t90-L3.*omega1.*t2.*t3.*t10-L1.*omega1.*t2.*t10.*t20;
t92 = L2.*omega2.*t7.*t9.*t20;
t93 = L3.*omega2.*t5.*t7.*t9;
t94 = t92+t93-L3.*omega2.*t5.*t8.*t10-L2.*omega2.*t8.*t10.*t20;
b = [omega2.*(m2.*(L1.*omega2.*t2.*t9.*t17-L1.*omega2.*t4.*t8.*t17+d2.*omega2.*t2.*t3.*t9-d2.*omega2.*t3.*t4.*t8).*(1.0./2.0)+m3.*t65.*(1.0./2.0)+m4.*t65.*(1.0./2.0)+m5.*t65.*(1.0./2.0))-omega1.*(m3.*(L1.*t2.*t34-L1.*t4.*t32+t2.*t3.*t38-t3.*t4.*t36).*(1.0./2.0)+m2.*(L1.*t2.*t46-L1.*t4.*t44+t2.*t3.*t42-t3.*t4.*t40).*(1.0./2.0)-m5.*(L1.*t4.*t58-L1.*t2.*t62+t3.*t4.*t50-t2.*t3.*t54).*(1.0./2.0)+m4.*(L1.*t2.*t69-L1.*t4.*t73-t3.*t4.*t71+t2.*t3.*t75).*(1.0./2.0))+omega4.*(m5.*(L4.*omega4.*t2.*t3.*t12-L4.*omega4.*t3.*t4.*t13+L1.*omega4.*t2.*t11.*t12-L1.*omega4.*t4.*t11.*t13).*(1.0./2.0)+m4.*(L1.*omega4.*t2.*t12.*t14-L1.*omega4.*t4.*t13.*t14+d4.*omega4.*t2.*t3.*t12-d4.*omega4.*t3.*t4.*t13).*(1.0./2.0))+t78.*(L1.*omega1.*t2.*t34-L1.*omega1.*t4.*t32+omega1.*t2.*t3.*t38-omega1.*t3.*t4.*t36).*(1.0./2.0)+t76.*(L1.*omega1.*t2.*t46-L1.*omega1.*t4.*t44+omega1.*t2.*t3.*t42-omega1.*t3.*t4.*t40).*(1.0./2.0)-t80.*(L1.*omega1.*t4.*t58-L1.*omega1.*t2.*t62+omega1.*t3.*t4.*t50-omega1.*t2.*t3.*t54).*(1.0./2.0)+t79.*(L1.*omega1.*t2.*t69-L1.*omega1.*t4.*t73-omega1.*t3.*t4.*t71+omega1.*t2.*t3.*t75).*(1.0./2.0)-omega3.*(m3.*(L1.*omega3.*t4.*t6.*t7-L1.*omega3.*t2.*t6.*t10+d3.*omega3.*t3.*t4.*t7-d3.*omega3.*t2.*t3.*t10).*(1.0./2.0)+m4.*t23.*(1.0./2.0)+m5.*t23.*(1.0./2.0))-m5.*omega5.*(L1.*omega5.*t4.*t26.*t27-L1.*omega5.*t2.*t26.*t30+d5.*omega5.*t3.*t4.*t27-d5.*omega5.*t2.*t3.*t30).*(1.0./2.0)-L1.*t4.*t76.*t77-L1.*t4.*t77.*t78-L1.*t4.*t77.*t79-L1.*t4.*t77.*t80-d1.*t4.*t77.*conj(m1);-omega1.*(m2.*(L1.*omega1.*t2.*t9.*t17-L1.*omega1.*t4.*t8.*t17+d2.*omega1.*t2.*t3.*t9-d2.*omega1.*t3.*t4.*t8).*(1.0./2.0)+m3.*t86.*(1.0./2.0)+m4.*t86.*(1.0./2.0)+m5.*t86.*(1.0./2.0))+omega4.*(m5.*(L4.*omega4.*t5.*t8.*t12-L4.*omega4.*t5.*t9.*t13+L2.*omega4.*t8.*t11.*t12-L2.*omega4.*t9.*t11.*t13).*(1.0./2.0)+m4.*(L2.*omega4.*t8.*t12.*t14-L2.*omega4.*t9.*t13.*t14+d4.*omega4.*t5.*t8.*t12-d4.*omega4.*t5.*t9.*t13).*(1.0./2.0))+omega2.*(m3.*(L2.*t9.*t32-L2.*t8.*t34+t5.*t9.*t36-t5.*t8.*t38).*(1.0./2.0)+m5.*(L2.*t9.*t58-L2.*t8.*t62+t5.*t9.*t50-t5.*t8.*t54).*(1.0./2.0)-m4.*(L2.*t8.*t69-L2.*t9.*t73-t5.*t9.*t71+t5.*t8.*t75).*(1.0./2.0)+m2.*(d2.*t9.*t44-d2.*t8.*t46+t9.*t17.*t40-t8.*t17.*t42).*(1.0./2.0))-t78.*(L2.*omega2.*t9.*t32-L2.*omega2.*t8.*t34+omega2.*t5.*t9.*t36-omega2.*t5.*t8.*t38).*(1.0./2.0)-t80.*(L2.*omega2.*t9.*t58-L2.*omega2.*t8.*t62+omega2.*t5.*t9.*t50-omega2.*t5.*t8.*t54).*(1.0./2.0)+t79.*(L2.*omega2.*t8.*t69-L2.*omega2.*t9.*t73-omega2.*t5.*t9.*t71+omega2.*t5.*t8.*t75).*(1.0./2.0)-t76.*(d2.*omega2.*t9.*t44-d2.*omega2.*t8.*t46+omega2.*t9.*t17.*t40-omega2.*t8.*t17.*t42).*(1.0./2.0)-omega3.*(m3.*(L2.*omega3.*t6.*t7.*t9-L2.*omega3.*t6.*t8.*t10+d3.*omega3.*t5.*t7.*t9-d3.*omega3.*t5.*t8.*t10).*(1.0./2.0)+m4.*t83.*(1.0./2.0)+m5.*t83.*(1.0./2.0))-m5.*omega5.*(L2.*omega5.*t9.*t26.*t27-L2.*omega5.*t8.*t26.*t30+d5.*omega5.*t5.*t9.*t27-d5.*omega5.*t5.*t8.*t30).*(1.0./2.0)-L2.*t9.*t77.*t78-L2.*t9.*t77.*t79-L2.*t9.*t77.*t80-d2.*t9.*t76.*t77;omega4.*(m5.*(L3.*omega4.*t7.*t11.*t12-L3.*omega4.*t10.*t11.*t13+L4.*omega4.*t7.*t12.*t20-L4.*omega4.*t10.*t13.*t20).*(1.0./2.0)+m4.*(L3.*omega4.*t7.*t12.*t14-L3.*omega4.*t10.*t13.*t14+d4.*omega4.*t7.*t12.*t20-d4.*omega4.*t10.*t13.*t20).*(1.0./2.0))-t80.*(L3.*omega3.*t10.*t58-L3.*omega3.*t7.*t62+omega3.*t10.*t20.*t50-omega3.*t7.*t20.*t54).*(1.0./2.0)+t79.*(L3.*omega3.*t7.*t69-L3.*omega3.*t10.*t73-omega3.*t10.*t20.*t71+omega3.*t7.*t20.*t75).*(1.0./2.0)-omega3.*(m5.*(L3.*t10.*t58-L3.*t7.*t62+t10.*t20.*t50-t7.*t20.*t54).*(-1.0./2.0)+m4.*(L3.*t7.*t69-L3.*t10.*t73-t10.*t20.*t71+t7.*t20.*t75).*(1.0./2.0)+m3.*(d3.*t7.*t34-d3.*t10.*t32+t6.*t7.*t38-t6.*t10.*t36).*(1.0./2.0))+t78.*(d3.*omega3.*t7.*t34-d3.*omega3.*t10.*t32+omega3.*t6.*t7.*t38-omega3.*t6.*t10.*t36).*(1.0./2.0)+omega1.*(m3.*(L1.*omega1.*t4.*t6.*t7-L1.*omega1.*t2.*t6.*t10+d3.*omega1.*t3.*t4.*t7-d3.*omega1.*t2.*t3.*t10).*(1.0./2.0)+m4.*t91.*(1.0./2.0)+m5.*t91.*(1.0./2.0))+omega2.*(m3.*(L2.*omega2.*t6.*t7.*t9-L2.*omega2.*t6.*t8.*t10+d3.*omega2.*t5.*t7.*t9-d3.*omega2.*t5.*t8.*t10).*(1.0./2.0)+m4.*t94.*(1.0./2.0)+m5.*t94.*(1.0./2.0))+m5.*omega5.*(L3.*omega5.*t7.*t26.*t30-L3.*omega5.*t10.*t26.*t27+d5.*omega5.*t7.*t20.*t30-d5.*omega5.*t10.*t20.*t27).*(1.0./2.0)-L3.*t10.*t77.*t79-L3.*t10.*t77.*t80-d3.*t10.*t77.*t78;-omega1.*(m5.*(L4.*omega1.*t2.*t3.*t12-L4.*omega1.*t3.*t4.*t13+L1.*omega1.*t2.*t11.*t12-L1.*omega1.*t4.*t11.*t13).*(1.0./2.0)+m4.*(L1.*omega1.*t2.*t12.*t14-L1.*omega1.*t4.*t13.*t14+d4.*omega1.*t2.*t3.*t12-d4.*omega1.*t3.*t4.*t13).*(1.0./2.0))-omega2.*(m5.*(L4.*omega2.*t5.*t8.*t12-L4.*omega2.*t5.*t9.*t13+L2.*omega2.*t8.*t11.*t12-L2.*omega2.*t9.*t11.*t13).*(1.0./2.0)+m4.*(L2.*omega2.*t8.*t12.*t14-L2.*omega2.*t9.*t13.*t14+d4.*omega2.*t5.*t8.*t12-d4.*omega2.*t5.*t9.*t13).*(1.0./2.0))-omega3.*(m5.*(L3.*omega3.*t7.*t11.*t12-L3.*omega3.*t10.*t11.*t13+L4.*omega3.*t7.*t12.*t20-L4.*omega3.*t10.*t13.*t20).*(1.0./2.0)+m4.*(L3.*omega3.*t7.*t12.*t14-L3.*omega3.*t10.*t13.*t14+d4.*omega3.*t7.*t12.*t20-d4.*omega3.*t10.*t13.*t20).*(1.0./2.0))-t80.*(L4.*omega4.*t12.*t58-L4.*omega4.*t13.*t62+omega4.*t11.*t12.*t50-omega4.*t11.*t13.*t54).*(1.0./2.0)+t79.*(d4.*omega4.*t13.*t69-d4.*omega4.*t12.*t73-omega4.*t12.*t14.*t71+omega4.*t13.*t14.*t75).*(1.0./2.0)+omega4.*(m5.*(L4.*t12.*t58-L4.*t13.*t62+t11.*t12.*t50-t11.*t13.*t54).*(1.0./2.0)-m4.*(d4.*t13.*t69-d4.*t12.*t73-t12.*t14.*t71+t13.*t14.*t75).*(1.0./2.0))-m5.*omega5.*(L4.*omega5.*t12.*t26.*t27-L4.*omega5.*t13.*t26.*t30+d5.*omega5.*t11.*t12.*t27-d5.*omega5.*t11.*t13.*t30).*(1.0./2.0)-L4.*t12.*t77.*t80-d4.*t12.*t77.*t79;t80.*(d5.*omega5.*t30.*t58-d5.*omega5.*t27.*t62+omega5.*t26.*t30.*t50-omega5.*t26.*t27.*t54).*(-1.0./2.0)+m5.*omega5.*(d5.*t30.*t58-d5.*t27.*t62+t26.*t30.*t50-t26.*t27.*t54).*(1.0./2.0)+m5.*omega1.*(L1.*omega1.*t4.*t26.*t27-L1.*omega1.*t2.*t26.*t30+d5.*omega1.*t3.*t4.*t27-d5.*omega1.*t2.*t3.*t30).*(1.0./2.0)+m5.*omega2.*(L2.*omega2.*t9.*t26.*t27-L2.*omega2.*t8.*t26.*t30+d5.*omega2.*t5.*t9.*t27-d5.*omega2.*t5.*t8.*t30).*(1.0./2.0)-m5.*omega3.*(L3.*omega3.*t7.*t26.*t30-L3.*omega3.*t10.*t26.*t27+d5.*omega3.*t7.*t20.*t30-d5.*omega3.*t10.*t20.*t27).*(1.0./2.0)+m5.*omega4.*(L4.*omega4.*t12.*t26.*t27-L4.*omega4.*t13.*t26.*t30+d5.*omega4.*t11.*t12.*t27-d5.*omega4.*t11.*t13.*t30).*(1.0./2.0)-d5.*t30.*t77.*t80];
