function K = nmodeK_3(in1,in2,in3,g,in5,in6,in7)
%NMODEK_3
%    K = NMODEK_3(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    11-Dec-2018 13:25:27

L1 = in2(1,:);
L2 = in2(2,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
omega1 = in6(1,:);
omega2 = in6(2,:);
omega3 = in6(3,:);
t2 = conj(L1);
t3 = omega1.*t2;
t4 = conj(L2);
t5 = omega2.*t4;
t6 = conj(d3);
t7 = omega3.*t6;
t8 = t3+t5+t7;
t9 = L1.*omega1;
t10 = L2.*omega2;
t11 = d3.*omega3;
t12 = t9+t10+t11;
t13 = omega1.^2;
t14 = d2.*omega2;
t15 = t9+t14;
t16 = conj(d2);
t17 = omega2.*t16;
t18 = t3+t17;
t19 = conj(m2);
t20 = conj(g);
t21 = conj(m3);
t22 = L1.*omega2.*t4;
t23 = L2.*omega2.*t2;
t24 = t22+t23;
t25 = m3.*t24.*(1.0./2.0);
t26 = L1.*omega2.*t16;
t27 = d2.*omega2.*t2;
t28 = t26+t27;
t29 = m2.*t28.*(1.0./2.0);
t30 = t25+t29;
t31 = omega2.*t30;
t32 = L1.*omega3.*t6;
t33 = d3.*omega3.*t2;
t34 = t32+t33;
t35 = m3.*omega3.*t34.*(1.0./2.0);
t36 = L1.*omega1.*omega2.*t4;
t37 = L2.*omega1.*omega2.*t2;
t38 = t36+t37;
t39 = L1.*omega1.*omega2.*t16;
t40 = d2.*omega1.*omega2.*t2;
t41 = t39+t40;
t42 = L1.*omega1.*t4;
t43 = L2.*omega1.*t2;
t44 = t42+t43;
t45 = m3.*t44.*(1.0./2.0);
t46 = L1.*omega1.*t16;
t47 = d2.*omega1.*t2;
t48 = t46+t47;
t49 = m2.*t48.*(1.0./2.0);
t50 = t45+t49;
t51 = omega2.^2;
t52 = L2.*omega3.*t6;
t53 = d3.*omega3.*t4;
t54 = t52+t53;
t55 = m3.*omega3.*t54.*(1.0./2.0);
t56 = L1.*omega1.*omega3.*t6;
t57 = d3.*omega1.*omega3.*t2;
t58 = t56+t57;
t59 = L1.*omega1.*t6;
t60 = d3.*omega1.*t2;
t61 = t59+t60;
t62 = L2.*omega2.*omega3.*t6;
t63 = d3.*omega2.*omega3.*t4;
t64 = t62+t63;
t65 = L2.*omega2.*t6;
t66 = d3.*omega2.*t4;
t67 = t65+t66;
K = reshape([t31+t35+t21.*(L1.*omega1.*t8-L1.*t2.*t13.*2.0+omega1.*t2.*t12).*(1.0./2.0)+t19.*(L1.*omega1.*t18-L1.*t2.*t13.*2.0+omega1.*t2.*t15).*(1.0./2.0)-omega1.*(m3.*(L1.*t8+t2.*t12-L1.*omega1.*t2.*2.0).*(1.0./2.0)+m2.*(L1.*t18+t2.*t15-L1.*omega1.*t2.*2.0).*(1.0./2.0))+d1.*t20.*conj(m1)+L1.*t19.*t20+L1.*t20.*t21,-omega1.*t50+omega2.*t50-t21.*t38.*(1.0./2.0)-t19.*t41.*(1.0./2.0),t21.*t58.*(-1.0./2.0)-m3.*omega1.*t61.*(1.0./2.0)+m3.*omega3.*t61.*(1.0./2.0),-t31+omega1.*t30-t21.*t38.*(1.0./2.0)-t19.*t41.*(1.0./2.0),t55+t21.*(L2.*omega2.*t8-L2.*t4.*t51.*2.0+omega2.*t4.*t12).*(1.0./2.0)+t19.*(d2.*omega2.*t18-d2.*t16.*t51.*2.0+omega2.*t15.*t16).*(1.0./2.0)+omega1.*t50-omega2.*(m3.*(L2.*t8+t4.*t12-L2.*omega2.*t4.*2.0).*(1.0./2.0)+m2.*(d2.*t18+t15.*t16-d2.*omega2.*t16.*2.0).*(1.0./2.0))+L2.*t20.*t21+d2.*t19.*t20,t21.*t64.*(-1.0./2.0)-m3.*omega2.*t67.*(1.0./2.0)+m3.*omega3.*t67.*(1.0./2.0),-t35-t21.*t58.*(1.0./2.0)+m3.*omega1.*t34.*(1.0./2.0),-t55-t21.*t64.*(1.0./2.0)+m3.*omega2.*t54.*(1.0./2.0),t21.*(d3.*omega3.^2.*t6.*-2.0+d3.*omega3.*t8+omega3.*t6.*t12).*(1.0./2.0)+d3.*t20.*t21+m3.*omega1.*t61.*(1.0./2.0)+m3.*omega2.*t67.*(1.0./2.0)-m3.*omega3.*(d3.*t8+t6.*t12-d3.*omega3.*t6.*2.0).*(1.0./2.0)],[3,3]);
