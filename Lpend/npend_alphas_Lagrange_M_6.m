function M = npend_alphas_Lagrange_M_6(in1,in2,in3,g,in5,in6,in7)
%NPEND_ALPHAS_LAGRANGE_M_6
%    M = NPEND_ALPHAS_LAGRANGE_M_6(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:32:44

I1 = in1(1,:);
I2 = in1(2,:);
I3 = in1(3,:);
I4 = in1(4,:);
I5 = in1(5,:);
I6 = in1(6,:);
L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
d5 = in3(5,:);
d6 = in3(6,:);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
m4 = in5(4,:);
m5 = in5(5,:);
m6 = in5(6,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
theta3 = in7(3,:);
theta4 = in7(4,:);
theta5 = in7(5,:);
theta6 = in7(6,:);
t2 = cos(theta1);
t3 = conj(L1);
t4 = sin(theta1);
t5 = t2.^2;
t6 = L1.*t3.*t5.*2.0;
t7 = t4.^2;
t8 = L1.*t3.*t7.*2.0;
t9 = t6+t8;
t10 = conj(d1);
t11 = cos(theta2);
t12 = conj(L2);
t13 = sin(theta2);
t14 = L1.*t2.*t11.*t12;
t15 = L2.*t2.*t3.*t11;
t16 = L1.*t4.*t12.*t13;
t17 = L2.*t3.*t4.*t13;
t18 = t14+t15+t16+t17;
t19 = conj(d2);
t20 = cos(theta3);
t21 = conj(L3);
t22 = sin(theta3);
t23 = L1.*t2.*t20.*t21;
t24 = L3.*t2.*t3.*t20;
t25 = L1.*t4.*t21.*t22;
t26 = L3.*t3.*t4.*t22;
t27 = t23+t24+t25+t26;
t28 = conj(d3);
t29 = cos(theta4);
t30 = conj(L4);
t31 = sin(theta4);
t32 = L1.*t2.*t29.*t30;
t33 = L4.*t2.*t3.*t29;
t34 = L1.*t4.*t30.*t31;
t35 = L4.*t3.*t4.*t31;
t36 = t32+t33+t34+t35;
t37 = conj(d4);
t38 = cos(theta5);
t39 = conj(L5);
t40 = sin(theta5);
t41 = conj(d5);
t42 = cos(theta6);
t43 = conj(d6);
t44 = sin(theta6);
t45 = m3.*t18.*(1.0./2.0);
t46 = m4.*t18.*(1.0./2.0);
t47 = m5.*t18.*(1.0./2.0);
t48 = m6.*t18.*(1.0./2.0);
t49 = L1.*t2.*t11.*t19;
t50 = d2.*t2.*t3.*t11;
t51 = L1.*t4.*t13.*t19;
t52 = d2.*t3.*t4.*t13;
t53 = t49+t50+t51+t52;
t54 = m2.*t53.*(1.0./2.0);
t55 = t45+t46+t47+t48+t54;
t56 = t11.^2;
t57 = L2.*t12.*t56.*2.0;
t58 = t13.^2;
t59 = L2.*t12.*t58.*2.0;
t60 = t57+t59;
t61 = L2.*t11.*t20.*t21;
t62 = L3.*t11.*t12.*t20;
t63 = L2.*t13.*t21.*t22;
t64 = L3.*t12.*t13.*t22;
t65 = t61+t62+t63+t64;
t66 = L2.*t11.*t29.*t30;
t67 = L4.*t11.*t12.*t29;
t68 = L2.*t13.*t30.*t31;
t69 = L4.*t12.*t13.*t31;
t70 = t66+t67+t68+t69;
t71 = m4.*t27.*(1.0./2.0);
t72 = m5.*t27.*(1.0./2.0);
t73 = m6.*t27.*(1.0./2.0);
t74 = L1.*t2.*t20.*t28;
t75 = d3.*t2.*t3.*t20;
t76 = L1.*t4.*t22.*t28;
t77 = d3.*t3.*t4.*t22;
t78 = t74+t75+t76+t77;
t79 = m3.*t78.*(1.0./2.0);
t80 = t71+t72+t73+t79;
t81 = m4.*t65.*(1.0./2.0);
t82 = m5.*t65.*(1.0./2.0);
t83 = m6.*t65.*(1.0./2.0);
t84 = L2.*t11.*t20.*t28;
t85 = d3.*t11.*t12.*t20;
t86 = L2.*t13.*t22.*t28;
t87 = d3.*t12.*t13.*t22;
t88 = t84+t85+t86+t87;
t89 = m3.*t88.*(1.0./2.0);
t90 = t81+t82+t83+t89;
t91 = t20.^2;
t92 = L3.*t21.*t91.*2.0;
t93 = t22.^2;
t94 = L3.*t21.*t93.*2.0;
t95 = t92+t94;
t96 = L3.*t20.*t29.*t30;
t97 = L4.*t20.*t21.*t29;
t98 = L3.*t22.*t30.*t31;
t99 = L4.*t21.*t22.*t31;
t100 = t96+t97+t98+t99;
t101 = m5.*t36.*(1.0./2.0);
t102 = m6.*t36.*(1.0./2.0);
t103 = L1.*t2.*t29.*t37;
t104 = d4.*t2.*t3.*t29;
t105 = L1.*t4.*t31.*t37;
t106 = d4.*t3.*t4.*t31;
t107 = t103+t104+t105+t106;
t108 = m4.*t107.*(1.0./2.0);
t109 = t101+t102+t108;
t110 = m5.*t70.*(1.0./2.0);
t111 = m6.*t70.*(1.0./2.0);
t112 = L2.*t11.*t29.*t37;
t113 = d4.*t11.*t12.*t29;
t114 = L2.*t13.*t31.*t37;
t115 = d4.*t12.*t13.*t31;
t116 = t112+t113+t114+t115;
t117 = m4.*t116.*(1.0./2.0);
t118 = t110+t111+t117;
t119 = m5.*t100.*(1.0./2.0);
t120 = m6.*t100.*(1.0./2.0);
t121 = L3.*t20.*t29.*t37;
t122 = d4.*t20.*t21.*t29;
t123 = L3.*t22.*t31.*t37;
t124 = d4.*t21.*t22.*t31;
t125 = t121+t122+t123+t124;
t126 = m4.*t125.*(1.0./2.0);
t127 = t119+t120+t126;
t128 = t29.^2;
t129 = L4.*t30.*t128.*2.0;
t130 = t31.^2;
t131 = L4.*t30.*t130.*2.0;
t132 = t129+t131;
t133 = L1.*t2.*t38.*t39;
t134 = L5.*t2.*t3.*t38;
t135 = L1.*t4.*t39.*t40;
t136 = L5.*t3.*t4.*t40;
t137 = t133+t134+t135+t136;
t138 = m6.*t137.*(1.0./2.0);
t139 = L1.*t2.*t38.*t41;
t140 = d5.*t2.*t3.*t38;
t141 = L1.*t4.*t40.*t41;
t142 = d5.*t3.*t4.*t40;
t143 = t139+t140+t141+t142;
t144 = m5.*t143.*(1.0./2.0);
t145 = t138+t144;
t146 = L2.*t11.*t38.*t39;
t147 = L5.*t11.*t12.*t38;
t148 = L2.*t13.*t39.*t40;
t149 = L5.*t12.*t13.*t40;
t150 = t146+t147+t148+t149;
t151 = m6.*t150.*(1.0./2.0);
t152 = L2.*t11.*t38.*t41;
t153 = d5.*t11.*t12.*t38;
t154 = L2.*t13.*t40.*t41;
t155 = d5.*t12.*t13.*t40;
t156 = t152+t153+t154+t155;
t157 = m5.*t156.*(1.0./2.0);
t158 = t151+t157;
t159 = L3.*t20.*t38.*t39;
t160 = L5.*t20.*t21.*t38;
t161 = L3.*t22.*t39.*t40;
t162 = L5.*t21.*t22.*t40;
t163 = t159+t160+t161+t162;
t164 = m6.*t163.*(1.0./2.0);
t165 = L3.*t20.*t38.*t41;
t166 = d5.*t20.*t21.*t38;
t167 = L3.*t22.*t40.*t41;
t168 = d5.*t21.*t22.*t40;
t169 = t165+t166+t167+t168;
t170 = m5.*t169.*(1.0./2.0);
t171 = t164+t170;
t172 = L4.*t29.*t38.*t39;
t173 = L5.*t29.*t30.*t38;
t174 = L4.*t31.*t39.*t40;
t175 = L5.*t30.*t31.*t40;
t176 = t172+t173+t174+t175;
t177 = m6.*t176.*(1.0./2.0);
t178 = L4.*t29.*t38.*t41;
t179 = d5.*t29.*t30.*t38;
t180 = L4.*t31.*t40.*t41;
t181 = d5.*t30.*t31.*t40;
t182 = t178+t179+t180+t181;
t183 = m5.*t182.*(1.0./2.0);
t184 = t177+t183;
t185 = t38.^2;
t186 = t40.^2;
t187 = L1.*t2.*t42.*t43;
t188 = d6.*t2.*t3.*t42;
t189 = L1.*t4.*t43.*t44;
t190 = d6.*t3.*t4.*t44;
t191 = t187+t188+t189+t190;
t192 = m6.*t191.*(1.0./2.0);
t193 = L2.*t11.*t42.*t43;
t194 = d6.*t11.*t12.*t42;
t195 = L2.*t13.*t43.*t44;
t196 = d6.*t12.*t13.*t44;
t197 = t193+t194+t195+t196;
t198 = m6.*t197.*(1.0./2.0);
t199 = L3.*t20.*t42.*t43;
t200 = d6.*t20.*t21.*t42;
t201 = L3.*t22.*t43.*t44;
t202 = d6.*t21.*t22.*t44;
t203 = t199+t200+t201+t202;
t204 = m6.*t203.*(1.0./2.0);
t205 = L4.*t29.*t42.*t43;
t206 = d6.*t29.*t30.*t42;
t207 = L4.*t31.*t43.*t44;
t208 = d6.*t30.*t31.*t44;
t209 = t205+t206+t207+t208;
t210 = m6.*t209.*(1.0./2.0);
t211 = L5.*t38.*t42.*t43;
t212 = d6.*t38.*t39.*t42;
t213 = L5.*t40.*t43.*t44;
t214 = d6.*t39.*t40.*t44;
t215 = t211+t212+t213+t214;
t216 = m6.*t215.*(1.0./2.0);
M = reshape([I1+m2.*t9.*(1.0./2.0)+m3.*t9.*(1.0./2.0)+m4.*t9.*(1.0./2.0)+m5.*t9.*(1.0./2.0)+m6.*t9.*(1.0./2.0)+m1.*(d1.*t5.*t10.*2.0+d1.*t7.*t10.*2.0).*(1.0./2.0),t55,t80,t109,t145,t192,t55,I2+m3.*t60.*(1.0./2.0)+m4.*t60.*(1.0./2.0)+m5.*t60.*(1.0./2.0)+m6.*t60.*(1.0./2.0)+m2.*(d2.*t19.*t56.*2.0+d2.*t19.*t58.*2.0).*(1.0./2.0),t90,t118,t158,t198,t80,t90,I3+m4.*t95.*(1.0./2.0)+m5.*t95.*(1.0./2.0)+m6.*t95.*(1.0./2.0)+m3.*(d3.*t28.*t91.*2.0+d3.*t28.*t93.*2.0).*(1.0./2.0),t127,t171,t204,t109,t118,t127,I4+m5.*t132.*(1.0./2.0)+m6.*t132.*(1.0./2.0)+m4.*(d4.*t37.*t128.*2.0+d4.*t37.*t130.*2.0).*(1.0./2.0),t184,t210,t145,t158,t171,t184,I5+m6.*(L5.*t39.*t185.*2.0+L5.*t39.*t186.*2.0).*(1.0./2.0)+m5.*(d5.*t41.*t185.*2.0+d5.*t41.*t186.*2.0).*(1.0./2.0),t216,t192,t198,t204,t210,t216,I6+m6.*(d6.*t42.^2.*t43.*2.0+d6.*t43.*t44.^2.*2.0).*(1.0./2.0)],[6,6]);