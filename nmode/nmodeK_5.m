function K = nmodeK_5(in1,in2,in3,g,in5,in6,in7)
%NMODEK_5
%    K = NMODEK_5(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    11-Dec-2018 19:37:19

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
t2 = conj(L1);
t3 = omega1.*t2;
t4 = conj(L2);
t5 = omega2.*t4;
t6 = L1.*omega1;
t7 = L2.*omega2;
t8 = L3.*omega3;
t9 = conj(L3);
t10 = omega3.*t9;
t11 = conj(d3);
t12 = omega3.*t11;
t13 = t3+t5+t12;
t14 = d3.*omega3;
t15 = t6+t7+t14;
t16 = conj(L4);
t17 = omega4.*t16;
t18 = conj(d5);
t19 = omega5.*t18;
t20 = t3+t5+t10+t17+t19;
t21 = omega1.^2;
t22 = L4.*omega4;
t23 = d5.*omega5;
t24 = t6+t7+t8+t22+t23;
t25 = L1.*omega2.*t4;
t26 = L2.*omega2.*t2;
t27 = t25+t26;
t28 = conj(d2);
t29 = d2.*omega2;
t30 = t6+t29;
t31 = omega2.*t28;
t32 = t3+t31;
t33 = L1.*omega3.*t9;
t34 = L3.*omega3.*t2;
t35 = t33+t34;
t36 = d4.*omega4;
t37 = t6+t7+t8+t36;
t38 = conj(d4);
t39 = omega4.*t38;
t40 = t3+t5+t10+t39;
t41 = conj(m2);
t42 = conj(g);
t43 = conj(m3);
t44 = conj(m4);
t45 = conj(m5);
t46 = L1.*omega1.*omega2.*t4;
t47 = L2.*omega1.*omega2.*t2;
t48 = t46+t47;
t49 = m3.*t27.*(1.0./2.0);
t50 = m4.*t27.*(1.0./2.0);
t51 = m5.*t27.*(1.0./2.0);
t52 = L1.*omega2.*t28;
t53 = d2.*omega2.*t2;
t54 = t52+t53;
t55 = m2.*t54.*(1.0./2.0);
t56 = t49+t50+t51+t55;
t57 = omega2.*t56;
t58 = L1.*omega1.*omega3.*t9;
t59 = L3.*omega1.*omega3.*t2;
t60 = t58+t59;
t61 = m4.*t35.*(1.0./2.0);
t62 = m5.*t35.*(1.0./2.0);
t63 = L1.*omega3.*t11;
t64 = d3.*omega3.*t2;
t65 = t63+t64;
t66 = m3.*t65.*(1.0./2.0);
t67 = t61+t62+t66;
t68 = omega3.*t67;
t69 = L1.*omega4.*t16;
t70 = L4.*omega4.*t2;
t71 = t69+t70;
t72 = m5.*t71.*(1.0./2.0);
t73 = L1.*omega4.*t38;
t74 = d4.*omega4.*t2;
t75 = t73+t74;
t76 = m4.*t75.*(1.0./2.0);
t77 = t72+t76;
t78 = omega4.*t77;
t79 = L1.*omega5.*t18;
t80 = d5.*omega5.*t2;
t81 = t79+t80;
t82 = m5.*omega5.*t81.*(1.0./2.0);
t83 = L1.*omega1.*omega2.*t28;
t84 = d2.*omega1.*omega2.*t2;
t85 = t83+t84;
t86 = L1.*omega1.*t4;
t87 = L2.*omega1.*t2;
t88 = t86+t87;
t89 = m3.*t88.*(1.0./2.0);
t90 = m4.*t88.*(1.0./2.0);
t91 = m5.*t88.*(1.0./2.0);
t92 = L1.*omega1.*t28;
t93 = d2.*omega1.*t2;
t94 = t92+t93;
t95 = m2.*t94.*(1.0./2.0);
t96 = t89+t90+t91+t95;
t97 = omega2.^2;
t98 = L2.*omega3.*t9;
t99 = L3.*omega3.*t4;
t100 = t98+t99;
t101 = L2.*omega2.*omega3.*t9;
t102 = L3.*omega2.*omega3.*t4;
t103 = t101+t102;
t104 = m4.*t100.*(1.0./2.0);
t105 = m5.*t100.*(1.0./2.0);
t106 = L2.*omega3.*t11;
t107 = d3.*omega3.*t4;
t108 = t106+t107;
t109 = m3.*t108.*(1.0./2.0);
t110 = t104+t105+t109;
t111 = omega3.*t110;
t112 = L2.*omega4.*t16;
t113 = L4.*omega4.*t4;
t114 = t112+t113;
t115 = m5.*t114.*(1.0./2.0);
t116 = L2.*omega4.*t38;
t117 = d4.*omega4.*t4;
t118 = t116+t117;
t119 = m4.*t118.*(1.0./2.0);
t120 = t115+t119;
t121 = omega4.*t120;
t122 = L2.*omega5.*t18;
t123 = d5.*omega5.*t4;
t124 = t122+t123;
t125 = m5.*omega5.*t124.*(1.0./2.0);
t126 = L1.*omega1.*omega3.*t11;
t127 = d3.*omega1.*omega3.*t2;
t128 = t126+t127;
t129 = L1.*omega1.*t9;
t130 = L3.*omega1.*t2;
t131 = t129+t130;
t132 = m4.*t131.*(1.0./2.0);
t133 = m5.*t131.*(1.0./2.0);
t134 = L1.*omega1.*t11;
t135 = d3.*omega1.*t2;
t136 = t134+t135;
t137 = m3.*t136.*(1.0./2.0);
t138 = t132+t133+t137;
t139 = L2.*omega2.*omega3.*t11;
t140 = d3.*omega2.*omega3.*t4;
t141 = t139+t140;
t142 = L2.*omega2.*t9;
t143 = L3.*omega2.*t4;
t144 = t142+t143;
t145 = m4.*t144.*(1.0./2.0);
t146 = m5.*t144.*(1.0./2.0);
t147 = L2.*omega2.*t11;
t148 = d3.*omega2.*t4;
t149 = t147+t148;
t150 = m3.*t149.*(1.0./2.0);
t151 = t145+t146+t150;
t152 = omega3.^2;
t153 = L3.*omega4.*t16;
t154 = L4.*omega4.*t9;
t155 = t153+t154;
t156 = m5.*t155.*(1.0./2.0);
t157 = L3.*omega4.*t38;
t158 = d4.*omega4.*t9;
t159 = t157+t158;
t160 = m4.*t159.*(1.0./2.0);
t161 = t156+t160;
t162 = omega4.*t161;
t163 = L3.*omega5.*t18;
t164 = d5.*omega5.*t9;
t165 = t163+t164;
t166 = m5.*omega5.*t165.*(1.0./2.0);
t167 = L1.*omega1.*omega4.*t16;
t168 = L4.*omega1.*omega4.*t2;
t169 = t167+t168;
t170 = L1.*omega1.*omega4.*t38;
t171 = d4.*omega1.*omega4.*t2;
t172 = t170+t171;
t173 = L1.*omega1.*t16;
t174 = L4.*omega1.*t2;
t175 = t173+t174;
t176 = m5.*t175.*(1.0./2.0);
t177 = L1.*omega1.*t38;
t178 = d4.*omega1.*t2;
t179 = t177+t178;
t180 = m4.*t179.*(1.0./2.0);
t181 = t176+t180;
t182 = L2.*omega2.*omega4.*t16;
t183 = L4.*omega2.*omega4.*t4;
t184 = t182+t183;
t185 = L2.*omega2.*omega4.*t38;
t186 = d4.*omega2.*omega4.*t4;
t187 = t185+t186;
t188 = L2.*omega2.*t16;
t189 = L4.*omega2.*t4;
t190 = t188+t189;
t191 = m5.*t190.*(1.0./2.0);
t192 = L2.*omega2.*t38;
t193 = d4.*omega2.*t4;
t194 = t192+t193;
t195 = m4.*t194.*(1.0./2.0);
t196 = t191+t195;
t197 = L3.*omega3.*omega4.*t16;
t198 = L4.*omega3.*omega4.*t9;
t199 = t197+t198;
t200 = L3.*omega3.*omega4.*t38;
t201 = d4.*omega3.*omega4.*t9;
t202 = t200+t201;
t203 = L3.*omega3.*t16;
t204 = L4.*omega3.*t9;
t205 = t203+t204;
t206 = m5.*t205.*(1.0./2.0);
t207 = L3.*omega3.*t38;
t208 = d4.*omega3.*t9;
t209 = t207+t208;
t210 = m4.*t209.*(1.0./2.0);
t211 = t206+t210;
t212 = omega4.^2;
t213 = L4.*omega5.*t18;
t214 = d5.*omega5.*t16;
t215 = t213+t214;
t216 = m5.*omega5.*t215.*(1.0./2.0);
t217 = L1.*omega1.*omega5.*t18;
t218 = d5.*omega1.*omega5.*t2;
t219 = t217+t218;
t220 = L1.*omega1.*t18;
t221 = d5.*omega1.*t2;
t222 = t220+t221;
t223 = L2.*omega2.*omega5.*t18;
t224 = d5.*omega2.*omega5.*t4;
t225 = t223+t224;
t226 = L2.*omega2.*t18;
t227 = d5.*omega2.*t4;
t228 = t226+t227;
t229 = L3.*omega3.*omega5.*t18;
t230 = d5.*omega3.*omega5.*t9;
t231 = t229+t230;
t232 = L3.*omega3.*t18;
t233 = d5.*omega3.*t9;
t234 = t232+t233;
t235 = L4.*omega4.*omega5.*t18;
t236 = d5.*omega4.*omega5.*t16;
t237 = t235+t236;
t238 = L4.*omega4.*t18;
t239 = d5.*omega4.*t16;
t240 = t238+t239;
K = reshape([t57+t68+t78+t82+t43.*(L1.*omega1.*t13-L1.*t2.*t21.*2.0+omega1.*t2.*t15).*(1.0./2.0)+t45.*(L1.*omega1.*t20-L1.*t2.*t21.*2.0+omega1.*t2.*t24).*(1.0./2.0)+t41.*(L1.*omega1.*t32-L1.*t2.*t21.*2.0+omega1.*t2.*t30).*(1.0./2.0)+t44.*(L1.*omega1.*t40-L1.*t2.*t21.*2.0+omega1.*t2.*t37).*(1.0./2.0)-omega1.*(m3.*(L1.*t13+t2.*t15-L1.*omega1.*t2.*2.0).*(1.0./2.0)+m5.*(L1.*t20+t2.*t24-L1.*omega1.*t2.*2.0).*(1.0./2.0)+m2.*(L1.*t32+t2.*t30-L1.*omega1.*t2.*2.0).*(1.0./2.0)+m4.*(L1.*t40+t2.*t37-L1.*omega1.*t2.*2.0).*(1.0./2.0))+d1.*t42.*conj(m1)+L1.*t41.*t42+L1.*t42.*t43+L1.*t42.*t44+L1.*t42.*t45,-omega1.*t96+omega2.*t96-t43.*t48.*(1.0./2.0)-t44.*t48.*(1.0./2.0)-t45.*t48.*(1.0./2.0)-t41.*t85.*(1.0./2.0),-omega1.*t138+omega3.*t138-t44.*t60.*(1.0./2.0)-t45.*t60.*(1.0./2.0)-t43.*t128.*(1.0./2.0),-omega1.*t181+omega4.*t181-t45.*t169.*(1.0./2.0)-t44.*t172.*(1.0./2.0),t45.*t219.*(-1.0./2.0)-m5.*omega1.*t222.*(1.0./2.0)+m5.*omega5.*t222.*(1.0./2.0),-t57+omega1.*t56-t43.*t48.*(1.0./2.0)-t44.*t48.*(1.0./2.0)-t45.*t48.*(1.0./2.0)-t41.*t85.*(1.0./2.0),t111+t121+t125+t43.*(L2.*omega2.*t13-L2.*t4.*t97.*2.0+omega2.*t4.*t15).*(1.0./2.0)+t45.*(L2.*omega2.*t20-L2.*t4.*t97.*2.0+omega2.*t4.*t24).*(1.0./2.0)+t44.*(L2.*omega2.*t40-L2.*t4.*t97.*2.0+omega2.*t4.*t37).*(1.0./2.0)+t41.*(d2.*omega2.*t32-d2.*t28.*t97.*2.0+omega2.*t28.*t30).*(1.0./2.0)+omega1.*t96-omega2.*(m3.*(L2.*t13+t4.*t15-L2.*omega2.*t4.*2.0).*(1.0./2.0)+m5.*(L2.*t20+t4.*t24-L2.*omega2.*t4.*2.0).*(1.0./2.0)+m4.*(L2.*t40+t4.*t37-L2.*omega2.*t4.*2.0).*(1.0./2.0)+m2.*(d2.*t32+t28.*t30-d2.*omega2.*t28.*2.0).*(1.0./2.0))+L2.*t42.*t43+L2.*t42.*t44+L2.*t42.*t45+d2.*t41.*t42,-omega2.*t151+omega3.*t151-t44.*t103.*(1.0./2.0)-t45.*t103.*(1.0./2.0)-t43.*t141.*(1.0./2.0),-omega2.*t196+omega4.*t196-t45.*t184.*(1.0./2.0)-t44.*t187.*(1.0./2.0),t45.*t225.*(-1.0./2.0)-m5.*omega2.*t228.*(1.0./2.0)+m5.*omega5.*t228.*(1.0./2.0),-t68+omega1.*t67-t44.*t60.*(1.0./2.0)-t45.*t60.*(1.0./2.0)-t43.*t128.*(1.0./2.0),-t111+omega2.*t110-t44.*t103.*(1.0./2.0)-t45.*t103.*(1.0./2.0)-t43.*t141.*(1.0./2.0),t162+t166+t45.*(L3.*omega3.*t20-L3.*t9.*t152.*2.0+omega3.*t9.*t24).*(1.0./2.0)+t44.*(L3.*omega3.*t40-L3.*t9.*t152.*2.0+omega3.*t9.*t37).*(1.0./2.0)+t43.*(d3.*omega3.*t13-d3.*t11.*t152.*2.0+omega3.*t11.*t15).*(1.0./2.0)+omega1.*t138+omega2.*t151-omega3.*(m5.*(L3.*t20+t9.*t24-L3.*omega3.*t9.*2.0).*(1.0./2.0)+m4.*(L3.*t40+t9.*t37-L3.*omega3.*t9.*2.0).*(1.0./2.0)+m3.*(d3.*t13+t11.*t15-d3.*omega3.*t11.*2.0).*(1.0./2.0))+L3.*t42.*t44+L3.*t42.*t45+d3.*t42.*t43,-omega3.*t211+omega4.*t211-t45.*t199.*(1.0./2.0)-t44.*t202.*(1.0./2.0),t45.*t231.*(-1.0./2.0)-m5.*omega3.*t234.*(1.0./2.0)+m5.*omega5.*t234.*(1.0./2.0),-t78+omega1.*t77-t45.*t169.*(1.0./2.0)-t44.*t172.*(1.0./2.0),-t121+omega2.*t120-t45.*t184.*(1.0./2.0)-t44.*t187.*(1.0./2.0),-t162+omega3.*t161-t45.*t199.*(1.0./2.0)-t44.*t202.*(1.0./2.0),t216+t45.*(L4.*omega4.*t20-L4.*t16.*t212.*2.0+omega4.*t16.*t24).*(1.0./2.0)+t44.*(d4.*omega4.*t40-d4.*t38.*t212.*2.0+omega4.*t37.*t38).*(1.0./2.0)+omega1.*t181+omega2.*t196+omega3.*t211-omega4.*(m5.*(L4.*t20+t16.*t24-L4.*omega4.*t16.*2.0).*(1.0./2.0)+m4.*(d4.*t40+t37.*t38-d4.*omega4.*t38.*2.0).*(1.0./2.0))+L4.*t42.*t45+d4.*t42.*t44,t45.*t237.*(-1.0./2.0)-m5.*omega4.*t240.*(1.0./2.0)+m5.*omega5.*t240.*(1.0./2.0),-t82-t45.*t219.*(1.0./2.0)+m5.*omega1.*t81.*(1.0./2.0),-t125-t45.*t225.*(1.0./2.0)+m5.*omega2.*t124.*(1.0./2.0),-t166-t45.*t231.*(1.0./2.0)+m5.*omega3.*t165.*(1.0./2.0),-t216-t45.*t237.*(1.0./2.0)+m5.*omega4.*t215.*(1.0./2.0),t45.*(d5.*omega5.^2.*t18.*-2.0+d5.*omega5.*t20+omega5.*t18.*t24).*(1.0./2.0)+d5.*t42.*t45+m5.*omega1.*t222.*(1.0./2.0)+m5.*omega2.*t228.*(1.0./2.0)+m5.*omega3.*t234.*(1.0./2.0)+m5.*omega4.*t240.*(1.0./2.0)-m5.*omega5.*(d5.*t20+t18.*t24-d5.*omega5.*t18.*2.0).*(1.0./2.0)],[5,5]);
