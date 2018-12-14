function M = npend_alphas_Lagrange_M_8(in1,in2,in3,g,in5,in6,in7)
%NPEND_ALPHAS_LAGRANGE_M_8
%    M = NPEND_ALPHAS_LAGRANGE_M_8(IN1,IN2,IN3,G,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    03-Dec-2018 23:33:31

I1 = in1(1,:);
I2 = in1(2,:);
I3 = in1(3,:);
I4 = in1(4,:);
I5 = in1(5,:);
I6 = in1(6,:);
I7 = in1(7,:);
I8 = in1(8,:);
L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
L6 = in2(6,:);
L7 = in2(7,:);
d1 = in3(1,:);
d2 = in3(2,:);
d3 = in3(3,:);
d4 = in3(4,:);
d5 = in3(5,:);
d6 = in3(6,:);
d7 = in3(7,:);
d8 = in3(8,:);
m1 = in5(1,:);
m2 = in5(2,:);
m3 = in5(3,:);
m4 = in5(4,:);
m5 = in5(5,:);
m6 = in5(6,:);
m7 = in5(7,:);
m8 = in5(8,:);
theta1 = in7(1,:);
theta2 = in7(2,:);
theta3 = in7(3,:);
theta4 = in7(4,:);
theta5 = in7(5,:);
theta6 = in7(6,:);
theta7 = in7(7,:);
theta8 = in7(8,:);
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
t41 = L1.*t2.*t38.*t39;
t42 = L5.*t2.*t3.*t38;
t43 = L1.*t4.*t39.*t40;
t44 = L5.*t3.*t4.*t40;
t45 = t41+t42+t43+t44;
t46 = conj(d5);
t47 = cos(theta6);
t48 = conj(L6);
t49 = sin(theta6);
t50 = L1.*t2.*t47.*t48;
t51 = L6.*t2.*t3.*t47;
t52 = L1.*t4.*t48.*t49;
t53 = L6.*t3.*t4.*t49;
t54 = t50+t51+t52+t53;
t55 = conj(d6);
t56 = cos(theta7);
t57 = conj(L7);
t58 = sin(theta7);
t59 = conj(d7);
t60 = cos(theta8);
t61 = conj(d8);
t62 = sin(theta8);
t63 = m3.*t18.*(1.0./2.0);
t64 = m4.*t18.*(1.0./2.0);
t65 = m5.*t18.*(1.0./2.0);
t66 = m6.*t18.*(1.0./2.0);
t67 = m7.*t18.*(1.0./2.0);
t68 = m8.*t18.*(1.0./2.0);
t69 = L1.*t2.*t11.*t19;
t70 = d2.*t2.*t3.*t11;
t71 = L1.*t4.*t13.*t19;
t72 = d2.*t3.*t4.*t13;
t73 = t69+t70+t71+t72;
t74 = m2.*t73.*(1.0./2.0);
t75 = t63+t64+t65+t66+t67+t68+t74;
t76 = t11.^2;
t77 = L2.*t12.*t76.*2.0;
t78 = t13.^2;
t79 = L2.*t12.*t78.*2.0;
t80 = t77+t79;
t81 = L2.*t11.*t20.*t21;
t82 = L3.*t11.*t12.*t20;
t83 = L2.*t13.*t21.*t22;
t84 = L3.*t12.*t13.*t22;
t85 = t81+t82+t83+t84;
t86 = L2.*t11.*t29.*t30;
t87 = L4.*t11.*t12.*t29;
t88 = L2.*t13.*t30.*t31;
t89 = L4.*t12.*t13.*t31;
t90 = t86+t87+t88+t89;
t91 = L2.*t11.*t38.*t39;
t92 = L5.*t11.*t12.*t38;
t93 = L2.*t13.*t39.*t40;
t94 = L5.*t12.*t13.*t40;
t95 = t91+t92+t93+t94;
t96 = L2.*t11.*t47.*t48;
t97 = L6.*t11.*t12.*t47;
t98 = L2.*t13.*t48.*t49;
t99 = L6.*t12.*t13.*t49;
t100 = t96+t97+t98+t99;
t101 = m4.*t27.*(1.0./2.0);
t102 = m5.*t27.*(1.0./2.0);
t103 = m6.*t27.*(1.0./2.0);
t104 = m7.*t27.*(1.0./2.0);
t105 = m8.*t27.*(1.0./2.0);
t106 = L1.*t2.*t20.*t28;
t107 = d3.*t2.*t3.*t20;
t108 = L1.*t4.*t22.*t28;
t109 = d3.*t3.*t4.*t22;
t110 = t106+t107+t108+t109;
t111 = m3.*t110.*(1.0./2.0);
t112 = t101+t102+t103+t104+t105+t111;
t113 = m4.*t85.*(1.0./2.0);
t114 = m5.*t85.*(1.0./2.0);
t115 = m6.*t85.*(1.0./2.0);
t116 = m7.*t85.*(1.0./2.0);
t117 = m8.*t85.*(1.0./2.0);
t118 = L2.*t11.*t20.*t28;
t119 = d3.*t11.*t12.*t20;
t120 = L2.*t13.*t22.*t28;
t121 = d3.*t12.*t13.*t22;
t122 = t118+t119+t120+t121;
t123 = m3.*t122.*(1.0./2.0);
t124 = t113+t114+t115+t116+t117+t123;
t125 = t20.^2;
t126 = L3.*t21.*t125.*2.0;
t127 = t22.^2;
t128 = L3.*t21.*t127.*2.0;
t129 = t126+t128;
t130 = L3.*t20.*t29.*t30;
t131 = L4.*t20.*t21.*t29;
t132 = L3.*t22.*t30.*t31;
t133 = L4.*t21.*t22.*t31;
t134 = t130+t131+t132+t133;
t135 = L3.*t20.*t38.*t39;
t136 = L5.*t20.*t21.*t38;
t137 = L3.*t22.*t39.*t40;
t138 = L5.*t21.*t22.*t40;
t139 = t135+t136+t137+t138;
t140 = L3.*t20.*t47.*t48;
t141 = L6.*t20.*t21.*t47;
t142 = L3.*t22.*t48.*t49;
t143 = L6.*t21.*t22.*t49;
t144 = t140+t141+t142+t143;
t145 = m5.*t36.*(1.0./2.0);
t146 = m6.*t36.*(1.0./2.0);
t147 = m7.*t36.*(1.0./2.0);
t148 = m8.*t36.*(1.0./2.0);
t149 = L1.*t2.*t29.*t37;
t150 = d4.*t2.*t3.*t29;
t151 = L1.*t4.*t31.*t37;
t152 = d4.*t3.*t4.*t31;
t153 = t149+t150+t151+t152;
t154 = m4.*t153.*(1.0./2.0);
t155 = t145+t146+t147+t148+t154;
t156 = m5.*t90.*(1.0./2.0);
t157 = m6.*t90.*(1.0./2.0);
t158 = m7.*t90.*(1.0./2.0);
t159 = m8.*t90.*(1.0./2.0);
t160 = L2.*t11.*t29.*t37;
t161 = d4.*t11.*t12.*t29;
t162 = L2.*t13.*t31.*t37;
t163 = d4.*t12.*t13.*t31;
t164 = t160+t161+t162+t163;
t165 = m4.*t164.*(1.0./2.0);
t166 = t156+t157+t158+t159+t165;
t167 = m5.*t134.*(1.0./2.0);
t168 = m6.*t134.*(1.0./2.0);
t169 = m7.*t134.*(1.0./2.0);
t170 = m8.*t134.*(1.0./2.0);
t171 = L3.*t20.*t29.*t37;
t172 = d4.*t20.*t21.*t29;
t173 = L3.*t22.*t31.*t37;
t174 = d4.*t21.*t22.*t31;
t175 = t171+t172+t173+t174;
t176 = m4.*t175.*(1.0./2.0);
t177 = t167+t168+t169+t170+t176;
t178 = t29.^2;
t179 = L4.*t30.*t178.*2.0;
t180 = t31.^2;
t181 = L4.*t30.*t180.*2.0;
t182 = t179+t181;
t183 = L4.*t29.*t38.*t39;
t184 = L5.*t29.*t30.*t38;
t185 = L4.*t31.*t39.*t40;
t186 = L5.*t30.*t31.*t40;
t187 = t183+t184+t185+t186;
t188 = L4.*t29.*t47.*t48;
t189 = L6.*t29.*t30.*t47;
t190 = L4.*t31.*t48.*t49;
t191 = L6.*t30.*t31.*t49;
t192 = t188+t189+t190+t191;
t193 = m6.*t45.*(1.0./2.0);
t194 = m7.*t45.*(1.0./2.0);
t195 = m8.*t45.*(1.0./2.0);
t196 = L1.*t2.*t38.*t46;
t197 = d5.*t2.*t3.*t38;
t198 = L1.*t4.*t40.*t46;
t199 = d5.*t3.*t4.*t40;
t200 = t196+t197+t198+t199;
t201 = m5.*t200.*(1.0./2.0);
t202 = t193+t194+t195+t201;
t203 = m6.*t95.*(1.0./2.0);
t204 = m7.*t95.*(1.0./2.0);
t205 = m8.*t95.*(1.0./2.0);
t206 = L2.*t11.*t38.*t46;
t207 = d5.*t11.*t12.*t38;
t208 = L2.*t13.*t40.*t46;
t209 = d5.*t12.*t13.*t40;
t210 = t206+t207+t208+t209;
t211 = m5.*t210.*(1.0./2.0);
t212 = t203+t204+t205+t211;
t213 = m6.*t139.*(1.0./2.0);
t214 = m7.*t139.*(1.0./2.0);
t215 = m8.*t139.*(1.0./2.0);
t216 = L3.*t20.*t38.*t46;
t217 = d5.*t20.*t21.*t38;
t218 = L3.*t22.*t40.*t46;
t219 = d5.*t21.*t22.*t40;
t220 = t216+t217+t218+t219;
t221 = m5.*t220.*(1.0./2.0);
t222 = t213+t214+t215+t221;
t223 = m6.*t187.*(1.0./2.0);
t224 = m7.*t187.*(1.0./2.0);
t225 = m8.*t187.*(1.0./2.0);
t226 = L4.*t29.*t38.*t46;
t227 = d5.*t29.*t30.*t38;
t228 = L4.*t31.*t40.*t46;
t229 = d5.*t30.*t31.*t40;
t230 = t226+t227+t228+t229;
t231 = m5.*t230.*(1.0./2.0);
t232 = t223+t224+t225+t231;
t233 = t38.^2;
t234 = L5.*t39.*t233.*2.0;
t235 = t40.^2;
t236 = L5.*t39.*t235.*2.0;
t237 = t234+t236;
t238 = L5.*t38.*t47.*t48;
t239 = L6.*t38.*t39.*t47;
t240 = L5.*t40.*t48.*t49;
t241 = L6.*t39.*t40.*t49;
t242 = t238+t239+t240+t241;
t243 = m7.*t54.*(1.0./2.0);
t244 = m8.*t54.*(1.0./2.0);
t245 = L1.*t2.*t47.*t55;
t246 = d6.*t2.*t3.*t47;
t247 = L1.*t4.*t49.*t55;
t248 = d6.*t3.*t4.*t49;
t249 = t245+t246+t247+t248;
t250 = m6.*t249.*(1.0./2.0);
t251 = t243+t244+t250;
t252 = m7.*t100.*(1.0./2.0);
t253 = m8.*t100.*(1.0./2.0);
t254 = L2.*t11.*t47.*t55;
t255 = d6.*t11.*t12.*t47;
t256 = L2.*t13.*t49.*t55;
t257 = d6.*t12.*t13.*t49;
t258 = t254+t255+t256+t257;
t259 = m6.*t258.*(1.0./2.0);
t260 = t252+t253+t259;
t261 = m7.*t144.*(1.0./2.0);
t262 = m8.*t144.*(1.0./2.0);
t263 = L3.*t20.*t47.*t55;
t264 = d6.*t20.*t21.*t47;
t265 = L3.*t22.*t49.*t55;
t266 = d6.*t21.*t22.*t49;
t267 = t263+t264+t265+t266;
t268 = m6.*t267.*(1.0./2.0);
t269 = t261+t262+t268;
t270 = m7.*t192.*(1.0./2.0);
t271 = m8.*t192.*(1.0./2.0);
t272 = L4.*t29.*t47.*t55;
t273 = d6.*t29.*t30.*t47;
t274 = L4.*t31.*t49.*t55;
t275 = d6.*t30.*t31.*t49;
t276 = t272+t273+t274+t275;
t277 = m6.*t276.*(1.0./2.0);
t278 = t270+t271+t277;
t279 = m7.*t242.*(1.0./2.0);
t280 = m8.*t242.*(1.0./2.0);
t281 = L5.*t38.*t47.*t55;
t282 = d6.*t38.*t39.*t47;
t283 = L5.*t40.*t49.*t55;
t284 = d6.*t39.*t40.*t49;
t285 = t281+t282+t283+t284;
t286 = m6.*t285.*(1.0./2.0);
t287 = t279+t280+t286;
t288 = t47.^2;
t289 = L6.*t48.*t288.*2.0;
t290 = t49.^2;
t291 = L6.*t48.*t290.*2.0;
t292 = t289+t291;
t293 = L1.*t2.*t56.*t57;
t294 = L7.*t2.*t3.*t56;
t295 = L1.*t4.*t57.*t58;
t296 = L7.*t3.*t4.*t58;
t297 = t293+t294+t295+t296;
t298 = m8.*t297.*(1.0./2.0);
t299 = L1.*t2.*t56.*t59;
t300 = d7.*t2.*t3.*t56;
t301 = L1.*t4.*t58.*t59;
t302 = d7.*t3.*t4.*t58;
t303 = t299+t300+t301+t302;
t304 = m7.*t303.*(1.0./2.0);
t305 = t298+t304;
t306 = L2.*t11.*t56.*t57;
t307 = L7.*t11.*t12.*t56;
t308 = L2.*t13.*t57.*t58;
t309 = L7.*t12.*t13.*t58;
t310 = t306+t307+t308+t309;
t311 = m8.*t310.*(1.0./2.0);
t312 = L2.*t11.*t56.*t59;
t313 = d7.*t11.*t12.*t56;
t314 = L2.*t13.*t58.*t59;
t315 = d7.*t12.*t13.*t58;
t316 = t312+t313+t314+t315;
t317 = m7.*t316.*(1.0./2.0);
t318 = t311+t317;
t319 = L3.*t20.*t56.*t57;
t320 = L7.*t20.*t21.*t56;
t321 = L3.*t22.*t57.*t58;
t322 = L7.*t21.*t22.*t58;
t323 = t319+t320+t321+t322;
t324 = m8.*t323.*(1.0./2.0);
t325 = L3.*t20.*t56.*t59;
t326 = d7.*t20.*t21.*t56;
t327 = L3.*t22.*t58.*t59;
t328 = d7.*t21.*t22.*t58;
t329 = t325+t326+t327+t328;
t330 = m7.*t329.*(1.0./2.0);
t331 = t324+t330;
t332 = L4.*t29.*t56.*t57;
t333 = L7.*t29.*t30.*t56;
t334 = L4.*t31.*t57.*t58;
t335 = L7.*t30.*t31.*t58;
t336 = t332+t333+t334+t335;
t337 = m8.*t336.*(1.0./2.0);
t338 = L4.*t29.*t56.*t59;
t339 = d7.*t29.*t30.*t56;
t340 = L4.*t31.*t58.*t59;
t341 = d7.*t30.*t31.*t58;
t342 = t338+t339+t340+t341;
t343 = m7.*t342.*(1.0./2.0);
t344 = t337+t343;
t345 = L5.*t38.*t56.*t57;
t346 = L7.*t38.*t39.*t56;
t347 = L5.*t40.*t57.*t58;
t348 = L7.*t39.*t40.*t58;
t349 = t345+t346+t347+t348;
t350 = m8.*t349.*(1.0./2.0);
t351 = L5.*t38.*t56.*t59;
t352 = d7.*t38.*t39.*t56;
t353 = L5.*t40.*t58.*t59;
t354 = d7.*t39.*t40.*t58;
t355 = t351+t352+t353+t354;
t356 = m7.*t355.*(1.0./2.0);
t357 = t350+t356;
t358 = L6.*t47.*t56.*t57;
t359 = L7.*t47.*t48.*t56;
t360 = L6.*t49.*t57.*t58;
t361 = L7.*t48.*t49.*t58;
t362 = t358+t359+t360+t361;
t363 = m8.*t362.*(1.0./2.0);
t364 = L6.*t47.*t56.*t59;
t365 = d7.*t47.*t48.*t56;
t366 = L6.*t49.*t58.*t59;
t367 = d7.*t48.*t49.*t58;
t368 = t364+t365+t366+t367;
t369 = m7.*t368.*(1.0./2.0);
t370 = t363+t369;
t371 = t56.^2;
t372 = t58.^2;
t373 = L1.*t2.*t60.*t61;
t374 = d8.*t2.*t3.*t60;
t375 = L1.*t4.*t61.*t62;
t376 = d8.*t3.*t4.*t62;
t377 = t373+t374+t375+t376;
t378 = m8.*t377.*(1.0./2.0);
t379 = L2.*t11.*t60.*t61;
t380 = d8.*t11.*t12.*t60;
t381 = L2.*t13.*t61.*t62;
t382 = d8.*t12.*t13.*t62;
t383 = t379+t380+t381+t382;
t384 = m8.*t383.*(1.0./2.0);
t385 = L3.*t20.*t60.*t61;
t386 = d8.*t20.*t21.*t60;
t387 = L3.*t22.*t61.*t62;
t388 = d8.*t21.*t22.*t62;
t389 = t385+t386+t387+t388;
t390 = m8.*t389.*(1.0./2.0);
t391 = L4.*t29.*t60.*t61;
t392 = d8.*t29.*t30.*t60;
t393 = L4.*t31.*t61.*t62;
t394 = d8.*t30.*t31.*t62;
t395 = t391+t392+t393+t394;
t396 = m8.*t395.*(1.0./2.0);
t397 = L5.*t38.*t60.*t61;
t398 = d8.*t38.*t39.*t60;
t399 = L5.*t40.*t61.*t62;
t400 = d8.*t39.*t40.*t62;
t401 = t397+t398+t399+t400;
t402 = m8.*t401.*(1.0./2.0);
t403 = L6.*t47.*t60.*t61;
t404 = d8.*t47.*t48.*t60;
t405 = L6.*t49.*t61.*t62;
t406 = d8.*t48.*t49.*t62;
t407 = t403+t404+t405+t406;
t408 = m8.*t407.*(1.0./2.0);
t409 = L7.*t56.*t60.*t61;
t410 = d8.*t56.*t57.*t60;
t411 = L7.*t58.*t61.*t62;
t412 = d8.*t57.*t58.*t62;
t413 = t409+t410+t411+t412;
t414 = m8.*t413.*(1.0./2.0);
M = reshape([I1+m2.*t9.*(1.0./2.0)+m3.*t9.*(1.0./2.0)+m4.*t9.*(1.0./2.0)+m5.*t9.*(1.0./2.0)+m6.*t9.*(1.0./2.0)+m7.*t9.*(1.0./2.0)+m8.*t9.*(1.0./2.0)+m1.*(d1.*t5.*t10.*2.0+d1.*t7.*t10.*2.0).*(1.0./2.0),t75,t112,t155,t202,t251,t305,t378,t75,I2+m3.*t80.*(1.0./2.0)+m4.*t80.*(1.0./2.0)+m5.*t80.*(1.0./2.0)+m6.*t80.*(1.0./2.0)+m7.*t80.*(1.0./2.0)+m8.*t80.*(1.0./2.0)+m2.*(d2.*t19.*t76.*2.0+d2.*t19.*t78.*2.0).*(1.0./2.0),t124,t166,t212,t260,t318,t384,t112,t124,I3+m4.*t129.*(1.0./2.0)+m5.*t129.*(1.0./2.0)+m6.*t129.*(1.0./2.0)+m7.*t129.*(1.0./2.0)+m8.*t129.*(1.0./2.0)+m3.*(d3.*t28.*t125.*2.0+d3.*t28.*t127.*2.0).*(1.0./2.0),t177,t222,t269,t331,t390,t155,t166,t177,I4+m5.*t182.*(1.0./2.0)+m6.*t182.*(1.0./2.0)+m7.*t182.*(1.0./2.0)+m8.*t182.*(1.0./2.0)+m4.*(d4.*t37.*t178.*2.0+d4.*t37.*t180.*2.0).*(1.0./2.0),t232,t278,t344,t396,t202,t212,t222,t232,I5+m6.*t237.*(1.0./2.0)+m7.*t237.*(1.0./2.0)+m8.*t237.*(1.0./2.0)+m5.*(d5.*t46.*t233.*2.0+d5.*t46.*t235.*2.0).*(1.0./2.0),t287,t357,t402,t251,t260,t269,t278,t287,I6+m7.*t292.*(1.0./2.0)+m8.*t292.*(1.0./2.0)+m6.*(d6.*t55.*t288.*2.0+d6.*t55.*t290.*2.0).*(1.0./2.0),t370,t408,t305,t318,t331,t344,t357,t370,I7+m8.*(L7.*t57.*t371.*2.0+L7.*t57.*t372.*2.0).*(1.0./2.0)+m7.*(d7.*t59.*t371.*2.0+d7.*t59.*t372.*2.0).*(1.0./2.0),t414,t378,t384,t390,t396,t402,t408,t414,I8+m8.*(d8.*t60.^2.*t61.*2.0+d8.*t61.*t62.^2.*2.0).*(1.0./2.0)],[8,8]);
