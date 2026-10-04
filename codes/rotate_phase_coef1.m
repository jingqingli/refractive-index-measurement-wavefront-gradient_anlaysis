function coef_new = rotate_phase_coef1(theta, coef)
% ---------------------------------------------------------
% rotate_phase_coef
% 将二维二次相位多项式系数在坐标旋转 theta 后进行变换
%
% 输入:
% theta  - 旋转角 (rad)
% coef   - 原始系数 [a b c d e f]
%          phi = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f
%
% 输出:
% coef_new - 旋转后的系数 [a' b' c' d' e' f']
%
% ---------------------------------------------------------

% a = coef.p20;
% b = coef.p11;
% c = coef.p02;
d = coef.p10;
e = coef.p01;
f = coef.p00;

ct = cos(theta);
st = sin(theta);

%% 二次项
% a2 = a*ct^2 + b*ct*st + c*st^2;
% b2 = (c-a)*sin(2*theta) + b*cos(2*theta);
% c2 = a*st^2 - b*ct*st + c*ct^2;

%% 一次项
d2 = d*ct + e*st;
e2 = -d*st + e*ct;

%% 常数项
f2 = f;
a2=0;b2=0;c2=0;
coef_new = [a2 b2 c2 d2 e2 f2];

end