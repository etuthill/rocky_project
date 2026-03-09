% Rocky_5_closed_loop_poles_auto.m
clear
close all

syms s a b l g Kp Ki Jp Ji Ci

Hvtheta = -s/l/(s^2-g/l);
K = Kp + Ki/s;
M = a*b/(s+a);
J = Jp + Ji/s + Ci/s^2;
Mfb = M/(1+M*J);
Hcloop = 1/(1-Hvtheta*Mfb*K);

% system parameters
g = 9.81;
l = 0.4873;
a = 23.2135;
b = 0.003;

Hcloop_sub = subs(Hcloop);

best_cost = inf;
best_poles = [];
best_gains = [];
best_TF = [];

sigma_range = 3:0.5:5;
omega_range = 3:0.5:5;

p3_range = -8:-2:-12;
p4_range = -20:-5:-35;
p5_range = -25:-5:-45;

for sigma = sigma_range
for omega = omega_range
for p3 = p3_range
for p4 = p4_range
for p5 = p5_range

try

p1 = -sigma + 1i*omega;
p2 = -sigma - 1i*omega;

tgt_poly = expand((s-p1)*(s-p2)*(s-p3)*(s-p4)*(s-p5));

[n,d] = numden(Hcloop_sub);

coeffs_d = coeffs(d,s,'All');
coeffs_t = coeffs(tgt_poly,s,'All');

coeffs_d = coeffs_d/coeffs_d(1);
coeffs_t = coeffs_t/coeffs_t(1);

sol = solve(coeffs_d(2:6)==coeffs_t(2:6),Jp,Ji,Kp,Ki,Ci);

if isempty(sol)
continue
end

Kp_val = double(sol.Kp);
Ki_val = double(sol.Ki);
Ji_val = double(sol.Ji);
Jp_val = double(sol.Jp);
Ci_val = double(sol.Ci);

if ~isreal([Kp_val Ki_val Ji_val Jp_val Ci_val])
continue
end

TFexpr = subs(Hcloop,{Kp,Ki,Jp,Ji,Ci},{Kp_val,Ki_val,Jp_val,Ji_val,Ci_val});

TFstring = char(TFexpr);

s = tf('s');
TFH = eval(TFstring);

info = stepinfo(TFH);

Ts = info.SettlingTime;
Mp = info.Overshoot;

cost = Ts + 0.005*Mp;

if cost < best_cost

best_cost = cost;
best_poles = [p1 p2 p3 p4 p5];
best_gains = [Kp_val Ki_val Ji_val Jp_val Ci_val];
best_TF = TFH;

end

catch
end

end
end
end
end
end

disp('Best poles:')
disp(best_poles)

disp('Controller gains [Kp Ki Ji Jp Ci]:')
disp(best_gains)

figure
step(best_TF)
title('Best Step Response')

figure
impulse(best_TF)
title('Best Impulse Response')