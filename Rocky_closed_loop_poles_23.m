% Rocky_closed_loop_poles_23.m
%
% 1) Symbolically calculates closed loop transfer function of a disturbannce
% rejection PI control system for Rocky.
% No motor model (M =1). With motor model (1st order TF)
%
% 2) Specify location of (target)poles based on desired reponse. The number of
% poles = denominator polynomial of closed loop TF
%
% 3) Extract the closed loop denomiator poly and set = polynomial of target
% poles
%
% 4) Solve for Ki, Kp to match coefficients of polynomials. In general,
% this will be underdefined and will not be able to place poles in exact
% locations. In this, the control constants can be found exactly
%
% 5) Plot impulse and step response to see closed-loop behavior.
%
% based on code by SG. last modified 2/25/23 CL
clear all;
close all;
syms s a b l g Kp Ki % define symbolic variables
Hvtheta = -s/l/(s^2-g/l); % TF from velocity to angle of pendulum
K = Kp + Ki/s; % TF of the PI angle controller
M = a*b/(s+a) % TF of motor (1st order model)
% M = 1; % TF without motor
%
%
% closed loop transfer function from disturbance d(t)totheta(t)
Hcloop = 1/(1-Hvtheta*M*K) % use this for no motor feedback
pretty(simplify(Hcloop)) % to display the total transfer function
% Substitute parameters and solve
% system parameters
g = 9.81;
l = 22*2.54/100 %effective length
a = 14; %nominal motor parameters
b = 1/400; %nominal motor parameters
Hcloop_sub = subs(Hcloop) % sub parameter values into Hcloop
% specify locations of the target poles,
% choose # based on order of Htot denominator
% e.g., want some oscillations, want fast decay, etc.
p1 = -1+ i % dominant pole pair
p2 = -1 -i % dominant pole pair
p3 = -20
% target characteristic polynomial
% if motor model (TF) is added, order of polynomial will increases
tgt_char_poly = (s-p1)*(s-p2)*(s-p3)
npoly = 3
% get the denominator from Hcloop_sub
[n d] = numden(Hcloop_sub)
% find the coefficients of the denominator polynomial TF
coeffs_denom = coeffs(d, s)
% divide though the coefficient of the highest power term
coeffs_denom = coeffs(d, s)/(coeffs_denom(end))
% find coefficients of the target charecteristic polynomial
coeffs_tgt = coeffs(tgt_char_poly, s)
% solve the system of equations setting the coefficients of the
% polynomial in the target to the actual polynomials
solutions = solve(coeffs_denom(1:npoly-1) == coeffs_tgt(1:npoly-1), Kp, Ki)
% display the solutions as double precision numbers
Kp = double(solutions.Kp)
Ki = double(solutions.Ki)
% reorder coefficients for the check polynomial
for ii = 1:length(coeffs_denom)
chk_coeffs_denom(ii) = coeffs_denom(length(coeffs_denom) + 1 - ii);
end
closed_loop_poles = vpa (roots(subs(chk_coeffs_denom)), npoly )
% Plot impulse response of closed-loop system
TFstring = char(subs(Hcloop));
% Define 's' as transfer function variable
s = tf('s');
% Evaluate the expression
eval(['TFH = ',TFstring]);
figure (1)
impulse(TFH); %plot the impulse reponse
figure(2)
step(TFH) %plot the step response
syms s Kp Ki
Hcloop_sym = Hcloop_sub;
a = 23.2135;
pole_sets = [
    0.48 0.48 0.04
    0.47 0.47 0.06
    0.49 0.49 0.02
    0.48 0.48 0.05
    0.48 0.48 0.03
    0.47 0.47 0.07
    0.49 0.49 0.01
    0.48 0.48 0.06
];
npoly = 3;
Hcloop_tf = cell(size(pole_sets,1),1);

legend_labels = cell(size(pole_sets,1),1);

for k = 1:size(pole_sets,1)

    p1 = -pole_sets(k,1)*a;
    p2 = -pole_sets(k,2)*a;
    p3 = -pole_sets(k,3)*a;

    tgt_char_poly = expand((s - p1)*(s - p2)*(s - p3));
    coeffs_tgt = coeffs(tgt_char_poly, s);
    coeffs_tgt = coeffs_tgt / coeffs_tgt(end);

    [num_sym, den_sym] = numden(Hcloop_sym);
    coeffs_denom = coeffs(den_sym, s);
    coeffs_denom = coeffs_denom / coeffs_denom(end);

    sol = solve(coeffs_denom(1:npoly-1) == coeffs_tgt(1:npoly-1), [Kp, Ki]);
    Kp_val = double(sol.Kp);
    Ki_val = double(sol.Ki);

    Hcloop_num = subs(Hcloop_sym, {Kp, Ki}, {Kp_val, Ki_val});
    [num_sym2, den_sym2] = numden(Hcloop_num);

    num = double(coeffs(num_sym2, s));
    den = double(coeffs(den_sym2, s));

    Hcloop_tf{k} = tf(num, den);

    legend_labels{k} = sprintf("[%.3f %.3f %.3f]", pole_sets(k,1), pole_sets(k,2), pole_sets(k,3));

end

figure;
hold on;
for k = 1:size(pole_sets,1)
    step(Hcloop_tf{k});
end
hold off;
legend(legend_labels, "Location","best");
title("Step Responses for All Pole Sets");

figure;
hold on;
for k = 1:size(pole_sets,1)
    impulse(Hcloop_tf{k});
end
hold off;
legend(legend_labels, "Location","best");
title("Impulse Responses for All Pole Sets");

%%

alpha   = 23.2135;
beta_eff = 0.0030; 
ell_eff = 9.81/0.4873;
omega_n = sqrt(ell_eff);


p1 = -alpha/3;
p2 = -alpha/3;
p3 = -alpha/3;

Kp = (ell_eff/(alpha*beta_eff))*(alpha^2/3 + omega_n^2);
Ki = (ell_eff/(alpha*beta_eff))*(alpha^3/27 + alpha*omega_n^2);

Kp_val = double(Kp);
Ki_val = double(Ki);


tgt_char_poly = expand((s - p1)*(s - p2)*(s - p3));
den_coeffs = double(sym2poly(tgt_char_poly));
num_coeff = 1;

TFH = tf(num_coeff, den_coeffs);


t = 0:0.01:10;
figure(1); impulse(TFH, t);
figure(2); step(TFH, t);

Kp_val
Ki_val

disp("Desired poles:")
disp([p1 p2 p3])

disp("Actual poles from TFH:")
disp(pole(TFH))