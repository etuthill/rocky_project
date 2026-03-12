% Rocky_5_closed_loop_poles_auto.m
clear
close all
syms s a b l g Kp Ki Jp Ji Ci

% Pendulum and motor TFs
Hvtheta = -s/l/(s^2-g/l);       % Pendulum transfer function
K = Kp + Ki/s;                  % PI angle controller
M = a*b/(s+a);                  % Motor first-order model
J = Jp + Ji/s + Ci/s^2;         % Controller around motor
Mfb = M/(1 + M*J);              % Closed-loop motor TF
Hcloop = 1/(1 - Hvtheta*Mfb*K); % Closed-loop TF from disturbance

g = 9.81;
l = 0.4873;
a = 23.2135;
b = 0.003;

Hcloop_sub = subs(Hcloop);

sigma_range = 1.5:0.5:3;    % Real part of complex pair (moderate speed)
omega_range = 2:0.5:4;      % Imag part of complex pair

p3_range = -3:-1:-7;        % Real pole 1
p4_range = -5:-2:-12;       % Real pole 2
p5_range = -6:-2:-15;       % Real pole 3

%% Initialize best solution
best_cost = inf;
best_poles = [];
best_gains = [];
best_TF = [];

max_Ts_allowed = 2;   % maximum allowed settling time (seconds)
penalty_factor = 50;  % penalty if Ts > max_Ts_allowed

%% Loop over pole combinations
for sigma = sigma_range
    for omega = omega_range
        for p3 = p3_range
            for p4 = p4_range
                for p5 = p5_range
                    try
                        % Target poles
                        p1 = -sigma + 1i*omega;
                        p2 = -sigma - 1i*omega;

                        tgt_poly = expand((s-p1)*(s-p2)*(s-p3)*(s-p4)*(s-p5));

                        % Extract denominator of Hcloop
                        [n,d] = numden(Hcloop_sub);
                        coeffs_d = coeffs(d,s,'All');
                        coeffs_t = coeffs(tgt_poly,s,'All');

                        % Normalize coefficients
                        coeffs_d = coeffs_d/coeffs_d(1);
                        coeffs_t = coeffs_t/coeffs_t(1);

                        % Solve for controller gains
                        sol = solve(coeffs_d(2:6) == coeffs_t(2:6), Jp, Ji, Kp, Ki, Ci);

                        if isempty(sol)
                            continue
                        end

                        Kp_val = double(sol.Kp);
                        Ki_val = double(sol.Ki);
                        Ji_val = double(sol.Ji);
                        Jp_val = double(sol.Jp);
                        Ci_val = double(sol.Ci);

                        % Skip non-real solutions
                        if ~isreal([Kp_val Ki_val Ji_val Jp_val Ci_val])
                            continue
                        end

                        % Build numeric TF
                        TFexpr = subs(Hcloop, {Kp,Ki,Jp,Ji,Ci}, {Kp_val,Ki_val,Jp_val,Ji_val,Ci_val});
                        TFstring = char(TFexpr);
                        s = tf('s');
                        TFH = eval(TFstring);

                        % Step response info
                        info = stepinfo(TFH);
                        Ts = info.SettlingTime;
                        Mp = info.Overshoot;

                        % Gain penalty
                        gain_penalty = 1e-3*sum(abs([Kp_val Ki_val Ji_val Jp_val Ci_val]));

                        % Total cost with penalties
                        if Ts > max_Ts_allowed
                            cost = Ts + 0.005*Mp + penalty_factor*(Ts - max_Ts_allowed) + gain_penalty;
                        else
                            cost = Ts + 0.005*Mp + gain_penalty;
                        end

                        % Keep best solution
                        if cost < best_cost
                            best_cost = cost;
                            best_poles = [p1 p2 p3 p4 p5];
                            best_gains = [Kp_val Ki_val Ji_val Jp_val Ci_val];
                            best_TF = TFH;
                        end

                    catch
                        continue
                    end
                end
            end
        end
    end
end

%% Display results
disp('Best poles:')
disp(best_poles)

disp('Controller gains [Kp Ki Ji Jp Ci]:')
disp(best_gains)

%% Plot responses
figure
step(best_TF)
title('Best Step Response')

figure
impulse(best_TF)
title('Best Impulse Response')