s = tf('s');

g_val = 9.81;
l_val = 0.4873;
omega_n = sqrt(g_val/l_val);

alpha = 23.2135;
beta  = 0.003;

frac_set = 0.05:0.05:0.95;

best_ts = inf;
best_TF = [];
best_Kp = 0;
best_Ki = 0;
best_fracs = [];

for f1 = frac_set
    for f2 = frac_set

        f3 = 1 - f1 - f2;
        if f3 <= 0
            continue
        end

        p1 = -alpha*f1;
        p2 = -alpha*f2;
        p3 = -alpha*f3;

        desired_poles = [p1 p2 p3];
        desired_char_poly = poly(desired_poles);

        A1 = desired_char_poly(3);
        A0 = desired_char_poly(4);

        Kp_val = (l_val/(alpha*beta))*(A1 + omega_n^2);
        Ki_val = (l_val/(alpha*beta))*(A0 + alpha*omega_n^2);

        Den = s^3 + alpha*s^2 + ...
              (Kp_val*alpha*beta/l_val - omega_n^2)*s + ...
              (Ki_val*alpha*beta/l_val - alpha*omega_n^2);

        Num = (s + alpha)*(s^2 - omega_n^2);

        TFH = Num/Den;

        ts = stepinfo(TFH).SettlingTime;

        if ts < best_ts
            best_ts = ts;
            best_TF = TFH;
            best_Kp = Kp_val;
            best_Ki = Ki_val;
            best_fracs = [f1 f2 f3];
        end
    end
end

disp("Best fractions:")
disp(best_fracs)

disp("Best settling time:")
disp(best_ts)

disp("Best Kp:")
disp(best_Kp)

disp("Best Ki:")
disp(best_Ki)

figure
step(best_TF)
title("Best Design - Step Response")

figure
impulse(best_TF)
title("Best Design - Impulse Response")

%%
disp("Closed-loop poles:")
disp(pole(best_TF))