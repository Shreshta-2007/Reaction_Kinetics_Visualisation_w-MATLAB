%% Constants and reactor basis
R  = 8.314;           % gas constant, J/(mol K)
V  = 2.0;             % reactor volume, L
n0 = 0.025;           % initial amount of N2O5, mol
C0 = n0 / V;          % initial concentration, mol/L

%% REAL DATA: first-order rate constants of N2O5 decomposition
T_data = [349 368 378 396];               % temperature, K
k_data = [0.0068 0.041 0.13 0.50];      % first-order rate constant, 1/s

t_half_25C_lit = 4.03e4;                  % half-life at 25 C, s
k_25C_lit      = log(2) / t_half_25C_lit;


%% PART 1: Decay curves at the measured temperatures
t = linspace(0, 300, 1000);             % time, s

figure;
subplot(1, 2, 1); hold on; grid on;
for i = 1:numel(T_data)
    plot(t, C0 * exp(-k_data(i) * t), 'LineWidth', 1.5);
end
yline(C0/2, '--k', 'C_0/2');
xlabel('Time (s)'); ylabel('Concentration of N_2O_5 (mol/L)');
title('(a) N_2O_5 decay, linear axes');
legend('349 K', '368 K', '378 K', '396 K', 'Location', 'northeast');

subplot(1, 2, 2); hold on; grid on;
for i = 1:numel(T_data)
    plot(t, C0 * exp(-k_data(i) * t), 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log'); ylim([1e-3 * C0, 1.5 * C0]);
xlabel('Time (s)'); ylabel('Concentration (log scale)');
title('(b) N_2O_5 decay, semi-log axes');
legend('349 K', '368 K', '378 K', '396 K', 'Location', 'southwest');

exportgraphics(gcf, 'concentration_vs_time.png', 'Resolution', 200)

fprintf('PART 1: measured rate constants and half-lives\n');
fprintf('  T (K)     k (1/s)    t_half (s)\n');
fprintf('  %5d    %8.4f    %9.2f\n', [T_data; k_data; log(2) ./ k_data]);


%% PART 2: Solve the ODE numerically and compare with the formula
k = k_data(1);                          % use the 349 K rate constant
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-12);
[t_num, C_num] = ode45(@(t, C) -k * C, [0 600], C0, opts);
C_exact = C0 * exp(-k * t_num);

figure;
subplot(1, 2, 1); hold on; grid on;
plot(t_num, C_exact, 'b-', 'LineWidth', 2);
plot(t_num, C_num, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Concentration of N_2O_5 (mol/L)');
title(sprintf('(a) ode45 vs analytical, T = %d K', T_data(1)));
legend('Analytical', 'Numerical (ode45)');

subplot(1, 2, 2);
semilogy(t_num, abs(C_num - C_exact), 'k', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Absolute Error (mol/L)');
title('(b) Absolute error of the numerical solution');

exportgraphics(gcf, 'concentration_vs_time.png', 'Resolution', 200)

fprintf('\nPART 2: max error of ode45 = %.2e mol/L\n', max(abs(C_num - C_exact)));


%% PART 3: Arrhenius analysis of the real data
p     = polyfit(1 ./ T_data, log(k_data), 1);   
Ea    = -p(1) * R;                              % activation energy, J/mol
A_pre = exp(p(2));                              % pre-exponential factor, 1/s
k_fun = @(T) A_pre * exp(-Ea ./ (R * T));       % fitted k(T)

ln_k_fit = polyval(p, 1 ./ T_data);
R2 = 1 - sum((log(k_data) - ln_k_fit).^2) / sum((log(k_data) - mean(log(k_data))).^2);

T_line = linspace(290, 400, 200);

figure;
subplot(1, 2, 1); hold on; grid on;
plot(1000 ./ T_data, log(k_data), 'ko', 'MarkerFaceColor', 'k');
plot(1000 ./ T_line, log(k_fun(T_line)), 'r-', 'LineWidth', 1.5);
plot(1000 / 298.15, log(k_25C_lit), 's', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
xlabel('1000 / T  (K^{-1})'); ylabel('ln k');
title(sprintf('(a) Arrhenius plot, E_a = %.1f kJ/mol, R^2 = %.4f', Ea / 1e3, R2));
legend('Measured data', 'Arrhenius fit', 'Literature at 25 C', 'Location', 'northeast');

subplot(1, 2, 2);
semilogy(T_line, log(2) ./ k_fun(T_line) / 3600, 'r-', 'LineWidth', 1.5); hold on; grid on;
semilogy(T_data, log(2) ./ k_data / 3600, 'ko', 'MarkerFaceColor', 'k');
xlabel('Temperature (K)'); ylabel('Half-life (hours, log scale)');
title('(b) Half-life from the Arrhenius fit');
legend('Fitted model', 'Measured data', 'Location', 'northeast');

exportgraphics(gcf, 'concentration_vs_time.png', 'Resolution', 200)

k_val_25 = k_fun(298.15);
if k_val_25 < k_25C_lit
    comp_text = 'below';
else
    comp_text = 'above';
end

fprintf('\nPART 3: Arrhenius fit of the real data\n');
fprintf('  Ea = %.1f kJ/mol,  A = %.2e 1/s,  R^2 = %.4f\n', Ea / 1e3, A_pre, R2);
fprintf('  Check at 25 C (298.15 K):\n');
fprintf('    model:      k = %.2e 1/s, t_half = %.2f h\n', k_val_25, log(2) / k_val_25 / 3600);
fprintf('    literature: k = %.2e 1/s, t_half = %.2f h\n', k_25C_lit, t_half_25C_lit / 3600);
fprintf('    model is %.0f%% %s the literature rate constant\n', ...
        100 * abs(k_val_25 / k_25C_lit - 1), comp_text);


%% PART 4: Thermodynamics of N2O5(g) -> 2 NO2(g) + 1/2 O2(g)
dHf_N2O5 = 11.3;     S_N2O5 = 355.7;     % kJ/mol, J/(mol K)
dHf_NO2  = 33.18;    S_NO2  = 240.06;    
dHf_O2   = 0;        S_O2   = 205.14;    

dH = (2 * dHf_NO2 + 0.5 * dHf_O2 - dHf_N2O5) * 1e3;            % J/mol
dS =  2 * S_NO2 + 0.5 * S_O2 - S_N2O5;                        % J/(mol K)

T_th = linspace(290, 400, 120);
dG   = dH - T_th * dS;                      
K    = exp(-dG ./ (R * T_th));              

% Use a standard for-loop instead of arrayfun for equilibrium conversion
X_1  = zeros(size(K));
X_01 = zeros(size(K));
for i = 1:numel(K)
    X_1(i)  = alpha_eq(K(i), 1.0);
    X_01(i) = alpha_eq(K(i), 0.1);
end

t4 = linspace(0, 300, 600);

figure;
subplot(1, 3, 1); hold on; grid on;
plot(T_th, dG / 1e3, 'b-', 'LineWidth', 1.5);
yline(0, '--k');
xlabel('Temperature (K)'); ylabel('\DeltaG^\circ (kJ/mol)');
title('(a) Gibbs energy');

subplot(1, 3, 2); hold on; grid on;
semilogy(T_th, 1 - X_1,  'r-', 'LineWidth', 1.5);
semilogy(T_th, 1 - X_01, 'b--', 'LineWidth', 1.5);
xlabel('Temperature (K)'); ylabel('N_2O_5 left at equilibrium');
title('(b) Equilibrium conversion');
legend('P = 1 bar', 'P = 0.1 bar', 'Location', 'northeast');

subplot(1, 3, 3); hold on; grid on;
for i = 1:numel(T_data)
    plot(t4, n0 * dH * (1 - exp(-k_data(i) * t4)), 'LineWidth', 1.5);
end
xlabel('Time (s)'); ylabel('Heat absorbed (J)');
title('(c) Heat absorbed');
legend('349 K', '368 K', '378 K', '396 K', 'Location', 'southeast');

exportgraphics(gcf, 'concentration_vs_time.png', 'Resolution', 200)

fprintf('\nPART 4: thermodynamics (per mol N2O5)\n');
fprintf('  dH = %.2f kJ/mol,  dS = %.2f J/(mol K)\n', dH / 1e3, dS);
fprintf('  T (K)    dG (kJ/mol)        K     X_eq(1 bar)  X_eq(0.1 bar)\n');
for T_val = [298.15 323.15 348.15 373.15 398.15]
    dGi = dH - T_val * dS;  
    Ki = exp(-dGi / (R * T_val));
    fprintf('  %6.2f   %9.2f   %10.1f   %9.4f   %10.4f\n', T_val, dGi / 1e3, Ki, ...
            alpha_eq(Ki, 1.0), alpha_eq(Ki, 0.1));
end


%% Local function for equilibrium extent
function a = alpha_eq(K_val, P)
    f = @(a) (2*a).^2 .* (a/2).^0.5 ./ ((1 - a) .* (1 + 1.5*a).^0.5) * sqrt(P) - K_val;
    a = fzero(f, [1e-12, 1 - 1e-12]);
end