%% ==============================================================
%%                        Main script-CPSO                       
%% ==============================================================

clear; clc; close all;

%% 1) Carica parametri
params = params();
fprintf("=== Parametri caricati ===\n");
disp(params);

%% 2) Generazione soluzione target sintetica
fprintf("\nGenerazione target sintetico...\n");
p_true = params.p_true;
[t, u_target] = forward_solver(p_true, params);
params.target_u = u_target;

fprintf("Soluzione target generata. Min/Max u(t): %.3f / %.3f\n", ...
    min(u_target), max(u_target));

figure; 
plot(t, u_target, 'LineWidth', 2);
title('Soluzione target sintetica'); xlabel('t'); ylabel('u_{target}(t)');
grid on;

%% 3) Setup CPSO
nvars = params.M;
lb = params.lb * ones(1, nvars);
ub = params.ub * ones(1, nvars);

% Funzione obiettivo
objfun = @(p_mat) objective_cpso(p_mat, params);

options = struct();
options.particles = 100;
options.max_velocity = 0.5 * (ub - lb);
options.verbose = 1;
options.early_stopping_delta = 1e-6;
options.early_stopping_patience = 20;
options.Cognitive_constant = 2.05;
options.Social_constant = 2.05;
options.n_jobs = 1;

%% 4) Creazione CPSO
cps = CPSO(objfun, nvars, lb, ub, options);

%% 5) Esecuzione CPSO
[max_iters] = 200;
[J_best, p_best, J_history] = cps.run(max_iters);

fprintf("\n=== CPSO concluso ===\n");
fprintf("Miglior valore J = %.6f\n", J_best);
disp("p_best ottenuto:");
disp(p_best);

%% 6) Simulazione con controllo ottimale
[t, u_best] = forward_solver(p_best, params);

%% 7) Postprocessing centralizzato
postprocess(t, u_best, u_target, p_best, params);

%% 8) Salvataggio risultati principali
save('results/results_cpso.mat', ...
    'p_best','J_best','u_best','u_target','t','params','J_history');

fprintf("\nRisultati salvati in results_cpso.mat\n");
