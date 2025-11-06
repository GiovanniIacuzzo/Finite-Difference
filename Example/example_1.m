%% ESEMPIO 1: Operatore di Laplace 1D con analisi di convergenza
% -------------------------------------------------------------------------
% -u_xx = f(x),  x in (0,1)
% u(0) = u(1) = 0
%
% f(x) = sin(pi*x)  ->  u(x) = (1/pi^2)*sin(pi*x)
% -------------------------------------------------------------------------

clear; close all; clc;

N_values = [25, 50, 100, 200, 400, 800];
err_L2 = zeros(size(N_values));

for k = 1:length(N_values)
    N = N_values(k);
    h = 1/(N+1);
    x = linspace(0,1,N+2)';

    % Laplaciano discreto con segno coerente: -u_xx
    e = ones(N,1);
    A = spdiags([-e 2*e -e], -1:1, N, N) / h^2;

    % Sorgente e soluzione esatta
    f = sin(pi*x(2:end-1));
    u_exact = (1/pi^2)*sin(pi*x);

    % Risoluzione sistema -u_xx = f -> A u = f
    u = A \ f;
    u_full = [0; u; 0];

    % Errore L2
    err_L2(k) = sqrt(h * sum((u_full - u_exact).^2));

    fprintf('N = %4d | h = %.3e | Errore L2 = %.3e\n', N, h, err_L2(k));
end

% Ordine di convergenza
rate = log(err_L2(1:end-1)./err_L2(2:end)) ./ log(2);
fprintf('\nTABELLA CONVERGENZA:\n');
fprintf('---------------------------------------------\n');
fprintf('   N       h          Errore L2      Ordine\n');
fprintf('---------------------------------------------\n');
for k = 1:length(N_values)
    if k == 1
        fprintf('%4d   %.3e   %.3e     ----\n', N_values(k), 1/(N_values(k)+1), err_L2(k));
    else
        fprintf('%4d   %.3e   %.3e     %.2f\n', N_values(k), 1/(N_values(k)+1), err_L2(k), rate(k-1));
    end
end
fprintf('---------------------------------------------\n');
fprintf('Ordine medio ≈ %.2f\n', mean(rate));

% Grafico log-log
figure;
loglog(1./(N_values+1), err_L2, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('h'); ylabel('Errore L_2');
title('Convergenza dell''errore L_2 per -u'''' = sin(\pi x)');