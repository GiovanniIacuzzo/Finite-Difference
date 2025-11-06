%% ESEMPIO 3: Equazione di Poisson 2D con analisi di convergenza
% -------------------------------------------------------------------------
% - (u_xx + u_yy) = f(x,y),  (x,y) in (0,1)x(0,1)
% u = 0 sui bordi
%
% f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
% => u(x,y) = sin(pi*x)*sin(pi*y)
% -------------------------------------------------------------------------

clear; close all; clc;

N_values = [10, 20, 40, 80];
err_L2 = zeros(size(N_values));
time_solve = zeros(size(N_values));

for k = 1:length(N_values)
    N = N_values(k);
    h = 1/(N+1);
    x = linspace(0,1,N+2)';
    [X,Y] = meshgrid(x,x);

    % Matrice Laplaciana 2D
    e = ones(N,1);
    L1 = spdiags([e -2*e e], -1:1, N, N) / h^2;
    I = speye(N);
    L2 = kron(I,L1) + kron(L1,I);

    % Sorgente e soluzione esatta
    f = 2*pi^2 * sin(pi*X(2:end-1,2:end-1)) .* sin(pi*Y(2:end-1,2:end-1));
    f_vec = reshape(f', N^2, 1);

    u_exact = sin(pi*X).*sin(pi*Y);

    % Risoluzione
    tic;
    u_vec = L2 \ (-f_vec);
    time_solve(k) = toc;

    % Ricostruzione matrice
    u_num = zeros(N+2,N+2);
    u_num(2:end-1,2:end-1) = reshape(u_vec, N, N)';

    % Errore L2
    err_L2(k) = sqrt(h^2 * sum((u_num(:) - u_exact(:)).^2));

    fprintf('N = %3d | h = %.3e | Errore L2 = %.3e | Tempo = %.3f s\n', ...
        N, h, err_L2(k), time_solve(k));
end

% Ordine di convergenza
rate = log(err_L2(1:end-1)./err_L2(2:end)) ./ log(2);
fprintf('\nTABELLA CONVERGENZA:\n');
fprintf('-------------------------------------------------------------\n');
fprintf('   N       h          Errore L2       Ordine      Tempo (s)\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:length(N_values)
    if k == 1
        fprintf('%4d   %.3e   %.3e     ----        %.3f\n', N_values(k), 1/(N_values(k)+1), err_L2(k), time_solve(k));
    else
        fprintf('%4d   %.3e   %.3e     %.2f        %.3f\n', N_values(k), 1/(N_values(k)+1), err_L2(k), rate(k-1), time_solve(k));
    end
end
fprintf('-------------------------------------------------------------\n');
fprintf('Ordine medio ≈ %.2f\n', mean(rate));

% Grafico log-log errore
figure;
loglog(1./(N_values+1), err_L2, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('h'); ylabel('Errore L_2');
title('Convergenza dell''errore L_2 per Poisson 2D');
