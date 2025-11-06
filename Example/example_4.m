%% ESEMPIO 4: Equazione dell'onda 1D con differenze finite
% -------------------------------------------------------------------------
% u_tt = c^2 * u_xx
% Dominio: x ∈ [0,1], t > 0
% Condizioni: u(0,t)=u(1,t)=0
% Dati iniziali: u(x,0)=sin(pi*x), u_t(x,0)=0
% Soluzione esatta: u(x,t)=sin(pi*x)*cos(pi*c*t)
% -------------------------------------------------------------------------

clear; close all; clc;

% PARAMETRI
Nx = 200; Lx = 1; c = 1;
dx = Lx/(Nx-1);
x = linspace(0,Lx,Nx);
dt = 0.9 * dx / c;   % CFL < 1
Tmax = 2;
Nt = ceil(Tmax/dt);

fprintf('Griglia: Nx=%d, dx=%.3e, dt=%.3e, CFL=%.2f\n', Nx, dx, dt, c*dt/dx);

% CONDIZIONI INIZIALI
u0 = sin(pi*x);               % u(x,0)
u1 = u0;                      % u_t=0 -> prima iterazione identica
u_prev = u0;                  % u^{n-1}
u = u1;                       % u^{n}

% EVOLUZIONE TEMPORALE
for n = 1:Nt
    t = n*dt;
    % Aggiornamento con schema a due passi (centrato in tempo e spazio)
    u_next = 2*u - u_prev + (c*dt/dx)^2 * ([u(2:end) 0] - 2*u + [0 u(1:end-1)]);
    
    % Condizioni al contorno
    u_next(1) = 0; u_next(end) = 0;

    % Aggiorna
    u_prev = u;
    u = u_next;

    % Plot
    if mod(n,40)==0
        plot(x, u, 'b', 'LineWidth', 1.5); hold on;
        plot(x, sin(pi*x).*cos(pi*c*t), '--r', 'LineWidth', 1.2);
        title(sprintf('Equazione dell''onda: t=%.2f', t));
        xlabel('x'); ylabel('u(x,t)');
        legend('Numerica','Analitica','Location','best');
        grid on; drawnow; hold off;
    end
end

% ERRORE finale
u_exact = sin(pi*x).*cos(pi*c*Tmax);
err_L2 = sqrt(dx*sum((u - u_exact).^2));
fprintf('Errore L2 finale = %.3e\n', err_L2);
