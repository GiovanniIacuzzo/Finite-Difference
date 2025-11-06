%% ESEMPIO 2: Equazione del calore 2D con differenze finite
% -------------------------------------------------------------------------
% u_t = alpha * (u_xx + u_yy)
% Dominio: (x,y) ∈ [0,1]^2
% Condizioni: Dirichlet omogenee (u = 0 sui bordi)
%
% Interpretazione: diffusione di una "macchia di calore" iniziale (gaussiana)
% che si espande e si attenua nel tempo.
% -------------------------------------------------------------------------

clear; close all; clc;

% PARAMETRI
Nx = 60; Ny = 60;
Lx = 1;  Ly = 1;
dx = Lx / (Nx-1);
dy = Ly / (Ny-1);
alpha = 0.01;

% Stabilità CFL per schema esplicito
dt_max = (min(dx,dy)^2) / (4*alpha);
dt = 0.25 * (min(dx,dy)^2) / alpha;
Tmax = 0.02;
Nt = ceil(Tmax / dt);

fprintf('Parametri numerici:\n');
fprintf('  dx = %.3e, dy = %.3e, dt = %.3e (dt_max = %.3e)\n', dx, dy, dt, dt_max);
fprintf('  Numero di passi temporali = %d\n\n', Nt);

% GRIGLIA
[x,y] = meshgrid(linspace(0,Lx,Nx), linspace(0,Ly,Ny));

% CONDIZIONE INIZIALE: picco gaussiano centrato in (0.5,0.5)
u = exp(-50*((x-0.5).^2 + (y-0.5).^2));

% SALVO ENERGIA INIZIALE (norma L2)
E0 = sum(u(:).^2) * dx * dy;

% LOOP TEMPORALE
for n = 1:Nt
    % Aggiornamento con schema esplicito
    u = u + alpha * dt * del2(u, dx, dy);
    
    % Condizioni al contorno (Dirichlet omogenee)
    u(1,:) = 0; u(end,:) = 0;
    u(:,1) = 0; u(:,end) = 0;

    % Ogni 50 passi: visualizzazione e statistiche
    if mod(n,50)==0 || n==Nt
        % Energia totale (misura di quanto "calore" rimane)
        E = sum(u(:).^2) * dx * dy;

        % Valore massimo (temperatura massima nel dominio)
        umax = max(u(:));

        fprintf('t = %.4f  -->  max(u) = %.4f,  Energia = %.4e\n', n*dt, umax, E);

        % Plot 2D
        surf(x,y,u);
        shading interp; view(2);
        title(sprintf('Tempo t = %.4f', n*dt));
        xlabel('x'); ylabel('y');
        colorbar;
        drawnow;
    end
end

fprintf('\nRiduzione energia totale: %.2f %%\n', 100*(1 - E/E0));
