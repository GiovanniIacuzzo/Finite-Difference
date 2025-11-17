function [t, u] = forward_solver(p, params)
% ===============================================================
% FORWARD_SOLVER
%
% Integra la EDO:
%     u'' + alpha u' + beta u = g(t)
%
% tramite differenze finite centrali,
% con inizializzazione di ordine 2.
%
% INPUT:
%   p      = parametri del controllo (dimensione M)
%   params = struttura con parametri del sistema e discretizzazione
%
% OUTPUT:
%   t      = tempi discreti (1 x N)
%   u      = soluzione numerica (1 x N)
%
% Dipendenze: build_control.m
% ===============================================================

    % --- Lettura parametri ---
    alpha = params.alpha;
    beta  = params.beta;
    dt    = params.dt;
    T     = params.T;
    u0    = params.u0;
    v0    = params.v0;

    % --- Griglia temporale ---
    t = 0:dt:T;
    N = length(t);

    % ===============================================================
    % Costruzione controllo g(t)
    % ===============================================================
    g = build_control(p, params);

    % Se il controllo non coincide in lunghezza, fai un reshape robusto
    if length(g) ~= N
        warning("build_control ha restituito g di lunghezza %d, atteso %d. Viene riadattato.", length(g), N);
        g = interp1(linspace(0,1,length(g)), g, linspace(0,1,N), 'linear');
    end

    % --- Array soluzione ---
    u = zeros(1, N);

    % ===============================================================
    % Inizializzazione:
    %
    % u(1) = u0
    % u(2) ≈ u0 + dt v0 + 0.5 dt^2 u''(0)
    %
    % dove u''(0) = g(0) - α u'(0) - β u(0)
    % ===============================================================

    u(1) = u0;

    % Derivata seconda iniziale dal modello
    udd0 = g(1) - alpha * v0 - beta * u0;

    % Passo iniziale di ordine 2
    u(2) = u0 + dt * v0 + 0.5 * dt^2 * udd0;

    % ===============================================================
    % Precomputazione coefficienti per velocizzare CPSO
    % ===============================================================
    A = (1/dt^2) + (alpha/(2*dt));     % coeff moltiplicativo di u(n+1)

    invA = 1 / A;  % per evitare divisioni ripetute

    c1 = 2/dt^2;
    c2 = 1/dt^2;
    c3 = alpha/(2*dt);

    % ===============================================================
    % Ciclo temporale esplicito centrato
    % ===============================================================
    for n = 2 : N-1
        % Costruzione termine destro (B)
        B = g(n) ...
            + c1*u(n) ...
            - c2*u(n-1) ...
            - c3*u(n-1) ...
            - beta*u(n);

        % Aggiornamento esplicito
        u(n+1) = B * invA;
    end

end
