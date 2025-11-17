function g_out = build_control(p, params)
%BUILD_CONTROL Costruisce il controllo g(t) a partire dal vettore p

    %% --- Controlli preliminari ---
    if ~isfield(params,'tspan') || ~isfield(params,'Nt')
        error('params deve contenere almeno tspan e Nt.');
    end
    if params.Nt < 2
        error('Nt deve essere almeno 2.');
    end

    % Tipo di controllo
    if isfield(params,'control_type')
        type = params.control_type;
    else
        type = 'piecewise';
    end

    % Griglia temporale
    t0 = params.tspan(1);
    tf = params.tspan(2);
    Nt = params.Nt;
    t = linspace(t0, tf, Nt);

    % Assicura limiti sui parametri
    if isfield(params,'p_min')
        p = max(p, params.p_min);
    end
    if isfield(params,'p_max')
        p = min(p, params.p_max);
    end

    %% --- Controllo singolo o popolazione ---
    if isvector(p)
        Np = length(p);
        % singolo vettore
        g_out = build_single_control(p(:), type, t, Np);
    else
        % matrice: ogni riga -> un controllo
        N_particles = size(p,1);
        g_out = zeros(N_particles, Nt);
        for k = 1:N_particles
            g_out(k,:) = build_single_control(p(k,:)', type, t, size(p,2));
        end
    end
end

%% --- Funzione interna per costruire g(t) da singolo vettore p ---
function g = build_single_control(p, type, t, Np)
    switch lower(type)
        case 'piecewise'
            g = zeros(1, length(t));
            idx = floor((t - t(1)) / (t(end) - t(1)) * Np) + 1;
            idx(idx > Np) = Np;
            g = p(idx)';

        case 'spline'
            if Np < 2
                g = p(1) * ones(1, length(t));
            else
                ti = linspace(t(1), t(end), Np);
                g = ppval(spline(ti, p), t);
            end

        otherwise
            error('Tipo di controllo non riconosciuto. Usa ''piecewise'' o ''spline''.');
    end
end
