function J = objective_cpso(p_mat, params)
    % Accetta matrice p_mat (N_particles x Np)
    N_particles = size(p_mat,1);
    J = zeros(N_particles,1);

    for k = 1:N_particles
        p = p_mat(k,:);              % singolo vettore di controllo
        g = build_control(p, params); % costruisce controllo (piecewise o spline)
        [t, u] = forward_solver(p, params);
        diff = u - params.target_u;
        J_data = trapz(t, diff.^2);
        J_reg  = params.lambda * sum(p.^2);
        J(k) = J_data + J_reg;
    end
end
