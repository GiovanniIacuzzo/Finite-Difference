function postprocess(t, u_best, u_target, p_best, params)
fprintf("\n=== POSTPROCESS: analisi risultati ===\n");

%% === 0) Ricostruisci controllo ottimale e controllo vero =======
g_best = build_control(p_best, params);
g_best = ensure_length(g_best, length(t));

if isfield(params, 'p_true')
    g_true = build_control(params.p_true, params);
    g_true = ensure_length(g_true, length(t));
else
    g_true = [];
end

%% === 1) Confronto soluzioni: target vs ottimizzata ========================
figure('Name','Confronto Soluzioni','NumberTitle','off');
subplot(2,1,1);
plot(t, u_target, 'k--', 'LineWidth', 2); hold on;
plot(t, u_best, 'r', 'LineWidth', 2);
ylabel('u(t)');
title('Soluzione Target vs Ottimizzata');
legend('Target','Ottimizzata');
grid on;

% Errore
subplot(2,1,2);
plot(t, u_best - u_target, 'm', 'LineWidth', 1.5);
xlabel('t'); ylabel('Errore');
title('Errore Assoluto nel Tempo');
grid on;

%% === 2) Errore relativo (percentuale) =====================================
err_rel = abs((u_best - u_target) ./ max(1e-12, abs(u_target)));

figure('Name','Errore Relativo','NumberTitle','off');
plot(t, err_rel*100, 'LineWidth', 1.7);
xlabel('t'); ylabel('Errore [%]');
title('Errore Relativo Percentuale');
grid on;

%% === 3) Controllo ottimale vs target ======================================
figure('Name','Controllo Recuperato','NumberTitle','off');

if ~isempty(g_true)
    plot(t, g_true, 'k--', 'LineWidth', 2); hold on;
end

plot(t, g_best, 'r', 'LineWidth', 1.8);

xlabel('t'); ylabel('g(t)');
legend({'Controllo Vero','Controllo Ottimale'}, ...
        'Location','best');
title('Confronto Controllo Target vs Ottimizzato');
grid on;

%% === 4) Spettro di Fourier dell'errore ====================================
err = u_best - u_target;
N = length(t);
dt = t(2) - t(1);
f = (0:N-1) / (N*dt);
E = abs(fft(err));

figure('Name','Spettro errore','NumberTitle','off');
plot(f(1:floor(N/2)), E(1:floor(N/2)), 'LineWidth', 1.5);
xlabel('Frequenza [Hz]'); ylabel('|FFT(errore)|');
title('Spettro in Frequenza dell Errore');
grid on;

%% === 5) Heatmap errore (opzionale ma molto utile) ==========================
figure('Name','Heatmap errore normalizzato','NumberTitle','off');
imagesc(t, 1, err); 
colormap hot; colorbar;
title('Errore nel tempo (heatmap)');
xlabel('t'); set(gca,'YTick',[]);

%% === 6) Tabella sintetica risultati =======================================
fprintf("\n--- Statistiche Errore ---\n");
fprintf("Min errore   = %.4e\n", min(err));
fprintf("Max errore   = %.4e\n", max(err));
fprintf("Mean abs err = %.4e\n", mean(abs(err)));
fprintf("RMSE         = %.4e\n", sqrt(mean(err.^2)));

if ~isempty(g_true)
    fprintf("\n--- Controllo ---\n");
    fprintf("Norma differenza controllo: %.4e\n", norm(g_best - g_true));
end

%% === 7) Salvataggio grafici ===============================================
savefig('img/post_confronto_soluzioni.fig');
savefig('img/post_controllo.fig');
savefig('img/post_errore_relativo.fig');
savefig('img/post_fft_errore.fig');

fprintf("\nGrafici salvati...\n");

end

%% --- Funzione ausiliaria per adattare lunghezza vettore controllo -------
function g_out = ensure_length(g_in, N)
    if length(g_in) ~= N
        g_out = interp1(linspace(0,1,length(g_in)), g_in, linspace(0,1,N), 'linear');
    else
        g_out = g_in;
    end
end
