function params = params()
%PARAMS  Definisce tutti i parametri del progetto di controllo ottimo.
%
% OUTPUT:
%   params : struttura contenente tutti i parametri del modello,
%            discretizzazione, controllo, target.
%
% STRUTTURA PRINCIPALE:
%   Sistema dinamico:
%       alpha, beta, u0, v0
%
%   Discretizzazione temporale:
%       T, dt, Nt, tspan
%
%   Controllo:
%       M                  -> numero intervalli del controllo (dimensione p)
%       control_type       -> 'piecewise' oppure 'spline'
%       lb, ub             -> limiti ammessi per p
%       lambda             -> peso regolarizzazione
%       smooth_spline      -> parametro smoothing
%
%   Target:
%       p_true             -> controllo "vero"
%       target_noise       -> aggiunta rumore al target
%
% ===============================================================

%% =======================
%   1) PARAMETRI SISTEMA
% ========================
params.alpha = 0.5;      % smorzamento
params.beta  = 2.0;      % rigidezza
params.u0    = 0.0;      % posizione iniziale
params.v0    = 0.0;      % velocità iniziale


%% ============================
%   2) DISCRETIZZAZIONE TEMPO
% ============================
params.T  = 10;          % tempo finale
params.dt = 0.01;        % passo temporale
params.Nt = floor(params.T/params.dt) + 1;
params.tspan = [0 params.T];
params.t = linspace(0, params.T, params.Nt);   % utile in molte parti


%% ==================
%   3) CONTROLLO
% ==================
params.M  = 20;                   % dimensione vettore p
params.control_type = 'spline';   % 'piecewise' o 'spline'

params.lb = -5;
params.ub = 5;
params.p_min = params.lb;
params.p_max = params.ub;

params.lambda = 0.01;    % regolarizzazione L2 sul controllo

params.smooth_spline = 1e-3;

%% =========================
%   5) CONTROLLO
% =========================
tt = linspace(0,1,params.M);
params.p_true = 2*sin(2*pi*tt) - cos(3*pi*tt);

params.p_true = min(max(params.p_true, params.lb), params.ub);

params.target_noise = 0.00;


%% ============================
%   6) VERIFICHE AUTOMATICHE
% ============================
if params.M < 5
    warning('M è molto piccolo: il controllo potrebbe risultare troppo rigido.');
end

if params.dt > 0.05
    warning('Il passo temporale dt è grande: la soluzione dell''ODE potrebbe perdere accuratezza.');
end

if params.lambda < 1e-4
    warning('Regolarizzazione lambda molto piccola: rischio oscillazioni di p.');
end

end
