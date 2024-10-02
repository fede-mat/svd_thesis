function [x, t, S, E, I, R, D] = sistema_pde_ND(L, N, T, K, sa1, sb, ea1, eb, ia1, ib, ra1, rb, da1, db, ...
                                            u0_s, u0_e, u0_i, u0_r, u0_d, ...
                                            beta_i, beta_e, alpha, gamma_e, gamma_i, delta, ...
                                            a_e, n_pop, v_s, v_e, v_i, v_r);
% ---- Risoluzione di un sistema di PDE ----
% Sistema di PDE con condizioni al bordo di Dirichlet
% e condizioni iniziali
% -----------------------------------------------
% Sintassi:
% [x, t, S, E, I, R, D] = sistema_pde_ND(L, N, T, K, sa1, sb, ea1, eb, ia1, ib, ra1, rb, da1, db, ...
%                                             u0_s, u0_e, u0_i, u0_r, u0_d, ...
%                                             beta_i, beta_e, alpha, gamma_e, gamma_i, delta, ...
%                                             a_e, n_pop, v_s, v_e, v_i, v_r)
%
% Input:
%   L semiampiezza intervallo spaziale (0,L)
%   N numero di sottointervalli in (0,L)
%   T estremo finale intervallo temporale (0,T)
%   K numero di sottointervalli in (0,T)
%   sa1, ea1, ia1, ra1, da1 sono le funzioni che descrivono la condizione
%   di Neumann nel primo estremo
%   sb, eb, ib, rb, db sono le funzioni che descrivono la condizione di Dirichlet al 
%   secondo estremo
%   u0_s, u0_e, u0_i, u0_r, u0_d funzioni che descrivono le condizioni iniziali
%   beta_i, beta_e, alpha, gamma_e, gamma_i, delta, a_e, n_pop, v_s, v_e, v_i, v_r
%   sono tutti parametri delle PDE. 
% Output:
%   x vettore dei nodi spaziali
%   t vettore dei nodi temporali
%   S, E, I, R, D soluzioni numeriche del sistema

% Calcolo passo di discretizzazione in spazio e tempo
h = L / N;
tau = T / K;

% Inizializzazione dei vettori t e x
t = linspace(0, T, K+1)';
x = linspace(0, L, N+1)';

% Inizializzazione delle matrici soluzione
S = zeros(N+1, K+1);
E = zeros(N+1, K+1);
I = zeros(N+1, K+1);
R = zeros(N+1, K+1);
D = zeros(N+1, K+1);

% Condizioni iniziali
S(:, 1) = u0_s(x);
E(:, 1) = u0_e(x);
I(:, 1) = u0_i(x);
R(:, 1) = u0_r(x);
D(:, 1) = u0_d(x);

% Condizioni al contorno
S(end, :) = sb(t);
E(end, :) = eb(t);
I(end, :) = ib(t);
R(end, :) = rb(t);
D(end, :) = db(t);

% Costruzione delle matrici A
e = ones(N, 1);
A_s = spdiags([-e 2*e -e], [-1 0 1], N, N) * (tau / h^2) * v_s*n_pop;
A_e = spdiags([-e 2*e -e], [-1 0 1], N, N) * (tau / h^2) * v_e*n_pop;
A_i = spdiags([-e 2*e -e], [-1 0 1], N, N) * (tau / h^2) * v_i*n_pop;
A_r = spdiags([-e 2*e -e], [-1 0 1], N, N) * (tau / h^2) * v_r*n_pop;
A_s(1,1)=2*tau/h^2*v_s*n_pop;
A_s(1,2)=-2*tau/h^2*v_s*n_pop;
A_e(1,1)=2*tau/h^2*v_e*n_pop;
A_e(1,2)=-2*tau/h^2*v_e*n_pop;
A_r(1,1)=2*tau/h^2*v_r*n_pop;
A_r(1,2)=-2*tau/h^2*v_r*n_pop;
A_i(1,1)=2*tau/h^2*v_i*n_pop;
A_i(1,2)=-2*tau/h^2*v_i*n_pop;

% Ciclo iterativo
for k = 1:K
    % Termini noti e correzioni per le condizioni al bordo
    s = S(1:end-1, k);
    e = E(1:end-1, k);
    i = I(1:end-1, k);
    r = R(1:end-1, k);
    F_s = -beta_i * (1 - a_e / n_pop) .* s .* i - beta_e * (1 - a_e / n_pop) .* s .* e;
    F_e = beta_i * (1 - a_e / n_pop) .* s .* i + beta_e * (1 - a_e / n_pop) .* s .* e - (alpha + gamma_e) * e;
    F_i = alpha * e - (gamma_i + delta) * i;
    F_r = gamma_e * e + gamma_i*i;
    F_d = delta * i;
    F_s(1)=F_s(1)+2*n_pop*v_s/h*sa1(t(k+1)); %N
    F_s(end)=F_s(end)+1*n_pop*v_s/h^2*sb(t(k+1)); %D
    F_e(1)=F_e(1)+2*n_pop*v_e/h*ea1(t(k+1));
    F_e(end)=F_e(end)+1*n_pop*v_e/h^2*eb(t(k+1));
    F_i(1)=F_i(1)+2*n_pop*v_i/h*ia1(t(k+1));
    F_i(end)=F_i(end)+1*n_pop*v_i/h^2*ib(t(k+1));
    F_r(1)=F_r(1)+2*n_pop*v_r/h*ra1(t(k+1));
    F_r(end)=F_r(end)+1*n_pop*v_r/h^2*rb(t(k+1));
    F_d(1)=F_d(1);
    F_d(end)=F_d(end);

    % Aggiornamento delle soluzioni
    S(1:end-1, k+1) = (eye(N) + A_s) \ (s + tau * F_s);
    E(1:end-1, k+1) = (eye(N) + A_e) \ (e + tau * F_e);
    I(1:end-1, k+1) = (eye(N) + A_i) \ (i + tau * F_i);
    R(1:end-1, k+1) = (eye(N) + A_r) \ (r + tau * F_r);
    D(1:end-1, k+1) = D(1:end-1, k) + tau * F_d;
end
