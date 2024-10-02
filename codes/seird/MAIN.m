%% Dati
clc
clear all

% Parametri del dominio spaziale e temporale
L = 1; % Semiampiezza dell'intervallo spaziale (0, L)
N = 500; % Numero di sottointervalli spaziali
T = 30; % Estremo finale dell'intervallo temporale (0, T)
K = 120; % Numero di sottointervalli temporali

x=linspace(0,1,N+1);
x=x';
t=linspace(0,T,K+1);
t=t';

% Condizioni al contorno nulle
sa1 = @(t) 0.0 .* t; % Condizione al contorno per S a x=0
sa = @(t) 0.0 .* t;
sb = @(t) 0.0 .* t; % Condizione al contorno per S a x=L
ea1 = @(t) 0.0 .* t; % Condizione al contorno per E a x=0
ea = @(t) 0.0 .* t;
eb = @(t) 0.0 .* t; % Condizione al contorno per E a x=L
ia1 = @(t) 0.0 .* t; % Condizione al contorno per I a x=0
ia = @(t) 0.0 .* t;
ib = @(t) 0.0 .* t; % Condizione al contorno per I a x=L
ra1 = @(t) 0.0 .* t; % Condizione al contorno per R a x=0
ra = @(t) 0.0 .* t;
rb = @(t) 0.0 .* t; % Condizione al contorno per R a x=L
da1 = @(t) 0.0 .* t; % Condizione al contorno per D a x=0
da = @(t) 0.0 .* t;
db = @(t) 0.0 .* t; % Condizione al contorno per D a x=L

% Condizioni iniziali
u0_s = @(x) exp(-(x + 1).^4) + exp(-((x - 0.35).^2) / 10^-2) + 1/8 * (exp(-((x - 0.62).^4) / 10^-5) + exp(-((x - 0.52).^4) / 10^-5) + exp(-((x - 0.42).^4) / 10^-5)) + 1/4 * exp(-((x - 0.375).^4) / 10^-5); % Popolazione inizialmente sana
u0_e = @(x) 1/20.*exp(-((x-0.75).^4)./10^-5); % Nessuno esposto inizialmente
u0_i = @(x) 0.0.*x; % Condizione iniziale per I, infetti nel centro
u0_r = @(x) 0.0 .*x; % Nessuno recuperato inizialmente
u0_d = @(x) 0.0 .*x; % Nessun morto inizialmente

% Parametri del modello
beta_i = 0.375;
beta_e = 0.375;
alpha = 0.09375;
gamma_e = 0.125;
gamma_i = 0.03125;
delta = 0.0046875;
a_e = 0; 
n_pop = 1; % Popolazione normale, costante per semplicità
v_s = 3.75*10^-5;
v_e = 0.75*10^-3;
v_i = 0.75*10^-10;
v_r = 3.75*10^-5;

%% Problema con condizioni Neumann-Dirichet

[x_nd, t_nd, S_nd, E_nd, I_nd, R_nd, D_nd] = sistema_pde_ND(L, N, T, K, sa1, sb, ea1, eb, ia1, ib, ra1, rb, da1, db, ...
                                            u0_s, u0_e, u0_i, u0_r, u0_d, ...
                                            beta_i, beta_e, alpha, gamma_e, gamma_i, delta, ...
                                            a_e, n_pop, v_s, v_e, v_i, v_r);
% Visualizzazione delle soluzioni
figure;
subplot(3, 2, 1);
surf(t_nd, x_nd, S_nd);
title('S_{nd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('S_{nd}');

subplot(3, 2, 2);
surf(t_nd, x_nd, E_nd);
title('E_{nd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('E_{nd}');

subplot(3, 2, 3);
surf(t_nd, x_nd, I_nd);
title('I_{nd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('I_{nd}');

subplot(3, 2, 4);
surf(t_nd, x_nd, R_nd);
title('R_{nd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('R_{nd}');

subplot(3, 2, 5);
surf(t_nd, x_nd, D_nd);
title('D_{nd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('D_{nd}');

%% Problema con condizioni Dirichlet-Dirichlet
[x_dd, t_dd, S_dd, E_dd, I_dd, R_dd, D_dd] = sistema_pde_DD(L, N, T, K, sa, sb, ea, eb, ia, ib, ra, rb, da, db, ...
                                            u0_s, u0_e, u0_i, u0_r, u0_d, ...
                                            beta_i, beta_e, alpha, gamma_e, gamma_i, delta, ...
                                            a_e, n_pop, v_s, v_e, v_i, v_r);
% Visualizzazione delle soluzioni
figure;
subplot(2, 3, 1);
surf(t_dd, x_dd, S_dd);
title('S_{dd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('S_dd');

subplot(2, 3, 2);
surf(t_dd, x_dd, E_dd);
title('E_{dd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('E_dd');

subplot(2, 3, 3);
surf(t_dd, x_dd, I_dd);
title('I_{dd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('I_dd');

subplot(2, 3, 4);
surf(t_dd, x_dd, R_dd);
title('R_{dd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('R_dd');

subplot(2, 3, 5);
surf(t_dd, x_dd, D_dd);
title('D_{dd}');
xlabel('Tempo');
ylabel('Spazio');
zlabel('D_dd');

%% Confonto fra le due soluzioni

%% Consideriamo Ora la DMD
X=[S_nd;E_nd;I_nd;R_nd;D_nd];
X1=X(:,1:end-1);
X2=X(:,2:end);
r=10;
dt=t_nd(2)-t_nd(1);
[Phi,omega,lambda,b,Xdmd] = DMD(X1,X2,r,dt);
figure;
surf(Xdmd)
figure;
surf(abs(Xdmd-X(:,1:end-1)))

%% Confronto Grafico fra la Soluzione esatta e la DMD
figure;
%% Primo plot grafico 
% Rappresentazione della soluzione numerica confrontata con quella esatta
for k = 1:K-1
    % Subplot per S_nd
    subplot(2, 3, 1);
    plot(x_nd, S_nd(:,k), '-b')   % Linea blu
    hold on
    plot(x_nd, Xdmd(1:N+1, k), '--r')  % Linea rossa tratteggiata
    hold off
    legend('S_{nd}', 'S_{dmd}')
    xlabel('x')
    ylabel(sprintf('S_{nd}(:, t = %.2f)', t(k)))
    title('Confronto S_{nd} vs S_{dmd}')

    % Subplot per E_nd
    subplot(2, 3, 2);
    plot(x_nd, E_nd(:,k), '-b')   % Linea blu
    hold on
    plot(x_nd, Xdmd(N+2:2*(N+1), k), '--r')  % Linea rossa tratteggiata
    hold off
    legend('E_{nd}', 'E_{dmd}')
    xlabel('x')
    ylabel(sprintf('E_{nd}(:, t = %.2f)', t(k)))
    title('Confronto E_{nd} vs E_{dmd}')

    % Subplot per I_nd
    subplot(2, 3, 3);
    plot(x_nd, I_nd(:,k), '-b')   % Linea blu
    hold on
    plot(x_nd, Xdmd(2*(N+1)+1:3*(N+1), k), '--r')  % Linea rossa tratteggiata
    hold off
    legend('I_{nd}', 'I_{dmd}')
    xlabel('x')
    ylabel(sprintf('I_{nd}(:, t = %.2f)', t(k)))
    title('Confronto I_{nd} vs I_{dmd}')

    % Subplot per R_nd
    subplot(2, 3, 4);
    plot(x_nd, R_nd(:,k), '-b')   % Linea blu
    hold on
    plot(x_nd, Xdmd(3*(N+1)+1:4*(N+1), k), '--r')  % Linea rossa tratteggiata
    hold off
    legend('R_{nd}', 'R_{dmd}')
    xlabel('x')
    ylabel(sprintf('R_{nd}(:, t = %.2f)', t(k)))
    title('Confronto R_{nd} vs R_{dmd}')

    % Subplot per D_nd
    subplot(2, 3, 5);
    plot(x_nd, D_nd(:,k), '-b')   % Linea blu
    hold on
    plot(x_nd, Xdmd(4*(N+1)+1:5*(N+1), k), '--r')  % Linea rossa tratteggiata
    hold off
    legend('D_{nd}', 'D_{dmd}')
    xlabel('x')
    ylabel(sprintf('D_{nd}(:, t = %.2f)', t(k)))
    title('Confronto D_{nd} vs D_{dmd}')

    % Aggiornamento della figura
    drawnow
    pause(0.01);  % Aggiunge un breve ritardo per permettere l'aggiornamento dinamico
end


%%
figure(1)
%% secondo plot grafico

for k = 1:K-1
    % Prima immagine
    figure(1);
    
    % Subplot per S_nd
    subplot(1, 3, 1);
    plot(x_nd, S_nd(:,k), '-b')
    hold on
    plot(x_nd, Xdmd(1:N+1, k), '--r')  % Sostituisce gli * con una linea tratteggiata rossa
    hold off
    legend('S_{nd}', 'S_{dmd}')
    title('Confronto S_{nd} e S_{dmd}')
    xlabel('x')
    ylabel(sprintf('S_{nd}(:, t = %.2f)', t(k)))

    % Subplot per E_nd
    subplot(1, 3, 2);
    plot(x_nd, E_nd(:,k), '-b')
    hold on
    plot(x_nd, Xdmd(N+2:2*(N+1), k), '--r')  % Sostituisce gli * con una linea tratteggiata rossa
    hold off
    legend('E_{nd}', 'E_{dmd}')
    title('Confronto E_{nd} e E_{dmd}')
    xlabel('x')
    ylabel(sprintf('E_{nd}(:, t = %.2f)', t(k)))

    % Subplot per I_nd
    subplot(1, 3, 3);
    plot(x_nd, I_nd(:,k), '-b')
    hold on
    plot(x_nd, Xdmd(2*(N+1)+1:3*(N+1), k), '--r')  % Sostituisce gli * con una linea tratteggiata rossa
    hold off
    legend('I_{nd}', 'I_{dmd}')
    title('Confronto I_{nd} e I_{dmd}')
    xlabel('x')
    ylabel(sprintf('I_{nd}(:, t = %.2f)', t(k)))

    % Aggiornamento della figura
    drawnow
    pause(0.01);  % Aggiunge un breve ritardo per permettere l'aggiornamento dinamico
end

%% 
figure(2);
%% terzo plot garfico

for k = 1:K-1
    % Seconda immagine
    figure(2);
    
    % Subplot per R_nd
    subplot(1, 2, 1);
    plot(x_nd, R_nd(:,k), '-b')
    hold on
    plot(x_nd, Xdmd(3*(N+1)+1:4*(N+1), k), '--r')  % Sostituisce gli * con una linea tratteggiata rossa
    hold off
    legend('R_{nd}', 'R_{dmd}')
    title('Confronto R_{nd} e R_{dmd}')
    xlabel('x')
    ylabel(sprintf('R_{nd}(:, t = %.2f)', t(k)))

    % Subplot per D_nd
    subplot(1, 2, 2);
    plot(x_nd, D_nd(:,k), '-b')
    hold on
    plot(x_nd, Xdmd(4*(N+1)+1:5*(N+1), k), '--r')  % Sostituisce gli * con una linea tratteggiata rossa
    hold off
    legend('D_{nd}', 'D_{dmd}')
    title('Confronto D_{nd} e D_{dmd}')
    xlabel('x')
    ylabel(sprintf('D_{nd}(:, t = %.2f)', t(k)))

    % Aggiornamento della figura
    drawnow
    pause(0.01);  % Aggiunge un breve ritardo per permettere l'aggiornamento dinamico
end
