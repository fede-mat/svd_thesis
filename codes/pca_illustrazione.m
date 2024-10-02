% Parametri della distribuzione gaussiana
mu = [0 0 0]; % Media
sigma = [1 0.5 0.2; 0.5 1 0.3; 0.2 0.3 1]; % Matrice di covarianza

% Generazione del dataset casuale distribuito gaussiano in R^3
rng(42); % Impostiamo il seed per la riproducibilità
n = 400; % Numero di punti nel dataset
data = mvnrnd(mu, sigma, n); % Generiamo dati casuali

% Calcolo della PCA
coeff = pca(data); % Calcoliamo i coefficienti della PCA

% Plot del dataset e dei vettori della PCA

scatter3(data(:,1), data(:,2), data(:,3), 'MarkerEdgeColor', 'k'); % Plot dei punti
hold on;
quiver3(mean(data(:,1)), mean(data(:,2)), mean(data(:,3)), ...
        coeff(1,1), coeff(2,1), coeff(3,1), 'r', 'LineWidth', 2); % Primo vettore della PCA
quiver3(mean(data(:,1)), mean(data(:,2)), mean(data(:,3)), ...
        coeff(1,2), coeff(2,2), coeff(3,2), 'g', 'LineWidth', 2); % Secondo vettore della PCA
quiver3(mean(data(:,1)), mean(data(:,2)), mean(data(:,3)), ...
        coeff(1,3), coeff(2,3), coeff(3,3), 'b', 'LineWidth', 2); % Terzo vettore della PCA
xlabel('X');
ylabel('Y');
zlabel('Z');
title('PCA con vettori');
grid on;
axis equal;
legend({'Dati', 'Componente 1', 'Componente 2', 'Componente 3'}, 'Location', 'best');
