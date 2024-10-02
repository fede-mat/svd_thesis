function compress_image(image_path, k_vector)
    img = imread(image_path);
    
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    img_double = double(img);
    [U, S, V] = svd(img_double);

    num_plots = length(k_vector) + 1;  % numero totale di subplot (immagine originale + immagini compresse)
    rows = 2;  % numero di righe
    cols = ceil(num_plots / rows);  % numero di colonne

    figure;
    set(gcf, 'Position', [100, 100, 1500, 600]);  % dimensione della finestra della figura
    subplot(rows, cols, 1);
    imshow(img);
    title('Immagine Originale', 'FontSize', 14);

    for i = 1:length(k_vector)
        k = k_vector(i);
        U_k = U(:, 1:k);
        S_k = S(1:k, 1:k);
        V_k = V(:, 1:k);
        img_compressed = U_k * S_k * V_k';
        img_compressed = uint8(img_compressed);
        
        subplot(rows, cols, i + 1);
        imshow(img_compressed);
        title(['k = ', num2str(k)], 'FontSize', 14);
    end

    % Aggiusta gli spazi tra i plot
    set(gcf, 'Color', 'w');  % Sfondo bianco per la figur
end