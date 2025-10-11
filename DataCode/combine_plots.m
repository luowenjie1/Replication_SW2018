%% ===================== Combine All Pollutant Plots =====================

% Create a single large figure (2 rows × 4 columns)
combinedFig = figure('Visible', 'on', 'Position', [100, 100, 2000, 1000]); % ← 加大画布尺寸

% Define layout
tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact'); % ← 减少空白间距

% Loop through pollutant plots
for i = 1:length(pollutants)
    nexttile(i);
    img = imread(fullfile(figDir, [pollutants{i}, '.png']));
    imshow(img, 'Border', 'tight'); % ← 移除边框空白
    title(upper(pollutants{i}), 'FontWeight', 'bold', 'FontSize', 16);
end

% Add a main title
sgtitle('Pollution Emission Counterfactuals (1990–2008)', ...
    'FontSize', 22, 'FontWeight', 'bold');

% Save the combined figure
saveas(combinedFig, fullfile(figDir, 'all_pollutants_combined.png'));
