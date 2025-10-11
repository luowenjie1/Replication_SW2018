%% ===================== Plot Implicit Pollution Taxes (1990–2008) =====================

% Define pollutants and years
pollutants = {'CO','NOX','PM10','PM25','SO2','VOC','CO2'};
years = 1990:2008;

% Create a new figure
figure('Position', [100, 100, 1000, 600]);
hold on; grid on; box on;

% Plot each pollutant’s implicit tax trend
for i = 1:length(pollutants)
    plot(years, implicit_tax(i, :), 'LineWidth', 2, 'DisplayName', pollutants{i});
end

% Add labels and legend
xlabel('Year', 'FontSize', 14);
ylabel('Implicit Pollution Tax', 'FontSize', 14);
title('Trends in Implicit Pollution Taxes (1990–2008)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside');
set(gca, 'FontSize', 12);

% Save figure
saveas(gcf, fullfile(figDir, 'implicit_tax_trends.png'));
