%% ===================== Plot Competition (1990–2008) =====================

% Define years
years = 1990:2008;

Gamma_hat_foreign_2D = squeeze(Gamma_hat_foreign(1, :, :, :));
Gamma_hat_foreign_2D = squeeze(mean(Gamma_hat_foreign_2D, 1));
Gamma_hat_US_2D = squeeze(mean(Gamma_hat_US, 1));
Gamma_hat_US_2D = squeeze(mean(Gamma_hat_US_2D, 1));

% Plot each competitiveness trend
clf;
figure('Position', [100, 100, 1000, 600]);
hold on; grid on; box on;
plot(years, Gamma_hat_foreign_2D(1,:), 'LineWidth', 2, 'DisplayName', 'Changes in Foreign to Foreign Competitiveness');
plot(years, Gamma_hat_foreign_2D(2,:), 'LineWidth', 2, 'DisplayName', 'Changes in Foreign to US Competitiveness');
xlabel('Year', 'FontSize', 14);
ylabel('Changes', 'FontSize', 14);
title('Foreign Competitiveness (1990–2008)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest');
set(gca, 'FontSize', 14);
saveas(gcf, fullfile(figDir, 'foreign_competitive.png'));

clf;
figure('Position', [100, 100, 1000, 600]);
hold on; grid on; box on;
plot(years, Gamma_hat_US_2D(1,:), 'LineWidth', 2, 'DisplayName', 'Changes in US to Foreign Competitiveness');
plot(years, Gamma_hat_US_2D(2,:), 'LineWidth', 2, 'DisplayName', 'Changes in US to US Competitiveness');
xlabel('Year', 'FontSize', 14);
ylabel('Changes', 'FontSize', 14);
title('US Competitiveness (1990–2008)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest');
set(gca, 'FontSize', 14);
saveas(gcf, fullfile(figDir, 'US_competitive.png'));

close all;