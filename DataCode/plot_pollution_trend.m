% 初始化图
figure;
hold on;
years = [1990,1996,1999,2002,2005,2008];
% 循环绘制每种污染物
for k = 1:length(pollutants)-1
    % 对第1维（行业）求和，得到该污染物在每个年份的总排放量
    total_pollutant = squeeze(sum(poll(:, :, k), 1));  % 尺寸为 [6×1]
    
    % 绘图
    plot(years, total_pollutant, 'LineWidth', 1.8);
end

% 添加图例、标题、标签
legend(pollutants, 'Location', 'northeast');
xlabel('Year');
ylabel('Total Emissions');
title('Total Emissions by Pollutant Over Time');
grid on;
hold off;