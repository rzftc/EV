clc; clear; close all;

%% 1. 加载数据
% 请确保此处的文件名与您实际生成的文件名一致
dataFileName = 'main_potential_5min_1000_8am.mat'; 

if exist(dataFileName, 'file')
    load(dataFileName);
    fprintf('数据文件 %s 加载成功。\n', dataFileName);
else
    error('错误：未找到文件 %s，请检查文件名或路径。', dataFileName);
end

%% 2. 数据提取与时间轴构建
% 提取核心变量
S_agg_data = results.S_agg;
lambda_data = results.lambda;

% 处理 NaN (将其视为0，防止线条断裂，视情况可选)
% S_agg_data(isnan(S_agg_data)) = 0; 

% 获取数据长度
num_steps = length(S_agg_data);

% 构建时间轴 (基于 8:00 开始，5分钟步长)
start_hour = 8;
dt_minutes = 5;
time_hours = start_hour + (0:num_steps-1) * (dt_minutes / 60);

%% 3. 绘图核心逻辑
% 创建图形窗口，设置背景为白色，尺寸适中 (影响输出图片的长宽比)
fig = figure('Color', 'w', 'Position', [200, 200, 1000, 500]);

hold on;

% 绘制 Lambda (控制信号) - 蓝色虚线
plot(time_hours, lambda_data, 'b--', 'LineWidth', 1.2); 

% 绘制 S_agg (聚合状态) - 红色实线
plot(time_hours, S_agg_data, 'r-', 'LineWidth', 1.5);

% 绘制 0 刻度参考线
yline(0, 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'Alpha', 0.5);

hold off;

%% 4. 图形修饰 (中文设置)
% 设置坐标轴范围
ylim([-1.2, 1.2]);
xlim([start_hour, time_hours(end)]);

% 设置图例 (字体放大)
lgd = legend('控制信号 \lambda', '电量偏差度 S_{agg}', ...
    'Location', 'best', 'FontSize', 16);
set(lgd, 'Interpreter', 'tex');

% 设置坐标轴标签 (字体放大)
xlabel('时间 (小时)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('归一化状态数值', 'FontSize', 20, 'FontWeight', 'bold');

% 开启网格
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
% 设置坐标轴刻度字体 (字体放大)
set(gca, 'FontSize', 16);

% 确保无标题
title('');

%% 5. 保存为高分辨率图片
outputFileName = '电量偏差度与控制信号对比图.png';
fprintf('正在保存图片到: %s ...\n', outputFileName);

% 方法 A: 使用 exportgraphics (MATLAB R2020a 及以上版本推荐，自动裁剪白边)
try
    exportgraphics(fig, outputFileName, 'Resolution', 300);
catch
    % 方法 B: 使用 print (兼容旧版本 MATLAB)
    % -dpng: PNG格式, -r300: 300 DPI
    print(fig, outputFileName, '-dpng', '-r300');
end

fprintf('图片保存成功！\n');