%% plot_future_potential.m
% 功能：读取各DER基线潜力，根据2030/2035/2040年倍率推演并绘图
% 输出：高DPI的PNG和EMF文件（无标题，中文文件名）

clear; close all; clc;

%% 1. 基础设置
% 绘图字体设置 (防止中文乱码)
set(0, 'DefaultAxesFontName', 'Microsoft YaHei');
set(0, 'DefaultTextFontName', 'Microsoft YaHei');
set(0, 'DefaultLegendFontName', 'Microsoft YaHei');

% 资源文件列表 (与 verify_Joint_Potential_Accuracy.m 保持一致)
der_files = {
    % 'DER1.mat', 'EV';
    'DER2.mat', 'AC';
    'DER3.mat', 'PV';
    'DER4.mat', 'P2G';
    'DER5.mat', 'IDC';
    % 'DER6.mat', 'Industry'
};

%% 2. 读取并聚合当前基线数据
fprintf('正在读取并聚合各资源基线数据...\n');

% 初始化
min_len = inf;
data_list = struct();

% 第一次遍历：确定最小长度并加载数据
for i = 1:size(der_files, 1)
    filename = der_files{i, 1};
    label = der_files{i, 2};
    
    if ~exist(filename, 'file')
        error('未找到文件: %s。请先运行 test_DER_x.m 系列脚本生成数据。', filename);
    end
    
    tmp = load(filename);
    if ~isfield(tmp, 'results')
        error('%s 中未找到 results 结构体', filename);
    end
    res = tmp.results;
    
    % 提取真实潜力 (True/Individual Sum) 作为基准
    switch label
        case 'EV'
            u_t = res.EV_Up_Individual_Sum(:);
            d_t = res.EV_Down_Individual_Sum(:);
        case 'AC'
            if isfield(res, 'Agg_P_Potential_Up_History')
                u_t = res.Agg_P_Potential_Up_History(:);
                d_t = res.Agg_P_Potential_Down_History(:);
            else
                error('DER2 (AC) 数据缺失 Agg_P_Potential_Up_History');
            end
        case 'PV'
            u_t = res.PV_Up(:);
            d_t = res.PV_Down(:);
        case 'P2G'
            u_t = res.P2G_Up(:);
            d_t = res.P2G_Down(:);
        case 'IDC'
            u_t = res.IDC_Up(:);
            d_t = res.IDC_Down(:);
        case 'Industry'
            u_t = res.Ind_Up(:);
            d_t = res.Ind_Down(:);
        otherwise
            u_t = []; d_t = [];
    end
    
    data_list(i).up = u_t;
    data_list(i).down = d_t;
    min_len = min(min_len, length(u_t));
end

% 截断对齐并求和
dt_short = 5; % 分钟
Total_Up_Base = zeros(min_len, 1);
Total_Down_Base = zeros(min_len, 1);

for i = 1:length(data_list)
    Total_Up_Base = Total_Up_Base + data_list(i).up(1:min_len);
    Total_Down_Base = Total_Down_Base + data_list(i).down(1:min_len);
end

% 生成时间轴 (小时)
time_axis = (0:min_len-1)' * dt_short / 60; 
% 假设从6点开始 (与 DER1/EV 保持一致，若不同DER起始时间不同需额外对齐，此处假设一致)
start_hour = 6; 
time_axis = time_axis + start_hour;

fprintf('基线数据聚合完成。时间步数: %d\n', min_len);

%% 3. 未来场景推演
% 2030年: 100倍
scale_2030 = 300;
Up_2030 = Total_Up_Base * scale_2030;
Down_2030 = Total_Down_Base * scale_2030;

% 2035年: 120倍
scale_2035 = 360;
Up_2035 = Total_Up_Base * scale_2035;
Down_2035 = Total_Down_Base * scale_2035;

% 2040年: 135倍
scale_2040 = 405;
Up_2040 = Total_Up_Base * scale_2040;
Down_2040 = Total_Down_Base * scale_2040;

%% 4. 绘图与保存

% 定义绘图通用参数
x_ticks = start_hour : 6 : (start_hour + 24); % 每6小时一个刻度
line_width = 1.5;
% 颜色配置 (红/绿/蓝 或 自定义)
c1 = [0.85, 0.33, 0.10]; % 2030
c2 = [0.93, 0.69, 0.13]; % 2035
c3 = [0, 0.45, 0.74];    % 2040

% === 图1：联合上调潜力 ===
fig1 = figure('Color', 'w', 'Position', [100, 100, 1000, 600], 'Visible', 'on');
hold on;
plot(time_axis, Up_2030, 'LineStyle', '-', 'LineWidth', line_width, 'Color', c1, 'DisplayName', '2030年');
plot(time_axis, Up_2035, 'LineStyle', '--', 'LineWidth', line_width, 'Color', c2, 'DisplayName', '2035年');
plot(time_axis, Up_2040, 'LineStyle', '-.', 'LineWidth', line_width, 'Color', c3, 'DisplayName', '2040年');

% 设置坐标轴
xlabel('时间 (小时)', 'FontSize', 14);
ylabel('上调潜力 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'XTick', x_ticks);
xlim([time_axis(1), time_axis(end)]);
% 不设置 Title

% 保存图1
filename_up = '联合上调潜力预测_2030_2040';
fprintf('正在保存: %s ...\n', filename_up);
print(fig1, [filename_up '.png'], '-dpng', '-r600'); % 600 DPI PNG
print(fig1, [filename_up '.emf'], '-dmeta', '-r600'); % EMF

% === 图2：联合下调潜力 ===
fig2 = figure('Color', 'w', 'Position', [150, 150, 1000, 600], 'Visible', 'on');
hold on;
plot(time_axis, Down_2030, 'LineStyle', '-', 'LineWidth', line_width, 'Color', c1, 'DisplayName', '2030年');
plot(time_axis, Down_2035, 'LineStyle', '--', 'LineWidth', line_width, 'Color', c2, 'DisplayName', '2035年');
plot(time_axis, Down_2040, 'LineStyle', '-.', 'LineWidth', line_width, 'Color', c3, 'DisplayName', '2040年');

% 设置坐标轴
xlabel('时间 (小时)', 'FontSize', 14);
ylabel('下调潜力 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'XTick', x_ticks);
xlim([time_axis(1), time_axis(end)]);
% 不设置 Title

% 保存图2
filename_down = '联合下调潜力预测_2030_2040';
fprintf('正在保存: %s ...\n', filename_down);
print(fig2, [filename_down '.png'], '-dpng', '-r600'); % 600 DPI PNG
print(fig2, [filename_down '.emf'], '-dmeta', '-r600'); % EMF

fprintf('所有绘图已完成并保存。\n');