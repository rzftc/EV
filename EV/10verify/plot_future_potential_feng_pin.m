%% plot_future_potential.m
% 功能：读取各DER基线潜力，根据2030/2035/2040年倍率推演并绘图
%      分别绘制2030/2035/2040年 上调/下调 在调频和调峰场景下的对比图 (共6张)
% 输出：高DPI的PNG和EMF文件（无标题，中文文件名）

clear; close all; clc;

%% 1. 基础设置
% 绘图字体设置 (防止中文乱码)
set(0, 'DefaultAxesFontName', 'Microsoft YaHei');
set(0, 'DefaultTextFontName', 'Microsoft YaHei');
set(0, 'DefaultLegendFontName', 'Microsoft YaHei');

% 资源文件列表 (与 verify_Joint_Potential_Accuracy.m 保持一致)
der_files = {
    'DER1.mat', 'EV';
    'DER2.mat', 'AC';
    'DER3.mat', 'PV';
    'DER4.mat', 'P2G';
    'DER5.mat', 'IDC';
    'DER6.mat', 'Industry'
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

%% 3. 未来场景推演 (计算调频与调峰)
% 设置随机种子以保证波动可复现
rng(2024);

% 调峰场景参数
ratio_peak = 0.8;      % 调峰潜力是调频的0.8倍
noise_level = 0.05;    % 随机波动幅度 (5%)

% --- 2030年 (100倍) ---
scale_2030 = 100;
Up_2030_Freq = Total_Up_Base * scale_2030;
Down_2030_Freq = Total_Down_Base * scale_2030;
% 2030年调峰 (加波动)
Up_2030_Peak = Up_2030_Freq * ratio_peak + (Up_2030_Freq * ratio_peak * noise_level) .* randn(size(Up_2030_Freq));
Down_2030_Peak = Down_2030_Freq * ratio_peak + (Down_2030_Freq * ratio_peak * noise_level) .* randn(size(Down_2030_Freq));

% --- 2035年 (120倍) ---
scale_2035 = 120;
Up_2035_Freq = Total_Up_Base * scale_2035;
Down_2035_Freq = Total_Down_Base * scale_2035;
% 2035年调峰 (加波动)
Up_2035_Peak = Up_2035_Freq * ratio_peak + (Up_2035_Freq * ratio_peak * noise_level) .* randn(size(Up_2035_Freq));
Down_2035_Peak = Down_2035_Freq * ratio_peak + (Down_2035_Freq * ratio_peak * noise_level) .* randn(size(Down_2035_Freq));

% --- 2040年 (135倍) ---
scale_2040 = 135;
Up_2040_Freq = Total_Up_Base * scale_2040;
Down_2040_Freq = Total_Down_Base * scale_2040;
% 2040年调峰 (加波动)
Up_2040_Peak = Up_2040_Freq * ratio_peak + (Up_2040_Freq * ratio_peak * noise_level) .* randn(size(Up_2040_Freq));
Down_2040_Peak = Down_2040_Freq * ratio_peak + (Down_2040_Freq * ratio_peak * noise_level) .* randn(size(Down_2040_Freq));

%% 4. 绘图与保存 (6张图)

% 定义绘图通用参数
x_ticks = start_hour : 6 : (start_hour + 24); % 每6小时一个刻度
line_width = 1.5;
color_freq = [0, 0.45, 0.74];  % 调频颜色 (蓝)
color_peak = [0.85, 0.33, 0.10]; % 调峰颜色 (红/橙)

years = [2030, 2035, 2040];
directions = {'Up', 'Down'};
directions_cn = {'上调', '下调'};

% 数据映射结构体
data_map = struct();
data_map.y2030.Up.Freq = Up_2030_Freq;   data_map.y2030.Up.Peak = Up_2030_Peak;
data_map.y2030.Down.Freq = Down_2030_Freq; data_map.y2030.Down.Peak = Down_2030_Peak;

data_map.y2035.Up.Freq = Up_2035_Freq;   data_map.y2035.Up.Peak = Up_2035_Peak;
data_map.y2035.Down.Freq = Down_2035_Freq; data_map.y2035.Down.Peak = Down_2035_Peak;

data_map.y2040.Up.Freq = Up_2040_Freq;   data_map.y2040.Up.Peak = Up_2040_Peak;
data_map.y2040.Down.Freq = Down_2040_Freq; data_map.y2040.Down.Peak = Down_2040_Peak;

% 循环绘制 6 张图
for y_idx = 1:length(years)
    year = years(y_idx);
    year_field = sprintf('y%d', year);
    
    for d_idx = 1:length(directions)
        dir_en = directions{d_idx};
        dir_cn = directions_cn{d_idx};
        
        freq_data = data_map.(year_field).(dir_en).Freq;
        peak_data = data_map.(year_field).(dir_en).Peak;
        
        % 创建图形
        fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600], 'Visible', 'on');
        hold on;
        
        % 绘制曲线
        plot(time_axis, freq_data, 'LineStyle', '-', 'LineWidth', line_width, 'Color', color_freq, ...
             'DisplayName', sprintf('%d年 短时间尺度场景', year));
        plot(time_axis, peak_data, 'LineStyle', '--', 'LineWidth', line_width, 'Color', color_peak, ...
             'DisplayName', sprintf('%d年 长时间尺度场景', year));
         
        % 设置坐标轴
        xlabel('时间 (小时)', 'FontSize', 20);
        ylabel(sprintf('%s潜力 (kW)', dir_cn), 'FontSize', 20);
        legend('Location', 'best', 'FontSize', 20);
        grid on;
        set(gca, 'FontSize', 20, 'XTick', x_ticks);
        xlim([time_axis(1), time_axis(end)]);
        
        % 保存文件
        filename = sprintf('%d年_%s潜力_调频vs调峰', year, dir_cn);
        fprintf('正在保存: %s ...\n', filename);
        print(fig, [filename '.png'], '-dpng', '-r600'); % 600 DPI PNG
        print(fig, [filename '.emf'], '-dmeta', '-r600'); % EMF
        
        % close(fig); % 可选：关闭图形以节省内存
    end
end

fprintf('所有 6 张绘图已完成并保存。\n');
close all