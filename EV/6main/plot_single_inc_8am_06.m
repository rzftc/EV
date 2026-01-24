%% plot_single_inc.m
clc;
clear;
close all;

%% 1. 加载数据
resultsFile = 'main_potential_5min_1000_8am_bound_06.mat';
fprintf('正在加载结果文件: %s\n', resultsFile);

if ~exist(resultsFile, 'file')
    warning(['未找到基础结果文件 "%s"，前序图表可能无法绘制。\n' ...
        '请确保文件存在。'], resultsFile);
else
    data = load(resultsFile); % 加载 'results' 结构体
    results = data.results;
    fprintf('结果加载完毕。\n');
end

%% 2. 准备绘图参数
% 仿真参数 (必须与 main_potential_agg_ind.m 匹配)
dt_short = 5; % 默认短步长为 5 分钟
simulation_start_hour = 8; % 仿真开始时间
selected_ev = 945; % 选择绘制的EV编号
% selected_ev = 945; % 选择绘制的EV编号

% 计算时间轴 (小时) - 这里的 total_steps 基于 results 计算
if exist('results', 'var')
    total_steps = length(results.lambda);
    time_hours = ((0:total_steps-1) * dt_short / 60) + simulation_start_hour;
else
    % 如果基础文件不存在，为后续图表预定义时间轴（假设24小时，dt=5min）
    total_steps = 24 * 60 / 5; 
    time_hours = ((0:total_steps-1) * 5 / 60) + simulation_start_hour;
end

% 定义新的坐标轴刻度
x_ticks = [8, 14, 20, 26, 32];
x_tick_labels = {'D1 08:00', 'D1 14:00', 'D1 20:00', 'D2 02:00', 'D2 08:00'};

fprintf('坐标轴已设置为 %d:00 (D1) 到 %d:00 (D2)。\n', simulation_start_hour, simulation_start_hour + 24);

%% --------------------------------------------------
% 图 1: 功率跟踪效果分析
% 保存为: 功率跟踪效果分析.png
% --------------------------------------------------
if exist('results', 'var')
    fprintf('正在绘制图 1 (功率对比)...\n');
    fig1 = figure('Name', '功率跟踪效果分析', 'Position', [100 100 1000 400], 'NumberTitle', 'off');

    % 绘制目标功率（橙色半透明区域）
    area(time_hours, results.P_tar, ...
        'FaceColor', [1 0.6 0], ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', ...
        'DisplayName', '目标功率');
    hold on;

    % 绘制实际功率（深绿色实线）
    plot(time_hours, results.P_agg, ...
        'LineWidth', 1.5, ...
        'Color', [0 0.5 0], ...
        'DisplayName', '实际功率');
    hold off;

    % 坐标轴和标签设置 (字体放大)
    xlabel('时间 (小时)', 'FontSize', 20);
    ylabel('功率 (kW)', 'FontSize', 20);
    set(gca, 'FontSize', 16); % 坐标轴刻度字体
    xlim([simulation_start_hour, simulation_start_hour + 24]); % [6, 30]
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    ylim_max = max([results.P_agg, results.P_tar]) * 1.1;
    ylim([0 ylim_max]);
    grid on;
    legend('Location', 'northwest', 'FontSize', 16);

    % 保存图像 (无标题，高DPI，中文名)
    print(fig1, '功率跟踪效果分析.png', '-dpng', '-r600');
end

%% --------------------------------------------------
% [修改] 图 2: 聚合模型实时功率与实际统计功率对比 (拆分为全景图与局部放大图)
% 保存为: 聚合功率计算验证_全景.png, 聚合功率计算验证_放大.emf
% --------------------------------------------------
if exist('results', 'var') && isfield(results, 'EV_Power')
    fprintf('正在绘制图 2 (聚合功率计算验证: 拆分保存)...\n');

    % --- 1. 计算误差与放大区域参数 ---
    % 计算两者的绝对误差
    diff_power = abs(results.P_agg - results.EV_Power);
    [~, max_diff_idx] = max(diff_power); % 找到误差最大的索引
    
    % 定义放大窗口 (最大误差时刻的前后 1.5 小时)
    t_center = time_hours(max_diff_idx);
    zoom_span = 0.2; % 单侧跨度 (小时)
    t_zoom_start = max(simulation_start_hour, t_center - zoom_span);
    t_zoom_end = min(simulation_start_hour + 24, t_center + zoom_span);
    
    % 获取放大区域的数据范围 (用于设置子图Y轴)
    zoom_mask = (time_hours >= t_zoom_start) & (time_hours <= t_zoom_end);
    % 增加一点余量防止曲线贴边
    y_zoom_max = max([results.P_agg(zoom_mask), results.EV_Power(zoom_mask)]) * 1.05;
    y_zoom_min = min([results.P_agg(zoom_mask), results.EV_Power(zoom_mask)]) * 0.95;

    % --- 2. 绘制并保存全景图 (Mother Plot) ---
    fig2_main = figure('Name', '聚合功率计算验证_全景', 'Position', [150 150 1000 400], 'NumberTitle', 'off');
    hold on;
    
    % 绘制 P_agg (黑色实线)
    plot(time_hours, results.P_agg, ...
        'LineWidth', 2.5, ...
        'Color', 'k', ... 
        'DisplayName', '单体功率之和');

    % 绘制 EV_Power (红色虚线)
    plot(time_hours, results.EV_Power, ...
        'LineStyle', '--', ...
        'LineWidth', 2.0, ...
        'Color', 'r', ... 
        'DisplayName', '聚合模型功率');
    
    % 在全景图中绘制蓝色虚线框，指示被放大的区域
    rectangle('Position', [t_zoom_start, y_zoom_min, (t_zoom_end - t_zoom_start), (y_zoom_max - y_zoom_min)], ...
              'EdgeColor', 'b', 'LineStyle', '-.', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    hold off;

    % 全景图坐标轴设置 (字体放大)
    xlabel('时间 (小时)', 'FontSize', 20);
    ylabel('功率 (kW)', 'FontSize', 20);
    set(gca, 'FontSize', 16);
    xlim([simulation_start_hour, simulation_start_hour + 24]); 
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    legend('Location', 'best', 'FontSize', 16);

    % 保存全景图
    print(fig2_main, '聚合功率计算验证_全景.png', '-dpng', '-r600');

    % --- 3. 绘制并保存局部放大图 (Child Plot) ---
    % 【修改】保存为 EMF 矢量格式，并去除网格线
    fig2_zoom = figure('Name', '聚合功率计算验证_放大', 'Position', [200 200 600 400], 'NumberTitle', 'off');
    hold on;
    
    % 重绘曲线 (保持样式一致)
    plot(time_hours, results.P_agg, 'k-', 'LineWidth', 2.5, 'DisplayName', '单体功率之和');
    plot(time_hours, results.EV_Power, 'r--', 'LineWidth', 2.0, 'DisplayName', '聚合模型功率');
    
    hold off;
    
    % 设置局部坐标轴范围 (聚焦误差最大处)
    xlim([t_zoom_start, t_zoom_end]);
    ylim([y_zoom_min, y_zoom_max]);
    
    % 局部图坐标轴设置 (字体特大)
    set(gca, 'FontSize', 24); % 保持较大的刻度字体
    grid off; % <--- 【修改】去除网格线
 
    % === 关键修改：保存为 EMF 矢量图 ===
    % 设置属性以确保 EMF 转换时不带背景色
    set(gcf, 'Color', 'none'); 
    set(gca, 'Color', 'none');
    set(gcf, 'InvertHardcopy', 'off'); 

    % 保存放大图
    print(fig2_zoom, '聚合功率计算验证_放大.emf', '-dmeta');
else
    if exist('results', 'var')
        warning('results 结构体中缺少 EV_Power 字段，跳过图 2 绘制。');
    end
end
%% --------------------------------------------------
% 图 3: Lambda 与 单台EV SOC 协同分析
% 保存为: 单体SOC与Lambda.png
% --------------------------------------------------
if exist('results', 'var')
    fprintf('正在绘制图 3 (单车SOC vs Lambda)...\n');

    if selected_ev > size(results.EV_S_original, 1)
        warning('selected_ev (%d) 大于EV总数 (%d)。跳过图 3 绘制。', selected_ev, size(results.EV_S_original, 1));
    else
        fig3 = figure('Name', sprintf('EV%d-Lambda&SOC协同', selected_ev), 'Position', [200 200 1000 400], 'NumberTitle', 'off');

        % 左侧坐标轴（SOC）
        yyaxis left;
        main_soc_ind = plot(time_hours, results.EV_S_original(selected_ev, :), ...
            'LineWidth', 1.2, ...
            'Color', [0.8 0.2 0.2], ...
            'DisplayName', '电量偏差度');
        ylabel('电量偏差度 (-1~1)', 'FontSize', 20, 'Color', [0.8 0.2 0.2]);
        ylim([-2 ,2]);
        set(gca, 'YColor', [0.8 0.2 0.2]);

        % 右侧坐标轴（Lambda）
        yyaxis right;
        main_lambda_ind = plot(time_hours, results.lambda, ...
            'LineWidth', 1.2, ...
            'Color', [0.2 0.4 0.8], ...
            'DisplayName', '\lambda^*');
        ylabel('\lambda^*', 'FontSize', 20, 'Color', [0.2 0.4 0.8]);
        % 使用与图2相同的动态Y轴范围
        ylim([-2 ,2]);
        set(gca, 'YColor', [0.2 0.4 0.8]);

        % 公共设置 (字体放大)
        xlabel('时间 (小时)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        xlim([simulation_start_hour, simulation_start_hour + 24]); % [6, 30]
        set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
        grid on;
        legend([main_soc_ind, main_lambda_ind], 'Location', 'northwest', 'FontSize', 16);
        
        % 保存图像 (无标题，高DPI，中文名)
        print(fig3, '单体SOC与Lambda.png', '-dpng', '-r600');
    end
end

%% --------------------------------------------------
% 图 4: 单台EV 电量对比 (基线 vs 实际)
% 保存为: 单体电量对比.png
% --------------------------------------------------
if exist('results', 'var')
    fprintf('正在绘制图 4 (单车电量: 基线 vs 实际)...\n');

    if selected_ev > size(results.EV_E_actual, 1)
        warning('selected_ev (%d) 大于EV总数 (%d)。跳过图 4 绘制。', selected_ev, size(results.EV_E_actual, 1));
    elseif ~isfield(results, 'EV_E_baseline') || ~isfield(results, 'EV_E_actual')
        warning('结果文件中缺少电量数据 (EV_E_baseline 或 EV_E_actual)。跳过图 4 绘制。');
    else
        fig4 = figure('Name', sprintf('EV%d-电量对比', selected_ev), 'Position', [250 250 1000 400], 'NumberTitle', 'off');

        % 绘制基线电量
        plot(time_hours, results.EV_E_baseline(selected_ev, :), ...
            '--', 'LineWidth', 2, ...
            'Color', [0.5 0.5 0.5], ...
            'DisplayName', '基线电量 (Baseline)');
        hold on;

        % 绘制实际电量
        plot(time_hours, results.EV_E_actual(selected_ev, :), ...
            '-', 'LineWidth', 2, ...
            'Color', [0 0.6 0.8], ...
            'DisplayName', '实际电量 (Actual)');
        hold off;

        % 坐标轴和标签设置 (字体放大)
        xlabel('时间 (小时)', 'FontSize', 20);
        ylabel('电量 (kWh)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        xlim([simulation_start_hour, simulation_start_hour + 24]); % [6, 30]
        set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
        grid on;
        legend('Location', 'best', 'FontSize', 16);
        
        % 保存图像 (无标题，高DPI，中文名)
        print(fig4, '单体电量对比.png', '-dpng', '-r600');
    end
end

%% --------------------------------------------------
% 算例 A (图 5): 用户行为不确定性分布 (入网/离网时间)
% 保存为: 用户行为分布.png
% --------------------------------------------------
if exist('results', 'var')
    fprintf('正在绘制图 5 (用户行为分布)...\n');

    if ~isfield(results, 'EV_t_in') || ~isfield(results, 'EV_t_dep')
        warning('结果文件中缺少时间分布数据 (EV_t_in 或 EV_t_dep)。跳过图 5 绘制。');
    else
        fig5 = figure('Name', '用户行为不确定性分布', 'Position', [300 300 1000 400], 'NumberTitle', 'off');
        
        % 转换时间单位到小时
        t_in_h = results.EV_t_in / 60;
        t_dep_h = results.EV_t_dep / 60;
        
        % 子图1: 入网时间分布 (字体放大)
        subplot(1, 2, 1);
        histogram(t_in_h, 24, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
        hold on;
        xlabel('入网时间 (小时)', 'FontSize', 20);
        ylabel('概率密度', 'FontSize', 20);
        grid on;
        set(gca, 'FontSize', 16);
        
        % 子图2: 离网时间分布 (字体放大)
        subplot(1, 2, 2);
        histogram(t_dep_h, 24, 'Normalization', 'pdf', 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
        xlabel('离网时间 (小时)', 'FontSize', 20);
        ylabel('概率密度', 'FontSize', 20);
        grid on;
        set(gca, 'FontSize', 16);
        
        % 保存图像
        print(fig5, '用户行为分布.png', '-dpng', '-r600');
    end
end

%% --------------------------------------------------
% 算例 B (图 6): 聚合体调节能力边界 (可行域)
% 保存为: 聚合调节边界.png
% --------------------------------------------------
if exist('results', 'var')
    fprintf('正在绘制图 6 (聚合调节边界)...\n');

    if ~isfield(results, 'P_base_agg') || ~isfield(results, 'EV_Up')
        warning('结果文件中缺少聚合功率数据。跳过图 6 绘制。');
    else
        fig6 = figure('Name', '聚合体调节能力边界', 'Position', [350 350 1000 400], 'NumberTitle', 'off');
        
        % 计算绝对边界
        P_upper_bound = results.P_base_agg + results.EV_Up;
        P_lower_bound = results.P_base_agg + results.EV_Down;
        
        hold on;
        
        % 绘制实际聚合功率
        plot(time_hours, results.P_agg, '-', 'LineWidth', 1.5, 'Color', [0.1 0.6 0.3], 'DisplayName', '实际功率 (Actual)');
        
        hold off;
        
        % 坐标轴设置 (字体放大)
        xlabel('时间 (小时)', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        xlim([simulation_start_hour, simulation_start_hour + 24]); 
        set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
        grid on;
        legend('Location', 'best', 'FontSize', 16);
        
        % 保存图像
        print(fig6, '聚合调节边界.png', '-dpng', '-r600');
    end
end

%% --------------------------------------------------
% 图 7: 不同调节时长上调节能力对比 (5min vs 15min vs 60min)
% 保存为: 不同调节时长上调节能力对比.png
% --------------------------------------------------
fprintf('正在绘制图 7 (不同调节时长上调节能力对比)...\n');

% 定义需加载的文件名、步长与标签
file_list = {'main_potential_5min_1000_8am_bound_06.mat', 'main_potential_15min_1000_8am_bound_06.mat', 'main_potential_60min_1000_8am_bound_06.mat'};
dt_list = [5, 15, 60]; % 对应每个文件的步长
legend_labels = {'t_{adj}=5min', 't_{adj}=15min', 't_{adj}=60min'};
line_colors = {'r', 'g', 'b'}; 

% 创建新画布
fig7 = figure('Name', '不同调节时长上调节能力对比', 'Position', [400 400 1000 400], 'NumberTitle', 'off');
hold on;

% 绘制零线
yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

file_found_count = 0;
for i = 1:length(file_list)
    target_file = file_list{i};
    current_dt = dt_list(i); % 获取当前文件的步长
    
    if exist(target_file, 'file')
        fprintf('正在加载对比数据: %s (dt=%d min)\n', target_file, current_dt);
        tmp_data = load(target_file);
        file_found_count = file_found_count + 1;
        
        if isfield(tmp_data, 'results')
            tmp_res = tmp_data.results;
            
            % 为当前文件单独计算时间轴
            tmp_total_steps = length(tmp_res.EV_Up);
            tmp_time_hours = ((0:tmp_total_steps-1) * current_dt / 60) + simulation_start_hour;

            % 绘制上调潜力 (实线)
            plot(tmp_time_hours, tmp_res.EV_Up, ...
                'LineStyle', '-', ...
                'Color', line_colors{i}, ...
                'LineWidth', 1.5, ...
                'DisplayName', [legend_labels{i} ' 上调 (Up)']);
        else
            warning('文件 %s 中未找到 results 结构体。', target_file);
        end
    else
        fprintf('提示: 未找到文件 %s，跳过绘制。\n', target_file);
    end
end

if file_found_count > 0
    hold off;

    % 坐标轴设置 (字体放大)
    xlabel('时间 (小时)', 'FontSize', 20);
    ylabel('上调节潜力 (kW)', 'FontSize', 20);
    set(gca, 'FontSize', 16);
    xlim([simulation_start_hour, simulation_start_hour + 24]); % [6, 30]
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    
    % 图例设置 (字体放大)
    legend('Location', 'west', 'FontSize', 16);

    % 保存图像
    print(fig7, '不同调节时长上调节能力对比.png', '-dpng', '-r600');
else
    fprintf('未找到任何对比文件，图 7 绘制跳过。\n');
    close(fig7);
end

%% --------------------------------------------------
% 图 8: 不同调节时长下调节能力对比 (5min vs 15min vs 60min)
% 保存为: 不同调节时长下调节能力对比.png
% --------------------------------------------------
fprintf('正在绘制图 8 (不同调节时长下调节能力对比)...\n');

% 创建新画布
fig8 = figure('Name', '不同调节时长下调节能力对比', 'Position', [400 100 1000 400], 'NumberTitle', 'off');
hold on;

% 绘制零线
yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

file_found_count = 0;
for i = 1:length(file_list)
    target_file = file_list{i};
    current_dt = dt_list(i); % 获取当前文件的步长
    
    if exist(target_file, 'file')
        % 已经加载过文件，理论上可以直接用，但为了逻辑清晰重新加载或使用缓存
        % 这里选择重新加载以保证独立性
        tmp_data = load(target_file);
        file_found_count = file_found_count + 1;
        
        if isfield(tmp_data, 'results')
            tmp_res = tmp_data.results;
            
            % 为当前文件单独计算时间轴
            tmp_total_steps = length(tmp_res.EV_Down);
            tmp_time_hours = ((0:tmp_total_steps-1) * current_dt / 60) + simulation_start_hour;

            % 绘制下调潜力 (实线以清晰展示)
            plot(tmp_time_hours, tmp_res.EV_Down, ...
                'LineStyle', '-', ... 
                'Color', line_colors{i}, ...
                'LineWidth', 1.5, ...
                'DisplayName', [legend_labels{i} ' 下调 (Down)']);
        end
    end
end

if file_found_count > 0
    hold off;

    % 坐标轴设置 (字体放大)
    xlabel('时间 (小时)', 'FontSize', 20);
    ylabel('下调节潜力 (kW)', 'FontSize', 20);
    set(gca, 'FontSize', 16);
    xlim([simulation_start_hour, simulation_start_hour + 24]); % [6, 30]
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    
    % 图例设置 (字体放大)
    legend('Location', 'west', 'FontSize', 16);

    % 保存图像
    print(fig8, '不同调节时长下调节能力对比.png', '-dpng', '-r600');
else
    fprintf('未找到任何对比文件，图 8 绘制跳过。\n');
    close(fig8);
end

fprintf('所有图像绘制完成。\n');