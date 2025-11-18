%% --------------------------------------------------
% 绘制仿真结果分析图
%
% 依赖:
% 1. 结果文件: 'main_potential_agg_vs_individual_sum_results.mat'
%    (由 main_potential_agg_ind.m 生成)
% 2. 绘图函数: 
%    - plotPowerComparison.m
%    - plotLambdaAndAggSOC.m
%    - plotLambdaAndSOC.m
%% --------------------------------------------------

clc; 
clear; 
close all;

%% 1. 加载结果文件
fprintf('正在加载结果文件...\n');
resultsFile = 'result_1000_inc.mat';

if ~exist(resultsFile, 'file')
    error(['未找到结果文件: %s\n' ...
           '请首先运行 main_potential_agg_ind.m 来生成该文件。'], resultsFile);
end
load(resultsFile); % 加载 'results' 结构体
fprintf('结果加载完毕。\n');

%% 2. 定义绘图所需的关键参数
% 这些参数必须与 main_potential_agg_ind.m 运行时使用的参数一致

% 从 initializeFromExcel.m 可知, dt_short 默认为 5 分钟
dt_short = 5; 

% 选择一个EV用于绘制单车SOC图 (例如第10辆)
% 这与 main_parfor.m 和 plot_test.m 中的选择一致
selected_ev =55; 

fprintf('使用 dt_short = %d 分钟, 绘制 EV %d 的个体图像。\n', dt_short, selected_ev);

%% 3. 调用绘图函数

% --- 绘制 图 1: 功率对比 (对应 image_b5c4e8.png) ---
% 调用 plotPowerComparison.m
fprintf('正在绘制图 1 (功率对比)...\n');
try
    plotPowerComparison(results, dt_short);
catch ME
    warning('绘制图 1 (功率对比) 失败: %s', ME.message);
end

% --- 绘制 图 2: 聚合SOC vs Lambda (对应 image_b5c4e5.png) ---
% 调用 plotLambdaAndAggSOC.m
fprintf('正在绘制图 2 (聚合SOC vs Lambda)...\n');
try
    plotLambdaAndAggSOC(results, dt_short);
catch ME
    warning('绘制图 2 (聚合SOC) 失败: %s', ME.message);
end

% --- 绘制 图 3: 单车SOC vs Lambda (对应 image_b5c507.png) ---
% 调用 plotLambdaAndSOC.m
fprintf('正在绘制图 3 (单车SOC vs Lambda)...\n');
try
    plotLambdaAndSOC(results, dt_short, selected_ev);
catch ME
    warning('绘制图 3 (单车SOC) 失败: %s', ME.message);
end

fprintf('所有图像绘制完成。\n');