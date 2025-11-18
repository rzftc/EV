% --- 绘制 图 1: 单峰分布 (正态分布) ---
figure;

% 1. 生成模拟数据
% 我们根据图像估计参数
N = 30000;       % 数据点数量 (越多图像越平滑)
mu = 9;          % 均值 (峰值位置)
sigma = 3;       % 标准差 (分布的宽度)

% 生成正态分布随机数
data1 = normrnd(mu, sigma, [N, 1]);

% 2. 绘制直方图
h1 = histogram(data1);

% 3. 设置标准化模式 (使用 'probability' 来匹配 Y 轴刻度)
% 'probability' 模式: 所有条形的高度总和为 1
% 'pdf' 模式: 所有条形的面积总和为 1 (峰值会高很多)
h1.Normalization = 'probability';

% 4. 调整外观以匹配图像
h1.BinWidth = 0.2;              % 调整条形的宽度，越小越精细
h1.FaceColor = [0.2 0.2 0.2];  % 设置为深灰色
h1.EdgeColor = 'none';          % 去除条形边缘线

% 5. 设置坐标轴和标签
ax1 = gca; % 获取当前坐标轴
ax1.FontSize = 12;
ylabel('概率密度');
xlabel('时间/小时');

xlim([0 22]);
ylim([0 0.035]);
ax1.YTick = 0:0.005:0.035; % 设置Y轴刻度 (使用 0.005 间隔更接近原图)
box on; % 显示边框

% --- 绘制 图 2: 双峰分布 (混合分布) ---
figure;

% 1. 生成模拟数据
N_total = 30000;   % 总数据点

% 峰 1 (左侧) 的参数
prop1 = 0.15;      % 占比 15%
N1 = round(N_total * prop1);
mu1 = 1.5;
sigma1 = 0.7;
data_c1 = normrnd(mu1, sigma1, [N1, 1]);

% 峰 2 (右侧) 的参数
N2 = N_total - N1; % 剩余 85%
mu2 = 18;
sigma2 = 2.5;
data_c2 = normrnd(mu2, sigma2, [N2, 1]);

% 合并两组数据
data2 = [data_c1; data_c2];

% 清理数据：原始图像数据都大于0
data2 = max(data2, 0);

% 2. 绘制直方图
h2 = histogram(data2);

% 3. 设置标准化模式 (同样使用 'probability')
h2.Normalization = 'probability';

% 4. 调整外观以匹配图像
h2.BinWidth = 0.2;              % 保持和图1一致的条形宽度
h2.FaceColor = [0.2 0.2 0.2];  % 设置为深灰色
h2.EdgeColor = 'none';          % 去除条形边缘线

% 5. 设置坐标轴和标签
ax2 = gca; % 获取当前坐标轴
ax2.FontSize = 12;
ylabel('概率密度');
xlabel('时间/小时');

xlim([0 25]);
ylim([0 0.035]);
ax2.YTick = 0:0.005:0.035;
box on; % 显示边框