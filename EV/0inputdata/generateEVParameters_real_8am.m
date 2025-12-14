function generateEVParameters_real_8am(filePath, numEV, areaRatio, varargin)
% 增强版EV参数生成函数（支持车型筛选与区域类型指定）
% 输入参数：
%   filePath    - 输出文件路径（必选）
%   numEV       - EV总数（默认1000）
%   areaRatio   - 居民区比例（默认0.5，当areaType='混合'时生效）
%   varargin    - 键值对参数：
%       'ModelNames'  : 指定车型名称（默认全部车型）
%       'AreaType'    : 区域类型['居民区','工作区','混合']（默认'混合'）
%% 使用范例
% generateEVParameters_real_8am('resi_inc_1000.xlsx', 1000, 0, 'AreaType', '居民区');

%% 参数解析系统
p = inputParser;
addRequired(p, 'filePath', @ischar);
addOptional(p, 'numEV', 1000, @(x) x>0 && mod(x,1)==0);
addOptional(p, 'areaRatio', 0.5, @(x) x>=0 && x<=1);
addParameter(p, 'ModelNames', {}, @iscellstr);
addParameter(p, 'AreaType', '混合', @(x) any(strcmp(x,{'居民区','工作区','混合'})));
parse(p, filePath, numEV, areaRatio, varargin{:});

%% 精确车型参数库（增加快充功率参数）
models = [
    struct('Name','比亚迪秦PLUS DM-i', 'C',[18.32,26.86], 'SlowCharge',6.6, 'FastCharge',60,  'Ratio',0.28),...
    struct('Name','比亚迪海鸥',        'C',[30.08,38.88], 'SlowCharge',6.6, 'FastCharge',40,  'Ratio',0.25),...
    struct('Name','比亚迪宋PLUS DM-i','C',[18.32,26.86], 'SlowCharge',6.6, 'FastCharge',70,  'Ratio',0.20),...
    struct('Name','特斯拉Model Y',   'C',[60.0,78.4],   'SlowCharge',11,  'FastCharge',250, 'Ratio',0.15),...
    struct('Name','比亚迪元PLUS',     'C',[50.12,60.48],'SlowCharge',7,   'FastCharge',90,  'Ratio',0.12)...
    ];

%% 车型筛选逻辑
if ~isempty(p.Results.ModelNames)
    % 验证输入车型名称
    validModels = ismember({models.Name}, p.Results.ModelNames);
    if ~any(validModels)
        error('无效车型名称，可用车型: %s', strjoin({models.Name}, ', '));
    end
    models = models(validModels);
    % 重新归一化比例
    totalRatio = sum([models.Ratio]);
    for i = 1:length(models)
        models(i).Ratio = models(i).Ratio / totalRatio;
    end
end

%% 时间参数定义
% 居民区时间参数
residential_timeSlots = {
    struct('range',0:4,  'dist',[70,20,10], 'max_dur',[24,8,3]),
    struct('range',5:8,  'dist',[20,70,10], 'max_dur',[19,8,3]),
    struct('range',9:15, 'dist',[10,20,70], 'max_dur',[15,8,3]),
    struct('range',16:17,'dist',[20,70,10], 'max_dur',[8,8,3]),
    struct('range',18:23,'dist',[70,20,10], 'max_dur',[24,8,3])
    };
residential_mean_arrival_h = 19.0;
residential_std_dev_h = 2.0;

% 工作区时间参数
workplace_timeSlots = {
    struct('range',0:4,  'dist',[10,20,70], 'max_dur',[3,8,3]),
    struct('range',5:8,  'dist',[20,10,70], 'max_dur',[3,8,3]),
    struct('range',9:15, 'dist',[20,70,10], 'max_dur',[6,8,3]),
    struct('range',16:17,'dist',[10,20,70], 'max_dur',[3,8,3]),
    struct('range',18:23,'dist',[10,20,70], 'max_dur',[3,8,3])
    };
% 【修改点1】工作区入网时间：10.0 (即第一天早上10点)
workplace_mean_arrival_h = 10.0; 
workplace_std_dev_h = 1.0;

%% 区域分配增强系统
switch p.Results.AreaType
    case '居民区'
        isResidential = true(numEV,1);
    case '工作区'
        isResidential = false(numEV,1);
    case '混合'
        isResidential = rand(numEV,1) < p.Results.areaRatio;
    otherwise
        error('未知区域类型: %s', p.Results.AreaType);
end

%% 车型分配系统
model_probs = cumsum([models.Ratio]);
model_idx = arrayfun(@(x) find(x <= model_probs,1), rand(numEV,1));

%% 初始化数据表（完整字段）
data = table();
data.EV_ID = (1:numEV)';
data.Area = repmat("", numEV,1);
data.t_in = zeros(numEV,1);
data.t_dep = zeros(numEV,1);

%% 严格参数绑定与时间生成
for i = 1:numEV
    m = models(model_idx(i));
    data.C(i) = m.C(randi(numel(m.C))); % 电池容量二选一
    
    % 根据区域选择充电功率、时间参数
    if isResidential(i)
        data.P_N(i) = m.SlowCharge;
        data.Area(i) = "居民区";
        mean_arrival_h = residential_mean_arrival_h;
        std_dev_h = residential_std_dev_h;
        timeSlots = residential_timeSlots;
    else
        data.P_N(i) = m.FastCharge;
        data.Area(i) = "工作区";
        mean_arrival_h = workplace_mean_arrival_h;
        std_dev_h = workplace_std_dev_h;
        timeSlots = workplace_timeSlots;
    end
    
    % -------------------- (修改开始) --------------------
    
    % (步骤 1: 预先生成 E, eta)
    data.eta(i) = 0.85 + 0.1*rand();  % 充电效率
    data.E_ini(i) = data.C(i) * (0.1 + 0.3*rand());       % 初始电量10%-40%
    data.E_tar_set(i) = data.C(i) * (0.8 + 0.1999*rand()); % 目标电量80%-99.99%

    % (步骤 2: 计算真正的最小充电时间)
    delta_E = data.E_tar_set(i) - data.E_ini(i);
    if delta_E <= 0 || data.P_N(i) == 0 || data.eta(i) == 0
        t_min_charge_minutes = 0;
    else
        t_min_charge_minutes = (delta_E / (data.eta(i) * data.P_N(i))) * 60;
    end
    % 最小持续时间 (确保 t_dep > t_in)
    t_min_duration_minutes = max(1, ceil(t_min_charge_minutes)); % 至少停1分钟

    % 【修改点2】定义时间窗口：8:00 (480) 到 次日08:00 (1920)
    simulation_start_minutes = 8 * 60; % 480
    max_departure_minutes = 32 * 60;   % 1920 (次日08:00)
    
    % (步骤 3: 循环验证，直到找到满足所有约束的时间)
    valid = false;
    loop_count = 0; % 防止死循环
    
    while ~valid
        % 3a. 生成入网时间 (基于正态分布, 0-1439 min)
        generated_hour_float = mean_arrival_h + std_dev_h * randn();
        t_in_h_float_daily = mod(generated_hour_float, 24);
        t_in_minute = round(t_in_h_float_daily * 60);
        
        % 3b. 确定统计停车时长
        t_in_h_for_slot_lookup = floor(t_in_h_float_daily);
        slot_idx = find(cellfun(@(x) ismember(t_in_h_for_slot_lookup, x.range), timeSlots),1);
        if isempty(slot_idx)
             currentSlot = timeSlots{end}; 
        else
             currentSlot = timeSlots{slot_idx};
        end
        
        randVal = rand() * 100;
        if randVal <= currentSlot.dist(1) % 长时
            minDur_h = 7; maxDur_h = currentSlot.max_dur(1);
        elseif randVal <= sum(currentSlot.dist(1:2)) % 中时
            minDur_h = 3; maxDur_h = currentSlot.max_dur(2);
        else % 短时
            minDur_h = 0; maxDur_h = currentSlot.max_dur(3);
        end
        
        if maxDur_h < minDur_h
            duration_h = max(0, maxDur_h);
        else
            duration_h = minDur_h + (maxDur_h - minDur_h) * rand();
        end
        dep_add_minute_statistical = round(duration_h * 60);
        
        % 3c. 确定最终所需时长
        final_duration_minutes = max(dep_add_minute_statistical, t_min_duration_minutes);

        % 3d. 应用 8:00 (480) 平移逻辑
        % 如果生成的入网时间 < 8:00 (例如 2:00)，则认为是次日的 2:00
        if t_in_minute < simulation_start_minutes
            t_in_shifted = t_in_minute + (24 * 60); % (1440-1919)
        else
            t_in_shifted = t_in_minute; % (480-1439)
        end
        
        % 3e. 计算期望的 t_dep
        t_dep_calculated = t_in_shifted + final_duration_minutes;
        
        % 3f. 验证是否满足 32:00 (1920) 约束
        if (t_dep_calculated <= max_departure_minutes)
            valid = true;
            data.t_in(i) = t_in_shifted;
            data.t_dep(i) = t_dep_calculated;
        else
            valid = false;
            loop_count = loop_count + 1;
            if loop_count > 100 % 安全退出
                error('EV%d (P_N=%.1f) 无法在 8:00-32:00 窗口内找到满足要求的时间。', i, data.P_N(i));
            end
        end
    end % 结束 while ~valid
    
    % -------------------- (修改结束) --------------------
end

%% 电价参数系统（保持不变）
basePrice = 0.5; % 基准电价元/kWh
priceRange = 0.3; % 浮动范围
data.p_incentive = round(60*rand(numEV,1), 1); 
data.P_0 = basePrice + priceRange*randn(numEV,1)*0.2;
data.P_h_max = data.P_0 + priceRange*rand(numEV,1);
data.P_l_min = data.P_0 - priceRange*rand(numEV,1);

data.P_l_min = max(data.P_l_min, 0.2);           % 最低电价保护
data.P_h_max = max(data.P_h_max, data.P_0 + 0.1); % 合理价差保障

data.Delta_E_h_max = 0.1 * data.C;  % 高功率区能量变化
data.Delta_E_q_max = 0.05 * data.C; % 低功率区能量变化

price_std = 0.1;
data.p_real = data.P_0 + price_std * randn(numEV,1);

%% 其他参数系统
data.r = 0.05 * ones(numEV,1);        % SOC调节系数
data.state = repmat("LockOFF", numEV,1);  % 初始状态

%% 数据完整性验证
valid_C = [18.32,26.86,30.08,38.88,50.12,60.0,60.48,78.4];
valid_PN = [6.6,7,11,40,60,70,90,250];
assert(all(ismember(data.C,valid_C)),'电池容量校验失败');
assert(all(ismember(data.P_N,valid_PN)),'充电功率校验失败');
assert(all(data.t_dep > data.t_in),'时间顺序校验失败 (t_dep > t_in)');

% 功率约束最终验证
m3_all = (data.E_tar_set - data.E_ini)./(data.eta.*(data.t_dep - data.t_in)/60);
m3_all(isnan(m3_all)) = 0; % 处理 (E_tar==E_ini) / (duration > 0) 导致的 0/t = 0
m3_all(m3_all < 0) = 0; % 处理 (E_tar < E_ini) 的情况
assert(all(data.P_N >= m3_all),'存在不满足P_N >= m3 (平均所需功率) 的情况');

%% 数据存储
writetable(data, filePath);
fprintf('EV数据生成完成：%s\n居民区：%d辆（%.1f%%）\n工作区：%d辆（%.1f%%）\n',...
    filePath,...
    sum(isResidential), 100*mean(isResidential),...
    sum(~isResidential), 100*mean(~isResidential));

% (新增) 打印时间范围进行验证
fprintf('时间范围校验 (应在 8.00H - 32.00H 之间):\nMin t_in: %.2f (%.2f H)\nMax t_in: %.2f (%.2f H)\n',...
    min(data.t_in), min(data.t_in)/60, max(data.t_in), max(data.t_in)/60);
fprintf('Min t_dep: %.2f (%.2f H)\nMax t_dep: %.2f (%.2f H)\n',...
    min(data.t_dep), min(data.t_dep)/60, max(data.t_dep), max(data.t_dep)/60);
end