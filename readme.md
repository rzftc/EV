# 电动汽车(EV)集群调节潜力评估与多资源验证平台

## 项目简介

本项目专注于 **电动汽车（Electric Vehicles, EV）** 集群的灵活性调节能力评估。它不仅包含了精细化的单车电池模型和用户行为模型，还集成了集群聚合算法、用户不确定性分析模块，以及一个专门的 **联合仿真验证模块**，用于对比分析 EV 与光伏（PV）、储能（ESS）、数据中心（Data Center）等其他分布式能源（DER）的调节特性差异。

主要应用场景包括：

- 评估大规模 EV 车队在不同充电策略下的聚合调节功率（上调/下调）。
- 分析用户参与度、响应不确定性对调节潜力的影响。
- 验证 EV 模型在虚拟电厂（VPP）仿真中的准确性。

## 目录结构说明

代码库主要分为以下几个功能模块：

### 1. 数据生成与预处理 (`0inputdata/`)

- **功能**: 生成仿真所需的输入参数，处理真实世界数据，并模拟随机性。
- **核心脚本**:
   - `generateEVParameters.m`: 生成 EV 基础参数（电池容量、最大功率等）。
   - `randomize_ev_times_only.m`: 随机生成 EV 的到达和离开时间，模拟用户出行行为。
   - `predictPowerLimits.m`: 预测功率限制。

### 2. 初始化模块 (`1initialize/`)

- **功能**: 系统启动时的对象初始化。
- **核心脚本**:
   - `initializeFromExcel.m`: 从 Excel 文件读取数据并初始化 EV 对象数组。
   - `initializeParameters.m`: 设置全局仿真参数（时间步长、仿真天数等）。

### 3. 基线功率计算 (`2basepower/`)

- **功能**: 计算 EV 在不进行调节时的“基线”充电负荷。
- **核心脚本**:
   - `EVbaseP_ChargeUntilFull.m`: 模拟“即插即充”模式（一到站就满功率充直到充满）。
   - `EVbaseP_aggregate.m`: 计算集群的聚合基线功率。
   - `distributeBasePower.m`: 分配基线功率。

### 4. 用户不确定性分析 (`3useruncertainty/`)

- **功能**: 模拟用户行为的随机性和对激励信号的响应偏差。
- **核心脚本**:
   - `calculateParticipation.m`: 计算用户参与调节计划的概率。
   - `incentiveTempEV.m`: 基于价格或温度激励的用户响应模型。
   - `calculateDeltaE.m`: 计算能量偏差。

### 5. 核心状态更新与潜力计算 (`4EVupdate/`)

- **功能**: 仿真引擎的核心，负责每个时间步的状态更新和调节能力计算。
- **核心脚本**:
   - `calculateEVAdjustmentPotentia.m`: 计算单车或集群在当前时刻的最大上调和下调功率。
   - `calculateVirtualSOC.m`: 计算虚拟荷电状态（用于聚合模型）。
   - `aggregateEVs.m`: 执行 EV 集群的状态聚合。
   - `updateLockState.m`: 更新电池锁定状态（模拟电池保护机制）。

### 6. 多资源联合验证 (`10verify/`)

- **功能**: 对比 EV 与其他类型 DER 的调节特性，验证模型准确性。
- **子模块**:
   - `PV/`, `ESS/`, `DataCenter/`, `LargeLoad/`: 各类资源的独立仿真模型。

- **核心脚本**:
   - `joint_simulation_for_test.m`: 联合仿真测试脚本。
   - `verify_Joint_Potential_Accuracy.m`: 验证聚合潜力计算的准确性。
   - `test_DER_*.m`: 针对不同 DER 的测试脚本。

### 7. 主程序 (`6main/`)

- **功能**: 项目的主要入口，用于执行大规模潜力评估仿真。
- **核心脚本**:
   - `main_potential_agg_all_bound.m`: 计算并评估整体聚合调节能力的边界。
   - `main_potential_agg_ind_bound.m`: 基于个体行为的调节能力边界评估。
   - `calculate_accuracy_metrics.m`: 计算仿真结果的误差指标。

### 8. 可视化 (`5plot/`)

- **功能**: 绘制仿真结果图表。
- **核心脚本**:
   - `plotResults.m`: 通用结果绘图。
   - `plotPowerComparison.m`: 功率对比图。
   - `plotEVStatusEnhanced.m`: 增强版 EV 状态演变图（SOC、功率等）。

## 快速开始

### 环境要求

- MATLAB R2020b 或更高版本。

### 运行步骤

1. **数据准备**:
   - 确保 `0inputdata` 中的参数生成脚本已运行，或已有基础数据文件。

2. **运行潜力评估**:
   - 打开 `6main/main_potential_agg_all_bound.m`。
   - 运行脚本。程序将初始化 EV 群体，模拟其充放电过程，并输出聚合调节能力曲线。

3. **运行联合验证**:
   - 打开 `10verify/joint_simulation_for_test.m`。
   - 该脚本将同时加载 EV、PV、ESS 等模型，进行联合运行测试，以验证 EV 模型在多资源环境下的表现。

4. **查看结果**:
   - 仿真结束后，通过 `5plot` 文件夹中的脚本（如 `plotResults.m`）查看生成的功率曲线和 SOC 变化图。

## 核心算法特点

- **虚拟 SOC (Virtual SOC)**: 引入虚拟 SOC 概念来描述 EV 集群在聚合层面的能量状态，便于统一调度。
- **双向调节**: 支持评估 V2G (Vehicle-to-Grid) 放电潜力和充电功率调节潜力。
- **物理约束**: 严格考虑了电池容量、最大充放电功率、用户出行需求（到达/离开时间、目标 SOC）等物理约束。

## 注意事项

- 请注意检查 `1initialize` 中读取 Excel 文件的路径，确保与本地文件位置一致。
- `copyMFiles.m` 是用于文件管理的辅助脚本，非核心业务逻辑。

*本文档由自动化工具根据代码库内容生成。*