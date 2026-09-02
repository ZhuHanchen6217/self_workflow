# FcRH5 D9 单/双VHH表位设计报告 v3修订版

完成时间：2026-09-01 17:55:48（Asia/Shanghai）

## 修订结论

- 单VHH首选改为跨轻/重链足迹的 **N795/R796/L805/S808/E812（NRLSE）**，跨度23.14 Å；28个非D9比较域中无完整复现、无4/5近似复现。
- 原 **NRLS** 降为紧凑回退方案；它保留重链细节图，但没有读取轻链细节图。
- v2的L781/R785/D822/R829/T832与两张细节图的五个直接残基零重叠，因此只能视为D9特异性探索表面，不能视为Cevostamab图像支持表位。
- 4个VHH模板的CDR3（9/14/17/17 aa）都能在某个刚性取向下覆盖NRLSE 5/5锚点；原子碰撞差异为41/19/50/10，说明长度和环构象必须同时优化。
- 单VHH主库推荐CDR3 14–17 aa，覆盖14/15/16/17四档；9 aa仅作短环对照。
- 双VHH方案仍为A=NRL特异门控、B=DVAES协同臂，与单VHH NRLSE是不同设计路线。

## 关键文件

1. `report.html`：完整图文报告。
2. `structures/01_chainwise_anchor_design.pse`：细节图锚点、NRLSE及双位点。
3. `structures/05_multi_CDR3_template_coverage.pse`：4个CDR3模板对3套锚点的12个场景。
4. `tables/07_multi_CDR3_template_coverage.csv`：CDR3覆盖与碰撞明细。
5. `tables/08_detail_image_anchor_reconciliation.csv`：v2/v3与细节图逐残基核对。

## 证据边界

细节图是公司预测分析而非实验结构；多模板结果是未松弛刚性放置。最终必须以糖基化D9/ECD突变、全域反筛、细胞结合和多模型复合物预测验证。
