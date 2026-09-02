# FcRH5 D9 细节图约束表位与 VHH 设计报告 v4

完成时间：2026-09-01 18:20:22 中国标准时间（Asia/Shanghai）

## 核心结论

- v4 主推荐：**N795/R796/L805/N806/S808（NRLNS）**。它完整保留紫色细节图 R796/N806/S808，跨度 17.53 Å；28 个非 D9 比较域中无 5/5 或 4/5 复现。
- 蓝色细节图候选 **D790/V791/A811/E812/S814（DVAES）** 在 FcRH5-D7、FcRH2-Ig4、FcRH3-Ig5 完整复现，且局部骨架高度相似，不适合作为唯一特异性表位。
- NRLNS 与 v3 NRLSE 共享 4/5 锚点。两者五位点序列特异性相当，但 NRLNS 更紧凑，四模板平均 CDR3 锚点距离更优（1.10 vs 1.99 Å），因此作为第一优先；NRLSE 是 N806 糖链风险的平行备选。
- N806-L807-S808 是潜在 N-糖基化 sequon；当前 UniProt 本地记录未标注 N806 已糖基化，但必须实验确认其占位和可及性。
- 主库建议 CDR3 14–17 aa 覆盖多拓扑，9–12 aa 作为短环对照；必须加入 FcRH5-D7/D8、FcRH1-Ig3、FcRH2-Ig4、FcRH3-Ig5/Ig6 反筛。

## 如何阅读

先打开 `report.html`。需要旋转、拆分表面或检查具体残基时，再打开 `structures/` 中的 PSE。数值复核与后续筛选清单位于 `tables/`。

## 关键文件

- `report.html`：完整可迁移网页报告，使用相对路径。
- `structures/01_detail_constrained_epitopes.pse`：D9 上两套表位。
- `structures/02_all_domains_purple_NRLNS.pse`：NRLNS 全域比较。
- `structures/03_all_domains_blue_DVAES.pse`：DVAES 全域比较。
- `structures/04_multi_CDR3_epitope_coverage.pse`：3 表位 × 4 CDR3 模板。
- `tables/01_detail_constrained_candidate_summary.csv`：候选摘要。
- `tables/02_all_domain_sequence_surface_comparison.csv`：29 域序列和表面比较。
- `tables/05_NRLNS_vs_NRLSE.csv`：v4/v3 直接比较。
- `tables/06_multi_CDR3_epitope_coverage.csv`：CDR3 覆盖与碰撞明细。

## 证据边界

这是预测结构与用户细节图约束下的设计报告，不是实验表位鉴定。必须用糖基化 D9/ECD、点突变、完整 FcRH 家族细胞反筛和复合物预测验证。
