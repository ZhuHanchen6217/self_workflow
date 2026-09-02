# Cevostamab–FcRH5 D9表位、特异性与VHH设计分析报告

**报告版本：** result_report_v1  
**完成日期：** 2026-09-01  
**整理完成时刻：** 2026-09-01 11:57:44（UTC+08:00）  
**报告性质：** 基于用户提供结合图的图像约束表位预测、序列特异性分析、结构表面比较及专利实验一致性核对。

## 1. 本报告做了什么

本分析以用户上传的Cevostamab–FcRH5 D9整体结合图为主要依据，以两张局部细节图标出的残基作为空间锚点。文件夹中已有预测复合物只用于提供D9骨架和编号参考，没有把其原始结合姿态直接当作最终表位。

完成的分析包括：

1. 从整体图和细节图分别识别蓝链、紫链在绿色FcRH5 D9表面的接触区域。
2. 将D9局部编号换算到完整FcRH5蛋白编号。
3. 比较预测表位在FcRH5 D1–D8中的序列保守性。
4. 比较FcRH1、2、3、4、6所有UniProt注释Ig样ECD结构域中的同源表位。
5. 使用AlphaFold DB v6预测结构，将29个Ig样结构域对齐到FcRH5 D9，比较表面形状、映射残基、局部Cα RMSD、表面暴露面积SASA、pLDDT和化学组成。
6. 核对US11192950B2和US12030947B2中的D9删除、FcRH家族交叉反应、Biacore及膜近端功能实验。
7. 根据序列与表面结构结果，提出蓝表位、紫表位及跨区VHH的筛选与反筛策略。

## 2. 最终表位结论

### 2.1 蓝链表位

- 图中直接支持的核心残基：**E812、S814**。
- 建议用于VHH设计和突变验证的连续表位：**FcRH5 810–814，TAEHS**。
- 表面性质：约426 Å²，紧凑、偏极性/酸性，近似净电荷−1。

### 2.2 紫链表位

- 图中直接支持的核心残基：**R796、N806、S808**。
- 建议采用的扩展构象表位：**FcRH5 792–796 + 805–812**。
- 拼接序列：`TLGNR-LNLSLTAE`。
- 表面性质：约892 Å²，包含极性环、R796/E812正负电边缘，以及L805/L807/L809疏水台面。

### 2.3 联合表位

两条链共同覆盖的D9表面主要集中在 **FcRH5 792–814**：

`TLGNRSSPSGGASLNLSLTAEHS`

这是不连续的构象表位，而不是可以完全由一条短肽代表的线性表位。

## 3. 特异性结果如何解读

### 3.1 蓝表位不是独立的FcRH5特异表面

蓝表位`TAEHS`在下列结构域中完全相同：

| 结构域 | 蓝表位序列 | 局部Cα RMSD | SASA |
|---|---|---:|---:|
| FcRH5 D9 | TAEHS | 0.00 Å | 426 Å² |
| FcRH5 D7 | TAEHS | 0.33 Å | 371 Å² |
| FcRH2 Ig4 | TAEHS | 0.22 Å | 367 Å² |
| FcRH3 Ig5 | TAEHS | 0.36 Å | 347 Å² |

这说明它们不仅序列相同，局部表面形状也高度相似。FcRH5 D5同样为`TAEHS`，但局部RMSD约2.03 Å，提示其差异主要来自结构构象。

**解读：** 若VHH只包埋`TAEHS`，预计存在FcRH5 D7及FcRH2/3交叉反应风险。蓝表位路线必须把结合足迹扩展到D9特异的邻近表面。

### 3.2 紫表位具有更可利用的D9判别残基

| 结构域 | 紫表位映射序列 | 同一性 | 局部Cα RMSD | SASA |
|---|---|---:|---:|---:|
| FcRH5 D9 | TLGNR-LNLSLTAE | 13/13 | 0.00 Å | 892 Å² |
| FcRH2 Ig4 | TLGNS-FNLSLTAE | 11/13 | 0.45 Å | 790 Å² |
| FcRH3 Ig5 | TLGNS-FNLSLTAE | 11/13 | 0.37 Å | 795 Å² |
| FcRH5 D7 | TLGSS-FNLSLTAE | 10/13 | 0.31 Å | 761 Å² |
| FcRH1 Ig3 | TLGSR-FNLSLTEE | 10/13 | 0.30 Å | 879 Å² |
| FcRH3 Ig6 | TLGNI-FNLSLTTE | 10/13 | 0.18 Å | 813 Å² |
| FcRH5 D8 | TLGKI-FNLSLTTE | 9/13 | 0.58 Å | 861 Å² |

主要判别残基是：

- **R796**：FcRH2 Ig4、FcRH3 Ig5和FcRH5 D7对应为S；FcRH3 Ig6和D8对应为I。
- **L805**：高同源结构域多数为F。
- **A811**：可继续区分FcRH1 Ig3、FcRH3 Ig6及FcRH5 D8。

**解读：** 紫表位比蓝表位更适合开发FcRH5特异VHH，但FcRH2 Ig4和FcRH3 Ig5仍必须作为首要反筛抗原。

## 4. 与专利实验是否一致

总体上高度一致，但专利也证明选择性不是绝对零交叉反应。

- FIG.25B：1G7.v85结合全长FcRH5，删除D9后失去结合，支持预测表位位于D9。
- FIG.25A：在该细胞表面FACS条件下，1G7.v85结合FcRH5而未检出FcRH1–4结合。
- Biacore：1G7.v85对FcRH5约2.4 nM，对FcRH3约90 nM，即FcRH3约弱37.5倍。该结果与FcRH3 Ig5表面高度相似完全吻合。
- Table 9：多个亲和成熟变体对FcRH2或FcRH3出现弱阳性，而对FcRH1和FcRH4阴性；与结构风险排序一致。
- 膜近端实验：D9附近的1G7表位产生更强TCR信号和杀伤，支持其功能上的膜近端优势。
- 两份专利未测试FcRH6，因此FcRH6结论仍是结构预测，不能称为实验验证。

详细核对见[tables/06_patent_validation_summary.csv](tables/06_patent_validation_summary.csv)。

## 5. VHH设计建议

### 蓝表位路线

1. 不要只用`TAEHS`短肽进行正筛。
2. 使用完整折叠D9或全长ECD，促使CDR3跨到R796/L805附近的D9特异边界。
3. 强制反筛：FcRH5 D7、FcRH2 Ig4、FcRH3 Ig5。
4. 以FcRH5 D5作为“相同序列、不同构象”探针，判断VHH依赖的是序列还是表面形状。

### 紫表位路线

1. 优先让VHH同时依赖R796和L805。
2. 可让酸性/氢键受体环境识别R796，并用大小受限的疏水口袋偏好L805、排斥更大的F。
3. 尽可能同时覆盖A811，以进一步排除FcRH1 Ig3、FcRH3 Ig6和D8。
4. 第一梯队反筛：FcRH2 Ig4、FcRH3 Ig5和FcRH5 D7。
5. N806-L807-S808是潜在N-糖基化序列子，应确认天然蛋白上的糖基化占有率。

完整优先级见[tables/07_VHH_design_priority.csv](tables/07_VHH_design_priority.csv)。

## 6. 图片阅读方法

- 蓝色或紫色：与FcRH5 D9相同的映射残基。
- 橙色：同源结构域中发生氨基酸替换的映射位置。
- 白色/灰色：未纳入当前表位定义的蛋白表面。
- 化学表面图中：红色为酸性残基，蓝色为碱性残基，黄色为极性残基，灰色为疏水残基。
- `identity n/N`表示映射表位中与D9完全相同的残基数。
- `RMSD`表示全Ig域对齐后映射表位Cα的均方根偏差。越小表示局部形状越相似。
- `SASA`表示预测表面暴露面积。它描述可接近表面的大小，不直接等于结合能。

## 7. PSE使用方法

### D9-VHH设计PSE

[structures/01_blue_D9_VHH_design.pse](structures/01_blue_D9_VHH_design.pse)和[structures/02_purple_D9_VHH_design.pse](structures/02_purple_D9_VHH_design.pse)包含以下场景：

1. `01_surface_only`：无抗体遮挡的表位表面。
2. `02_with_antibody_context`：显示图像约束恢复的对应抗体链。
3. `03_cartoon_context`：cartoon及表位侧链。
4. `04_chemical_surface`：按残基化学类型着色。
5. `05_residue_labels`：完整FcRH5编号标签。

### 全结构域比较PSE

[structures/03_blue_all_domains_comparison.pse](structures/03_blue_all_domains_comparison.pse)和[structures/04_purple_all_domains_comparison.pse](structures/04_purple_all_domains_comparison.pse)包含29个已对齐的Ig样域。默认只显示FcRH5 D9；可在PyMOL对象面板逐个启用其他结构域进行同角度比较。

## 8. 文件目录说明

### images

| 文件 | 说明 |
|---|---|
| `01_blue_D9_surface.png` | 蓝表位无抗体遮挡的D9表面 |
| `02_purple_D9_surface.png` | 紫表位无抗体遮挡的D9表面 |
| `03_blue_D9_chemical_surface.png` | 蓝表位化学性质着色 |
| `04_purple_D9_chemical_surface.png` | 紫表位化学性质着色 |
| `05_blue_high_risk_comparison.png` | 蓝表位最高风险同源表面九宫格 |
| `06_purple_high_risk_comparison.png` | 紫表位最高风险同源表面九宫格 |
| `07_blue_all_domains_atlas.png` | 蓝表位在全部29个结构域上的表面图谱 |
| `08_purple_all_domains_atlas.png` | 紫表位在全部29个结构域上的表面图谱 |

### tables

| 文件 | 说明 |
|---|---|
| `01_epitope_residue_summary.csv` | 蓝/紫表位的核心、扩展和完整编号 |
| `02_domain_sequence_comparison.csv` | FcRH5 D1–D8及家族各Ig域的接触位点同一性 |
| `03_mapped_contact_residues.csv` | 每个D9表位位置向各同源域的逐残基映射 |
| `04_best_ECD_windows.csv` | 各FcRH蛋白完整ECD中最相似的连续窗口 |
| `05_structural_surface_metrics.csv` | 两表位在29个Ig域上的RMSD、SASA、pLDDT及化学组成 |
| `06_patent_validation_summary.csv` | 专利实验与本分析结论的逐项一致性 |
| `07_VHH_design_priority.csv` | 三条VHH设计路线、反筛抗原和判别残基 |

### structures

| 文件 | 说明 |
|---|---|
| `01_blue_D9_VHH_design.pse` | 蓝表位独立检查及VHH设计场景 |
| `02_purple_D9_VHH_design.pse` | 紫表位独立检查及VHH设计场景 |
| `03_blue_all_domains_comparison.pse` | 蓝表位29个同源结构域对齐比较 |
| `04_purple_all_domains_comparison.pse` | 紫表位29个同源结构域对齐比较 |

## 9. 关键限制

1. 核心残基来自截图直接标注；扩展表位来自图像约束姿态和表面邻域推断。
2. 图像约束复合物仍有少量原子冲突，PSE中的抗体场景用于表位定位而非能量评价。
3. AlphaFold结构用于同源表面比较，不代表所有环在天然糖基化、膜环境和结合状态下完全固定。
4. 专利验证的是D9依赖和家族选择性，没有提供逐残基丙氨酸扫描。
5. 最终热点需通过定点突变、SPR/BLI及细胞交叉反应实验确认。

## 10. 一句话结论

**蓝表位`TAEHS`在FcRH5 D7和FcRH2/3中序列与表面形状均高度保守，不适合单独追求FcRH5特异性；紫表位可利用R796/L805/A811的组合差异，更适合作为FcRH5选择性VHH的主要设计表面。**
