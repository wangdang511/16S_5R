# 自建 UDP（DN）生成说明

目标：生成 10 nt 内联 UDP 的候选池，用于 5′ 端加在 16S 引物池上（i5 加在正向引物、i7 加在反向引物）。序列全部由程序从头生成，不借用 Illumina 序列。

## 规则
**单条序列**：长度 10；GC 4–6/10；无同聚物 ≥3；无二核苷酸重复 ≥3 次；无 GGGG；前两位不同时为 G；不含 TruSeq/Nextera/ME/P5/P7 接头序列的 6-mer（含反向互补）；非回文（与自身反向互补的编辑距离 ≥4）；与 Illumina UDP0001–0384 的 i5/i7（含反向互补）不相同（汉明距离 ≥1）。

**集合（池级）**：同类序列（i5 之间、i7 之间）Levenshtein ≥4、与反向互补 ≥3。前 384 对（阶段 1–2）里 768 条序列两两编辑距离都 ≥4；扩展的 96 对只保证同类内部满足。

**配色**：每个 48 组、每个周期的 GC、A+C、A+T 占比 40–60%，A/C/G/T 各 ≥15%，i5、i7 各自满足；相邻 2 个 48 组组成 96 组块，同样满足。

## 为什么反向互补放宽到 3
严格要求 768 条序列两两编辑距离 ≥4 且反向互补 ≥4 时，字典序贪心最多凑出 464 条，不够 768；反向互补 ≥3 时有 836–867 条。480 对（960 条）时，全部序列两两 ≥4 也超过贪心的上限（约 935 条），所以扩展部分只要求同类内部。

## 文件
- `reproduce.sh`：从头复现的命令和参数（随机种子：阶段 2 为 4，扩展为 1；阶段 1 无随机性）。
- `data/stage1_indep_GCAT_D4_RC3.pkl`：阶段 1 独立集（836 条）；已核对重跑结果一致。
- `data/stage2_blocks_seed4.json`：阶段 2 的 16 个 48 组分配；`data/dn384_final.json`：验证后的 384 对。
- `data/ext_i5_*.pkl`、`data/ext_i7_*.pkl`：扩展阶段的候选；`data/dn480.json`：最终 480 对（i5、i7 各 480 条，前 384 对不变）。
- 代码都在 `python_5R/explore/`：`udp_denovo384.py`、`udp_denovo384_stage2.py`、`udp_denovo384_verify.py`、`udp_denovo_extend.py`、`udp_denovo_480.py`；管内评估 `udp_tail_eval.py`；选 48/96 组 `udp_dn_select384.py`。
- 序列表：`docs/primer_design/udp_tail/udp_denovo_384.tsv`、`udp_denovo_480.tsv`、`udp_denovo_new96.tsv`；Excel 和 csv 见 `docs/primer_design/udp_denovo_384.*`。

## 注意
- 阶段 2 的二聚体预筛指标只作软代价，不能预测管内分级（A 级比例约 21%），分级以 `udp_tail_eval.py` 的整管评估为准。
- 阶段 2 重跑需要约 15 分钟；评估整个池需要约 25 分钟/96 个管（4 核）。
