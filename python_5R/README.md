# python_5R — 5R / SMURF 流程的 Python 复现

| 文件 | 说明 |
|---|---|
| `smurf5r.py` | 与 `mFiles/*.m` 一一对应的 Python 实现（每个函数注释中标明对应的 .m 文件），另含引物覆盖度评估与从 FASTA 构建 5R 数据库的工具 |
| `primer_design.py` / `design_primers.py` | 引物重设计：把参考序列投影到 E. coli 坐标、扫描保守位点、按“同管引物位点不重叠”精确求解 1/2 管平铺、Tm 调整与二聚体检查 |
| `make_notebook.py` | 生成 `../notebooks/5R_SMURF_step_by_step.ipynb` |
| `../notebooks/5R_SMURF_step_by_step.ipynb` | 逐步复现 + 与官方 `example_results/` 比对 + 扩展分析（已执行，含输出） |
| `../docs/primer_redesign.html` | 引物重设计结果（加入 V4）与现有 5R 的坐标核对 |
| `../docs/5R_SMURF_pipeline.html` | 流程图解、算法讲解、代码问题清单与改进建议 |

## 环境

```bash
pip install numpy scipy pandas matplotlib jupyter
```

## 快速使用

```python
import smurf5r as s
cfg  = s.Config(kmer_len=126)
dbs  = s.load_bact_db('../GG_5R', cfg)                 # ~40 s
_, taxa = s.load_headers_taxonomy('../GG_5R')          # ~1.5 min（可缓存为 pickle）
for name, [(r1, r2)] in s.find_sample_pairs('fastq_dir').items():
    res = s.run_sample(r1, r2, dbs, cfg, taxa)         # ~10 s / 样本
    print(name, res['table'].sort_values('freq', ascending=False).head())
```

## 复现精度（示例数据）

* `ReadCountStats` 的 28 项计数与 MATLAB 结果逐项一致；
* 物种检出集合完全一致（67 / 286 个物种），丰度最大绝对差 3×10⁻⁶（EM 收敛容差量级）。

## 与 MATLAB 版的差异

* read 与 k-mer 的 ≤2 错配检索使用鸽巢原理分段哈希索引（结果与暴力比对相同，Notebook 中有验证），速度快约 3 个数量级；
* EM 不做依赖墙钟时间的 60 秒剪枝（MATLAB 版的行为与机器速度有关）；
* `Config(faithful=True)`（默认）保留原代码 R2 质量过滤误用 `q1` 的 bug 以便逐位复现，设为 `False` 即修正。
