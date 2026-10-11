# 在本地对猪基因组做引物脱靶扫描

`offtarget_primers.json` 含 36 条池的全部引物（28、20 条池是它的子集，按名称筛选即可）。

规则：引物 3′端 8 nt 必须完全匹配，整体错配 ≤5；展开所有简并序列逐条扫；正向/反向命中再配成 ≤2000 bp 的潜在产物。

## 步骤
1. 准备猪基因组的 2bit 文件（如 UCSC susScr11.2bit）和猪线粒体 fasta（如 NC_000845.1）。
2. 在仓库根目录运行（`<输出目录>` 就是本目录，里面已有 `offtarget_primers.json`）：

   PYTHONPATH=python_5R python3 python_5R/explore/offtarget.py <susScr11.2bit> <猪线粒体.fasta> docs/primer_design/offtarget_local

   全基因组约需几十分钟，需要 numpy。结果写到 `offtarget_hits_all.pkl`。
3. 汇总：

   PYTHONPATH=python_5R python3 docs/primer_design/offtarget_local/summarize_offtarget.py docs/primer_design/offtarget_local/offtarget_hits_all.pkl docs/primer_design/offtarget_local/offtarget_primers.json pig_offtarget.csv

4. 把 `pig_offtarget.csv`（或 pkl）发给我，我再和人基因组结果对比并改引物。

参照：人基因组 hg19 上，上一版 A5-R 有 1281 个 ≤2 错配位点，换成起点 1188 的新 A5-R 后降到 10 个。
