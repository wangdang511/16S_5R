"""由 iter5R_primers.csv / iter_scheme.csv 生成 docs/primer_5R_iteration.html"""
import re, html
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>5R 引物迭代</title>", old[:old.index("</style>") + 8], 1)
P = pd.read_csv(f"{HERE}/iter5R_primers.csv"); SCH = pd.read_csv(f"{HERE}/iter_scheme.csv")
pct = lambda x: f"{x * 100:.0f}%"
CFG = ["原引物", "1条 ≤4简并", "1条 ≤16简并", "2条 各≤4", "3条 各≤4"]

NUMRE = re.compile(r"^[\d.%/ –→-]+$")


def cell(i, v, seqcol):
    cls = "seq" if i == seqcol else ("n" if NUMRE.match(str(v)) else "")
    return '<td class="' + cls + '">' + html.escape(str(v)) + "</td>"


def tbl(rows, heads, seqcol=None):
    h = "".join(f"<th>{c}</th>" for c in heads)
    b = "".join("<tr>" + "".join(cell(i, v, seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


avg = P.groupby("config")[["n_oligos", "total_fold", "SILVA(48k)_key_mean", "SILVA(48k)_key_min", "GG留出(12k)_key_mean", "GG留出(12k)_key_min"]].mean().reindex(CFG)
t1 = tbl([[c, f"{r.n_oligos:.1f}", f"{r.total_fold:.1f}", pct(r["SILVA(48k)_key_mean"]), pct(r["SILVA(48k)_key_min"]), pct(r["GG留出(12k)_key_mean"]), pct(r["GG留出(12k)_key_min"])] for c, r in avg.iterrows()],
         ["配置", "每个位点寡核苷酸数", "每个位点总简并数", "SILVA 验证 主要门平均", "SILVA 验证 最差门", "GG 留出 主要门平均", "GG 留出 最差门"])
rows = []
for _, r in SCH.iterrows():
    rows.append([r.db, r.config, r.per_region, f"{r.mean_regions:.2f}", pct(r.ge1), pct(r.ge3), pct(r.all5), f"{r.exp_bits:.0f}", f"{pct(r.ge1_phylum_min)}（{r.ge1_worst}）", f"{pct(r.ge3_phylum_min)}（{r.ge3_worst}）"])
t2 = tbl(rows, ["验证集", "配置", "各区域扩增子覆盖 R1/R2/R3/R4/R5 %", "平均扩出区域数", "≥1 个区域", "≥3 个区域", "5 个全部", "每条序列平均信息量（bit）", "≥1 个区域：最差门", "≥3 个区域：最差门"])
per = []
for name in P.primer.unique():
    sub = P[P.primer == name].set_index("config")
    per.append([name, sub.loc["原引物", "site"], pct(sub.loc["原引物", "SILVA(48k)_key_min"]) + " → " + pct(sub.loc["1条 ≤16简并", "SILVA(48k)_key_min"]) + " → " + pct(sub.loc["3条 各≤4", "SILVA(48k)_key_min"]),
                pct(sub.loc["原引物", "SILVA(48k)_key_mean"]) + " → " + pct(sub.loc["1条 ≤16简并", "SILVA(48k)_key_mean"]) + " → " + pct(sub.loc["3条 各≤4", "SILVA(48k)_key_mean"])])
t3 = tbl(per, ["引物", "位置", "最差主要门：原 → 1 条 ≤16 简并 → 3 条", "主要门平均：原 → 1 条 ≤16 简并 → 3 条"])
seqrows = []
for _, r in P[P.config.isin(["原引物", "1条 ≤16简并", "2条 各≤4", "3条 各≤4"])].iterrows():
    seqrows.append([r.primer, r.config, r.n_oligos, r.total_fold, r.seqs])
t4 = tbl(seqrows, ["引物", "配置", "寡核苷酸数", "总简并数", "序列 5′→3′（多条用 | 分隔）"], seqcol=4)

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物迭代</div>
  <h1>5R 现有引物：简并碱基迭代与同位点多引物</h1>
  <p class="lede">保持 10 个引物位点和长度不变，只改引物序列：①加简并碱基；②同一位点放 2–3 条引物。设计只用 Greengenes 85% OTU 集和 SILVA 随机 1.2 万条，验证用没参与设计的 SILVA 4.8 万条和 Greengenes 1.2 万条留出序列。</p>
  <p><b>更正说明</b>：本页数字已用修复后的算法重跑。之前追加引物时起始序列没有按“尚未覆盖的序列”加权，多条引物的收益被低估。方案层面的最新结论见 <code>primer_v3v4_design.html</code>。</p>
  <p>数据库是本环境能拿到的 Greengenes 13_8（2013）和 SILVA 128（2016）。SILVA 138 和 GTDB 本环境下载不了，没有评估。这些是计算机模拟，没有做实验，也没有检查对人基因组和线粒体的特异性。</p>
</header>

<h2>一、每个位点的结果</h2>
<p>位置和长度不变，所以扩增子定义、SMURF 切除引物的长度都不用改。覆盖率按“最多 1 个错配、3′ 端 3 个碱基必须匹配”，主要门指 12 个，各占相同权重。下表是 10 个位点的平均值（验证集，不是设计集）。</p>
{t1}
<ul>
<li>只加简并碱基（每个位点仍是 1 条寡核苷酸，简并度限制在 4 倍以内），SILVA 验证集上主要门平均覆盖率从 76% 升到 86%，“最差门”（每个位点最差门覆盖率的平均）从 31% 升到 52%。放宽到 16 倍简并，分别是 90% 和 60%。</li>
<li>总简并数相同（平均约 11 种寡核苷酸）时，每个位点放 3 条引物（各 ≤4 倍简并）比单条 16 倍简并更好：主要门平均 92.4% 对 90.3%，最差门 74% 对 60%（SILVA 验证）。放 2 条引物（平均 7.4 种）与单条 16 倍简并相近（90.1%，63%）。</li>
<li>总简并数越高，每种寡核苷酸在混合物里的有效浓度越低，非特异性扩增的风险也会上升，这里没有评估。</li>
</ul>
<h3>每个引物</h3>
{t3}
<p>R1 的两条引物、R3-F 和 R5-R 的最差门提升仍然有限，说明这几个位点本身的序列在某些门里差异大，加简并解决不了，需要换位置。</p>

<h2>二、整套方案层面的效果</h2>
<p>在近全长序列上，每个区域要求正反引物同时匹配。“信息量”是该序列被扩出的区域所测到的碱基位置的保守度熵之和（按 2×126 bp 读长，去掉引物）。</p>
{t2}
<ul>
<li>5 个区域全部扩出的序列比例从 17%（SILVA）/ 19%（GG 留出）提高到 45–50% / 49–55%（“1 条 ≤4 简并”只到 36% / 39%）。</li>
<li>能扩出至少 3 个区域的序列从 66% / 72% 提高到 83–85% / 93–96%；平均每条序列扩出的区域数从 3.0 增加到 3.8–4.0。</li>
<li>至少扩出 1 个区域的序列在 SILVA 上：只加简并（每位点 1 条）从 86.5% 升到 87–90%；每位点放 2–3 条引物升到 97%–98%，因为追加的引物可以为某些类群（如 Bacteroidetes、Verrucomicrobia）单独设计。</li>
</ul>

<h2>三、具体序列</h2>
<p>多条引物用 | 分隔，反向引物按 5′→3′ 给出。</p>
{t4}

<h2>四、落地时必须改的地方</h2>
<ul>
<li><b>k-mer 数据库必须用新引物重建。</b>现有数据库里“某条参考序列在某区域是否被扩增”是按旧引物（≤2 错配）判定的。新引物让更多序列能扩增，但数据库仍把它们标成“不扩增”，SMURF 就认不出这些序列。<code>smurf5r.build_region_db_from_fasta</code> 已支持每个位点多条引物（传 list）。</li>
<li><b>区域分配</b>：<code>smurf5r.split_to_regions</code> 同样接受每个位点多条引物。MATLAB 版的 <code>get_configs.m</code> / <code>read_unireads_save_split_to_regions.m</code> 目前一个位点只处理一条引物字符串，需要改成逐条展开再拼接。</li>
<li><b>所有变体必须等长</b>：这里每个位点的各条引物长度相同，才能保持“切除固定长度的引物”这一步不变。</li>
<li>引物浓度：多条引物和高简并引物需要按各变体在混合物里的浓度做配比，并用已知组成的标准菌群（mock community）验证。</li>
</ul>

<h2>五、换成更新的数据库</h2>
<p>在你本地拿到 SILVA 138.2 或 GTDB 的 16S FASTA 和分类文件后：</p>
<pre style="overflow-x:auto;background:var(--surface);border:1px solid var(--line);padding:12px;border-radius:8px"><code>import primer_design as pdz, iterate_5R as it
ref0 = pdz.load_reference(GG_DIR)                       # 取 E. coli 坐标参照 ref0.ref_seq
S = pdz.reference_from_sequences(seqs, domain, phylum, ref0.ref_seq)   # seqs/domain/phylum: {{id: ...}}
# 把 S 作为设计与验证集（或再留一部分做验证），调用 it.run(R, S, H, out_prefix)</code></pre>

<footer>脚本：<code>python_5R/iterate_5R.py</code>、<code>python_5R/primer_design.py</code>、<code>python_5R/smurf5r.py</code>；数据：<code>docs/primer_design/iter5R_primers.csv</code>、<code>iter_scheme.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_5R_iteration.html", "w", encoding="utf-8").write(head + body)
print("ok")
