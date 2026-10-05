"""生成 docs/primer_tenericutes.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>支原体各位点覆盖</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: "—" if pd.isna(x) else f"{x * 100:.0f}%"
S = pd.read_csv(f"{HERE}/ten_site_cov.csv"); A = pd.read_csv(f"{HERE}/ten_amp_cov.csv"); G = pd.read_csv(f"{HERE}/ten_genus_cov.csv"); AD = pd.read_csv(f"{HERE}/ten_addon.csv"); AA = pd.read_csv(f"{HERE}/ten_addon_amp.csv")
order = ["A1-F", "A1-R", "A2-F", "A2-R", "A3-F", "A3-R", "A4-F", "A4-F(原 906–927)", "A4-R", "A4-R(原)", "A5-F", "A5-R"]
t1 = tbl([[r.site, int(r.oligos), f"{pc(r.GG_cov)}（{int(r.GG_n_valid)}）", f"{pc(r.SILVA_cov)}（{int(r.SILVA_n_valid)}）"] for _, r in S.set_index("site").loc[order].reset_index().iterrows()],
         ["位点（当前 V6·V7 缩短版）", "寡核苷酸数", "Greengenes 覆盖（序列数）", "SILVA 覆盖（序列数）"])
t2 = tbl([[r.site, int(r.oligos), f"{pc(r.GG_cov)}", f"{pc(r.SILVA_cov)}"] for _, r in S[S.site.str.startswith("5R")].iterrows()], ["5R 原引物位点", "寡核苷酸数", "Greengenes 覆盖", "SILVA 覆盖"])
t3 = tbl([[r.amplicon, pc(r.GG), pc(r.SILVA)] for _, r in A.iterrows()], ["扩增子", "Greengenes（两端都匹配）", "SILVA（两端都匹配）"])
g = G[G.genus == "Mycoplasma"]
t4 = tbl([[r.site, int(r.n), f"{int(r.covered)}（{r.covered / r.n * 100:.0f}%）"] for _, r in g.iterrows()], ["位点", "Greengenes 里的支原体属（Mycoplasma）序列数", "被覆盖"])
t5 = tbl([[r.site, r.addon, int(r.nexp), r.Tm, f"{pc(r.GG_before)} → {pc(r.GG_after)}", f"{pc(r.SILVA_before)} → {pc(r.SILVA_after)}"] for _, r in AD.iterrows()], ["位点", "补充寡核苷酸 5′→3′", "展开数", "Tm °C", "Greengenes 覆盖", "SILVA 覆盖"], seqcol=1)
t6 = tbl([[r.amplicon, f"{pc(r.GG_现有)} → {pc(r.GG_加补充)}", f"{pc(r.SILVA_现有)} → {pc(r.SILVA_加补充)}"] for _, r in AA.iterrows()], ["扩增子", "Greengenes", "SILVA"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>支原体（Tenericutes）在各个引物位点的覆盖</h1>
  <p class="lede">当前设计对支原体已经比 5R 好得多（5R 有两个位点几乎完全匹配不上），但仍有 4–5 个位点覆盖率偏低，最弱的是 V3 的反向引物（A2-R）和 V6·V7 缩短版的正向引物（A4-F）。每个弱位点加 1 条补充寡核苷酸，能把支原体覆盖率提到 89–98%。这是序列层面的模拟，没有实验验证。</p>
</header>

<h2>一、方法</h2>
<ul>
<li>用 Greengenes 13_8 里全部 Tenericutes 序列（202 条）和 SILVA 128 里全部 Tenericutes 序列（355 条，门由朴素贝叶斯按 Greengenes 推断，不是 SILVA 自己的分类），每个位点只算该位点有数据的序列。覆盖按“最多 1 个错配、3′ 端 3 个碱基匹配”算。</li>
<li>Greengenes 的支原体几乎都是 Mycoplasma（92 条）、Candidatus Phytoplasma（47）、Asteroleplasma（19）、Acholeplasma（19）、Candidatus Hepatoplasma（8）、Ureaplasma（7）。这个门里序列数少，百分比的误差大。</li>
</ul>

<h2>二、当前设计各位点的支原体覆盖</h2>
{t1}
<ul>
<li><b>最弱的位点</b>：A2-R（Greengenes 73%，SILVA 86%）、A4-F 缩短版（74% / 88%）、A5-F（85% / 77%）、A1-R（90% / 79%）、A4-R（89% / 88%）。</li>
<li>原来的 V6·V7 正向引物（906–927）对支原体几乎全覆盖（96% / 96%），缩短版 967–985 这个新位点掉到 74% / 88%，所以你重视支原体时，保留原来的 906–927 是对的；原来的 A4-R 版本也很弱（73% / 80%），改了一个碱基后提高到 89% / 88%。</li>
<li>A1-F、A5-R 在 SILVA 里有数据的序列很少（125、93 条），这两个数字不稳。</li>
</ul>

<h2>三、扩增子级别（两端都匹配）</h2>
{t3}
<ul>
<li>当前设计的 5 个扩增子里，支原体 Greengenes 覆盖率 67%（V6·V7 缩短版）–87%（V1·V2）；用原来的 906 正向引物，V6·V7 提到 85% / 84%。</li>
<li>5 个扩增子里，Greengenes 的支原体有 77% 至少能扩出 3 个扩增子，99% 至少扩出 1 个；SILVA 分别是 72% 和 98%。</li>
<li>5R 原引物（下表）对支原体的 R4（约 944–1104）和 R5（约 1175–1391）扩增子基本扩不出来（Greengenes 6.5% 和 1.6%），整体覆盖弱得多：</li>
</ul>
{t2}

<h2>四、没有被覆盖的是哪些属</h2>
{t4}
<ul>
<li>各位点没被覆盖的主要是 Mycoplasma 本身（A2-R 漏 23/88，A4-F 漏 16/89，A5-F 漏 18/92），其次是 Candidatus Phytoplasma（植物病原，A4-F 漏 21/46）和 Asteroleplasma（A2-R 漏 19/19，A4-R 漏 12/19）。</li>
<li>Ureaplasma 在 A5-F 全漏（7/7）。Acholeplasma 基本都能覆盖。</li>
<li>也就是说，如果你关心的是人体或细胞培养污染里的支原体（Mycoplasma、Ureaplasma），Asteroleplasma 和植物病原的 Phytoplasma 漏掉影响较小，但 Mycoplasma 本身仍漏了 7–26%。</li>
</ul>

<h2>五、补一条针对性的寡核苷酸能提高多少</h2>
<p>对每个弱位点，只拿“当前引物没覆盖的支原体序列”再设计一条简并引物（≤8 个展开序列，≤3 个简并位点），加到该位点：</p>
{t5}
{t6}
<ul>
<li>各位点的支原体覆盖率提高到 89–98%，V3、V8·V9 的扩增子覆盖提高最多（V3：Greengenes 72% → 96%，SILVA 84% → 96%；V8·V9：83% → 96%，72% → 89%）。</li>
<li>代价：每个弱位点 +1 条寡核苷酸、+8 个展开序列（A2-F 补 4 个）；补充引物的 Tm 偏低（53–60 °C），加的是不同的序列，所以会改变混合比例。<b>这些补充寡核苷酸我还没有检查二聚体、人基因组脱靶和对其他菌门的影响</b>，只是对支原体序列的覆盖模拟，不能直接下单。</li>
<li>V4（A3-F / A3-R）已经比较好（90–95%），不需要补。</li>
</ul>

<h2>六、限制</h2>
<ul>
<li>覆盖率是序列匹配规则的结果，不是实际 PCR 效率；序列数少，数值有几个百分点的误差。</li>
<li>Greengenes 的支原体里 Phytoplasma 占比较高（47/202），可能和你的样本不一致；如果你只关心 Mycoplasma 和 Ureaplasma，需要单独看这两个属（表 4）。</li>
<li>只看了支原体，没有评估补充引物对其他菌门的影响。</li>
</ul>
<footer>脚本：<code>python_5R/explore/ten1.py</code>、<code>ten2.py</code>、<code>ten3.py</code>；数据：<code>docs/primer_design/ten_*.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_tenericutes.html", "w", encoding="utf-8").write(head + body)
print("ok")
