"""生成 docs/primer_order_list.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>5 扩增子引物订购清单</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
D = pd.read_csv(f"{HERE}/order_list_5amp.csv")
t1 = tbl([[r["name"], r.amplicon, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", f"{r.GC * 100:.0f}%", int(r.expansions)] for _, r in D.iterrows()],
         ["订购名称", "区域", "序列 5′→3′（IUPAC）", "nt", "E. coli 位置", "Tm °C（展开范围）", "GC", "展开数"], seqcol=2)
mix = D[["name", "site", "n_at_site", "share", "uL_100uM_per_100uL_10x"]]
t2 = tbl([[r["name"], f"{r.share * 100:.0f}%" if r.n_at_site > 1 else "单条", f"{r.uL_100uM_per_100uL_10x:.2f}"] for _, r in mix.iterrows()] + [["水", "", "75.00"], ["合计", "", "100.00"]],
         ["订购名称", "在本位点内的比例", "100 µM 母液（µL）→ 100 µL 的 10× 混合液"])
t3 = tbl([[a, s] for a, s in D.drop_duplicates("amplicon")[["amplicon", "amplicon_span"]].values], ["区域", "扩增子范围（含引物）"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>5 扩增子引物订购清单（均衡后）</h1>
  <p class="lede">14 条寡核苷酸、86 个展开序列、10 个引物位点，覆盖 V1·V2、V3、V4、V6·V7、V8·V9 五个区域。这是计算机设计，没有实验验证，下单前请看第四节的限制。</p>
</header>

<h2>一、订购清单</h2>
{t1}
<ul>
<li>序列用 IUPAC 简并码（R=A/G，Y=C/T，W=A/T，M=A/C，K=G/T，S=G/C，D=A/G/T，B=C/G/T，V=A/C/G），按“混合碱基”合成，不要拆成单个展开序列分别订购。</li>
<li>反向引物的序列已经是 5′→3′ 的引物序列，位置是它在 E. coli 16S（顶链）上占据的范围。</li>
<li>只含基因特异序列，没有接头或 UDI。如果你的建库方案要在引物上加接头，请自己加，并重新检查二聚体和 Tm（我没有评估带接头的版本）。</li>
<li>Tm 是 primer3 在 50 mM Na⁺、2 mM Mg²⁺、0.2 mM dNTP、250 nM 引物下算的各展开序列平均；括号里是各展开序列的范围。</li>
<li>纯化：普通脱盐即可；规模按你的用量，这里没有额外要求。</li>
</ul>
{t3}

<h2>二、建议的混合比例</h2>
<p>规则：<b>每个引物位点（共 10 个）的总量相同</b>，终浓度 250 nM（和 Tm 计算用的一样）；同一位点有多条寡核苷酸时，按它们各自单独能匹配到的序列比例分，取整到 5%。</p>
{t2}
<p>用法：配好 100 µL 的 10× 混合液，每 25 µL PCR 加 2.5 µL（终浓度每个位点 250 nM，每条寡核苷酸取其位点内比例的部分）。简并寡核苷酸里每个展开序列的浓度是该寡核苷酸浓度除以展开数，比如 A1-F 的 8 个展开序列各约 31 nM。</p>

<h2>三、同一位点内的比例是怎么来的</h2>
<ul>
<li>A1-R：45% / 35% / 20%（三条互补，第三条在 5′ 端多出 2 nt，覆盖最少）。</li>
<li>A3-F：90% / 10%（第二条只补充少量序列）。</li>
<li>A5-F：50% / 50%。</li>
<li>这些比例只是起点；位点之间也都用了相同总量，没有考虑各扩增子长度和效率差异。</li>
</ul>

<h2>四、限制，下单前请注意</h2>
<ul>
<li><b>没有实验验证</b>：覆盖率是按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算的，数据库是 Greengenes 13_8 和 SILVA 128。</li>
<li><b>池里的热力学检查</b>：各寡核苷酸平均 Tm 54.8–65.8 °C（含各展开序列的范围是 50.6–67.6 °C）；最热的是 A2-F（65.7）和 A4-R（65.8），最冷的是 A1-F（54.8）；退火温度建议先做梯度（大致 52–60 °C）。池内所有寡核苷酸（含简并展开）两两检查，没有严重二聚体或发夹（A5-F.1 原来有一个临界自二聚体，已经处理，见下一条）。</li>
<li><b>A5-F.1 已改</b>：原序列 <code>GCTRCACRCRTGCTACAAT</code> 的展开序列 <code>GCTGCACGCGTGCTACAAT</code> 含回文核心 GCACGCGTGC，自二聚体 ΔG −9.0 kcal/mol（3′ 端 −1.4）；现在改为 <code>GCKACACRCRTGCTACAAT</code>（第 4 位 R→K，去掉了这个展开序列）。位点覆盖率没有下降（SILVA / GG 88.3% / 93.7% → 90.2% / 94.1%），人基因组位点从 34 个降到 17 个。</li>
<li><b>内联 UDP（每端 10 bp）+ 2×150 读长下，V6·V7（290 bp）和 V8·V9（288 bp）各缺 10 / 8 bp 测不通</b>，缺口在保守区，对属鉴定影响可以忽略，但需要合并读段的流程会丢掉这两个扩增子；缩短版（V6·V7 229 bp 等）没有这个问题，见 <code>primer_shorter_amplicons.html</code>。</li>
<li><b>扩增子长度不同</b>：V3 只有 186 bp，其余 249–290 bp，同一管里短的通常更容易扩出。混合比例没有校正这个偏倚，需要用已知组成的样本（标准菌群）实测后调整各位点的量。</li>
<li><b>覆盖率</b>：5 个扩增子全部扩出的序列，SILVA 约 56%，Greengenes 约 68%；V1·V2 的位点覆盖是最弱的一环。</li>
<li><b>人基因组</b>：A1-F（17 nt）在人基因组里有很多近似位点（≤2 个错配 2,329 个），在 FFPE 等人 DNA 多的样本里会浪费读段；其余引物已延长或移位，4 个扩增子的部分在 ≤2 个错配下没有潜在产物，详见 <code>primer_offtarget.html</code>。</li>
<li>属准确率（算入扩增失败）在这套方案里和 4 个扩增子基本持平，加 V1·V2 的好处主要是多出一个区域的信息，详见 <code>primer_v1v2.html</code>。</li>
</ul>
<footer>脚本：<code>python_5R/explore/order_list.py</code>；数据：<code>docs/primer_design/order_list_5amp.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_order_list.html", "w", encoding="utf-8").write(head + body)
print("ok")
