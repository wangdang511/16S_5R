"""生成 docs/primer_order_list_full_v5.html"""
import re, html
import pandas as pd
import sys; HERE = __file__.rsplit("/", 1)[0]; sys.path.insert(0, HERE); from naming import rename_text, new_order, REGION; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>完整引物订购清单（V5 版）</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
D = pd.read_csv(f"{HERE}/order_list_full_v5.csv"); PCV = pd.read_csv(f"{HERE}/v5_pair_cov.csv"); DE = pd.read_csv(f"{HERE}/v5_design_eval.csv"); S = pd.read_csv(f"{HERE}/v5pool_offtarget_summary.csv").set_index("mm"); OR = pd.read_csv(f"{HERE}/par_base_order.csv")
core = D[D.group != "备选"].copy(); alt = D[D.group == "备选"]
# 精简顺序（来自帕累托页面：22 条基线上逐条删除使综合得分下降最小的顺序）
amap = {"V5-R.2": "SMURF5-V5-R.2", "A4-F.2": "SMURF5-A4-F.2", "A1-R.4(补)": "SMURF5-A1-R.4（补充）", "A4-R.2(补)": "SMURF5-A4-R.2（补充）", "A1-R.3": "SMURF5-A1-R.3", "A2-F.2(补)": "SMURF5-A2-F.2（补充）", "A5-F.2": "SMURF5-A5-F.2", "A1-R.2": "SMURF5-A1-R.2", "A3-F.2": "SMURF5-A3-F.2", "A4-F.3(补)": "SMURF5-A4-F.3（补充）", "A5-F.3(补)": "SMURF5-A5-F.3（补充）"}
step = {amap[r.removed]: int(r.step) for _, r in OR.iterrows() if r.removed in amap}
fm = lambda v: "—" if pd.isna(v) else f"{v:.2f}"
t1 = tbl([["@@N:" + r.order.replace("-", "~") + "@@", r.recommend, r.region, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}", step.get(r.order, "保留"), fm(r.uL_with_sup), fm(r.uL_no_sup)] for _, r in core.iterrows()],
         ["订购名称", "建议", "区域", "序列 5′→3′（IUPAC）", "nt", "E. coli 位置", "Tm °C（展开范围）", "展开数", "人基因组位点 / 展开序列", "精简顺序（1 = 最先可删）", "µL 母液（含补充）", "µL 母液（不含补充）"], seqcol=3)
t2 = tbl([["@@N:" + r.order.replace("-", "~") + "@@", r.recommend, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}"] for _, r in alt.iterrows()], ["订购名称（备选）", "用途", "序列 5′→3′", "nt", "位置", "Tm °C", "展开数", "人基因组位点 / 展开序列"], seqcol=2)
tup = lambda s: tuple(float(x) for x in re.findall(r"[\d.]+", s)); pc = lambda x: f"{x * 100:.0f}%"
t3 = tbl([[r.pair, f"{pc(tup(r.GG)[0])} / {pc(tup(r.GG)[1])}", f"{pc(tup(r.SILVA)[0])} / {pc(tup(r.SILVA)[1])}", f"{pc(r.Ten_GG)} / {pc(r.Ten_SILVA)}"] for _, r in PCV.iterrows()], ["引物对（产物）", "GG：平均 / 最差门", "SILVA：平均 / 最差门", "支原体 GG / SILVA"])
t4 = tbl([[r.design, f"{r.ideal * 100:.2f}%", f"{r.abs60 * 100:.1f} / {r.abs300 * 100:.1f} / {r['frac0.8'] * 100:.1f}%"] for _, r in DE.iterrows()], ["方案（留出集 2,196 条）", "理想属准确率", "属准确率（三种“扩增失败也算”规则）"])
n_oligo = len(core); n_exp = int(core.expansions.sum())
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:1100px;max-width:100%;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>完整引物订购清单（V5 版）</h1>
  <p class="lede">{n_oligo} 条寡核苷酸、{n_exp} 个展开序列、11 个引物位点，覆盖 V1·V2、V3、V4、V5、V6·V7、V8·V9。相比上一版完整清单，V6·V7 正向引物改成 967 位点（229 bp），新增 V5 反向引物 907R 型（2 条）。表里新增“精简顺序”一列：想减少引物数量时，按这个顺序从 1 号开始往下删，影响最小；完整的帕累托前沿见 <code>primer_pareto.html</code>。这是计算机设计，没有实验验证。</p>
</header>

<h2>一、订购清单</h2>
{t1}
<ul>
<li>序列用 IUPAC 简并码，按混合碱基合成（R=A/G，Y=C/T，W=A/T，M=A/C，K=G/T，S=G/C，D=A/G/T，B=C/G/T，V=A/C/G）；反向引物的序列已经是 5′→3′ 的引物序列，位置是它在 E. coli 16S（顶链）上占据的范围；只含基因特异序列，没有接头和 UDP，带接头后的 Tm 和二聚体没有评估。</li>
<li><b>建议</b>：“必订”是核心；“建议”是支原体补充（A2-F.2、A4-R.2、A4-F.3、A5-F.3），人基因组命中低、Tm 合适，对其他菌门没有害处；“可选”的 A1-R.4 只在 SILVA 里提高 V1·V2 的支原体覆盖。</li>
<li><b>混合比例</b>：每个引物位点（共 11 个）总量相同，终浓度 250 nM；位点内有多条时，核心按单独覆盖率分，补充各占 10%。表里是“100 µL 的 10× 混合液里要加多少 µL 的 100 µM 母液”，11 个位点 × 2.5 µL = 27.5 µL，加 72.5 µL 水；每 25 µL PCR 加 2.5 µL。<b>精简时</b>删掉一条后，同位点其余引物按比例重新分配，位点总量不变。</li>
<li><b>扩增子</b>：V1·V2 258、V3 198、V4 249、V5 371（A3-F × V5-R）、V6·V7 229、V8·V9 288 bp。内联 UDP（每端 10 bp）+ 2×150 读长下，V8·V9 缺 8 bp（1363–1370，不是保守区）；V5 产物 371 bp 读的是两端，中间不读，对属鉴定影响不大。</li>
</ul>

<h2>二、备选（选一个版本时才订）</h2>
{t2}
<ul>
<li><b>A4-F（906 版）</b>：不加 V5 时用的 V6·V7 正向引物（SMURF5-A4-F.alt906），和 V5-R 的 907 位点互补，二者不能同时订；订它就不要订 V5-R 和 967 版的 A4-F.1/.2/.3。</li>
<li><b>V3 短版反向引物</b>（A2-R 505–523，186 bp，Tm 63.8 °C）：长版 A2-R（516–535）更好，对支原体覆盖高得多，代价是 Tm 66.9 °C；选短版就不要订长版 A2-R。</li>
<li><b>借鉴 SNAP 的两条</b>：V1_f 可替代 A1-F，覆盖相当，展开序列少 6 个，对人基因组命中 25 个/展开序列（A1-F 291 个）；907R 单条可替代 V5-R.1 和 V5-R.2，少 1 条、少 10 个展开序列，覆盖 97% / 97%（我们 2 条 98% / 98%）。</li>
</ul>

<h2>三、评估（V5 版）</h2>
{t3}
{t4}
<ul>
<li>池内所有寡核苷酸（含展开）两两检查：没有严重二聚体或发夹；平均 Tm 54.8–66.9 °C，退火温度建议先做梯度（大致 52–60 °C）。</li>
<li>人基因组（hg19，整池 22 条）：≤2 个错配 {int(S.loc[2].loci):,} 个位点、{int(S.loc[2].products)} 个潜在产物；≤3 个错配 {int(S.loc[3].products)} 个潜在产物；线粒体命中都是 3–4 个错配（A1-R.1 在 ND5 基因 2 个错配，没有成对的反向命中），没有潜在产物。命中最多的仍是 A1-F（17 nt，291 个/展开序列），借鉴 SNAP 的 V1_f 可以缓解。</li>
<li>引物位点组合：18 种 F×R 组合里，设计的目标产物是 198、229、249、258、288 bp 和 V5 的 372 bp，所有非目标产物都 ≥467 bp，比目标产物都长，不会出现更短的副产物；V5 的 372 bp 产物和 A3 的 249 bp 产物共用正向引物 A3-F，竞争，V5 读段的比例要靠实测调整（比如降低 A3-R、提高 V5-R 的比例）。</li>
<li>属准确率在两种规则下比不加 V5 的版本高约 2 个百分点，一种规则略降 0.4 个百分点；这些规则对“扩增失败”的处理很敏感，1 个百分点以内的差别不要当真。</li>
</ul>

<h2>四、限制</h2>
<ul>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算，不是实际 PCR 效率；数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本；没有实验验证。</li>
<li>混合比例没有校正扩增子长度差异引起的扩增偏倚，需要用已知组成的样本实测后调整。</li>
<li>V6·V7 用 967 版对支原体 Greengenes 覆盖降了约 6 个百分点（已加 1 条补充）。</li>
<li>人基因组检查只用序列匹配规则，没有热力学；A1-F 在人 DNA 多的样本里会浪费读段。</li>
</ul>
<footer>脚本：<code>python_5R/explore/v5_*.py</code>、<code>par*.py</code>；数据：<code>docs/primer_design/order_list_full_v5.csv</code>、<code>par_*.csv</code>。</footer>
</main></div>
"""
body = rename_text(body)
body = re.sub(r"@@N:(.*?)@@", lambda m: new_order(m.group(1).replace("~", "-")), body)
if "@@MAP@@" in body: body = body.replace("@@MAP@@", MAPHTML)
open(DOCS + "/primer_order_list_full_v5.html", "w", encoding="utf-8").write(head + body)
print("ok")
