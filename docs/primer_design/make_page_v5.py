"""生成 docs/primer_add_v5.html"""
import re, html, json
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>加入 V5 的评估</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: f"{x * 100:.0f}%"
tup = lambda s: tuple(float(x) for x in re.findall(r"[\d.]+", s))
SC = pd.read_csv(f"{HERE}/v5_site_compare.csv"); PCV = pd.read_csv(f"{HERE}/v5_pair_cov.csv"); DE = pd.read_csv(f"{HERE}/v5_design_eval.csv"); OT = pd.read_csv(f"{HERE}/v5_offtarget.csv"); CB = json.load(open(f"{HERE}/v5_combos.json")); O = pd.read_csv(f"{HERE}/order_list_full_v5.csv")
t1 = tbl([[r.site, int(r.oligos), int(r.expansions), f"{pc(tup(r.GG)[0])} / {pc(tup(r.GG)[1])}", f"{pc(tup(r.SILVA)[0])} / {pc(tup(r.SILVA)[1])}", r.Tm.replace("[", "").replace("]", "")] for _, r in SC.iterrows()], ["位点", "寡核苷酸数", "展开数", "Greengenes：主要门平均 / 最差门", "SILVA：主要门平均 / 最差门", "Tm °C"])
t2 = tbl([[r.pair, f"{pc(tup(r.GG)[0])} / {pc(tup(r.GG)[1])}", f"{pc(tup(r.SILVA)[0])} / {pc(tup(r.SILVA)[1])}", f"{pc(r.Ten_GG)} / {pc(r.Ten_SILVA)}"] for _, r in PCV.iterrows()], ["引物对（产物）", "GG：平均 / 最差门", "SILVA：平均 / 最差门", "支原体 GG / SILVA"])
t3 = tbl([[r.design, f"{r.ideal * 100:.2f}%", f"{r.abs60 * 100:.1f} / {r.abs300 * 100:.1f} / {r['frac0.8'] * 100:.1f}%"] for _, r in DE.iterrows()], ["方案（留出集 2,196 条）", "理想属准确率", "属准确率（三种“扩增失败也算”规则）"])
chg = O[O.order.isin(["SMURF5-A4-F.1", "SMURF5-A4-F.2", "SMURF5-A4-F.3（补充）", "SMURF5-V5-R.1", "SMURF5-V5-R.2"]) & (O.group != "备选")]
chg = chg.drop_duplicates("order", keep="last")
t4 = tbl([[r.order, r.recommend, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}", f"{r.uL_with_sup:.2f}", f"{r.uL_no_sup:.2f}"] for _, r in chg.iterrows()], ["订购名称", "建议", "序列 5′→3′", "nt", "E. coli 位置", "Tm °C（展开范围）", "展开数", "人基因组位点 / 展开序列", "µL 母液（含补充）", "µL 母液（不含补充）"], seqcol=2)
n_oligo = len(O[(O.group != "备选") & (O.recommend != "备选")]); n_exp = int(O[O.group != "备选"].expansions.sum())
t5 = tbl([[f"{f}×{r}", l] for f, r, l in CB], ["引物组合（我们的 V5 版）", "产物长度 bp"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:1040px;max-width:100%;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>把 V5 加进我们的方案</h1>
  <p class="lede">可以加，用成熟的 907R 位点（907–927，2 条寡核苷酸，对数据库匹配 98%）做 V5 的反向引物，V5 的读段由它和 V4 的正向引物（A3-F）配成的 371 bp 产物提供。代价是 V6·V7 的正向引物不能再用 906 位点（和 907R 互补），要换成 967 位点（229 bp，对支原体弱一点，加 1 条补充）。理想属准确率 95.38% → 95.58%，属准确率在两种规则下提高 2 个百分点、一种规则下略降。这是序列层面的计算比较，没有实验验证。</p>
</header>

<h2>一、先核对成熟引物的区域有没有明显问题</h2>
<p>你要求优先参考成熟引物选的区域。我把 Swift SNAP 里的位点（515F、806R 型、907R、967F 型）和我们的位点按同一套规则比：</p>
{t1}
<ul>
<li><b>907R（V5 反向）和 515F（V4 正向）没有明显问题</b>：907R 位点 2 条寡核苷酸 98% / 98%、最差门 95%；SNAP 自己的 907R 只用 1 条就有 97% / 97%、最差门 94% / 95%。515F 只用 1 条（1 个展开序列）就有 98% / 98%、最差门 96% / 92%，比我们的 A3-F（2 条，16 个展开序列，94% / 95%，最差门 78% / 84%）好。</li>
<li><b>806R 型（V4 反向）有明显问题</b>：SNAP 的 V4_r 无简并，73% / 78%，最差门只有 5% / 6%；我们在同一个位点做了简并（A3-R，96% / 93%，最差门 92% / 78%），所以保留我们的版本。967F 型：SNAP 6 条变体合起来 85% / 90%，最差门 36% / 62%，我们的 2 条版本 93% / 94%，最差门 63% / 84%，也保留我们的。</li>
<li><b>结构上互相冲突，不是覆盖率问题</b>：这几个成熟位点在 E. coli 上两两互补，不能放在同一管里：
  <ul>
  <li>907R（907–926）和我们 V6·V7 的 906–927 正向引物（926F 型）互补；</li>
  <li>515F（517–533）和我们 V3 的反向引物 A2-R（516–535）互补，这就是我们 V4 正向引物用 556–576 而不是 515F 的原因；</li>
  <li>V5 的成熟正向引物（785F/799F 型，GGATTAGATACCC…）和我们 V4 的反向引物 A3-R（785–804）互补。</li>
  </ul></li>
</ul>

<h2>二、几种结构的比较</h2>
<ul>
<li><b>方案 B（推荐，改动最小）</b>：保留现有 5 个扩增子，加 907R 型反向引物（2 条）做 V5，V6·V7 正向换成 967 位点（2 条 + 1 条支原体补充）。V5 由 A3-F × V5-R（556–927，371 bp）提供：读段的正向 120 nt 重复 V4，反向读段读 V5（788–906）。</li>
<li><b>方案 C（按成熟引物重组，类似 SNAP）</b>：去掉 V3 的反向引物，用 515F 做 V4 正向，A2-F × 806R 型（338–804）读 V3 前半、515F × 806R 型（287 bp）读 V4、515F × 907R（410 bp）读 V4 前半和 V5。寡核苷酸数更少，但 V3 只读到 433–477，约 69%。</li>
<li><b>方案 D（把 V4 和 V5 合成一个扩增子，去掉 A3-R）</b>：理想准确率反而降到 95.31%（−0.07 个百分点），因为 V4 后半段（约 697–786）没人读了。</li>
</ul>
<p>按窗口法算（每个扩增子两端读段能读到的碱基，含 10 bp 内联 UDP、2×150）的理想准确率：当前推荐 95.35%，方案 B 95.55%，方案 C 95.50%。把完整 V5（822–879）加上去的上限是 +0.22 个百分点，实际增益比这个小，因为 V5 读段只来自较长的产物。</p>
{t3}
<p>方案 B 在两种规则下属准确率提高约 2 个百分点（88.0 / 88.4 → 90.3 / 90.0），第三种规则（≥80% 自身位点）略降 0.4 个百分点（91.0% → 90.6%）；这几种规则对“扩增失败”的处理很敏感，差 1 个百分点以内的不要当作真实差别。</p>

<h2>三、方案 B 的引物对覆盖</h2>
{t2}
<ul>
<li>V5 产物（A3-F × V5-R）：GG 92%、SILVA 93%，最差门 77% / 83%，支原体 89% / 93%，覆盖不错。</li>
<li>V6·V7 换成 967 版后：GG 91%、SILVA 91%（906 版：93% / 94%），支原体 83% / 88%（906 版：89% / 87%）；最差门 68% / 81%（906 版 79% / 87%），所以代价是这一对覆盖略弱，主要在支原体 Greengenes 和个别门上，已经加了 1 条支原体补充（`RTACMCGAARAACCTTACC`）。</li>
</ul>

<h2>四、方案 B 的订购变化</h2>
<p>相对完整订购清单（推荐 18 条），改动：去掉 A4-F（906 正向）1 条，加 A4-F 967 版 2 条 + 1 条支原体补充，加 V5-R 2 条。推荐池 {n_oligo} 条寡核苷酸，展开 {n_exp} 个。完整清单在 <code>docs/primer_design/order_list_full_v5.csv</code>。</p>
{t4}
<ul>
<li>池内所有寡核苷酸（含展开）两两检查，没有严重二聚体或发夹。V5-R 对人基因组的命中很少（2.8 / 0.4 个/展开序列）；线粒体命中 4 个错配（1159，V5-R，没有产物）。</li>
<li>有 V5 后每个位点仍是相同总量（11 个位点 × 2.5 µL，总母液 27.5 µL，加 72.5 µL 水），位点内比例按表。</li>
<li>如果不想加 V5：保持原来的 A4-F（906 版，备选行 <code>SMURF5-A4-F.alt906</code>），不订 V5-R 和 967 版，两者互相冲突。</li>
</ul>

<h2>五、会不会扩出很多小片段</h2>
<p>你提醒的 Swift SNAP 是单管，所有正向引物都能和下游的所有反向引物配对：我按引物位置数了一遍，SNAP 有 20 种可能的 F×R 组合，其中 V3_f（341–357）× V2_r（372–391）只有 51 bp，是引物二聚体大小的产物，被它们脚本里的最小长度 130 过滤掉，浪费读段；V4_f × V4_r 287 bp、V7_f × V8_r 309 bp、V6′_f × V8_r 353 bp、V1_f × V2_r 383 bp、V4_f × V5_r 410 bp 等中等长度的产物彼此竞争，哪些是设计的目标、哪些是副产物，要看 Swift 的说明书。</p>
<p>我们的方案 B 有 18 种组合：设计的目标产物 198 / 229 / 249 / 258 / 288 bp，加上 V5 的 372 bp；所有非目标产物都 ≥467 bp，比目标产物都长，PCR 倾向于优先扩增短的，所以不会出现更短的副产物，不过 V5 的 372 bp 产物会和 A3 的 249 bp 产物竞争同一条正向引物，V5 读段的比例要靠实测调（比如适当降低 A3-R 的比例，提高 V5-R 的比例）。</p>
{t5}

<h2>六、限制</h2>
<ul>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算，不是实际 PCR 效率；数据库是 Greengenes 13_8 和 SILVA 128；没有实验验证。</li>
<li>V5 的读段只来自 371 bp 的产物，且和 249 bp 的 A3 产物共用正向引物 A3-F；真实的读段比例我无法从序列预测，需要用已知组成的样本调整比例。</li>
<li>方案 C 只做了读段窗口的准确率估计，没有做引物设计和二聚体检查。</li>
<li>V6·V7 用 967 版对支原体 Greengenes 覆盖降了约 6 个百分点，已经加补充，但补充没有再单独检查延长后的人基因组脱靶（unextended 版 1.8 个/展开序列，低）。</li>
</ul>
<footer>脚本：<code>python_5R/explore/v5_*.py</code>；数据：<code>docs/primer_design/v5_*.csv</code>、<code>order_list_full_v5.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_add_v5.html", "w", encoding="utf-8").write(head + body)
print("ok")
