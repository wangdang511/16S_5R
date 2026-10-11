"""生成 docs/primer_vs_swift_snap.html"""
import re, html, json
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>与 Swift SNAP 的比较</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: "—" if pd.isna(x) else f"{x * 100:.0f}%"
P = pd.read_csv(f"{HERE}/snap_positions.csv"); C = pd.read_csv(f"{HERE}/snap_site_cov.csv"); A = pd.read_csv(f"{HERE}/snap_all_sites.csv"); SS = pd.read_csv(f"{HERE}/snap_scheme_summary.csv").set_index("scheme")
CMP = pd.read_csv(f"{HERE}/snap_compare.csv"); PC = pd.read_csv(f"{HERE}/snap_pair_cov.csv").set_index("scheme"); OS = pd.read_csv(f"{HERE}/snap_offtarget_summary.csv"); OP = pd.read_csv(f"{HERE}/snap_offtarget_per_oligo.csv")
W = json.load(open(f"{HERE}/snap_windows.json")); PR = pd.read_csv(f"{HERE}/snap_products.csv"); PE = pd.read_csv(f"{HERE}/snap_pe_compare.csv"); O = pd.read_csv(f"{HERE}/order_list_full.csv")
five = [("R1-F", "F", 103, 120), ("R1-R", "R", 314, 332), ("R2-F", "F", 338, 355), ("R2-R", "R", 519, 536), ("R3-F", "F", 685, 702), ("R3-R", "R", 908, 927), ("R4-F", "F", 944, 964), ("R4-R", "R", 1087, 1104), ("R5-F", "F", 1175, 1193), ("R5-R", "R", 1374, 1391)]
VR = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}
# ---- 示意图
X0, XW = 60, 900; sx = lambda p: X0 + p / 1542 * XW
def track(y, title, prim, wins, col_f, col_r):
    s = f'<text x="4" y="{y + 18}" font-size="11" fill="currentColor">{html.escape(title)}</text>'
    for a, b in wins: s += f'<rect x="{sx(a):.1f}" y="{y + 26}" width="{max(1.5, sx(b) - sx(a)):.1f}" height="7" fill="#8ab4e8" opacity="0.7"/>'
    for o, a, b in prim:
        if o == "F": s += f'<polygon points="{sx(a):.1f},{y + 8} {sx(b) + 3:.1f},{y + 19} {sx(a):.1f},{y + 30}" fill="{col_f}" opacity="0.9"/>' if False else f'<rect x="{sx(a):.1f}" y="{y + 12}" width="{max(2.5, sx(b) - sx(a)):.1f}" height="9" fill="{col_f}"/>'
        else: s += f'<rect x="{sx(a):.1f}" y="{y + 38}" width="{max(2.5, sx(b) - sx(a)):.1f}" height="9" fill="{col_r}"/>'
    return s
svg = f'<svg viewBox="0 0 980 270" width="100%" role="img" aria-label="三种方案的引物位置和读段窗口" style="max-width:980px">'
for k, (a, b) in VR.items():
    svg += f'<rect x="{sx(a):.1f}" y="14" width="{sx(b) - sx(a):.1f}" height="14" fill="#c8a24a" opacity="0.35"/><text x="{(sx(a) + sx(b)) / 2:.1f}" y="11" font-size="9" text-anchor="middle" fill="currentColor">{k}</text>'
svg += f'<line x1="{sx(1):.1f}" y1="30" x2="{sx(1542):.1f}" y2="30" stroke="currentColor" opacity="0.4"/>'
snap_prim = [("F" if (r["name"].endswith("_f") or "_f" in r["name"]) else "R", r.start, r.end) for _, r in P.drop_duplicates(["start", "end", "strand"]).iterrows()]
svg += track(40, "Swift SNAP", snap_prim, W["snap"], "#2a7f62", "#c0562d")
f5 = [(o, a, b) for _, o, a, b in five]; svg += track(105, "5R 原方案", f5, W["five"], "#2a7f62", "#c0562d")
ours = []
for _, r in O[O.group != "备选"].iterrows():
    a, b = [int(x) for x in r.pos.replace("–", "-").split("-")]; ours.append((r.orient, a, b))
svg += track(170, "我们的推荐设计", ours, W["ours"], "#2a7f62", "#c0562d")
svg += '<text x="60" y="258" font-size="10" fill="currentColor">绿色 = 正向引物位点（上），橙色 = 反向引物位点（下），蓝色条 = 每条读段在引物之后能读到的窗口（SNAP、5R 每读 130 nt；我们含 10 bp 内联 UDP 的 2×150）；金色 = 高变区（E. coli 常用划分）。</text></svg>'
t_pos = tbl([[r["name"], r.seq, r.nt, f"{r.start}–{r.end}", "正向" if r.strand == "+" else "反向", r.mismatch] for _, r in P.iterrows()], ["SNAP 引物", "序列 5′→3′（去掉前缀 ^）", "nt", "E. coli 位置", "方向", "与 E. coli 错配数"], seqcol=1)
t_ours = tbl([["5R 原方案", " · ".join(f"{n} {a}–{b}" for n, o, a, b in five)], ["Swift SNAP", "V1_f 9–27 · V2_r 372–391 · V3_f 341–357 · V4_f 517–533 · V4_r 785–803 · V5_r 907–926 · V6_f 967–985 · V6'_f 1055–1070 · V7_f 1099–1114 · V8_r 1390–1407 · V9_r 1492–1507"],
                ["我们的推荐设计", "A1 8–24 / 246–265 · A2 338–356 / 516–535 · A3 556–576 / 785–804 · A4 906–927 / 1177–1195 · A5 1223–1241 / 1492–1510（F / R）"]], ["方案", "引物位点（E. coli 位置）"])
t_win = tbl([[r.scheme, int(r.bases), f"{r.ideal_acc * 100:.2f}%"] + [f"{r[k] * 100:.0f}%" for k in VR] for _, r in pd.concat([CMP[CMP.scheme.str.startswith('5R')], PE]).iterrows()], ["方案", "读到的碱基数", "理想属准确率（所有读段都得到）"] + [f"{k}" for k in VR])
t_site = tbl([[s, int(r.sites), int(r.oligos), int(r.exp), f"{r.GGm * 100:.0f}% / {r.SILm * 100:.0f}%", f"{r.GGmin * 100:.0f}% / {r.SILmin * 100:.0f}%", f"{r.Ten_GG * 100:.0f}% / {r.Ten_SIL * 100:.0f}%", int(r.n_weak)] for s, r in SS.iterrows()],
             ["方案", "位点数", "寡核苷酸数", "展开数", "各位点平均覆盖（主要门平均）GG / SILVA", "最差位点的最差门 GG / SILVA", "支原体平均覆盖 GG / SILVA", "平均覆盖 <80% 的位点数（GG）"])
sn = A[A.scheme == "Swift SNAP"]
t_snap = tbl([[r.site, int(r.oligos), int(r.expansions), f"{pc(r.GG_mean)} / {pc(r.GG_min)}", f"{pc(r.SILVA_mean)} / {pc(r.SILVA_min)}", f"{pc(r.Ten_GG)} / {pc(r.Ten_SILVA)}"] for _, r in sn.iterrows()],
             ["SNAP 位点", "寡核苷酸数", "展开数", "GG：主要门平均 / 最差门", "SILVA：主要门平均 / 最差门", "支原体 GG / SILVA"])
pp = PC
t_prod = tbl([[r.F, r.R, int(r.length), f"{r.GG_all * 100:.0f}% / {r.SILVA_all * 100:.0f}%", f"{r.GG_phylum_mean * 100:.0f}% / {r.SILVA_phylum_mean * 100:.0f}%"] for _, r in PR.iterrows()], ["正向引物", "反向引物", "产物长度 bp", "两条引物都能结合的序列比例 GG / SILVA（全部序列）", "主要门平均 GG / SILVA"])
t_pair = tbl([[s, int(r.pairs), f"{r.GG_pair_mean * 100:.0f}% / {r.SILVA_pair_mean * 100:.0f}%", f"{r.GG_pair_worst * 100:.0f}% / {r.SILVA_pair_worst * 100:.0f}%", f"{r.GG_ge_half * 100:.0f}% / {r.SILVA_ge_half * 100:.0f}%", f"{r.GG_all * 100:.0f}% / {r.SILVA_all * 100:.0f}%"] for s, r in pp.iterrows()],
             ["方案", "引物对数", "各引物对平均覆盖 GG / SILVA", "最差引物对的最差门 GG / SILVA", "至少一半引物对都扩出 GG / SILVA", "全部引物对都扩出 GG / SILVA"])
s_snap = OS[OS.set == "SNAP"].set_index("mm"); s_rec = OS[OS.set == "REC2"].set_index("mm")
hp = OP[(OP.set == "SNAP") & (OP.per_exp >= 20)]
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:1040px;max-width:100%;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 方案比较</div>
  <h1>5R、我们的 5 扩增子方案、Swift 16S SNAP 的引物位置和优缺点</h1>
  <p class="lede">三个方案的引物位置很不一样：5R 的引物集中在 V2、V3、V5、V6、V8 附近，缺 V1、V4、V7、V9；Swift SNAP 用 11 组引物（23 条寡核苷酸），读段窗口在 PE250 下覆盖 V1 到 V9，在 PE150 下漏掉 V2 中段和 V8 前半；我们的方案 10 个位点（18 条寡核苷酸，118 个展开序列），覆盖 V1、V2、V3、V4、V6、V7、V8、V9，缺 V5。在引物对序列的匹配上，我们的方案最好、SNAP 居中、5R 最差；SNAP 的优点是区域覆盖面最广（PE250 下 V1 到 V9 全读到）、引物简单（展开序列只有 34 个）、有现成的商品化流程，主要短板是有几个引物（V2_r、V4_r、V3_f、V7_f、V6_f）对一些菌门匹配很差，并且需要 PE250。这些都是序列层面的计算比较，没有实验验证。</p>
</header>

<p style="border:1px solid var(--warn);background:var(--surface);padding:10px 14px;border-radius:8px"><b>更正（读长）：</b>本页第一版把 SNAP 的读段按 130 nt 算（它脚本里的默认 READLEN=130，对应 PE150）。你说 SNAP 实际用 PE250（读段约 230 nt），这样 SNAP 的两条读段几乎把 V1 到 V9 全读到（理想属准确率 95.86%），“V2 只读到 21%、V8 71%”只在 PE150 下成立。下面的读段覆盖表同时列出 SNAP 在 PE250 和 PE150 下的结果；而我们的方案只按 PE150（含 10 bp 内联 UDP）算。</p>
<h2>一、引物位置（E. coli 16S 编号）</h2>
{svg}
{t_ours}
<p>SNAP 的引物位置是我把引物序列比对到 E. coli 16S 得到的（前缀 ^ 去掉；V2_r 的 8 条是同一个位点上的变体，位置在 372–391，是 V3 起点之后的一段，<b>不是</b> V2 的末端；V6_f 的 6 条在 967–985）：</p>
{t_pos}

<h2>二、读段覆盖哪些区域</h2>
<p>两端各读一条，每条读段从引物之后读一段：SNAP 在 PE250 下约 230 nt（也列了 PE150 下的 130 nt 作对照），5R 假设 130 nt，我们的推荐设计按内联 UDP 10 bp、2×150 来算。这个窗口只取决于引物位置，不依赖引物是怎么配对的：</p>
{t_win}
<ul>
<li><b>SNAP，PE250</b>：读到约 1,385 个碱基，V1 到 V9 全部读到，理想属准确率 95.86%，是三个方案里最高的。<b>SNAP，PE150</b>：读到约 1,240 个碱基，V2 只有 21%（V1_f 往后读到 157，V2_r 往前读到 242，中间 158–241 没读到）、V8 71%、V3 85%，理想属准确率 95.30%。</li>
<li><b>5R</b>：约 860 个碱基，缺 V1、V4、V7、V9；理想属准确率 94.05%。</li>
<li><b>我们（PE150）</b>：约 1,070 个碱基，缺 V5，V4 读到 99%；理想属准确率 95.35%；加了 V5 的方案 B 约 1,140 个碱基，95.55%。同样是 PE150，我们比 SNAP 高 0.05–0.25 个百分点；SNAP 用 PE250 才能读全，比我们高 0.3 个百分点（方案 B）到 0.5 个百分点（原推荐）。</li>
<li>SNAP 的引物对（哪条正向配哪条反向）在你给的文件里没有，我只能按名字和位置推测（V1_f–V2_r、V3_f–V4_r、V4_f–V5_r、V6_f 和 V6′_f–V8_r、V7_f–V9_r）；真实的产物长度要看 Swift 的产品说明书。V3_f（341–357）和 V2_r（372–391）只相隔 14 个碱基，如果它们配成一对，只会得到 51 bp 的短产物。</li>
</ul>

<h2>三、引物序列对数据库的匹配（覆盖率）</h2>
<p>按“最多 1 个错配、3′ 端 3 个碱基匹配”的规则，用 Greengenes 13_8 和 SILVA 128（留出集，每个主要门等权）：</p>
{t_site}
<p>SNAP 各位点：</p>
{t_snap}
<ul>
<li><b>SNAP 的弱点</b>：V2_r（8 条变体合起来 GG 60% / SILVA 68%，最差门只有 24% / 36%，支原体 20% / 55%）、V4_r（73% / 78%，最差门 5% / 6%）、V6_f（6 条合起来 85% / 90%，最差门 36% / 62%）、V3_f（平均 87% / 89%，但最差门只有 3% / 28%，是浮霉菌门 Planctomycetes）和 V7_f（平均 86% / 87%，最差门 6% / 21%）。这些引物几乎没有简并，是经典的通用引物。</li>
<li><b>SNAP 的优点</b>：V1_f、V4_f、V5_r、V8_r、V9_r 都在 94–99%，最差门也在 88–96%。</li>
<li><b>我们</b>：平均覆盖最高（GG 95% / SILVA 94%），最差位点的最差门也在 69% / 78%，没有平均覆盖低于 80% 的位点；代价是 118 个展开序列（SNAP 34 个，5R 14 个）。</li>
<li><b>5R</b>：10 个位点里有 7 个平均覆盖低于 80%，几个位点的最差门是 0%。</li>
</ul>
<p>把引物配成对之后的覆盖（SNAP 的配对是我推测的）：</p>
{t_pair}
<ul>
<li>SNAP 6 对引物里“全部扩出”的序列 Greengenes 约 47%、SILVA 约 43%（我们 5 对：72% / 58%；5R：15% / 14%）；至少一半的引物对扩出（GG / SILVA）：SNAP 98% / 83%，我们 97% / 83%，5R 69% / 62%。</li>
<li>对支原体（Tenericutes）：SNAP 各位点平均 GG 82% / SILVA 87%，我们 95% / 94%，5R 59% / 66%。SNAP 的 V2_r 和 V6_f 对支原体最弱。</li>
</ul>

<h2>四、人基因组和线粒体脱靶</h2>
<ul>
<li>SNAP：≤2 个错配的位点 {int(s_snap.loc[2].loci)} 个，≤3 个错配的潜在产物 {int(s_snap.loc[3].products)} 个；我们的推荐池：≤2 个错配 {int(s_rec.loc[2].loci):,} 个位点，≤3 个错配 {int(s_rec.loc[3].products)} 个潜在产物。SNAP 的寡核苷酸少、展开序列少，所以人基因组命中明显少。</li>
<li>SNAP 里每个展开序列命中人基因组最多的是 V3_f（218 个）、V6′_f（64 个）、V9_r（47 个，16 nt）、V1_f（25 个）；我们最多的是 A1-F（291 个）。</li>
<li>线粒体：SNAP 的 V9_r 在 12S rRNA 的 1556 位只有 1 个错配（只有这一个，所以单独不会形成产物）；我们对应的位点是 3–4 个错配。</li>
</ul>

<h2>五、单管里 SNAP 会不会扩出短片段：会</h2>
<p>SNAP 所有引物在同一管（按你的说法），任何一条正向引物都能和下游的任何一条反向引物配成产物。我用 Greengenes 和 SILVA 的留出序列（主要门）检查每对引物在序列上能不能同时结合（≤1 个错配、3′ 端 3 个碱基匹配），产物长度按 E. coli 位置算：</p>
{t_prod}
<ul>
<li><b>确认会有短片段</b>：V3_f（341–357）和 V2_r（372–391）两条引物的 3′ 端面对面，中间只隔 14 个碱基，产物只有 <b>51 bp</b>（含两条引物，中间只剩 14 个碱基的模板）。它不是引物二聚体，是模板特异的产物，在 SILVA 58%、Greengenes 61% 的序列上两条引物都能结合，所以绝大多数细菌都会形成这个 51 bp 产物；产物这么短，扩增效率最高，会优先消耗 V3_f 和 V2_r。</li>
<li>这个 51 bp 产物在他们的流程里会被丢掉：cutadapt 去引物后只剩 14 bp，而 <code>--minimum-length</code> 设成 130（READLEN），所以这些读段被浪费，不会进入分析；浪费多少比例取决于实际扩增效率，我无法从序列预测。</li>
<li>其他产物（287–590 bp）也会同时形成：V4_f × V4_r（287 bp）、V7_f × V8_r（309 bp）、V6′_f × V8_r（353 bp）、V1_f × V2_r（383 bp）、V4_f × V5_r（410 bp）等，在 57–95% 的序列上都能形成，彼此竞争同一条引物。</li>
<li>SNAP 引物池两两检查没有严重的引物二聚体（ΔG ≤ −9 kcal/mol 或 3′ 端 ≤ −6 的），所以短片段主要来自上面这种“引物位点靠得近”的模板产物，不是引物二聚体。</li>
<li>和我们的对比：我们的方案所有非目标产物都 ≥467 bp，比目标产物（198–290 bp，加 V5 后还有 372 bp）都长，不会出现比目标更短的产物。</li>
<li><b>局限</b>：这是纯序列层面的推断，没有考虑 Swift 实际可能用的引物浓度、引物配对（比如把 V3_f 和 V2_r 分在不同的反应里）、退火条件或其他设计。你给的文件里没有这些信息，真实情况要看 Swift 的产品说明书或实测。</li>
</ul>

<h2>六、分析流程的差别</h2>
<ul>
<li><b>16S-SNAPP（zip 里的流程）</b>：cutadapt 按引物文件两端去引物，每条读段截到 130 nt，DADA2（池化样本推断，R1/R2 直接拼接 <code>justConcatenate=TRUE</code>，不要求重叠，去嵌合体），RDP 11.5 做 BLAST 找模板，把来自不同区域的读段关联到同一个参考序列上，再合并计数，用 RDP Classifier 分类到属及以上。</li>
<li><b>5R/SMURF</b>（我们这个仓库）：每个区域一个 k-mer 数据库，把拼接后的读段（R1 + 反向互补 R2）用期望最大化分配到物种，不依赖 DADA2。</li>
<li>两个流程的共同点是读段对直接拼接、不要求重叠，所以我们之前讨论的“V8·V9 缺 8 bp 缺口”对这两个流程都不是问题；需要合并读段的流程才会受影响。</li>
<li>SNAPP 依赖 RDP 11.5 数据库和 RDP Classifier，只能到属；SMURF 的结果取决于你建的区域数据库。引物换了，SMURF 的数据库需要重建，SNAPP 不需要但要用它自己的引物文件。</li>
</ul>

<h2>七、各自的优缺点</h2>
{tbl([["Swift SNAP", "区域覆盖面最广（PE250 下 V1 到 V9 全部读到）；引物简单，展开序列只有 34 个；人基因组脱靶少；有商品化试剂盒和现成流程（含 DADA2 + RDP）；引物是经典通用引物，有大量文献使用。", "需要 PE250 才能读全（PE150 下 V2 中段、V8 前半读不到）；V2_r、V4_r、V6_f、V3_f、V7_f 对一些菌门匹配很差（最差门 3–36%），对支原体弱；引物位点之间有重叠或靠得很近（V3_f 和 V2_r 只隔 14 个碱基）；V9_r Tm 只有 49 °C；流程只到属，依赖 RDP 数据库。"],
          ["5R", "引物数最少（10 个，展开序列 14 个）；已有协议和历史数据，结果可比；SMURF 流程能分到物种。", "引物是 6 年前的，对今天的数据库匹配最差（平均覆盖 74–77%，7 个位点低于 80%）；缺 V1、V4、V7、V9；对支原体 V6、V8 扩增子几乎扩不出。"],
          ["我们的 5 扩增子", "引物匹配最好（平均覆盖 94–95%，没有低于 80% 的位点），支原体覆盖高；位点不重叠、扩增子短，读段能覆盖；人基因组脱靶经过检查和处理；有二聚体和 Tm 检查；理想属准确率最高或相当（95.35%，SNAP 95.30%）。", "没有 V5；寡核苷酸多、简并度高（18 条、118 个展开序列），Tm 范围 54.8–66.9 °C；没有实验验证；需要重建 SMURF 数据库；A1-F 对人基因组命中多。"]], ["方案", "优点", "缺点"])}

<h2>八、限制</h2>
<ul>
<li>SNAP 的引物是文件里的序列，我按 E. coli 序列推位置；引物是怎么配对、产物多长、各引物浓度比例，文件里没有，我没有假设来源，只在表里标了“假设配对”。</li>
<li>覆盖率按我们统一的序列匹配规则计算，不是实际 PCR 效率；数据库是 Greengenes 13_8 和 SILVA 128，不是 SNAP 用的 RDP 11.5；没有实验验证。</li>
<li>PE250 下的读段长度是我按“读长 250 减去引物约 20 nt”估计的 230 nt（200 nt 的结果相同，因为窗口都覆盖到）；真实读长要看数据质量。</li>
<li>理想准确率是代理指标（属级最近邻，假设所有读段都得到），不代表 SNAPP 或 SMURF 的实际表现。</li>
<li>我们的 A2-F 引物在这次比较里发现并修正了一个问题（见完整订购清单页面顶部的更正）。</li>
</ul>
<footer>脚本：<code>python_5R/explore/snap*.py</code>；数据：<code>docs/primer_design/snap_*.csv</code>、<code>Swift_16S_SNAP_primers_v2.fasta</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_vs_swift_snap.html", "w", encoding="utf-8").write(head + body)
print("ok")
