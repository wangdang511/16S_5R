"""生成 docs/primer_order_list_levels.html"""
import re, html
import numpy as np, pandas as pd
import sys; HERE = __file__.rsplit("/", 1)[0]; sys.path.insert(0, HERE); from naming import rename_text, new_order, REGION; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>12、16、24 条引物订购清单</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: f"{x * 100:.1f}%"; pc0 = lambda x: f"{x * 100:.0f}%"
J = pd.read_csv(f"{HERE}/lvl4_J.csv"); J["level"] = J.level.astype(str); J = J.set_index("level"); AMP = pd.read_csv(f"{HERE}/lvl4_amp.csv"); AMP["level"] = AMP.level.astype(str); POOL = pd.read_csv(f"{HERE}/lvl4_pool.csv"); POOL["level"] = POOL.level.astype(str); POOL = POOL.set_index("level"); OT = pd.read_csv(f"{HERE}/lvl4_offtarget.csv"); OT["level"] = OT.level.astype(str); W = pd.read_csv(f"{HERE}/order_list_levels.csv"); TV = pd.read_csv(f"{HERE}/lvl4_ten_valid.csv"); TV["level"] = TV.level.astype(str)
ot = lambda lv, mm, c: int(OT[(OT.level == lv) & (OT.mm == mm)][c].iloc[0])
all6 = lambda lv, db: float(AMP[(AMP.level == lv) & (AMP.db == db) & (AMP.amplicon == '全部 6 个扩增子')]['mean'].iloc[0])
LB = {'24': '24 条', '16': '16 条', '12v': '12 条 + V1_f', '12': '12 条'}
t1 = tbl([[LB[lv], int(POOL.loc[lv].expansions), f"{J.loc[lv].J:.4f}（{J.loc[lv].J_lo:.4f}–{J.loc[lv].J_hi:.4f}）", pc(J.loc[lv].accf), pc(J.loc[lv]["cov"]), pc(J.loc[lv].ten), f"{all6(lv, 'GG') * 100:.0f}% / {all6(lv, 'SILVA') * 100:.0f}%", f"{POOL.loc[lv].Tm_min}–{POOL.loc[lv].Tm_max}", int(POOL.loc[lv].severe), f"{ot(lv, 2, 'loci'):,} / {ot(lv, 3, 'products')}"] for lv in ("24", "16", "12v", "12")],
         ["规模", "展开序列数", "综合得分 J（95% 区间）", "属准确率（完整数据库）", "扩增子覆盖（主要门平均）", "支原体≥3 个扩增子", "6 个扩增子全部扩出 GG / SILVA", "Tm °C", "严重二聚体", "人基因组 ≤2 错配位点 / ≤3 错配产物"])
t2 = tbl([["24 → 16", "−8 条，−54 个展开序列", "0.0091（0.0060–0.0124）", "0.2（0.0–0.4）", "0.9（0.7–1.0）", "1.7（0.7–2.6）"], ["16 → 12", "−4 条，−24 个展开序列", "0.0417（0.0340–0.0504）", "0.1（−0.3–0.6）", "4.9（4.5–5.3）", "7.5（5.2–10.0）"], ["12 → 12 + V1_f", "A1-F 换成 V1_f，展开序列 −6", "0.0005（0.0001–0.0008）", "0.03（0.00–0.08）", "0.10（−0.01–0.20）", "0.0（0.0–0.0）"], ["16 → 12 + V1_f", "−4 条，−30 个展开序列", "0.0421（0.0345–0.0508）", "0.1（−0.3–0.6）", "5.0（4.7–5.4）", "7.5（5.2–10.0）"]], ["从 → 到", "省下的引物", "J 下降（95% 区间）", "属准确率下降（百分点）", "扩增子覆盖下降（百分点）", "支原体≥3 个扩增子下降（百分点）"])
names = ["V1·V2", "V3", "V4", "V5 (A3-F×V5-R)", "V6·V7", "V8·V9"]
rows = []
for nmx in names:
    r = [nmx]
    for lv in ('24', '16', '12v', '12'):
        g = AMP[(AMP.level == lv) & (AMP.amplicon == nmx)]
        gg = g[g.db == "GG"].iloc[0]; sv = g[g.db == "SILVA"].iloc[0]
        tg = TV[(TV.level == lv) & (TV.amplicon == nmx) & (TV.db == "GG")].iloc[0]; ts_ = TV[(TV.level == lv) & (TV.amplicon == nmx) & (TV.db == "SILVA")].iloc[0]
        r.append(f"{pc0(gg['mean'])} / {pc0(sv['mean'])}（最差门 {pc0(gg.worst)} / {pc0(sv.worst)}）；支原体 {pc0(tg.cov_valid)} / {pc0(ts_.cov_valid)}")
    rows.append(r)
t3 = tbl(rows, ["扩增子（覆盖：GG / SILVA 主要门平均；支原体 GG / SILVA）", "24 条", "16 条", "12 条 + V1_f", "12 条"])
def f(v): return "—" if pd.isna(v) else f"{v:.2f}"
t4 = tbl([["@@N:" + r.order.replace("-", "~") + "@@", "补充" if r.is_sup else "", r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}", f(r.uL_12), f(r.uL_12v), f(r.uL_16), f(r.uL_24)] for _, r in W.iterrows()],
         ["订购名称", "类型", "序列 5′→3′（IUPAC）", "nt", "E. coli 位置", "Tm °C（展开范围）", "展开数", "人基因组位点 / 展开序列（≤2 错配）", "12 条池 µL", "12 条 + V1_f 池 µL", "16 条池 µL", "24 条池 µL"], seqcol=2)
MAPHTML = "<p>订购名称按 6 个扩增子编号：16S-A1（V1·V2）、A2（V3）、A3（V4）、A4（V5）、A5（V6·V7）、A6（V8·V9）；V5 扩增子的正向引物和 V4 共用（16S-A3-F），所以 A4 只有反向引物（16S-A4-R）。末尾 s = 支原体补充，.SNAP = 借鉴 Swift SNAP 的引物。分析里用过的旧名在下表：</p>" + tbl([[re.sub(r"^SMURF5-", "", r.order), new_order(r.order), r.slot] for _, r in W.iterrows()], ["分析里用过的旧名（页面和早期文件里）", "订购名称（Excel、本页）", "旧位点名"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:1120px;max-width:100%;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>12、12+V1_f、16、24 条引物的订购清单（含选 16 条的理由）</h1>
  <p class="lede">四个引物池（12 条、12 条 + V1_f、16 条、24 条）都是同一个 24 条候选池的子集，一共订 24 条寡核苷酸，按下表给的体积配四个池，就可以做对比实验。我推荐 16 条：比 24 条少 8 条、少 54 个展开序列（−35%），综合得分只下降 0.009；再少到 12 条，得分下降 0.042，主要是扩增子覆盖和支原体明显变差。这是序列层面的计算，没有实验验证。</p>
</header>

<h2>一、四个池的比较</h2>
{t1}
<ul>
<li><b>24 条</b>：22 条 V5 版推荐引物（含 4 条支原体补充）加上 2 条借鉴的 SNAP 引物（V1_f 和 907R 型 V5_r）。其中 V1_f 和 A1-F、V5_r 和 V5-R 两条在同一个位点，是互相冗余的，所以 24 条是“能力上限”，不是推荐版本。</li>
<li><b>16 条</b>：每个位点的冗余变体都删掉，留下对覆盖贡献大的；A1-F 换成 SNAP 的 V1_f，V5 反向只留 1 条（V5-R.1），保留 3 条支原体补充中的 2 条（A4-F.3、A5-F.3）。</li>
<li><b>12 条</b>：几乎每个位点只剩 1–2 条（A3-F 保留 2 条，A5-F 只留 1 条，不含任何支原体补充）。</li>
<li><b>12 条 + V1_f</b>：12 条池里把 A1-F 换成 SNAP 的 V1_f，其余完全相同，用来把“引物条数”和“V1_f 带来的人基因组脱靶下降”分开。</li>
<li>综合得分 J 是三项等权平均：属准确率（数据库按完整序列、样本只有扩出的扩增子）、6 个扩增子覆盖（主要门平均）、支原体至少扩出 3 个扩增子的比例；括号里是对样本序列做 400 次自助重抽样的 95% 区间。</li>
</ul>

<h2>二、为什么选 16 条</h2>
{t2}
<ul>
<li><b>24 → 16 条</b>：J 只下降 0.009，其中属准确率几乎没变（0.2 个百分点），主要是扩增子覆盖降 0.9 个百分点、支原体≥3 个扩增子降 1.7 个百分点；区间不含 0，是真实的但很小的差别。</li>
<li><b>16 → 12 条</b>：J 下降 0.042，是前一步的 4.6 倍：扩增子覆盖降 4.9 个百分点，支原体≥3 个扩增子降 7.5 个百分点，6 个扩增子全部扩出的序列在 Greengenes 里从 67% 降到 53%；属准确率这一项没有明显差别。这是在 16 条这一点之后损失明显加快，所以 16 条是拐点。</li>
<li><b>16 条的额外好处来自 V1_f，不是条数</b>：24 条池同时含 A1-F，12 条池用 A1-F，人基因组 ≤2 个错配的位点是 3,125 个和 2,849 个；16 条池用 V1_f，降到 702 个；12 条 + V1_f 池只要 570 个，≤3 个错配的潜在产物 2 个（12 条池是 50 个）。换成 V1_f 对覆盖几乎没有影响（12 条池换了之后 J 只下降 0.0005）。线粒体命中四个池都只有 3–4 个错配（没有成对的反向命中），没有潜在产物。</li>
<li>属准确率这一项对引物数量很不敏感，所以单看准确率三个规模没有区别；要看覆盖和支原体才拉得开。</li>
</ul>

<h2>三、各扩增子的覆盖（四个池）</h2>
{t3}
<ul>
<li>12 条池里支原体最弱的是 V6·V7（GG 50%、SILVA 68%）和 V8·V9（SILVA 61%）：这两个位点只剩 1 条正向引物、没有补充；16 条池靠 A4-F.3 和 A5-F.3 补充把它们拉回到 81% / 84% 和 94% / 88%；V1·V2 的支原体覆盖在 16 条池（GG 90%、SILVA 77%）比 24 条池（89% / 84%）的 SILVA 略低（少了 A1-R.3 和 A1-R.4 补充）。</li>
<li>每个扩增子的覆盖：GG / SILVA 是主要门平均，括号里是最差的门；最差门的数字对序列少的门不稳。支原体的百分数只算该扩增子两端位点都有数据的序列；综合得分里的“支原体至少扩出 3 个扩增子”则把缺数据的序列按“未扩出”算，所以数值偏保守。</li>
</ul>

<h2>四、订购清单（共 24 条，一次订齐）</h2>
{t4}
<ul>
<li>序列用 IUPAC 简并码，按混合碱基合成；反向引物的序列已经是 5′→3′ 的引物序列，位置是它在 E. coli 16S（顶链）上占据的范围；只含基因特异序列，没有接头和 UDP，带接头后的 Tm 和二聚体没有评估。</li>
<li><b>配池方法</b>：每个池都是 100 µL 的 10× 混合液；表里“µL”列是该池里每条引物要加的 100 µM 母液体积，空白（—）表示这个池不含这条。四个池每个位点的总量相同（2.5 µL，终浓度 250 nM），所以引物总浓度一致，池之间只差引物的多样性；11 个位点共 27.5 µL，加 72.5 µL 水。每 25 µL PCR 加 2.5 µL。</li>
<li>同位点内的比例：核心引物按各自单独能匹配到的序列比例分，补充引物各占 10%。</li>
<li>24 条池里同位点的冗余引物（A1-F 和 V1_f；V5-R.1、V5-R.2 和 907R 型）按覆盖比例平分。</li>
<li>池内所有寡核苷酸（含展开）两两检查：四个池都没有严重二聚体或发夹。</li>
</ul>

<h2>五、名称对照</h2>
@@MAP@@

<h2>六、对比实验的建议</h2>
<ul>
<li>用同一批样本（已知组成的标准菌群，最好含支原体和浮霉菌；加上你关心的真实样本）、同一份 DNA、同一套循环条件，分别用 12、16、24 条池扩增，每个池至少 3 个重复；测序读段数抽平后比较。</li>
<li>比较的指标：每个扩增子的读段占比（看 V5、V3 的读出比例是否够）、标准菌群各物种的检出情况和相对丰度偏差、支原体和浮霉菌是否掉线、通过引物去除的读段比例、人 DNA 污染样本里的无效读段比例。</li>
<li>预期（按计算）：16 条和 24 条差别很小，12 条在 V3 和 V8·V9 的覆盖、支原体检出上明显变差；如果实验里 16 条和 24 条也有可见差别，说明补充引物或冗余变体比我估计的更重要。</li>
<li>注意：16 条池和“12 条 + V1_f”池用 V1_f，12、24 条池用 A1-F。比较引物条数的效果时，用“12 条 + V1_f”对“16 条”，或“12 条”对“24 条”，这样 V1·V2 的引物选择是一致的；看 V1_f 本身的效果时用“12 条”对“12 条 + V1_f”。</li>
</ul>

<h2>七、限制</h2>
<ul>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算，不是实际 PCR 效率；数据库是 Greengenes 13_8 和 SILVA 128；没有实验验证。</li>
<li>区间只反映样本序列的抽样波动，不包含引物选择本身的不确定性；综合得分 J 是我定的三项等权平均，你按自己的侧重（比如只看支原体）可能选出别的点。</li>
<li>删掉一条引物后只重新分配了位点内的比例，没有评估扩增效率和扩增子占比的变化，需要实测调整。</li>
</ul>
<footer>脚本：<code>python_5R/explore/lvl1.py</code>、<code>lvl2.py</code>、<code>par*.py</code>；数据：<code>docs/primer_design/order_list_levels.csv</code>、<code>lvl_*.csv</code>。</footer>
</main></div>
"""
body = rename_text(body)
body = re.sub(r"@@N:(.*?)@@", lambda m: new_order(m.group(1).replace("~", "-")), body)
if "@@MAP@@" in body: body = body.replace("@@MAP@@", MAPHTML)
open(DOCS + "/primer_order_list_levels.html", "w", encoding="utf-8").write(head + body)
print("ok")
