"""生成 docs/primer_pareto.html"""
import re, html, pickle
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>引物精简帕累托前沿</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
A = pd.read_csv(f"{HERE}/par_all_states.csv"); F1 = pd.read_csv(f"{HERE}/par_frontier_oligos.csv"); F2 = pd.read_csv(f"{HERE}/par_frontier_exp.csv"); LV = pd.read_csv(f"{HERE}/par_levels.csv"); OR = pd.read_csv(f"{HERE}/par_base_order.csv"); LOO = pd.read_csv(f"{HERE}/par_loo.csv")
pc = lambda x: f"{x * 100:.1f}%"
# ---- 图
W, H = 940, 330
def sx(v, lo, hi, x0, x1): return x0 + (v - lo) / (hi - lo) * (x1 - x0)
def sy(v, lo=0.78, hi=0.96): return 270 - (v - lo) / (hi - lo) * 230
svg = f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="精简帕累托前沿" style="max-width:{W}px">'
# 左图：J 与分量 vs 寡核苷酸数
x0, x1 = 50, 450
for v in (0.80, 0.85, 0.90, 0.95): svg += f'<line x1="{x0}" y1="{sy(v):.1f}" x2="{x1}" y2="{sy(v):.1f}" stroke="currentColor" opacity="0.12"/><text x="{x0 - 4}" y="{sy(v) + 3:.1f}" font-size="9" text-anchor="end" fill="currentColor">{v:.2f}</text>'
for k in range(11, 25, 1): svg += f'<text x="{sx(k, 11, 24, x0, x1):.1f}" y="288" font-size="9" text-anchor="middle" fill="currentColor">{k}</text>'
svg += f'<text x="{(x0 + x1) / 2}" y="304" font-size="10" text-anchor="middle" fill="currentColor">寡核苷酸数</text>'
fr = F1[F1.oligos >= 11].sort_values("oligos")
def line(col, color, w=2, dash=""):
    pts = " ".join(f"{sx(r.oligos, 11, 24, x0, x1):.1f},{sy(r[col]):.1f}" for _, r in fr.iterrows())
    return f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{w}" {dash}/>'
svg += line("J", "#2a7f62", 3) + line("cov", "#4a7fd0", 1.6) + line("ten", "#c0562d", 1.6) + line("accf", "#8a8f98", 1.6, 'stroke-dasharray="4 3"')
for col, color, lab, yy in (("J", "#2a7f62", "综合得分 J", 20), ("cov", "#4a7fd0", "扩增子覆盖（主要门平均）", 32), ("ten", "#c0562d", "支原体≥3 个扩增子", 44), ("accf", "#8a8f98", "属准确率（完整数据库）", 56)):
    svg += f'<line x1="{x0 + 8}" y1="{yy - 3}" x2="{x0 + 28}" y2="{yy - 3}" stroke="{color}" stroke-width="2"/><text x="{x0 + 32}" y="{yy}" font-size="9" fill="currentColor">{lab}</text>'
# 右图：J vs 展开数（全部状态 + 前沿）
rx0, rx1 = 520, 920
for v in (0.80, 0.85, 0.90, 0.95): svg += f'<line x1="{rx0}" y1="{sy(v):.1f}" x2="{rx1}" y2="{sy(v):.1f}" stroke="currentColor" opacity="0.12"/>'
for e in (60, 80, 100, 120, 140, 160): svg += f'<text x="{sx(e, 55, 160, rx0, rx1):.1f}" y="288" font-size="9" text-anchor="middle" fill="currentColor">{e}</text>'
svg += f'<text x="{(rx0 + rx1) / 2}" y="304" font-size="10" text-anchor="middle" fill="currentColor">展开序列总数</text>'
for _, r in A[(A.oligos >= 10) & (A.J > 0.78)].iterrows(): svg += f'<circle cx="{sx(r.expansions, 55, 160, rx0, rx1):.1f}" cy="{sy(r.J):.1f}" r="2" fill="#8a8f98" opacity="0.5"/>'
f2 = F2[(F2.oligos >= 10) & (F2.J > 0.78)].sort_values("expansions")
svg += '<polyline points="' + " ".join(f"{sx(r.expansions, 55, 160, rx0, rx1):.1f},{sy(r.J):.1f}" for _, r in f2.iterrows()) + '" fill="none" stroke="#2a7f62" stroke-width="2"/>'
for _, r in f2.iterrows(): svg += f'<circle cx="{sx(r.expansions, 55, 160, rx0, rx1):.1f}" cy="{sy(r.J):.1f}" r="3" fill="#2a7f62"/><text x="{sx(r.expansions, 55, 160, rx0, rx1):.1f}" y="{sy(r.J) - 6:.1f}" font-size="8" text-anchor="middle" fill="currentColor">{int(r.oligos)}</text>'
svg += f'<text x="{rx0}" y="22" font-size="10" fill="currentColor">右：所有评估过的子集（灰点）与前沿（绿线，标注寡核苷酸数）</text><text x="{x0}" y="14" font-size="10" fill="currentColor">左：前沿上每个规模的得分</text></svg>'
t_lv = tbl([[int(r.oligos), int(r.expansions), f"{r.J:.4f}", pc(r.accf), f"{r.GG_amp_mean * 100:.1f}% / {r.SILVA_amp_mean * 100:.1f}%", f"{r.Ten_GG_ge3 * 100:.0f}% / {r.Ten_SILVA_ge3 * 100:.0f}%", f"{r.GG_all * 100:.0f}% / {r.SILVA_all * 100:.0f}%", f"{r.Tm_min}–{r.Tm_max}", int(r.severe)] for _, r in LV.iterrows()],
            ["寡核苷酸数", "展开数", "综合得分 J", "属准确率（完整数据库）", "扩增子覆盖（主要门平均）GG / SILVA", "支原体≥3 个扩增子 GG / SILVA", "6 个扩增子全部扩出 GG / SILVA", "Tm °C", "严重二聚体"])
meta = pickle.load(open(__file__.rsplit("/", 1)[0] + "/../../python_5R/explore/_nothing.pkl", "rb")) if False else None
t_or = tbl([[int(r.step), r.removed, r.slot, r.seq, int(r.nexp), int(r.left), f"{r.J:.4f}", int(r.expansions)] for _, r in OR.iterrows() if r.step <= 14], ["删除顺序", "删掉的寡核苷酸", "位点", "序列 5′→3′", "展开数", "剩余寡核苷酸", "综合得分 J", "剩余展开数"], seqcol=3)
t_loo = tbl([[r["name"], r.slot, int(r.nexp), f"{r.dJ:.4f}", f"{r.dten * 100:.1f}", f"{r.dcov * 100:.1f}", f"{r.dall_GG * 100:.1f} / {r.dall_SILVA * 100:.1f}"] for _, r in LOO.iterrows() if r.dJ < 0.03], ["寡核苷酸", "位点", "展开数", "单独删掉的 J 下降", "支原体≥3 个扩增子下降（百分点）", "扩增子覆盖下降（百分点）", "全部扩出下降（百分点）GG / SILVA"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:1040px;max-width:100%;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>引物精简：效能和引物数量的帕累托前沿</h1>
  <p class="lede">在加了 V5 的 22 条推荐引物（加上借鉴来的 2 条 SNAP 引物作为候选，共 24 条）里，先删掉对结果影响最小的。结论：从 24 条减到 19 条（展开序列 154 → 120，−22%）几乎没有损失，减到 16–17 条（展开 100–104，−33%）只损失一点，再往下损失明显加快，低于 12 条会丢掉整个扩增子。这是序列层面的计算，没有实验验证。</p>
</header>

<h2>一、怎么算“效能”</h2>
<ul>
<li>候选共 24 条寡核苷酸，分布在 11 个引物位点（A1 到 A5 的正反向，加 V5 的反向）；其中有 2 条是借鉴来的 Swift SNAP 引物（V1_f 可替代 A1-F，907R 型 V5_r 可替代我们的 V5-R 两条）。</li>
<li>每删一条引物，重新算 6 个扩增子里每个扩增子能扩出多少序列（只要这个位点还有一条引物匹配，就算扩出）。综合得分 J 是三项的平均：①属准确率（数据库按完整序列算，样本只有扩出的扩增子，留出集 2,196 条）；②6 个扩增子覆盖（Greengenes 和 SILVA 主要门的平均）；③支原体至少扩出 3 个扩增子的比例。</li>
<li>属准确率这一项对引物数量很不敏感（0.93–0.94 之间），差 0.5 个百分点以内是噪声（约 10 条序列）；真正拉开差距的是覆盖率和支原体，所以我把它们并入综合得分。</li>
<li>用向后逐条删除（每步删掉使 J 下降最小的那条）和向前逐条添加两种贪心方法，再对每个规模做同位点互换的局部搜索，没有找到更好的；前沿是所有评估过的子集里不被更好的子集压过的那些。它是近似最优，不是穷举。</li>
</ul>
{svg}

<h2>二、几个有代表性的规模</h2>
{t_lv}
<ul>
<li><b>24 → 19 条</b>：J 0.9319 → 0.9296，展开 154 → 120。删掉的是冗余变体（V5-R.2、借鉴的 V5_r、A4-F.2、A1-R.4、A4-R.2）。</li>
<li><b>16–17 条</b>：J 0.9228–0.9272，展开 100–104，是性价比的拐点：再多删一条，J 下降从 0.004 变成 0.007–0.012。</li>
<li><b>14 条</b>：J 0.9036，支原体覆盖开始明显下降（至少扩出 3 个扩增子：GG 92%，SILVA 89%）。</li>
<li><b>11–12 条</b>：每个位点只剩 1 条（最少 11 个位点），J 0.87–0.88，6 个扩增子全部扩出 GG 49–53%、SILVA 38–41%。</li>
<li>表里这些代表规模都做了二聚体检查，没有严重二聚体。</li>
</ul>

<h2>三、先删哪些（对结果影响最小的顺序）</h2>
<p>在 22 条基线上（不含借鉴来的 SNAP 引物），逐条删掉使 J 下降最小的：</p>
{t_or}
<p>单独删掉某一条引物的影响（只列下降小于 0.03 的）：</p>
{t_loo}
<ul>
<li><b>最不值得保留的</b>：V5-R.2、A4-F.2（几乎没有贡献，J 下降 0.0002）；A1-R.4（补充）、A4-R.2（补充）、V5-R.1、A1-R.3（J 下降 0.001）。</li>
<li><b>支原体补充引物</b>的价值在支原体指标上：A2-F.2（补）支原体下降 1.1 个百分点，A4-F.3（补）1.6，A5-F.3（补）2.7，其他覆盖指标几乎没影响；如果你的样本里支原体不重要，补充引物可以先删。</li>
<li><b>不能删</b>：只剩一条的位点（A1-F 或 SNAP V1_f、A2-R、A3-R、A5-R、V5-R 里至少 1 条）：删掉就丢整个扩增子（J 下降 0.06–0.09）；A3-F.1 删掉最伤（J 下降 0.22，说明 A3-F.2 单独不够）。</li>
<li>单独删一条的影响会随着其他条已删而变化，所以应该看第二节的前沿点，不要把单条删除的影响简单相加。</li>
</ul>

<h2>四、有没有值得借鉴的现有方案</h2>
{tbl([["Swift SNAP V1_f（9–27，19 nt，2 个展开序列）", "替代 A1-F（17 nt，8 个展开序列）", "覆盖相当（GG / SILVA 97% / 94% 对 97% / 95%），展开序列少 6 个；对人基因组脱靶 25 个/展开序列，比 A1-F 的 291 个少得多。换了以后 J 几乎不变（0.9306 对 0.9285）。建议采用。"], ["Swift SNAP V5_r（907R，20 nt，1 条）", "替代我们的 V5-R 两条（12 个展开序列）", "覆盖 97% / 97%（我们 2 条 98% / 98%），少 1 条、少 10 个展开序列，J 只下降 0.001。想精简时可以采用。"], ["515F（V4 正向，1 条）", "替代 A3-F（2 条，16 个展开序列）", "覆盖更好（98% / 98%，最差门 96% / 92%），但位点和我们的 V3 反向引物 A2-R 互补，不能同用；要换就得把 V3 改成只用正向读段（方案 C，V3 只读到 69%）。"], ["SNAP V4_r（806R 型）、V6_f 的 6 条变体", "对应 A3-R、A4-F", "覆盖明显差（V4_r 73% / 78%，最差门 5%；V6_f 6 条合起来 85% / 90%），不借鉴，保留我们做过简并的版本。"], ["5R 原引物", "—", "整体覆盖最差（平均 74–77%，7 个位点低于 80%），没有可借鉴的位点。"]], ["来源引物", "对应我们的", "评估和建议"])}
<p>同时换上两条借鉴引物（V1_f 和 907R 型 V5_r）：21 条寡核苷酸、展开 134 个，J 0.9293，比 22 条基线（150 个展开，J 0.9285）少 16 个展开序列，J 不降。</p>

<h2>五、限制</h2>
<ul>
<li>J 是我为权衡定的综合得分（三项等权），不同侧重（比如只看支原体，或只看属准确率）会选出不同的前沿；表里同时给出了各分量，你可以按自己的侧重选。</li>
<li>属准确率用留出集 2,196 条序列的代理指标，对 1 个百分点以内的差别不稳；覆盖率是序列匹配规则的结果，不是实际 PCR 效率；没有实验验证。</li>
<li>删除只按匹配覆盖来评估，没有评估对各位点混合比例、扩增效率的影响；删掉一条引物后，同位点其他引物的比例需要重新分配。</li>
<li>前沿是贪心 + 局部搜索得到的近似最优，不是穷举全部 2<sup>24</sup> 种组合。</li>
</ul>
<footer>脚本：<code>python_5R/explore/par1.py</code>—<code>par6.py</code>；数据：<code>docs/primer_design/par_*.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_pareto.html", "w", encoding="utf-8").write(head + body)
print("ok")
