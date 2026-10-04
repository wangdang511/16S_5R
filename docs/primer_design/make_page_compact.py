"""生成 docs/primer_compact_designs.html"""
import re, html, json
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>更少扩增子的设计</title>", old[:old.index("</style>") + 8], 1)
SV = pd.read_csv(f"{HERE}/subset_val.csv"); SJ = pd.read_csv(f"{HERE}/subset_joined.csv").set_index("combo")
BB = pd.read_csv(f"{HERE}/block_backward.csv"); GL = pd.read_csv(f"{HERE}/gs_loose290.csv"); GS = pd.read_csv(f"{HERE}/gs_strict290.csv"); S5 = pd.read_csv(f"{HERE}/seeded_loose290_K5.csv")
NUMRE = re.compile(r"^[\d.%/ –→-]+$")
def cell(i, v): return '<td class="' + ("n" if NUMRE.match(str(v)) else "") + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads):
    return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(i, v) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: f"{x * 100:.1f}%"
p0 = lambda x: f"{x * 100:.0f}%"

top = SV.sort_values(["n", "val_acc"], ascending=[True, False]).groupby("n").head(3)
rows = []
for _, r in top.iterrows():
    j = SJ.loc[r.combo]
    rows.append([int(r.n), r.combo, r.spans, int(r.bases), pc(r.val_acc), pc(j.acc_all), int(r.sites), int(r.oligos_multi), p0(j.p_all_SILVA), p0(j.p_all_GG), p0(j.worst_phylum_all_SILVA)])
t_sub = tbl(rows, ["扩增子数", "覆盖的区域", "扩增子（E. coli 位置）", "测序碱基数", "属准确率（验证）", "属准确率（全部 5,899 条）", "引物位点数", "每位点 2–3 条时的寡核苷酸数", "全部扩增子都扩出：SILVA", "同：GG", "同：最差门（SILVA）"])
bb = BB.set_index("n_blocks")
t_bb = tbl([[n, int(bb.loc[n, "bases"]), pc(bb.loc[n, "acc"])] for n in [75, 60, 50, 40, 30, 25, 20, 15, 12, 10]], ["保留的 20 nt 区块数", "保留的碱基数", "设计集准确率（逐块删除后）"])
sr = []
for nm, G in [("宽松（单引物主要门平均 ≥70%，最差门 ≥40%）", GL), ("严格（≥85%，≥60%）", GS)]:
    for _, r in G.iterrows():
        sr.append([nm, int(r.K), pc(r.val_acc), p0(r.minq), r.amplicons])
for _, r in S5.iterrows():
    sr.append(["宽松，K=5（分段初值 + 坐标下降）", 5, pc(r.val_acc), p0(r.minq), r.amplicons])
t_search = tbl(sr, ["位点池", "扩增子数", "属准确率（验证）", "池内最差引物的最差门覆盖", "搜索得到的扩增子"])

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>更少的扩增子能不能更好：以属准确率为目标的搜索</h1>
  <p class="lede">在单管、引物位于保守区、扩增子 ≤290 bp、同管不重叠的约束下，用属水平最近邻准确率直接做目标，找比上一轮 5 扩增子方案更好的策略。结论是准确率已经接近天花板，真正能改进的是在差不多的准确率下减少扩增子和引物。</p>
  <p>这是计算机模拟，没有实验验证。数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本。属准确率只测到属，不到种。</p>
</header>

<h2>一、结论</h2>
<ul>
<li><b>天花板很低。</b>用整条 16S（20–1510）做最近邻，验证集属准确率是 <b>95.7%</b>；5R 现有 94.0%；上一轮 5 扩增子方案 95.1%（距天花板 0.6 个百分点）。再怎么设计，5 个扩增子最多再多 0.5 个百分点。</li>
<li><b>补上 V5 的收益很小。</b>把上一轮方案里 V4 和 V6 之间的空隙（785–927，覆盖 V5）加进去，准确率从 95.1% 升到 95.4%，增加 143 nt，但要多一个扩增子，而且同管放不下（空隙只有 108 nt）。不值得。</li>
<li><b>真正的改进是减少扩增子。</b>在上一轮方案的扩增子里挑子集：<b>4 个扩增子（V3+V4+V6V7+V8V9）94.0%，3 个扩增子（V3+V4+V8V9）93.5%</b>，都与 5R 现有（94.0%）相当，但分别少 1 个和 2 个扩增子，位点数 8 和 6（5R 是 10）。</li>
<li><b>如果更看重准确率：</b>保持 5 个扩增子（95.1%）；<b>如果更看重稳健性和引物数量：</b>4 个扩增子，去掉 V1·V2（V1·V2 的引物最弱，全部扩增子都扩出的比例比其他区域低）。</li>
<li>我用自由搜索（不限定在上一轮的位点）没有找到比上一轮 5 扩增子更好的 5 扩增子方案（最好 94.8%），原因是搜索只能从保守性合格的位点池里选，池比较小（见第四节）。</li>
</ul>

<h2>二、评估方法</h2>
<ul>
<li>5,899 条带属标签的近全长（20–1510 位都有数据）细菌序列，590 个属，来自 Greengenes 97% 代表序列，用修正过的参照做 E. coli 坐标投影。随机分成设计集（1,100 条）和验证集（4,799 条，其中有同属参照的 4,741 条作为查询）。</li>
<li>属准确率：留一最近邻，用扩增子内部被测到的碱基（去掉引物）的错配数找最近的另一条序列，看它是否同属，并列按比例计分。读长 2×150，所以 ≤300 bp 的扩增子内部全部被测到。</li>
<li>搜索用设计集，最后的数字用验证集；第五节的子集表另外给出用全部 5,899 条（自参照）算的数字，样本更多，标准误约 0.3%。</li>
<li>差别在 0.3 个百分点以内不要当作有差别。</li>
</ul>

<h2>三、哪些区段有用：逐块删除</h2>
<p>从整条 16S 出发，每次删掉“删了之后准确率下降最少（或上升）”的 20 nt 区块（在设计集上选，不管引物位点）：</p>
{t_bb}
<ul>
<li>删到只剩约 500 nt，设计集准确率还高于全部碱基（92.8%）：有些位点主要增加属内的噪声，所以<b>碱基越多不一定越好</b>。这也是为什么 4–5 个扩增子就接近天花板。</li>
<li>保留 600–800 nt 时，验证集准确率是 95.0%–95.5%，保留的区段集中在 V1–V2（80–239）、V4（580–759）、V6（980–1059）、V7（1120–1159）、V8（1280–1359）和 V9（1420–1479）。这和上一轮方案选中的区域一致。</li>
<li>注意：这里是在设计集上贪心选块，验证集上的数字（95.6% / 95.5% / 95.0% / 94.1%，对应 1200 / 800 / 600 / 400 nt）才有意义。</li>
</ul>

<h2>四、自由搜索：不限定位点，直接找最优扩增子</h2>
<p>候选位点是前面扫描出的、两个库都能用的位点（去掉了有缺口的）；按单引物覆盖率分成宽松和严格两个位点池。每个槽位遍历全部合法扩增子，取使设计集准确率最高的，反复迭代。</p>
{t_search}
<ul>
<li>K=4 时 94.5%，K=5 时 94.8%，比上一轮方案的 95.1% 低；严格池里放不下第 4 个扩增子。结论：位点池太小，限制了搜索；上一轮方案用的是后来针对多条引物重新设计过的位点，不受这个池的限制。</li>
<li>宽松池在两端（V1·V2 和 V8·V9）选的引物最差门覆盖只有约 47%–49%，说明想要更高准确率时，会选用质量较差的引物位点。</li>
<li>早期的模拟退火在 K≥5 时没有找到合法的初始解，K=4 的结果还不如 K=3，说明没收敛；这里的结果来自贪心加坐标下降，仍然是局部最优，不是全局最优的保证。</li>
</ul>

<h2>五、在上一轮 5 个扩增子里挑子集</h2>
<p>5 个扩增子是 V1·V2（11–260）、V3（337–530）、V4（558–802）、V6·V7（910–1195）、V8·V9（1222–1509），引物已经设计好、做过验证。下表是每个扩增子数下准确率最高的 3 个组合。“全部扩增子都扩出”用上一轮每位点 2–3 条引物的设计，在没参与设计的 SILVA（10,443 条）和 Greengenes（3,151 条）上算。</p>
{t_sub}
<ul>
<li><b>V4 是最重要的单个扩增子</b>：只用 V4 就有 90%，加上 V8·V9 到 91.9%，再加 V3 到 93.5%。</li>
<li><b>V1·V2 是最不稳的扩增子</b>：它的引物对最差门的覆盖，单引物是 44% / 40%（SILVA / Greengenes），每位点多条引物后是 68% / 67%；单独一个 V1·V2 扩增子，SILVA 上只有 70% 的序列能扩出（其他扩增子 82%–93%）。含它的组合“全部扩增子都扩出”的比例明显更低（5 个扩增子只有 58% / 68%）。去掉它，准确率下降约 1.1 个百分点（95.1% → 94.0%），稳健性好很多。</li>
<li>4 个扩增子里，去掉 V3 的组合（V1V2+V4+V6V7+V8V9）准确率最高（94.2%），去掉 V1·V2 的组合（V3+V4+V6V7+V8V9）稳健性最好（94.0%，全部扩出 70% / 82%）。</li>
<li>每个位点 2–3 条引物时，3 个扩增子是 17–18 种寡核苷酸，4 个是 23–24 种，5 个是 29 种。单引物版本每个位点 1 条，但 V4 的正向引物对 Bacteroidetes 只有 2%–8% 覆盖，所以 V4 的正向位点必须放 2 条引物。</li>
</ul>

<h2>六、建议</h2>
<div class="tbl"><table><thead><tr><th>目标</th><th>方案</th><th>属准确率（验证）</th><th>扩增子 / 引物位点</th></tr></thead><tbody>
<tr><td>最高准确率</td><td>5 个扩增子（上一轮方案）</td><td class="n">95.1%</td><td class="n">5 / 10</td></tr>
<tr><td>准确率与稳健性折中</td><td>4 个：V3 + V4 + V6V7 + V8V9</td><td class="n">94.0%</td><td class="n">4 / 8</td></tr>
<tr><td>最少引物，与 5R 现有相当</td><td>3 个：V3 + V4 + V8V9</td><td class="n">93.5%</td><td class="n">3 / 6</td></tr>
<tr><td>5R 现有（对照）</td><td>5 个，R1–R5</td><td class="n">94.0%</td><td class="n">5 / 10</td></tr>
</tbody></table></div>
<p>三个方案的具体引物序列在 <code>primer_v3v4_design.html</code>（A2、A3、A4、A5 对应 V3、V4、V6·V7、V8·V9；A1 对应 V1·V2）。</p>

<h2>七、局限</h2>
<ul>
<li>属准确率是一个代理指标：只测属，没有种；Greengenes 的属标签偏向培养过的属；最近邻按错配数找，和 SMURF 的似然模型不同。</li>
<li>第五节的“全部扩增子都扩出”只是引物覆盖，没有考虑不同扩增子的扩增效率差异，也没有把“没扩出的扩增子”对准确率的影响算进去：真实样本里序列缺失某个扩增子时，准确率会比表中低。</li>
<li>不同扩增子数的差别在 0.3–1 个百分点，比标准误（约 0.3%）大不了多少；要确认需要更大的带标签数据和实验验证。</li>
<li>数据库不是最新版本，GTDB 本环境下载不了。k-mer 数据库需要用新引物重建，SMURF 里写死的 <code>nR = 5</code> 要改。</li>
</ul>
<footer>脚本：<code>python_5R/tiling_search.py</code>、<code>python_5R/explore/</code>；数据：<code>docs/primer_design/subset_*.csv</code>、<code>block_backward.csv</code>、<code>gs_*.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_compact_designs.html", "w", encoding="utf-8").write(head + body)
print("ok")
