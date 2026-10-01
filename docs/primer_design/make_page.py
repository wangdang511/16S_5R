"""由 designs_final.json / profile_*.csv / frontier.csv 生成 docs/primer_redesign.html"""
import json, re, html
import numpy as np, pandas as pd

HERE = __file__.rsplit("/", 1)[0]
OUT = HERE.rsplit("/", 1)[0] + "/primer_redesign.html"
old = open(HERE.rsplit("/", 1)[0] + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = old[:old.index("</style>") + 8]                       # 复用样式
head = re.sub(r"<title>.*?</title>", "<title>5R 引物重设计</title>", head, 1)
D = json.load(open(f"{HERE}/designs_final.json"))
prof = pd.read_csv(f"{HERE}/profile_bact.csv", index_col=0)
fr = pd.read_csv(f"{HERE}/frontier.csv")
fr.columns = ["key min", "pools", "n V", "complete V", "amplicons"]
HO = pd.read_csv(f"{HERE}/holdout.csv")
VR = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879),
      "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}
FIVE = [("R1", 103, 332), ("R2", 338, 536), ("R3", 685, 927), ("R4", 944, 1104), ("R5", 1175, 1391)]


def track_svg():
    W, X0, X1 = 900, 105, 885
    sx = lambda p: X0 + (p - 1) / 1541 * (X1 - X0)
    s = prof.top_freq.reindex(range(1, 1543)).fillna(0).rolling(9, center=True, min_periods=1).mean()
    pts = " ".join(f"{sx(p):.1f},{150 - 120 * (v - 0.4) / 0.6:.1f}" for p, v in s.items())
    rows = [("5R 现有", [(n, a, b, 0) for n, a, b in FIVE])]
    plan = D["A_two_pool"]["amplicons"]
    rows.append(("方案 A 管 1", [(a["amplicon"], a["start"], a["end"], 0) for a in plan if a["pool"] == "A"]))
    rows.append(("方案 A 管 2", [(a["amplicon"], a["start"], a["end"], 0) for a in plan if a["pool"] == "B"]))
    o = [f'<svg viewBox="0 0 {W} 300" width="100%" style="min-width:680px" role="img" aria-label="16S 保守性与引物平铺">']
    for v, (a, b) in VR.items():
        o.append(f'<rect x="{sx(a):.1f}" y="20" width="{sx(b) - sx(a):.1f}" height="140" fill="var(--accent-soft)"/>'
                 f'<text x="{(sx(a) + sx(b)) / 2:.1f}" y="14" font-size="11" text-anchor="middle" fill="var(--muted)">{v}</text>')
    for y, lab in [(150, "40%"), (90, "70%"), (30, "100%")]:
        o.append(f'<line x1="{X0}" x2="{X1}" y1="{y}" y2="{y}" stroke="var(--line)"/><text x="{X0 - 6}" y="{y + 4}" font-size="10" text-anchor="end" fill="var(--muted)">{lab}</text>')
    o.append(f'<polyline points="{pts}" fill="none" stroke="var(--ink)" stroke-width="1.2"/>')
    o.append(f'<text x="{X0}" y="178" font-size="11" fill="var(--muted)">每个位点最常见碱基的频率（细菌，9 nt 滑动平均）</text>')
    y = 196
    for lab, items in rows:
        o.append(f'<text x="{X0 - 6}" y="{y + 11}" font-size="11" text-anchor="end" fill="var(--ink)">{lab}</text>')
        for n, a, b, _ in items:
            col = "var(--accent)"
            o.append(f'<rect x="{sx(a):.1f}" y="{y}" width="{max(sx(b) - sx(a), 2):.1f}" height="14" rx="2" fill="{col}"/>'
                     f'<text x="{sx(a) + 3:.1f}" y="{y + 11}" font-size="10" fill="var(--surface)">{n}</text>')
        y += 24
    o.append(f'<text x="{X0}" y="{y + 18}" font-size="11" fill="var(--muted)">E. coli 16S 位置 1–1542；同一行内的扩增子互不相交，不同行之间可以重叠。</text>')
    o.append("</svg>")
    return "\n".join(o)


def tbl(df, cols, heads, num=()):
    h = "".join(f"<th>{c}</th>" for c in heads)
    b = ""
    for _, r in df.iterrows():
        b += "<tr>" + "".join(
            f'<td class="{"seq" if c == "primer" else ("n" if c in num else "")}">{html.escape(str(r[c]))}</td>' for c in cols) + "</tr>"
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def plan_tables(name):
    P = pd.DataFrame(D[name]["primers"])
    A = pd.DataFrame(D[name]["amplicons"])
    P["bacteria"] = (P.bacteria * 100).round(0).astype(int).astype(str) + "%"
    P["key_min"] = (P.key_min * 100).round(0).astype(int).astype(str) + "%（" + P.worst + "）"
    P["key_mean"] = (P.key_mean * 100).round(0).astype(int).astype(str) + "%"
    P["archaea"] = (P.archaea * 100).round(0).astype(int).astype(str) + "%"
    ta = tbl(A, ["amplicon", "pool", "start", "end", "length", "regions", "F_site", "R_site"],
             ["扩增子", "管", "起点", "终点", "长度 bp", "完整覆盖", "正向位点", "反向位点"], num=("start", "end", "length"))
    tp = tbl(P, ["name", "pool", "site", "primer", "length", "fold", "GC", "Tm", "key_mean", "key_min", "archaea"],
             ["引物", "管", "E. coli 位置", "序列 5′→3′", "nt", "简并数", "GC", "Tm °C", "主要门平均", "主要门最低", "古菌"],
             num=("length", "fold", "GC", "Tm", "key_mean", "key_min", "archaea"))
    return ta + tp


fr_show = fr.copy()
fr_show["amplicons"] = fr_show["amplicons"].fillna("")
frt = tbl(fr_show, ["key min", "pools", "n V", "complete V", "amplicons"],
          ["主要门最低覆盖率 ≥", "管数", "完整覆盖的可变区数", "完整覆盖的可变区", "扩增子（起点-终点，A/B 为管）"], num=("key min", "pools", "n V"))

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>5R 引物重设计：加入 V4，避免引物位点重叠</h1>
  <p class="lede">用 Greengenes 13_8 的 4,797 条 85% OTU 代表序列（4,536 条细菌、261 条古菌）统一投影到 <i>E. coli</i> 坐标，扫描每一个 18–20 nt 位点的保守性，再在“扩增子 180–250 bp、同管引物位点互不相交”的约束下求最优平铺。</p>
  <p>本页是计算机模拟的设计方案。覆盖率按“最多 1 个错配、3′ 端 3 个碱基必须匹配”计算，数据来自较小的 85% 聚类集合，没有做过任何实验验证，下单合成之前需要先在更大的数据库和实际样本上复核。</p>
</header>

<h2>一、现有 5R 引物的坐标核对</h2>
<p>把 10 条引物逐条定位到 <i>E. coli</i> 坐标（允许简并碱基）。结论：五个扩增子之间<b>没有重叠</b>，相邻扩增子的最近引物位点间隔分别为 R1/R2 5 nt、R3/R4 16 nt，其余相隔 100 nt 以上，所以不存在“短产物抢引物”的情况。扩增子内部（去掉引物）完整覆盖了 V2、V3、V5、V6、V8；V1、V4（576–682）、V7、V9 没有被覆盖。</p>
<div class="tbl"><table><thead><tr><th>区域</th><th>正向位点</th><th>反向位点</th><th>扩增子（含引物）</th><th>与上一个扩增子的间隔</th></tr></thead><tbody>
<tr><td>R1</td><td class="n">103–120</td><td class="n">314–332</td><td class="n">103–332（230 bp）</td><td>–</td></tr>
<tr><td>R2</td><td class="n">338–355</td><td class="n">519–536</td><td class="n">338–536（199 bp）</td><td class="n">5 nt</td></tr>
<tr><td>R3</td><td class="n">685–702</td><td class="n">908–927</td><td class="n">685–927（243 bp）</td><td class="n">148 nt</td></tr>
<tr><td>R4</td><td class="n">944–964</td><td class="n">1087–1104</td><td class="n">944–1104（161 bp）</td><td class="n">16 nt</td></tr>
<tr><td>R5</td><td class="n">1175–1193</td><td class="n">1374–1391</td><td class="n">1175–1391（217 bp）</td><td class="n">70 nt</td></tr>
</tbody></table></div>
<p>更值得关心的是引物本身的通用性。按 12 个主要细菌门各占相同权重计算，10 条现有引物里有 6 条在某个主要门上的覆盖率低于 30%（例如 R1 反向引物对 Fusobacteria 为 0%，R1 正向引物对 Bacteroidetes 为 12%，R4 反向引物对 Planctomycetes 为 29%）。</p>

<h2>二、保守性与可用位点</h2>
<figure class="figure">{track_svg()}
<figcaption>灰蓝色底纹为可变区。引物必须落在曲线接近 100% 的位置，扩增子才能同时覆盖相邻的可变区。</figcaption></figure>
<p>允许 18 nt 窗口里最多 2 个位置低于 85% 保守度（可以用简并碱基补上）时，全长 16S 上只有 14 个这样的“岛”：8–27、329–360、503–537、777–807、882–902、904–942、949–966、968–986、1046–1085、1175–1198、1226–1243、1367–1410、1490–1514、1521–1542。引物只能落在这些岛里，所以扩增子的边界基本由它们决定。几个直接的结果：</p>
<ul>
<li><b>V4（576–682）有两条路</b>：①正向引物在 503–537 岛（515F 附近）、反向引物在 777–807 岛，扩增子 290 bp，主要门最低覆盖率可以做到 85% 以上，但超过了 250 bp；②正向引物放在 556 附近、反向引物在 773–791，扩增子 236 bp，在 250 bp 以内，但 Bacteroidetes 的覆盖率只有 66% 左右。我把 V4 的上限单独放宽到 292 bp，选了第一条。</li>
<li><b>V1、V2 需要放宽覆盖率要求才能纳入</b>：当要求主要门最低覆盖率 ≥ 50% 时，优化器会选出 51–277 的扩增子（227 bp）同时覆盖 V1 和 V2；要求 ≥ 60% 以上时 V1 和 V2 都不会被选中。现有 5R 的 R1 引物覆盖 V2，但它的两条引物对 Bacteroidetes（12%）和 Fusobacteria（0%）几乎不扩增。</li>
<li><b>V9（1435–1465）靠近 3′ 末端</b>，Greengenes 代表序列大量被截短（1492 位处只有 41% 的序列有数据），覆盖率估计偏低、不确定性大。</li>
</ul>

<h2>三、设计约束与优化方法</h2>
<ul>
<li>扩增子含引物 180–250 bp；V4 单独放宽到 292 bp。</li>
<li><b>同一管内</b>任意两个扩增子的区间 [正向引物起点, 反向引物终点] 互不相交，且间隔 ≥ 20 nt。这样一管里不会有一条引物落在另一个扩增子内部，也就不会出现短的跨区产物；剩下的只有相隔很远的外侧引物配对，方案 A 中最短的是 F1+R3（341–986，646 bp），在 FFPE 等降解样本里很难扩出来。</li>
<li><b>不同管之间</b>允许重叠。把 V3、V4、V5 这类相邻的区域分到两管，是用实验步骤换覆盖度的办法。</li>
<li>引物最多 3 个简并位点、总简并数 ≤ 8；在候选位点上按 12 个主要细菌门等权（其余细菌门合计占 1 份）贪心加入简并碱基。</li>
<li>目标函数：Σ 可变区权重 × 被扩增子内部（去掉引物）覆盖的比例，V9 权重 0.6，其余 1。用动态规划精确求解 1 管和 2 管的最优方案。</li>
<li>候选位点筛选：12 个主要门的平均覆盖率与最低覆盖率同时满足阈值。</li>
<li>Tm：固定 3′ 端，把 5′ 端缩短或延长到 17–24 nt，让 Tm 靠近 62 °C（最近邻法，50 mM Na⁺、2 mM Mg²⁺、0.2 mM dNTP、250 nM 引物）。</li>
</ul>

<h2>四、覆盖度与引物质量的取舍</h2>
<p>下表给出主要门最低覆盖率要求逐步放宽时，单管和两管方案能完整覆盖的可变区。</p>
{frt}
<p>读法：要求每个主要门的最低覆盖率 ≥ 75% 时，单管只能完整覆盖 V4、V6、V8 三个可变区，两管可以覆盖 V3–V8 六个。放宽到 60% 才出现 V9，放宽到 50% 才能碰到 V1、V2，而此时引物对个别门的覆盖已经很差。</p>

<h2>五、与现有 5R 的比较</h2>
<div class="tbl"><table><thead><tr><th></th><th>现有 5R</th><th>方案 A</th><th>方案 B</th><th>方案 C</th></tr></thead><tbody>
<tr><td>PCR 管数</td><td class="n">1</td><td class="n">2</td><td class="n">1</td><td class="n">2</td></tr>
<tr><td>扩增子数</td><td class="n">5</td><td class="n">5</td><td class="n">3</td><td class="n">6</td></tr>
<tr><td>完整覆盖的可变区</td><td>V2 V3 V5 V6 V8</td><td>V3 V4 V5 V6 V7 V8</td><td>V4 V6 V8</td><td>V3–V9</td></tr>
<tr><td>同管内引物位点重叠</td><td>无</td><td>无</td><td>无</td><td>无</td></tr>
{"".join(f'<tr><td>{lab}（{n} 条引物平均）</td>' + "".join(f'<td class="n">{HO[HO.set == s][c].mean() * 100:.1f}%</td>' for s in ["5R current", "A_two_pool", "B_single_pool", "C_two_pool_full"]) + '</tr>' for lab, c, n in [("主要门平均覆盖率", "key_mean", "每组"), ("最差主要门覆盖率", "key_min", "每组"), ("古菌覆盖率", "arch", "每组")])}
</tbody></table></div>
<p>覆盖率是用 12,000 条不在设计集合里的 97% OTU 序列（11,734 条细菌、266 条古菌）评估的：最多 1 个错配、3′ 端 3 个碱基必须匹配。这批序列和设计集合属于同一个数据库，很多是设计集合序列的近缘种，所以数字仍然偏乐观，不能当作真实样本里的扩增效率。现有 5R 引物本来允许 2 个错配，在这个更严格的标准下看起来偏低。</p>
<p>方案 A 相对现有 5R 的<b>得失</b>：新增 V4 和 V7；失去 V2（V1、V9 仍未覆盖）；引物的覆盖更均匀；代价是需要两管 PCR，并重新建 k-mer 数据库。方案 A 只有 V4 这一条扩增子的两条引物对古菌都有较高覆盖（A2-F 68%、A2-R 98%），其余扩增子不能扩增古菌。</p>

<h2>六、推荐方案 A：两管，覆盖 V3–V8</h2>
<p>最低主要门覆盖率 ≥ 75%、平均 ≥ 85%。每个扩增子都与同管其他扩增子至少间隔 20 nt。</p>
{plan_tables("A_two_pool")}
<p>3′ 端互补检查（≥5 nt，已按简并展开所有组合）：{"；".join(f"{d['pool']} 管 {d['pair']}（{d['nt']} nt）" for d in D["A_two_pool"]["dimers"])}。其中多条是同一条简并引物的不同展开版本之间互补（自二聚体），这是简并引物的普遍现象，需要在合成前用 OligoAnalyzer 或 Primer3 核对实际最差组合。</p>

<h2>七、备选方案</h2>
<h3>方案 B：单管，覆盖 V4、V6、V8</h3>
<p>适合只做一次 PCR、降解严重的样本。三个扩增子是 V4（515–804，290 bp）、V6（905–1082，178 bp）、V8（1180–1406，227 bp）。V4 扩增子偏长，FFPE 样本里效率会下降。</p>
{plan_tables("B_single_pool")}
<h3>方案 C：两管，覆盖 V3–V9</h3>
<p>在方案 A 的基础上增加 V9（1330–1509，180 bp）。V9 的正向引物对 Chloroflexi 的覆盖率只有 61%，反向引物与 A2-R 之间有 6 nt 的 3′ 互补，适合在完整 DNA 样本上尝试。</p>
{plan_tables("C_two_pool_full")}

<h2>八、必须知道的局限</h2>
<ul>
<li><b>Tm 不够均匀</b>：方案 A 的引物 Tm 从 52 到 74 °C 都有，A1-R 和 A2-F 的 GC 含量高达 76%。因为 3′ 端必须落在保守位置，只能调整 5′ 端，调不齐。真实使用前应该用 Primer3 或 IDT 的工具在每个位点上重新做 Tm 平衡，必要时微调 3′ 端的位置（向左右移动 1–3 nt）。</li>
<li><b>参考集偏小且偏旧</b>：设计用的是 4,797 条 85% OTU 代表序列，对每个门的深度有限，OP11、OD1 这类候选门只有几十条，统计不稳；Greengenes 13_8 本身是 2013 年的数据。本环境的网络策略下载不了 SILVA 和 GTDB，换用这两个库重新评估是下一步。</li>
<li><b>古菌</b>：方案 A 里只有 V4 扩增子（A2）的两条引物能扩增古菌，其他区域几乎不能，所以古菌只能得到一个区域的信息。若需要古菌，要另外加入专门引物。</li>
<li><b>下游算法需要改</b>：SMURF 现在假定 5 个区域、数据库按 5 个区域建好。换成 6 个区域和两管需要用 <code>build_region_db_from_fasta</code> 重新生成 k-mer 数据库，并把 <code>scott_format_newer_func.m</code> 里写死的 <code>nR = 5</code> 改掉。两管还需要考虑两管各自的总读数不同，要在 A 矩阵里体现。</li>
<li><b>实验验证</b>：至少用一个已知组成的 mock community，把每个区域、每条引物的实际扩增效率测出来，再决定要不要调整。</li>
</ul>

<footer>脚本：<code>python_5R/design_primers.py</code> · <code>python_5R/primer_design.py</code>；数据：<code>docs/primer_design/</code>（全部候选位点 sites_key.csv、保守性剖面、设计结果 designs_final.json）。</footer>
</main></div>
"""
open(OUT, "w", encoding="utf-8").write(head + body)
print("written", OUT, len(head + body))
