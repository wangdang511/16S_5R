"""生成 docs/primer_offtarget.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>人基因组脱靶检查</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
S = pd.read_csv(f"{HERE}/offtarget_summary.csv"); O = pd.read_csv(f"{HERE}/offtarget_per_oligo.csv")
rows = []
for sn in ("5R 现有", "最终 4 扩增子", "V1V2 追加"):
    d = S[S.set == sn].set_index("mm")
    rows.append([sn, int(d.oligos.iloc[0])] + [f"{int(d.loc[m].loci):,}" for m in (1, 2, 3)] + [int(d.loc[m].products) for m in (2, 3)])
t1 = tbl(rows, ["引物集合", "寡核苷酸数", "≤1 错配的位点", "≤2 错配的位点", "≤3 错配的位点", "≤2 错配的潜在产物", "≤3 错配的潜在产物"])
O = O.sort_values("loci_per_expansion", ascending=False)
E = pd.read_csv(f"{HERE}/extend_candidates.csv"); D = pd.read_csv(f"{HERE}/extend_design_eval.csv"); EP = pd.read_csv(f"{HERE}/extend_offtarget_per_oligo.csv"); ES = pd.read_csv(f"{HERE}/extend_offtarget_summary.csv")
CH = {("A2-F", 0): 3, ("A2-R", 0): 2, ("A4-R", 0): 3, ("A5-F", 0): 3, ("A5-F", 1): 3, ("A1-F", 0): 0, ("A1-R", 0): 1, ("A1-R", 1): 1}
r3 = []
for (sn, i), e in CH.items():
    o = E[(E.site == sn) & (E.i == i) & (E.e == 0)].iloc[0]; n = E[(E.site == sn) & (E.i == i) & (E.e == e)].iloc[0]
    r3.append([f"{sn}_{i}", o.seq, n.seq if e else "（未延长）", f"+{e}", f"{o.Tm} → {n.Tm}", f"{o.SILVA * 100:.0f}% / {o.GG * 100:.0f}% → {n.SILVA * 100:.0f}% / {n.GG * 100:.0f}%", f"{o.loci2 / o.nexp:g} → {n.loci2 / n.nexp:g}"])
t3 = tbl(r3, ["寡核苷酸", "原序列", "延长后", "延长 nt", "Tm °C", "位点覆盖率 SILVA / GG", "每个展开序列的人基因组位点（≤2 错配）"], seqcol=2)
t4 = tbl([[r.design, int(r.oligos), int(r.expansions), f"{r.Tm_min}–{r.Tm_max}", int(r.severe), r.SILVA_amp, r.GG_amp, f"{r.SILVA_all * 100:.1f}% / {r.GG_all * 100:.1f}%", f"{r.ideal * 100:.1f}%", f"{r.abs60 * 100:.1f} / {r.abs300 * 100:.1f} / {r['frac0.8'] * 100:.1f}%"] for _, r in D.iterrows()],
         ["方案", "寡核苷酸", "展开数", "Tm °C", "严重二聚体", "扩增子覆盖 SILVA", "扩增子覆盖 GG", "全部扩出 SILVA / GG", "理想准确率", "属准确率（三种规则）"])
t5 = tbl([[r.set.replace("扩增子延长后", " 个扩增子（延长后）") if False else r.set, int(r.mm), f"{int(r.loci):,}", int(r.products)] for _, r in ES.iterrows()], ["延长后的集合", "错配上限", "人基因组位点", "潜在产物（80–1500 bp）"])
t2 = tbl([[r.set, r.oligo, r.seq, r.nt, r.expansions, f"{r.loci_le2mm:,}", f"{r.loci_per_expansion:g}"] for _, r in O.iterrows() if r.loci_per_expansion >= 2],
         ["集合", "寡核苷酸", "序列 5′→3′", "nt", "展开数", "≤2 错配的位点", "每个展开序列的位点数"], seqcol=2)
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>对人基因组和线粒体的特异性检查</h1>
  <p class="lede">线粒体基本安全，没有找到能形成产物的引物对；但对核基因组，新设计里 16–17 nt 的短寡核苷酸有大量近似匹配位点，在人 DNA 占比高的样本（比如 FFPE）里会浪费读段。5R 原引物（≥18 nt）几乎没有这个问题。</p>
  <p>这是计算机检查：只用序列匹配，没有做实验，也没有算热力学。</p>
</header>

<h2>一、方法</h2>
<ul>
<li>基因组：hg19（UCSC，含 chrM 和未定位的 contig，用的是 bioconda 里的 BSgenome 包），另外单独用 rCRS（NC_012920）检查线粒体。NCBI、UCSC 和 Ensembl 网站在这个环境里访问不了，所以用的是包里的版本，不是 hg38 或 T2T。</li>
<li>每条引物的所有简并展开序列，在两条链上找：3′ 端 8 nt 必须完全匹配，整体错配数 ≤1、2、3 的位点；再把方向相对、3′ 端相对、间距 80–1500 bp 的命中配成潜在产物（任意两个引物都可配对，因为是同一个管）。</li>
</ul>

<h2>二、结果</h2>
{t1}
<ul>
<li><b>线粒体（rCRS）</b>：最终 4 个扩增子的寡核苷酸有 4 个命中（≤4 错配），都在 12S rRNA 基因里：A2-R_1 在 883（2 个错配）、A3-R_4 在 1072（4 个）、A5-R_9 在 1556（3 个和 4 个）。V1·V2 追加的引物有 5 个命中，最好的是 A1-R_0 在 13851（1 个错配，ND5 基因），其余 4 个错配。5R 只有 R3-R 在 1160 一个命中（4 个错配）。所有集合在线粒体上都找不到相对方向的命中对，所以没有潜在产物。</li>
<li><b>核基因组</b>：5R 在 ≤2 个错配下有 114 个位点、没有潜在产物；最终 4 个扩增子有 9,788 个位点、15 个潜在产物；V1·V2 追加有 3,283 个位点、1 个潜在产物。放宽到 ≤3 个错配，最终 4 个扩增子有 490 个潜在产物。</li>
<li><b>风险集中在短寡核苷酸</b>：见下表。≥19 nt 的引物（A3-F、A3-R、A4-F、A5-R、V1·V2 的 A1-R_2）每个展开序列只有 0–2 个位点，和 5R 差不多；16–17 nt 的（A4-R_6、A2-F、A5-F、A2-R、A1-F）有几百到近 1,000 个。</li>
</ul>
{t2}

<h2>三、含义和建议（延长之前）</h2>
<ul>
<li>脱靶产物不是 16S，不会被 SMURF 的 16S 数据库匹配到，主要的后果是在人 DNA 很多的样本里浪费读段、拉低有效数据量；会不会真的扩出，取决于退火温度和人 DNA 量，这里没法判断。</li>
<li>最直接的修正：给这几条短寡核苷酸的 5′ 端各延长 2–4 个碱基（沿 16S 模板），这样既减少人基因组的匹配，也能同时抬高它们的 Tm（A5-F 的两条本来就偏低），代价是覆盖率会略降，需要重新评估。这一步我还没做。</li>
<li>不建议为此放弃新设计；也要承认 5R 在这点上更好。</li>
</ul>

<h2>四、按建议延长短寡核苷酸后重新评估</h2>
<p>对 16–17 nt 的寡核苷酸，沿 16S 模板在 5′ 端延长 1–6 个碱基（新增碱基取该寡核苷酸所覆盖序列中的多数碱基），再逐个比较 Tm、位点覆盖率、人基因组命中数。选择规则：每个展开序列在人基因组 ≤2 错配的位点数 ≲30，Tm ≲66 °C，覆盖率损失尽量小。选定的延长长度：</p>
{t3}
<p>整个池评估（A1 = V1·V2；覆盖率是 Greengenes 和 SILVA 留出集，“全部扩出”指所有扩增子都扩出）：</p>
{t4}
{t5}
<ul>
<li><b>人基因组脱靶明显降低</b>：最终 4 个扩增子在 ≤2 个错配下，位点从 9,788 降到 324，潜在产物从 15 降到 0；≤3 个错配的潜在产物从 490 降到 1。V1·V2 追加的位点从 3,283 降到 2,641，≤3 个错配的潜在产物从 67 降到 41（A1-F 没法延长，剩下的主要来自它）。线粒体的命中不变，仍然没有潜在产物。</li>
<li><b>代价</b>：扩增子覆盖率下降约 0–4 个百分点（SILVA 81/78/80/79 → 80/78/78/75，GG 95/91/93/95 → 92/91/91/93），“全部扩出”SILVA 66.8% → 62.7%，GG 80.6% → 76.1%；属准确率（算入扩增失败，三种规则）91.4/90.9/91.5% → 90.2/90.4/91.3%，理想准确率不变（94.2%）。主要原因是 A2-F：延长 1 nt 就让覆盖率从 95%/98% 降到 90%/88%（我的匹配规则是总错配 ≤1，与长度无关，所以更长的引物更难满足；真实 PCR 对长引物更宽容，所以这是偏保守的估计）。</li>
<li><b>Tm 范围变宽</b>：56.5–69.6 °C（原来 50.8–63.5 °C）。A2-R 是 GC 约 78% 的 18 nt，Tm 69.6 °C，是池里最热的；最低的 A5-F_8 约 56.6 °C，比原来的 50.8 °C 好。退火温度要按这个范围重新选。</li>
<li><b>一个临界的自二聚体</b>：A5-F_7 延长后（`GCTRCACRCRTGCTACAAT`）自身 ΔG −9.0 kcal/mol（3′ 端 −3.5），刚到我设的阈值；任何长度的延长都会这样，不延长（16 nt）则没有，但人基因组位点 142 个/展开序列。3′ 端稳定性不强，所以我认为风险不高，但需要实验验证。</li>
<li><b>A1-F 没有延长</b>：它位于 8–24，再往 5′ 延长会超出数据库里序列的数据范围（约 20% 序列在 8 位之前没有数据），评估不了，所以保持 17 nt，人基因组位点 2,329 个（291 个/展开序列）。如果要加 V1·V2，这是剩下的弱点。</li>
</ul>
<h2>五、局限</h2>
<ul>
<li>hg19 不是最新版本，也缺少 rDNA 阵列等重复区域；人 18S/28S rRNA 没有单独检查，核 rDNA 可能比表中更多。</li>
<li>“3′ 端 8 nt 完全匹配 + 总错配数”是粗略规则，没考虑错配位置、GC、二级结构，也没算热力学；位点数不等于会扩出。</li>
<li>只检查了人，没有检查宿主以外的污染或样本里的其他真核生物。</li>
</ul>
<footer>脚本：<code>python_5R/explore/offtarget.py</code>、<code>extend_ext*.py</code>；数据：<code>docs/primer_design/offtarget_*.csv</code>、<code>offtarget_primers.json</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_offtarget.html", "w", encoding="utf-8").write(head + body)
print("ok")
