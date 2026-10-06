"""引物命名：旧的分析内部名称 → 订购名称（16S-A1 … 16S-A6，对应 6 个扩增子）。

扩增子编号（沿基因从 5′ 到 3′）：
  A1 V1·V2 · A2 V3 · A3 V4 · A4 V5 · A5 V6·V7 · A6 V8·V9
V5 扩增子的正向引物与 V4 的正向引物（A3-F）共用，所以 A4 只有反向引物（A4-R）。
分析里的旧名称：A1–A3 不变，旧 V5-R → A4-R，旧 A4（V6·V7）→ A5，旧 A5（V8·V9）→ A6。
"""
import re
SLOT = {"A1": "A1", "A2": "A2", "A3": "A3", "V5": "A4", "A4": "A5", "A5": "A6"}
_TOK = re.compile(r"(?<![A-Za-z0-9_])(SMURF5-|S5-)?(A[1-5]|V5)-([FR])(?![A-Za-z0-9])")
def rename_text(s):
    """把文本里所有旧的引物位点名（A1-F、V5-R.1、SMURF5-A4-F.2 …）换成新名；SMURF5-/S5- 前缀换成 16S-。"""
    def f(m):
        pre = "16S-" if m.group(1) else ""
        return f"{pre}{SLOT[m.group(2)]}-{m.group(3)}"
    out = _TOK.sub(f, s)
    return re.sub(r"((?:16S-)?A[1-6]-[FR](?:\.\w+)?)(?:（补充）|\(补\))", r"\1s", out)
def new_order(old):
    """订购名称（ASCII）：旧订购名（如 SMURF5-A1-R.4（补充））→ 16S-A1-R.4s"""
    n = old.replace("SMURF5-", "").replace("（补充）", "s").replace("(补)", "s")
    n = n.replace("A1-F.V1f（SNAP）", "A1-F.SNAP").replace("V5-R.907R（SNAP）", "V5-R.SNAP")
    n = n.replace("A1-F.alt-V1f", "A1-F.SNAP").replace("V5-R.alt-907R单条", "V5-R.SNAP").replace("A2-R.alt短版", "A2-R.alt-short")
    n = rename_text(n)
    return "16S-" + n
def new_slot(old):
    a, b = old.split("-")
    return f"{SLOT[a]}-{b}"
REGION = {"A1": "V1·V2", "A2": "V3", "A3": "V4", "A4": "V5", "A5": "V6·V7", "A6": "V8·V9"}
if __name__ == "__main__":
    for t in ["SMURF5-A1-F.V1f（SNAP）", "SMURF5-V5-R.907R（SNAP）", "SMURF5-A1-R.4（补充）", "SMURF5-A4-F.alt906", "SMURF5-A5-R", "SMURF5-V5-R.1", "SMURF5-A3-F.1"]:
        print(t, "→", new_order(t))
    print(rename_text("A3-F × V5-R 和 A4-F.2、A5-R、SMURF5-A4-F.3（补充）；V1_f 与 V4_r 不变"))
