#!/bin/bash
# 复现自建 UDP（DN）480 对。在仓库根目录运行：bash docs/primer_design/udp_denovo/reproduce.sh <输出目录>
# 依赖：python3, numpy, rapidfuzz, primer3-py, openpyxl；评估部分（udp_tail_eval.py）还需要克隆 wangdang511/Olivar_primer 到 /home/user/wangdang511/olivar_primer。
# 全程用到的输入：docs/primer_design/primer_order_6amp_pool33_v3.xlsx（33 条引物池）和 docs/primer_design/udp_tail/udp_all_1_384.tsv（Illumina UDP0001–0384 的 i5/i7，仅用于排除相同序列）。
set -e
OUT=${1:-/tmp/udp_denovo_repro}; mkdir -p $OUT/stage1 $OUT/ext
export PYTHONPATH=python_5R
X=docs/primer_design/primer_order_6amp_pool33_v3.xlsx; ILL=docs/primer_design/udp_tail/udp_all_1_384.tsv; E=python_5R/explore

# 阶段 1：单条规则筛选 + 独立集（Levenshtein ≥4，与反向互补 ≥3，符号顺序 GCAT 的字典序贪心；确定性，无随机数）→ 836 条
SYM=GCAT ORDER=lex python3 $E/udp_denovo384.py $X $ILL $OUT/stage1 4 3 -4 -3 11
# 阶段 2：从 836 条里选 768 条，分 16 个 48 组（8 个 i5 + 8 个 i7），每组每个周期配色平衡；随机种子 4，迭代 900000
python3 $E/udp_denovo384_stage2.py $X $OUT/stage1/indep.pkl $OUT/stage2_blocks.json 4 900000
# 验证并生成 384 对（i5/i7 配对 = 同块同序号）
python3 $E/udp_denovo384_verify.py $OUT/stage2_blocks.json $ILL $OUT/dn384_final.json
# 扩展到 480 对：保留 384 对，分别为 i5、i7 追加同类兼容候选（编辑距离 ≥4、反向互补 ≥3），符号顺序 i5 用 GCAT、i7 用 CGTA
python3 $E/udp_denovo_extend.py $OUT/dn384_final.json $ILL 4 3 GCAT $OUT/ext/ext_i5.pkl i5
python3 $E/udp_denovo_extend.py $OUT/dn384_final.json $ILL 4 3 CGTA $OUT/ext/ext_i7.pkl i7
# 选出新增 96 个 i5 + 96 个 i7（各分 2 个 48 组，配色平衡），种子 1
python3 $E/udp_denovo_480.py $OUT/dn384_final.json $OUT/ext/ext_i5.pkl $OUT/ext/ext_i7.pkl $OUT/dn480.json 1
echo "完成：$OUT/dn480.json（i5、i7 各 480 条；前 384 对与 data/dn384_final.json 相同）"
# 之后的管内评估：
#   python3 $E/udp_tail_eval.py $X <对应的 tsv> <输出 pkl> 58
#   python3 $E/udp_dn_select384.py ...  （从评估后的候选里选 48/96 组）
