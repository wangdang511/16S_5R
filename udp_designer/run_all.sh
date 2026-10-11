#!/bin/bash
# 一键流程（可断点续算，掉线后原命令重跑即可）。
# 用法：./run_all.sh <工作目录> <panel.csv> [config.json] [最多几轮=3] [现有条形码.txt]
#   每一轮 = s2 热点 → s3 单条尾巴(增量) → s3b 相对分级校准（只在第一次，随机组合整管实测定阈值） → s5 FR(增量) → s6 配对(用累积实测覆盖) → s7 整管实测 → 把代理漏掉的引物对补进热点。
#   一轮里实测没有漏掉的引物对且没有 C 级就提前结束。最后 s8 互换修复、s9 导出。
set -e
W=$1; P=$2; C=${3:-}; R=${4:-3}; EXC=${5:-}
CFG=(); [ -n "$C" ] && CFG=(--cfg "$C")
EX=(); [ -n "$EXC" ] && EX=(--exclude "$EXC")
cd "$(dirname "$0")"; mkdir -p "$W"
python3 s0_universe.py --work "$W" "${CFG[@]}" "${EX[@]}"
[ -f "$W/baseline.npz" ] || python3 s1_baseline.py --work "$W" --panel "$P" "${CFG[@]}"
for r in $(seq 1 "$R"); do
  echo "================ 第 $r 轮 ================"
  if [ ! -f "$W/tubes_r$r.summary.json" ]; then     # 已完成的轮次（掉线重跑时）直接跳过，只做收敛判断
  python3 s2_hotpairs.py --work "$W" --panel "$P" "${CFG[@]}"
  python3 s3_single.py   --work "$W" --panel "$P" "${CFG[@]}"
  if python3 s3b_calibrate.py --work "$W" --panel "$P" "${CFG[@]}" | tee /dev/stderr | grep -q "相对阈值"; then   # 刚完成相对分级校准：按新阈值与补充热点重跑 s2/s3（增量）
    mkdir -p "$W/pre_calib"; mv -f "$W/pool.json" "$W"/fr_*of*.pkl "$W/pre_calib/" 2>/dev/null || true     # 阈值变了：候选池与 FR（提前终止依赖阈值）作废
    python3 s2_hotpairs.py --work "$W" --panel "$P" "${CFG[@]}"; python3 s3_single.py --work "$W" --panel "$P" "${CFG[@]}"; fi
  [ -f "$W/pool.json" ] || python3 s4_pool.py --work "$W" "${CFG[@]}"
  python3 s5_fr.py       --work "$W" --panel "$P" "${CFG[@]}"
  MEAS=(); [ -f "$W/meas.json" ] && MEAS=(--meas "$W/meas.json")
  python3 s6_pair.py     --work "$W" "${CFG[@]}" "${MEAS[@]}" --out "pairs_r$r.json"
  python3 s7_tubes.py    --work "$W" --panel "$P" "${CFG[@]}" --pairs "pairs_r$r.json" --tag "r$r" --update-hot
  else echo "第 $r 轮已完成，跳过"; fi
  LAST=$r
  python3 - "$W" "$r" <<'PY' && break || true
import json,sys
s=json.load(open(f"{sys.argv[1]}/tubes_r{sys.argv[2]}.summary.json"))
done = (not s["missed_pairs"]) and s["counts"]["C"] == 0
print("本轮：漏掉的引物对", len(s["missed_pairs"]), "C 管", s["counts"]["C"], "→", "收敛，结束循环" if done else "继续下一轮")
sys.exit(0 if done else 1)
PY
done
python3 s8_repair.py --work "$W" --panel "$P" "${CFG[@]}" --pairs "pairs_r${LAST:-1}.json" --tag "r${LAST:-1}" --target B
python3 s9_export.py --work "$W" --panel "$P" "${CFG[@]}" --pairs pairs_repaired.json --tubes "tubes_r${LAST:-1}_rep" --out "$W/out" "${EX[@]}"
