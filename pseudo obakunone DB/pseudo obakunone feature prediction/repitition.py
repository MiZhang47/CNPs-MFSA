import os
import pandas as pd
from collections import defaultdict

# 设置文件夹路径
folder_path = r"C:\Users\zhang\Desktop\limonoid\limonoid-MS2-ref\TypeIV"

# 获取所有CSV文件路径
csv_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".csv")]

# 用于统计Formula在不同文件中出现的次数和记录对应的m/z
formula_info = defaultdict(lambda: {"count": 0, "mz_values": set()})

# 遍历每个CSV文件
for csv_file in csv_files:
    try:
        df = pd.read_csv(csv_file)
        if "Formula" not in df.columns or "m/z" not in df.columns:
            continue

        seen_formulas = set()
        for _, row in df.iterrows():
            formula = str(row["Formula"]).strip()
            mz_value = row["m/z"]

            if formula not in seen_formulas:
                formula_info[formula]["count"] += 1
                seen_formulas.add(formula)

            formula_info[formula]["mz_values"].add(mz_value)
    except Exception as e:
        print(f"❌ 错误读取文件 {csv_file}：{e}")

# 构建结果 DataFrame
total_files = len(csv_files)
result_rows = []

for formula, data in formula_info.items():
    if data["count"] > 1:  # 只保留重复的
        mz_combined = "; ".join([str(mz) for mz in sorted(data["mz_values"])])
        repeat_rate = data["count"] / total_files
        result_rows.append({
            "Formula": formula,
            "m/z values": mz_combined,
            "Repeat Count": data["count"],
            "Repeat Rate": round(repeat_rate, 4)
        })

result_df = pd.DataFrame(result_rows)
result_df = result_df.sort_values(by="Repeat Count", ascending=False)

# 保存结果
output_path = os.path.join(folder_path, "Repeated_Formulas_Summary-Type4.csv")
result_df.to_csv(output_path, index=False, encoding='utf-8-sig')

print(f"✅ 分析完成，结果已保存至：{output_path}")
