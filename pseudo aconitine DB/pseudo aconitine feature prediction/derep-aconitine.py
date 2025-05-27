import pandas as pd
import re
from collections import defaultdict

def parse_formula(formula: str) -> dict:
    """将化学式字符串转为元素-数量字典。"""
    pattern = r'([A-Z][a-z]?)(\d*)'
    counts = {}
    for elem, num in re.findall(pattern, formula):
        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)
    return counts

def collect_formulas_from_excel(excel_path: str) -> list:
    """
    1. 只从Excel中读取所有包含 'ion formula' 的列
    2. 合并这些列、去重并按字母序排序
    3. 返回唯一分子式列表
    """
    df = pd.read_excel(excel_path, engine='openpyxl')
    # 找出所有包含 'ion formula' 的列
    ion_cols = [col for col in df.columns if 'ion formula' in col.lower()]
    if not ion_cols:
        print("没有找到包含 'ion formula' 的列")
        return []
    # 合并所有相关列并展开为一列
    all_formulas = df[ion_cols].values.ravel()
    # 去除NaN和空字符串，转为字符串
    all_formulas = [str(f).strip() for f in all_formulas if pd.notna(f) and str(f).strip()]
    # 去重并排序
    unique_formulas = sorted(set(all_formulas))
    return unique_formulas

def classify_by_carbons(formulas: list) -> dict:
    """将分子式列表按碳原子数分组。"""
    groups = defaultdict(list)
    for f in formulas:
        counts = parse_formula(f)
        cnum = counts.get('C', 0)
        groups[cnum].append(f)
    return dict(groups)

# =========================== 主程序 =============================
excel_path = r"C:\Users\zhang\Desktop\aconitine\pseudo-aconitine library_fragments2.xlsx"

formulas = collect_formulas_from_excel(excel_path)

print(",".join(formulas))
print(f"\nTotal unique formulas: {len(formulas)}\n")

# 如果还要按碳数分类
groups = classify_by_carbons(formulas)
print("Classification by carbon atom count:")
for c in sorted(groups):
    lst = groups[c]
    print(f"  C={c} ({len(lst)}): {','.join(lst)}")
