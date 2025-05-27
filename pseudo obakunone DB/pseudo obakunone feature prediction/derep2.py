import pandas as pd
import re
from collections import defaultdict

def parse_formula(formula: str) -> dict:
    """
    将一个化学式字符串解析成 {元素: 数量} 的字典。
    例如 "C10H16O" → {'C':10, 'H':16, 'O':1}
    """
    pattern = r'([A-Z][a-z]?)(\d*)'
    counts = {}
    for elem, num in re.findall(pattern, formula):
        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)
    return counts

def collect_formulas1(csv_path: str) -> list:
    """
    1. 从 Excel 中读取所有包含 'ion formula' 的列；
    2. 从 CSV 中读取 'Formula' 列；
    3. 合并、去重并按字母序排序，返回唯一分子式列表。
    """
    # 读取 CSV 的 Formula 列
    df_csv = pd.read_csv(csv_path, encoding='utf-8-sig')
    csv_formulas = df_csv['Formula'].dropna().astype(str).tolist()

    # 合并、去重、排序
    unique = sorted(set(csv_formulas))
    return unique

def collect_formulas(excel_path: str, csv_path: str) -> list:
    """
    1. 从 Excel 中读取所有包含 'ion formula' 的列；
    2. 从 CSV 中读取 'Formula' 列；
    3. 合并、去重并按字母序排序，返回唯一分子式列表。
    """
    # 读取 Excel 并挑选含 'ion formula' 的列
    df_xl = pd.read_excel(excel_path, engine='openpyxl')
    ion_cols = [col for col in df_xl.columns if 'ion formula' in col.lower()]
    ion_formulas = (
        df_xl[ion_cols]
        .stack()
        .dropna()
        .astype(str)
        .tolist()
    )

    # 读取 CSV 的 Formula 列
    df_csv = pd.read_csv(csv_path, encoding='utf-8-sig')
    csv_formulas = df_csv['Formula'].dropna().astype(str).tolist()

    # 合并、去重、排序
    unique = sorted(set(ion_formulas + csv_formulas))
    return unique

def classify_by_carbons(formulas: list) -> dict:
    """
    将分子式列表按碳原子数分组，返回 {碳数: [分子式,...], ...} 的字典。
    """
    groups = defaultdict(list)
    for f in formulas:
        counts = parse_formula(f)
        cnum = counts.get('C', 0)
        groups[cnum].append(f)
    return dict(groups)

# excel_path = r"C:\Users\zhang\Desktop\pDB-limonoid\pseudo-pDB-limonoid library\pseudo-pDB-limonoid library.xlsx"
csv_path  = r"C:\Users\zhang\Desktop\limonoid\limonoid-MS2-ref\unique_formulas.csv"

# 收集所有唯一分子式
formulas = collect_formulas1(csv_path)

# 打印所有分子式及总数
print(",".join(formulas))
print(f"\nTotal unique formulas: {len(formulas)}\n")

# 按碳数分类并打印
groups = classify_by_carbons(formulas)
print("Classification by carbon atom count:")
for c in sorted(groups):
    lst = groups[c]
    print(f"  C={c} ({len(lst)}): {','.join(lst)}")
