import os
import glob
import re
import pandas as pd

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

def build_formula(counts: dict) -> str:
    """
    从元素计数字典重新生成化学式字符串，遵循 Hill 排序（C, H, 其他元素按字母顺序）。
    数量为 1 时省略数字。
    """
    parts = []
    # C, H 优先
    if 'C' in counts:
        parts.append(('C', counts['C']))
    if 'H' in counts:
        parts.append(('H', counts['H']))
    # 其它元素
    for elem in sorted(e for e in counts if e not in ('C', 'H')):
        parts.append((elem, counts[elem]))
    # 组合
    fmt = ''
    for elem, num in parts:
        if num > 0:
            fmt += elem + (str(num) if num != 1 else '')
    return fmt

def subtract_one_h(formula: str) -> str:
    """
    对一个化学式减去一个 H。若原式中没有 H，则返回原式不变。
    若 H 减至 0，则在结果中去掉 H 项。
    """
    counts = parse_formula(formula)
    if counts.get('H', 0) >= 1:
        counts['H'] -= 1
        if counts['H'] == 0:
            counts.pop('H')
    # 否则不做修改
    return build_formula(counts)

def collect_and_subtract_h(input_folder: str, output_csv: str):
    # 1. 读取所有 CSV 并提取 Formula 列去重
    csv_files = glob.glob(os.path.join(input_folder, '*.csv'))
    unique_formulas = set()
    for fp in csv_files:
        try:
            df = pd.read_csv(fp)
        except Exception as e:
            print(f"跳过 {fp}：读取失败 ({e})")
            continue
        if 'Formula' not in df.columns:
            print(f"跳过 {fp}：缺少 'Formula' 列")
            continue
        unique_formulas.update(df['Formula'].dropna().astype(str).tolist())

    if not unique_formulas:
        print("未找到任何有效的 Formula。")
        return

    # 2. 对每个去重后的式子减去一个 H
    sorted_orig = sorted(unique_formulas)
    result = []
    for orig in sorted_orig:
        newf = subtract_one_h(orig)
        result.append({'Original': orig, 'Minus_H': newf})

    # 3. 输出到 CSV
    out_df = pd.DataFrame(result)
    out_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"处理完成：共 {len(result)} 条记录，已保存到 {output_csv}")

if __name__ == '__main__':
    # ———— 请根据实际修改下面两行路径 ————
    # 输入：包含多个 CSV 文件的文件夹路径
    input_folder = r"C:\Users\zhang\Desktop\limonoid\limonoid-MS2-ref\summary"
    # 输出：去重后分子式写入的 CSV 文件路径
    output_csv = r'C:\Users\zhang\Desktop\limonoid\limonoid-MS2-ref\unique_formulas.csv'
    collect_and_subtract_h(input_folder, output_csv)

