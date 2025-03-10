import os
import sqlite3
from pyteomics import mgf, mass
import pandas as pd
from typing import Dict
import logging

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ------------------ 文件路径配置（方便修改） ------------------
MGF_FILE_PATH = 'C:/Users/20220/Desktop/makomokonkei_tricin.mgf'  # 输入 MGF 文件路径
DB_PATH = 'Positive Ion Formula_NL.db'  # 输入 SQLite 数据库文件路径
OUTPUT_CSV_PATH = 'C:/Users/20220/Desktop/makomokonkei_tricin_neutralloss.csv'  # 输出 CSV 文件路径
# ---------------------------------------------------------

# 强度阈值定义
INTENSITY_THRESHOLD = 1e5  # 强度大于 1e5 时才处理


# 计算给定分子式、电荷状态和 ppm 公差的 m/z 范围
def mz_range(formula: str, charge: int, ppm: float) -> tuple:
    mz = mass.calculate_mass(formula=formula, charge=charge)
    mz_min = mz - mz * ppm / 1e6
    mz_max = mz + mz * ppm / 1e6
    return abs(mz_min), abs(mz_max)  # 以元组形式返回 m/z 范围


# 从 SQLite 数据库读取分子式并计算其 m/z 范围
def read_formulas_from_db(db_path: str) -> dict:
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT formula FROM compounds")
        formulas = cursor.fetchall()
        conn.close()
        return {formula[0]: mz_range(formula[0], 1, 5) for formula in formulas}
    except sqlite3.Error as e:
        logging.error(f"Error reading database: {e}")
        raise


# 根据 m/z 值和母体分子式找到最佳匹配的分子式
def find_best_formula(mz_value: float, formula_ranges: dict, precursor_formula: str) -> str:
    matched = [
        (formula, calculate_difference(precursor_formula, formula))
        for formula, (mz_min, mz_max) in formula_ranges.items()
        if mz_min <= mz_value <= mz_max
    ]
    if matched:
        best_formula = sorted(matched, key=lambda x: x[1])[0][0]  # 返回质量差异最小的最佳配方
        return best_formula
    return 'NaN'


# 计算两个公式之间的质量差
def calculate_difference(formula1: str, formula2: str) -> float:
    try:
        mass1 = mass.calculate_mass(formula=formula1)
        mass2 = mass.calculate_mass(formula=formula2)
        return abs(mass1 - mass2)
    except mass.PyteomicsError:
        return float('NaN')  # 如果发现无效公式，则返回 NaN 作为质量差异


# 处理 MGF 文件并找到前体和产物离子的匹配分子式
def process_mgf_file(mgf_path: str, db_path: str, output_csv: str):
    # 检查文件路径是否存在
    if not os.path.exists(mgf_path):
        raise FileNotFoundError(f"MGF file not found: {mgf_path}")
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file not found: {db_path}")

    # 读取数据库中的分子式并计算 m/z 范围
    logging.info("Reading formulas from database...")
    formula_ranges = read_formulas_from_db(db_path)

    # 初始化结果存储列表
    rows = []

    # 读取 MGF 文件并逐一处理谱图
    logging.info("Processing MGF file...")
    spectra = list(mgf.read(mgf_path))  # 将生成器转换为列表，确保可以多次访问

    for spectrum in spectra:
        row_id = spectrum['params']['title']
        precursor_mz = spectrum['params']['pepmass'][0]
        product_ions = spectrum['m/z array']
        intensities = spectrum['intensity array']

        # 筛选产物离子：根据强度阈值过滤
        filtered_ions = [
            (ion, intensity) for ion, intensity in zip(product_ions, intensities)
            if intensity >= INTENSITY_THRESHOLD and ion < precursor_mz and ion > 100
        ]

        # 如果筛选后的列表为空，跳过处理
        if not filtered_ions:
            continue

        # 新增：对筛选后的离子按整数部分分组，并保留每组中强度最大的离子
        filtered_ions_df = pd.DataFrame(filtered_ions, columns=['m/z', 'intensity'])
        filtered_ions_df['integer_m/z'] = filtered_ions_df['m/z'].astype(int)  # 添加整数部分列
        filtered_ions_df = filtered_ions_df.loc[
            filtered_ions_df.groupby('integer_m/z')['intensity'].idxmax()
        ]  # 按整数部分分组，并保留强度最大的

        # 查找前体离子的分子式
        precursor_formula = find_best_formula(precursor_mz, formula_ranges, precursor_formula='')

        # 遍历保留的产物离子并处理
        for _, row in filtered_ions_df.iterrows():
            product_mz = row['m/z']
            product_formula = find_best_formula(product_mz, formula_ranges, precursor_formula)

            # 计算质量差
            mass_diff = calculate_difference(precursor_formula, product_formula)

            # 添加到结果行
            rows.append({
                'row_id': row_id,
                'precursor_mz': precursor_mz,
                'precursor_formula': precursor_formula,
                'product_mz': product_mz,
                'product_formula': product_formula,
                'mass_diff': mass_diff
            })

    # 将结果保存为 DataFrame 并导出到 CSV
    logging.info("Saving results to CSV...")
    result_df = pd.DataFrame(rows)
    result_df.to_csv(output_csv, index=False, float_format='%.6f')
    logging.info(f"Output saved to {output_csv}")


# 主程序入口
if __name__ == "__main__":
    try:
        process_mgf_file(MGF_FILE_PATH, DB_PATH, OUTPUT_CSV_PATH)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
