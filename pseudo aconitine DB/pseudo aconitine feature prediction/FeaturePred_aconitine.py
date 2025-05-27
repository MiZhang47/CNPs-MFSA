import pandas as pd
import re
import os

# 定义shift集合（你已给出，省略复制）
# Define shift value of A parts
A_shift_set = {
    'A1-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A1-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A1-3': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
    'A2-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A2-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A3-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A3-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A3-3': {'C': -2, 'H': -4, 'O': -2, 'N': 0},
    'A4-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
    'A4-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'A5-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0}
}

# Define shift value of B parts
B_shift_set = {
    'B1-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'B1-2': {'C': -2, 'H': -4, 'O': -2, 'N': 0},

    'B2-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'B2-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},

    'B3-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
    'B3-2': {'C': 0, 'H': -3, 'O': 0, 'N': -1},

    'B4-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
    'B4-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},

    'B5-1': {'C': -7, 'H': -6, 'O': -2, 'N': 0},
    'B5-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},

    'B6-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
    'B6-2': {'C': 0, 'H': -2, 'O': -1, 'N': 0},

    'B7-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
}

# Define shift value of C parts
C_shift_set = {
    'C1-1': {'C': -7, 'H': -6, 'O': -2, 'N': 0},
    'C1-2': {'C': 0, 'H': -2, 'O': -1, 'N': 0},

    'C2-1': {'C': -8, 'H': -8, 'O': -3, 'N': 0},
    'C2-2': {'C': 0, 'H': -2, 'O': -1, 'N': 0},

    'C3-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
}

# Define shift value of C parts
D_shift_set = {
    'D1-1': {'C': 0, 'H': -2, 'O': -1, 'N': 0},
    'D1-2': {'C': -1, 'H': -4, 'O': -1, 'N': 0},

    'D2-1': {'C': -1, 'H': -4, 'O': -1, 'N': 0},
}

# Define shift value of C parts
E_shift_set = {
    'E1': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'E2': {'C': -1, 'H': -2, 'O': 0, 'N': 0},
    'E3': {'C': -2, 'H': -4, 'O': 0, 'N': 0},
}

shift_sets = {'A': A_shift_set, 'B': B_shift_set, 'C': C_shift_set, 'D': D_shift_set, 'E': E_shift_set}

# 分子式解析函数
def parse_formula(formula):
    elements = {'C':0, 'H':0, 'O':0, 'N':0}
    for elem, num in re.findall(r'([CHON])(\d*)', formula):
        elements[elem] = int(num) if num else 1
    print(f"Parsed {formula}: {elements}")
    return elements

# 分子式拼接函数
def build_formula(elem_counts):
    formula = ''
    for elem in ['C', 'H', 'N', 'O']:
        count = elem_counts.get(elem, 0)
        if count > 0:
            formula += f"{elem}{count}"
    return formula

# 主函数
def predict_fragments(input_csv):
    df = pd.read_csv(input_csv)

    all_feature_ions = []

    for idx, row in df.iterrows():
        print(f"\nProcessing row {idx+1}")
        current_formula = parse_formula(row['MolecularFormula'])
        feature_ions = []

        for part in ['A', 'B', 'C', 'D', 'E']:
            part_key = row[f'{part} parts']
            shifts = [shift for shift in shift_sets[part] if shift.startswith(part_key)]
            print(f"Applying shifts for {part} parts '{part_key}': {shifts}")

            for shift in shifts:
                shift_value = shift_sets[part][shift]
                for elem in ['C', 'H', 'O', 'N']:
                    current_formula[elem] += shift_value[elem]
                new_formula = build_formula(current_formula)
                print(f"  After applying shift '{shift}': {new_formula}")
                feature_ions.append(new_formula)

        all_feature_ions.append(feature_ions)

    # 生成新表格
    max_ions = max(len(ions) for ions in all_feature_ions)
    ion_columns = [f'ion formula{i+1}' for i in range(max_ions)]

    ion_df = pd.DataFrame(all_feature_ions, columns=ion_columns)
    output_df = pd.concat([df, ion_df], axis=1)

    output_path = os.path.splitext(input_csv)[0] + '_fragments.csv'
    output_df.to_csv(output_path, index=False)
    print(f"\nResult saved to {output_path}")

# 你的csv文件路径
input_csv = r"C:\Users\zhang\Desktop\aconitine\pseudo-aconitine library.csv"
predict_fragments(input_csv)