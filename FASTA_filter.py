# FASTA_filter.py
import streamlit as st
import os

# --- 混合塩基の一致判定マップ（IUPACコード対応） ---
IUPAC_TABLE = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'G', 'C'}, '-': set()
}

def load_fasta(file_content):
    sequences = {}
    current_id = None
    for line in file_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ''
        else:
            if current_id:
                sequences[current_id] += line.upper()
    return sequences

def is_match(base1, base2):
    return bool(IUPAC_TABLE.get(base1.upper(), set()) & IUPAC_TABLE.get(base2.upper(), set()))

def calc_identity(seq, ref):
    match = 0
    total = 0
    for a, b in zip(seq, ref):
        if a == '-' or b == '-':
            continue
        total += 1
        if is_match(a, b):
            match += 1
    return match / total if total > 0 else 0

def filter_sequences_partial(sequences, ref_seq, threshold):
    filtered = {}
    for id, seq in sequences.items():
        for i in range(len(seq) - len(ref_seq) + 1):
            window = seq[i:i+len(ref_seq)]
            identity = calc_identity(window, ref_seq)
            if identity >= threshold:
                filtered[id] = window  # 部分配列のみ保存
                break
    return filtered


st.title("FASTAフィルタリングツール")

uploaded_file = st.file_uploader("FASTA形式またはTXT形式のファイルをアップロード", type=["fasta", "fa", "txt"])
ref_seq = st.text_input("参照配列（AGCTなど）を入力", max_chars=10000).upper()
threshold = st.slider("一致率の閾値（%）", 50, 100, 90)

if uploaded_file and ref_seq:
    content = uploaded_file.read().decode('utf-8')
    sequences = load_fasta(content)
    st.write(f"読み込んだサンプル数: {len(sequences)}")

    filtered = filter_sequences_partial(sequences, ref_seq, threshold / 100)
    st.write(f"閾値以上一致したサンプル数: {len(filtered)}")

    if filtered:
        st.success(f"{len(filtered)} 件の一致部分配列が見つかりました")
        output = '\n'.join(f">{id}\n{seq}" for id, seq in filtered.items())
        st.download_button("結果をダウンロード", output, file_name="filtered_partial.fasta")
    else:
        st.warning("一致した部分配列が見つかりませんでした")
