import streamlit as st

def parse_fasta(file):
    sequences = {}
    current_name = None
    for line in file:
        line = line.strip()
        if line.startswith('>'):
            current_name = line[1:]
            sequences[current_name] = ''
        elif current_name:
            sequences[current_name] += line.upper()
    return sequences

# IUPAC対応
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
    '-': set()
}

def match(base1, base2):
    return bool(IUPAC.get(base1, {base1}) & IUPAC.get(base2, {base2}))

def calc_identity(seq1, seq2):
    matches = 0
    total = 0
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        total += 1
        if match(a, b):
            matches += 1
    return matches / total if total > 0 else 0

def find_match_window(full_seq, ref_seq, threshold):
    ref_len = len(ref_seq)
    best_score = 0
    best_result = None
    for i in range(len(full_seq) - ref_len + 1):
        window = full_seq[i:i+ref_len]
        score = calc_identity(window, ref_seq)
        if score >= threshold and score > best_score:
            best_result = (i, window, score)
            best_score = score
    return best_result

# Streamlit UI
st.title("FASTAフィルタリングツール（スライド比較 & IUPAC対応）")

uploaded = st.file_uploader("FASTAまたはTXTファイルをアップロード", type=["fasta", "fa", "txt"])
ref_seq = st.text_input("参照配列（例：AAAGTG）").upper()
threshold = st.slider("一致率の閾値（％）", 50, 100, 90) / 100

if uploaded and ref_seq:
    fasta_text = uploaded.read().decode('utf-8').splitlines()
    seqs = parse_fasta(fasta_text)

    matched = {}
    for name, seq in seqs.items():
        result = find_match_window(seq, ref_seq, threshold)
        if result:
            matched[name] = seq

    st.markdown(f"###  {len(matched)} 件の一致した配列が見つかりました")

    if matched:
        output = ''
        for name, seq in matched.items():
            output += f">{name}\n{seq}\n"
        st.download_button(" FASTAをダウンロード", output, file_name="filtered.fasta")
    else:
        st.warning("条件に合致する配列が見つかりませんでした。")
