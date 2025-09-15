import io
import numpy as np
import pandas as pd
import streamlit as st

# --- Utilities (inline for MVP; move to crispr/ modules later) ---

def find_pams(seq, pam="NGG"):
    """Return indices where a PAM occurs on + strand.
       'N' matches any base; PAM index refers to PAM start."""
    seq = seq.upper()
    hits = []
    for i in range(len(seq) - len(pam) + 1):
        window = seq[i:i+len(pam)]
        if all(p == 'N' or p == b for p, b in zip(pam, window)):
            hits.append(i)
    return hits

def guides_from_pams(seq, pam_idxs, guide_len=20, pam="NGG"):
    """For SpCas9, guide is the 20 nt immediately upstream of PAM on + strand."""
    seq = seq.upper()
    records = []
    for p in pam_idxs:
        start = p - guide_len
        end = p
        if start >= 0:
            guide = seq[start:end]
            pam_seq = seq[p:p+len(pam)]
            records.append({"start": start, "end": end, "guide": guide, "pam": pam_seq, "strand": "+"})
    return records

def gc_content(s):
    s = s.upper()
    gc = sum(1 for c in s if c in ("G","C"))
    return round(100.0 * gc / max(1, len(s)), 1)

def has_homopolymer(s, k=4):
    s = s.upper()
    return any(base * k in s for base in "ACGT")

def hamming(a, b):
    if len(a) != len(b): return np.inf
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def local_offtarget_proxy(seq, guide, max_mismatches=2):
    """Toy proxy: count windows in the PROVIDED sequence that are within <= mismatches.
       Educational only; NOT a genome-wide search."""
    L = len(guide)
    seq = seq.upper()
    guide = guide.upper()
    count = 0
    for i in range(0, len(seq) - L + 1):
        window = seq[i:i+L]
        if hamming(window, guide) <= max_mismatches and window != guide:
            count += 1
    return count

def composite_score(gc, homopolymer, off_proxy):
    """Simple heuristic to rank (higher is better)."""
    # Penalize extreme GC, homopolymers, and off-target proxy
    gc_pen = abs(gc - 50) / 50  # 0 at 50%, 1 at 0% or 100%
    score = 1.0 - 0.5*gc_pen - 0.3*(1 if homopolymer else 0) - 0.2*np.tanh(off_proxy/3)
    return round(max(0.0, score), 3)

# --- Streamlit UI ---

st.set_page_config(page_title="CRISPR Guide Explorer (Educational)", layout="wide")
st.title("🧬 CRISPR Guide Explorer (Educational)")

with st.sidebar:
    st.markdown("**Input sequence (educational only)**")
    seq_mode = st.radio("Choose input source:", ["Paste sequence", "Upload FASTA (toy)"], index=0)
    pam = st.selectbox("PAM (SpCas9 default)", ["NGG"], index=0)
    guide_len = st.number_input("Guide length", min_value=18, max_value=23, value=20, step=1)
    max_mm = st.slider("Off-target proxy mismatches (toy)", 0, 3, 2)
    run_btn = st.button("Find guides")

seq = ""

if seq_mode == "Paste sequence":
    seq = st.text_area("Paste a DNA sequence (A/C/G/T). Use synthetic or toy regions.", height=180).strip()
else:
    upload = st.file_uploader("Upload a small FASTA (toy/synthetic).", type=["fa","fasta","txt"])
    if upload:
        from Bio import SeqIO
        content = upload.getvalue()
        handle = io.StringIO(content.decode())
        recs = list(SeqIO.parse(handle, "fasta"))
        if recs:
            options = {f"{r.id} (len {len(r.seq)})": str(r.seq) for r in recs}
            pick = st.selectbox("Choose record", list(options.keys()))
            seq = options.get(pick, "")

if run_btn:
    if not seq:
        st.warning("Please provide a sequence first.")
        st.stop()

    seq = "".join([c for c in seq.upper() if c in "ACGT"])  # basic sanitize
    if len(seq) < guide_len + 3:
        st.warning("Sequence is too short for scanning.")
        st.stop()

    pam_idxs = find_pams(seq, pam=pam)
    cands = guides_from_pams(seq, pam_idxs, guide_len=guide_len, pam=pam)

    rows = []
    for r in cands:
        guide = r["guide"]
        gc = gc_content(guide)
        hom = has_homopolymer(guide, k=4)
        offp = local_offtarget_proxy(seq, guide, max_mismatches=max_mm)
        score = composite_score(gc, hom, offp)
        rows.append({
            "start": r["start"],
            "end": r["end"],
            "guide": guide,
            "PAM": r["pam"],
            "strand": r["strand"],
            "GC%": gc,
            "homopolymer>=4": hom,
            "offtarget_proxy": offp,
            "score": score
        })

    if not rows:
        st.info("No PAMs found with current settings.")
        st.stop()

    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    st.subheader("Candidate guides")
    st.dataframe(df, use_container_width=True)

    st.download_button(
        "Download CSV",
        data=df.to_csv(index=False).encode(),
        file_name="guides_educational.csv",
        mime="text/csv"
    )

    # Simple position plot (plotly)
    import plotly.express as px
    pos_df = df.copy()
    pos_df["pos"] = pos_df["start"]
    fig = px.scatter(pos_df, x="pos", y="score", hover_data=["guide","GC%","offtarget_proxy"], title="Guide positions vs score")
    st.plotly_chart(fig, use_container_width=True)

st.caption("⚠️ Educational use only. No clinical or experimental guidance.")
