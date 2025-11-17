#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline UNIFICADO — Middle-down proteomics (GLOBAL vs ISOFORMS)
================================================================

Fluxo:
  1) Carrega múltiplos reports Spectronaut (um por enzima) e normaliza por amostra (wide).
     - Suporta dois formatos:
       A) run-wise:  colunas 'PG.MS2Quantity'/'PEP.MS2Quantity' + 'R.Condition' = {1,2,3}
       B) wide:      várias colunas terminando em '.PG.MS2Quantity' / '.PEP.MS2Quantity'
                    (condição inferida do header)
  2) Minera / anota isoformas não-canônicas (curadoria por prefixos custom + regras).
  3) Constrói 3 VISÕES paralelas (sem recursão de diretórios!):
       - GLOBAL_ALL    : tudo
       - GLOBAL_NO_ISO : exclui grupos mantidos como isoformas não-canônicas
       - ISOFORMS_ONLY : somente isoformas mantidas (mapeia ProteinID -> isoform_kept)
  4) Para cada visão:
       - FCs por enzima (log2): N+ vs CTRL, N0 vs N+, N0 vs CTRL
       - FC GLOBAL = mediana entre enzimas
       - TOP markers (GLOBAL) orientados por condição
       - Ranks MS2 (par, combinado, facet) — paleta fixa por ENZIMA
       - UpSets por enzima e por condição
       - Coverage (usa 'PG.Coverage (Global)', fallback 'PG.Coverage'), filtrado por presença MS2
  5) Opcional: Heatmap de FC (vlag, centrado em 0)

Requisitos:
  pip install pandas numpy matplotlib seaborn upsetplot
"""

import os, re, argparse, warnings, logging
from typing import Optional, List, Dict, Tuple, Iterable
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec


warnings.filterwarnings("ignore", category=FutureWarning, module=r"upsetplot\.plotting")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ============================
# Paletas fixas
# ============================
ENZYME_COLORS = {
    "Tryp":  "#D62728",  # red
    "LysC":  "#FFD700",  # yellow
    "AspN":  "#1f77b4",  # blue
    "Chymo": "#9467bd",  # purple
    "GluC":  "#2ca02c",  # greenish
}
COND_COLORS = {"CTRL": "#7f7f7f", "N0": "#1f77b4", "N+": "#17becf"}
def col_enzyme(e): return ENZYME_COLORS.get(str(e), "#7f7f7f")
def col_cond(c):   return COND_COLORS.get(str(c), "#7f7f7f")

# ============================
# Utils & dirs
# ============================
def ensure_root_dirs(root: str):
    os.makedirs(root, exist_ok=True)
    for sub in [
        "debug",
        "isoforms",
        "views/GLOBAL_ALL/rank_tables","views/GLOBAL_ALL/protein_rank_plots","views/GLOBAL_ALL/upset","views/GLOBAL_ALL/coverage","views/GLOBAL_ALL/heatmaps",
        "views/GLOBAL_NO_ISO/rank_tables","views/GLOBAL_NO_ISO/protein_rank_plots","views/GLOBAL_NO_ISO/upset","views/GLOBAL_NO_ISO/coverage","views/GLOBAL_NO_ISO/heatmaps",
        "views/ISOFORMS_ONLY/rank_tables","views/ISOFORMS_ONLY/protein_rank_plots","views/ISOFORMS_ONLY/upset","views/ISOFORMS_ONLY/coverage","views/ISOFORMS_ONLY/heatmaps",
    ]:
        os.makedirs(os.path.join(root, sub), exist_ok=True)

def safe_minmax(y: pd.Series, pad=0.5):
    if len(y)==0 or not np.isfinite(y).any(): return (-1,1)
    return (np.nanmin(y)-pad, np.nanmax(y)+pad)

def _force_series_1d(obj):
    import pandas as _pd
    return obj.iloc[:, 0] if isinstance(obj, _pd.DataFrame) else obj

def sanitize_long_table(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty: return df
    df = df.loc[:, ~df.columns.duplicated()].copy()
    for k in ["ProteinID", "Enzyme", "Condition", "SampleCol"]:
        if k in df.columns:
            df[k] = _force_series_1d(df[k])
            df[k] = df[k].astype(str).str.strip()
    if "Intensity" in df.columns:
        df["Intensity"] = pd.to_numeric(_force_series_1d(df["Intensity"]), errors="coerce").fillna(0.0)
    return df

# ============================
# Naming helpers
# ============================
RUNWISE_COND_MAP = {1: "CTRL", 2: "N+", 3: "N0", "1": "CTRL", "2": "N+", "3": "N0"}

def infer_enzyme_from_path(path: str) -> Optional[str]:
    t = os.path.basename(path).lower()
    if "aspn" in t:  return "AspN"
    if "gluc" in t:  return "GluC"
    if "lysc" in t:  return "LysC"
    if "tryp" in t:  return "Tryp"
    if "chymo" in t: return "Chymo"
    return None

def parse_filename_from_ms2col(colname: str) -> str:
    m = re.search(r"\]\s*(.+?)\.raw\.", colname, re.I)
    if m: return m.group(1).strip()
    m2 = re.search(r"\]\s*([^\.]+)\.", colname)
    if m2: return m2.group(1).strip()
    m3 = re.search(r"(.+?)\.(?:PG|PEP)\.MS2Quantity$", colname, re.I)
    if m3: return m3.group(1).split("]")[-1].strip()
    return colname

def cond_from_token(tok: Optional[str]) -> Optional[str]:
    if tok is None: return None
    u = str(tok).upper().strip().replace(" ","")
    u = u.replace("N-PLUS","N+").replace("N_PLUS","N+").replace("NPLUS","N+")
    if u.endswith("_N0") or "_N0_" in u or u=="N0": return "N0"
    if u.endswith("_C")  or "_C_"  in u or "CTRL" in u or "CONTROLE" in u: return "CTRL"
    if "N+" in u or u.endswith("_N") or "_N_" in u or u=="N": return "N+"
    return None

# ============================
# Loader (MS2 long-table)
# ============================
def pick_protein_id_col(df: pd.DataFrame, forced: Optional[str]) -> str:
    if forced and forced in df.columns: return forced
    for c in ["PG.ProteinGroups","PG.ProteinAccessions","PG.Genes","PG.ProteinNames"]:
        if c in df.columns: return c
    for c in df.columns:
        if df[c].dtype==object: return c
    raise SystemExit("Não encontrei coluna de ID proteico; use --protein-id-col.")

def pick_ms2_cols_wide(df: pd.DataFrame, prefer_pg=True):
    if prefer_pg:
        cols = [c for c in df.columns if c.endswith(".PG.MS2Quantity")]
        if not cols:
            cols = [c for c in df.columns if c.endswith(".PEP.MS2Quantity")]
    else:
        cols = [c for c in df.columns if c.endswith(".PEP.MS2Quantity")]
        if not cols:
            cols = [c for c in df.columns if c.endswith(".PG.MS2Quantity")]
    return cols

def load_reports_unified(paths: List[str], protein_id_col: Optional[str], prefer="PG",
                         norm="median", outdir: Optional[str]=None) -> pd.DataFrame:
    rows, dbg = [], []
    use_pg = (prefer.upper()=="PG")
    for p in paths:
        enz = infer_enzyme_from_path(p) or "UNK"
        df  = pd.read_csv(p, sep="\t", low_memory=False)
        pidcol = pick_protein_id_col(df, protein_id_col)
        prot_ids = df[pidcol].astype(str).str.split(";").str[0].str.strip()

        # formato?
        has_run_pg  = "PG.MS2Quantity"  in df.columns
        has_run_pep = "PEP.MS2Quantity" in df.columns
        wide_pg     = [c for c in df.columns if c.endswith(".PG.MS2Quantity")]
        wide_pep    = [c for c in df.columns if c.endswith(".PEP.MS2Quantity")]

        if has_run_pg or has_run_pep:
            if "R.Condition" not in df.columns:
                raise SystemExit(f"{os.path.basename(p)} (run-wise) sem R.Condition.")
            ms2_col = "PG.MS2Quantity" if (use_pg and has_run_pg) else ("PEP.MS2Quantity" if has_run_pep else "PG.MS2Quantity")
            intens = pd.to_numeric(df[ms2_col], errors="coerce").fillna(0.0)
            cond   = df["R.Condition"].map(RUNWISE_COND_MAP).fillna("UNK")
            sample = df["R.FileName"] if "R.FileName" in df.columns else ms2_col
            rows.append(pd.DataFrame({
                "ProteinID": prot_ids, "Enzyme": enz, "Condition": cond,
                "SampleCol": sample, "Intensity": intens
            }))
            dbg.append({"file": os.path.basename(p), "format":"run-wise", "ms2_col":ms2_col})

        elif wide_pg or wide_pep:
            cols = pick_ms2_cols_wide(df, prefer_pg=use_pg)
            block = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
            if norm and norm.lower() in {"median","med","libsize"}:
                med = block.replace(0,np.nan).median(axis=0, skipna=True)
                scale = med/np.nanmedian(med) if np.isfinite(med).any() else 1.0
                block = block.div(scale, axis=1).fillna(0.0)
            elif norm and norm.lower() in {"sum","tic"}:
                s = block.sum(axis=0)
                scale = s/np.nanmedian(s) if np.isfinite(s).any() else 1.0
                block = block.div(scale, axis=1).fillna(0.0)
            col2cond = {c: cond_from_token(parse_filename_from_ms2col(c)) for c in cols}
            for c in cols:
                rows.append(pd.DataFrame({
                    "ProteinID": prot_ids, "Enzyme": enz,
                    "Condition": col2cond[c] if col2cond[c] is not None else "UNK",
                    "SampleCol": c, "Intensity": block[c].values
                }))
            dbg.append({"file": os.path.basename(p), "format":"wide", "n_cols":len(cols)})
        else:
            dbg.append({"file": os.path.basename(p), "format":"none"})
            logging.warning(f"{os.path.basename(p)}: sem MS2Quantity (run-wise nem wide).")

    long_ms2 = (pd.concat(rows, ignore_index=True)
                if rows else pd.DataFrame(columns=["ProteinID","Enzyme","Condition","SampleCol","Intensity"]))
    if outdir:
        pd.DataFrame(dbg).to_csv(os.path.join(outdir,"debug","ms2_loader_debug.csv"), index=False)
        long_ms2.to_csv(os.path.join(outdir,"debug","long_ms2_table.tsv"), sep="\t", index=False)
    return sanitize_long_table(long_ms2)

# ============================
# Mineração de Isoformas (curadoria)
# ============================
STRIP_PREFIXES = ("CON__", "REV__", "DECOY_")

def strip_prefixes(pid: str) -> str:
    for pre in STRIP_PREFIXES:
        if pid.startswith(pre): return pid[len(pre):]
    return pid

def build_custom_regex(prefixes: List[str]) -> Optional[re.Pattern]:
    prefixes = [p for p in prefixes if p]
    return re.compile(rf"^(?:{'|'.join(map(re.escape, prefixes))})", re.I) if prefixes else None

def is_custom_id(pid: str, custom_regex: Optional[re.Pattern]) -> bool:
    return bool(custom_regex and custom_regex.match(pid))

def is_isoform(pid: str) -> bool:
    return "-" in strip_prefixes(pid)

def base_of_isoform(pid: str) -> str:
    s = strip_prefixes(pid)
    return s.split("-",1)[0] if "-" in s else s

def process_proteingroup_ids(protein_id_field: str, custom_regex: Optional[re.Pattern]) -> Dict:
    tokens = [t.strip() for t in str(protein_id_field).split(";") if t.strip()]
    custom_ids, uniprot_isoforms, other_canonicals = [], [], []
    if len(tokens)==1 and is_isoform(tokens[0]) and strip_prefixes(tokens[0]).count("-")==1:
        return dict(custom_ids=[], uniprot_isoforms=[tokens[0]], other_canonicals=[],
                    has_custom=False, only_noncanonical_isoforms=True,
                    multiple_isoforms_same_base=False, valid_isoforms=[tokens[0]],
                    keep_group=True, group_type="isoform_only_single")
    for t in tokens:
        if is_custom_id(t, custom_regex):
            custom_ids.append(t)
        elif is_isoform(t):
            uniprot_isoforms.append(t)
        else:
            other_canonicals.append(t)
    custom_ids = sorted(set(custom_ids))
    uniprot_isoforms = sorted(set(uniprot_isoforms))
    other_canonicals = sorted(set(other_canonicals))
    has_custom = len(custom_ids)>0
    only_iso   = (len(uniprot_isoforms)>0) and (len(other_canonicals)==0)
    # múltiplas isoformas com mesmo "base"?
    base_map: Dict[str, List[str]] = {}
    for iso in uniprot_isoforms:
        base_map.setdefault(base_of_isoform(iso), []).append(iso)
    multiple_same_base = any(len(v)>1 for v in base_map.values())
    valid_isoforms: List[str] = [] if multiple_same_base else uniprot_isoforms.copy()
    keep_custom = has_custom and only_iso and (not multiple_same_base)
    keep_iso    = (not has_custom) and only_iso and (not multiple_same_base)
    group_type  = "custom" if keep_custom else ("isoform_only" if keep_iso else "")
    keep_group  = keep_custom or keep_iso
    return dict(custom_ids=custom_ids, uniprot_isoforms=uniprot_isoforms,
                other_canonicals=other_canonicals, has_custom=has_custom,
                only_noncanonical_isoforms=only_iso,
                multiple_isoforms_same_base=multiple_same_base,
                valid_isoforms=valid_isoforms, keep_group=keep_group, group_type=group_type)

def annotate_isoforms_table(paths: List[str], protein_accessions_col: str,
                            custom_prefixes: List[str], outdir: str) -> pd.DataFrame:
    creg = build_custom_regex(custom_prefixes)
    annots = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
        if protein_accessions_col not in df.columns:
            raise SystemExit(f"{os.path.basename(p)}: coluna '{protein_accessions_col}' não encontrada.")
        recs = df[protein_accessions_col].astype(str).apply(lambda s: process_proteingroup_ids(s, creg))
        out = df.copy()
        for k in ["custom_ids","uniprot_isoforms","other_canonicals","has_custom",
                  "only_noncanonical_isoforms","multiple_isoforms_same_base","valid_isoforms",
                  "keep_group","group_type"]:
            out[k] = recs.apply(lambda d: d[k])
        out["__Enzyme__"] = infer_enzyme_from_path(p) or "UNK"
        annots.append(out)
    aa = pd.concat(annots, ignore_index=True)
    aa.to_csv(os.path.join(outdir,"isoforms","isoforms_annotations.tsv"), sep="\t", index=False)
    # resumo
    total = len(aa)
    kept  = int(aa["keep_group"].sum())
    custom_count = int((aa["group_type"] == "custom").sum())
    isoform_only_count = int((aa["group_type"] == "isoform_only").sum())
    isoform_only_single_count = int((aa["group_type"] == "isoform_only_single").sum())
    with open(os.path.join(outdir,"isoforms","isoforms_summary.txt"),"w") as f:
        f.write(f"Total protein groups: {total}\n")
        f.write(f"Groups kept (keep_group=True): {kept}\n")
        f.write(f"  custom: {custom_count}\n")
        f.write(f"  isoform_only: {isoform_only_count}\n")
        f.write(f"  isoform_only_single: {isoform_only_single_count}\n")
        f.write(f"Fraction kept: {kept/total if total else 0:.3f}\n")
    return aa

def explode_isoforms(ann_df: pd.DataFrame) -> pd.DataFrame:
    keep_df = ann_df[ann_df["keep_group"]].copy()
    e = keep_df.explode("valid_isoforms", ignore_index=True)
    return e.rename(columns={"valid_isoforms":"isoform_kept"})

# ============================
# FCs, TOPs, Ranks, UpSets, Coverage
# ============================
def summarize_by_protein(long: pd.DataFrame, agg_func="sum", pseudocount=1.0) -> Tuple[pd.DataFrame,pd.DataFrame]:
    if long is None or long.empty:
        return (pd.DataFrame(), pd.DataFrame())
    long = sanitize_long_table(long)
    if agg_func=="median":
        wide = (long.groupby(["ProteinID","Enzyme","Condition"], dropna=False)["Intensity"]
                    .median().unstack("Condition").fillna(0.0))
    else:
        wide = (long.groupby(["ProteinID","Enzyme","Condition"], dropna=False)["Intensity"]
                    .sum().unstack("Condition").fillna(0.0))
    for c in ["CTRL","N+","N0"]:
        if c not in wide.columns: wide[c]=0.0
    log2fc = lambda a,b: np.log2((a+pseudocount)/(b+pseudocount))
    wide["log2FC_N+_vs_CTRL"] = log2fc(wide["N+"], wide["CTRL"])
    wide["log2FC_N0_vs_N+"]   = log2fc(wide["N0"], wide["N+"])
    wide["log2FC_N0_vs_CTRL"] = log2fc(wide["N0"], wide["CTRL"])
    per_enzyme = wide.reset_index()
    metrics = ["log2FC_N+_vs_CTRL","log2FC_N0_vs_N+","log2FC_N0_vs_CTRL"]
    g = (per_enzyme.set_index(["ProteinID","Enzyme"])[metrics]
                    .groupby(level=0).median().reset_index())
    g["Enzyme"]="GLOBAL"
    per_enzyme_global = pd.concat([per_enzyme,g], ignore_index=True)
    return per_enzyme, per_enzyme_global

def export_top_markers_global(per_enzyme_global: pd.DataFrame, outdir: str, top_k=50):
    cond_to_col = {"CTRL":"log2FC_N+_vs_CTRL","N+":"log2FC_N0_vs_N+","N0":"log2FC_N0_vs_CTRL"}
    cond_sign   = {"CTRL":-1,"N+":-1,"N0":1}
    g = per_enzyme_global[per_enzyme_global["Enzyme"]=="GLOBAL"].copy()
    folder = os.path.join(outdir,"rank_tables"); os.makedirs(folder, exist_ok=True)
    for cond,col in cond_to_col.items():
        if col not in g.columns: continue
        w = g[["ProteinID",col]].dropna().copy()
        if w.empty: continue
        w["effect_for_cond"] = cond_sign[cond]*w[col]
        up  = w.sort_values("effect_for_cond", ascending=False).head(top_k).assign(Regulation="Up", Condition=cond)
        down= w.sort_values("effect_for_cond", ascending=True ).head(top_k).assign(Regulation="Down", Condition=cond)
        out = pd.concat([up,down], ignore_index=True).rename(columns={col:"log2FC_pair"})
        out.to_csv(os.path.join(folder,f"TOP_GLOBAL_{cond}.csv"), index=False)

def ranks_ms2_pair(long: pd.DataFrame, outdir: str):
    plots = os.path.join(outdir,"protein_rank_plots"); tabs=os.path.join(outdir,"rank_tables")
    os.makedirs(plots, exist_ok=True); os.makedirs(tabs, exist_ok=True)
    agg = (long.groupby(["Enzyme","Condition","ProteinID"])["Intensity"].sum().reset_index())
    agg["log10_Abundance"] = np.log10(agg["Intensity"]+1.0)
    for (enz,cond), g in agg.groupby(["Enzyme","Condition"]):
        w = g.sort_values("log10_Abundance", ascending=False).reset_index(drop=True)
        w["Rank"]=np.arange(1,len(w)+1)
        w.to_csv(os.path.join(tabs,f"MS2_rank_{enz}_{cond}.csv"), index=False)
        plt.figure(figsize=(9,6))
        plt.scatter(w["Rank"], w["log10_Abundance"], s=6, color=col_enzyme(enz))
        plt.xlabel("Rank dentro da enzima (1 = mais abundante)")
        plt.ylabel("log10(MS2 abundance + 1)")
        plt.title(f"MS2 Protein Rank — {enz} | {cond}")
        y0,y1=safe_minmax(w["log10_Abundance"],0.5); plt.ylim(y0,y1)
        plt.tight_layout(); plt.savefig(os.path.join(plots,f"MS2_rank_{enz}_{cond}.svg"), dpi=300); plt.close()

def ranks_ms2_combined(long: pd.DataFrame, outdir: str):
    plots = os.path.join(outdir,"protein_rank_plots"); tabs=os.path.join(outdir,"rank_tables")
    os.makedirs(plots, exist_ok=True); os.makedirs(tabs, exist_ok=True)
    agg = (long.groupby(["Enzyme","Condition","ProteinID"])["Intensity"].sum().reset_index())
    agg["log10_Abundance"]=np.log10(agg["Intensity"]+1.0)
    for cond, sub in agg.groupby("Condition"):
        dfc=sub.sort_values(["Enzyme","log10_Abundance"],ascending=[True,False]).copy()
        dfc["Rank"]=dfc.groupby("Enzyme")["log10_Abundance"].rank(method="first",ascending=False).astype(int)
        dfc.to_csv(os.path.join(tabs,f"MS2_rank_combined_{cond}.csv"), index=False)
        plt.figure(figsize=(10,7))
        for enz,g in dfc.groupby("Enzyme"):
            plt.scatter(g["Rank"], g["log10_Abundance"], s=7, label=str(enz), color=col_enzyme(enz))
        plt.xlabel("Rank dentro da enzima (1 = mais abundante)")
        plt.ylabel("log10(MS2 abundance + 1)")
        plt.title(f"MS2 Protein Rank — TODAS ENZIMAS | Condição: {cond}")
        plt.legend(title="Enzyme", markerscale=2, fontsize=8)
        y0,y1=safe_minmax(dfc["log10_Abundance"],0.5); plt.ylim(y0,y1)
        plt.tight_layout(); plt.savefig(os.path.join(plots,f"MS2_rank_combined_{cond}.svg"), dpi=300); plt.close()

def ranks_ms2_facet(long: pd.DataFrame, outdir: str):
    plots = os.path.join(outdir,"protein_rank_plots")
    agg = (long.groupby(["Enzyme","Condition","ProteinID"])["Intensity"].sum().reset_index())
    agg["log10_Abundance"]=np.log10(agg["Intensity"]+1.0)
    conds=[c for c in ["CTRL","N0","N+"] if c in agg["Condition"].unique()]
    if not conds: return
    y0,y1=safe_minmax(agg["log10_Abundance"],0.5)
    fig,axes=plt.subplots(1,len(conds),figsize=(6*len(conds),6),sharey=True)
    if len(conds)==1: axes=[axes]
    for ax,cond in zip(axes,conds):
        sub=agg[agg["Condition"]==cond].copy()
        sub=sub.sort_values(["Enzyme","log10_Abundance"],ascending=[True,False])
        sub["Rank"]=sub.groupby("Enzyme")["log10_Abundance"].rank(method="first",ascending=False).astype(int)
        for enz,g in sub.groupby("Enzyme"):
            ax.scatter(g["Rank"], g["log10_Abundance"], s=7, label=str(enz), color=col_enzyme(enz))
        ax.set_title(f"Condição: {cond}"); ax.set_ylim(y0,y1)
    axes[0].set_ylabel("log10(MS2 abundance + 1)")
    h,l=axes[-1].get_legend_handles_labels()
    if h: fig.legend(h,l,title="Enzyme",loc="upper center", ncol=min(5,len(l)))
    fig.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(os.path.join(plots,"MS2_rank_facet_all_conditions.svg"), dpi=300); plt.close()

def upset_series(sets_dict: Dict[str,set], index_order: Optional[List[str]]=None):
    keys = index_order if index_order else sorted(sets_dict.keys())
    all_ids = set().union(*sets_dict.values()) if sets_dict else set()
    memberships=[]
    for pid in all_ids:
        present = tuple(sorted([k for k in keys if pid in sets_dict[k]]))
        memberships.append(present)
    counts = Counter(memberships)
    def tup2bool(t): return tuple(k in t for k in keys)
    mi = pd.MultiIndex.from_tuples([tup2bool(t) for t in counts], names=keys)
    ser = pd.Series(list(counts.values()), index=mi)
    return ser, keys, counts

def export_intersections(counts, sets_dict, out_csv, combo_label):
    rows=[]
    for combo in sorted(counts.keys(), key=lambda t: (-counts[t], t)):
        if not combo: continue
        inter = set.intersection(*[sets_dict[k] for k in combo])
        rows.append({f"{combo_label} Combination":", ".join(combo),
                     "Protein Count": len(inter),
                     "Proteins": "; ".join(sorted(map(str, inter)))})
    pd.DataFrame(rows).to_csv(out_csv, index=False)

def _harmonize_upset_label_colors(ax_list: List[plt.Axes], keys: List[str], palette: Dict[str, str]):
    keyset = set(keys)
    for ax in ax_list:
        try:
            for tick in ax.get_yticklabels():
                t = tick.get_text()
                if t in keyset:
                    tick.set_color(palette.get(t, "#000000"))
        except Exception:
            pass
        try:
            for tick in ax.get_xticklabels():
                t = tick.get_text()
                if t in keyset:
                    tick.set_color(palette.get(t, "#000000"))
        except Exception:
            pass

def upset_by_enzyme(long: pd.DataFrame, outdir: str, presence_thr: float = 0.0, view_tag: str = ""):
    sets = {}
    for enz, sub in long.groupby("Enzyme"):
        keep = (sub["Intensity"] > presence_thr)
        sets[str(enz)] = set(sub.loc[keep,"ProteinID"].unique())
    series, order, counts = upset_series(sets)
    plt.figure(figsize=(10,6))
    u = UpSet(series, show_counts=True, sort_by='degree'); u.plot()
    fig = plt.gcf()
    _harmonize_upset_label_colors(fig.axes, list(sets.keys()), ENZYME_COLORS)
    p = os.path.join(outdir,"upset",f"upset_proteins_by_enzyme{view_tag}.svg")
    plt.savefig(p,dpi=300,bbox_inches="tight"); plt.close()
    export_intersections(counts, sets, os.path.join(outdir,"upset",f"proteins_intersections_by_enzyme{view_tag}.csv"), "Enzyme")
    return p

def upset_by_condition(long: pd.DataFrame, outdir: str, presence_thr: float = 0.0, view_tag: str = ""):
    sets = {}
    for cond, sub in long.groupby("Condition"):
        keep = (sub["Intensity"] > presence_thr)
        sets[str(cond)] = set(sub.loc[keep,"ProteinID"].unique())
    series, order, counts = upset_series(sets, index_order=["CTRL","N0","N+"] if all(k in sets for k in ["CTRL","N0","N+"]) else None)
    plt.figure(figsize=(8,6))
    u = UpSet(series, show_counts=True, sort_by='degree'); u.plot()
    fig = plt.gcf()
    _harmonize_upset_label_colors(fig.axes, list(sets.keys()), COND_COLORS)
    p = os.path.join(outdir,"upset",f"upset_proteins_by_condition{view_tag}.svg")
    plt.savefig(p,dpi=300,bbox_inches="tight"); plt.close()
    export_intersections(counts, sets, os.path.join(outdir,"upset",f"proteins_intersections_by_condition{view_tag}.csv"), "Condition")
    return p

# -------- Coverage ----------
def parse_row_coverage(val) -> float:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.nan
    parts = [p.strip().replace("%", "") for p in str(val).split(";") if p.strip()]
    nums = []
    for p in parts:
        try: nums.append(float(p))
        except: pass
    return float(np.median(nums)) if nums else np.nan

def _collect_coverage_rows_from_file(path: str,
                                     protein_id_col: str = "PG.ProteinAccessions") -> pd.DataFrame:
    enz = infer_enzyme_from_path(path) or "UNK"
    df  = pd.read_csv(path, sep="\t", low_memory=False)
    cov_col = "PG.Coverage (Global)" if "PG.Coverage (Global)" in df.columns else ("PG.Coverage" if "PG.Coverage" in df.columns else None)
    if cov_col is None:
        return pd.DataFrame(columns=["ProteinID","Enzyme","Condition","CoveragePct"])
    if protein_id_col not in df.columns:
        for c in ["PG.ProteinAccessions","PG.Genes","PG.ProteinGroups","PG.ProteinNames"]:
            if c in df.columns: protein_id_col=c; break
    pid = df[protein_id_col].astype(str).str.split(";").str[0].str.strip()
    cov = df[cov_col].apply(parse_row_coverage)
    cond = None
    if "R.Condition" in df.columns:
        cond = df["R.Condition"].map(RUNWISE_COND_MAP)
    out = pd.DataFrame({
        "ProteinID": pid, "Enzyme": enz,
        "Condition": cond if cond is not None else pd.Series([None]*len(df)),
        "CoveragePct": cov
    })
    return out

def build_coverage_panel(cov_dir: str, view_tag: str = ""):
    """
    Monta um painel tipo figura de paper com:
      A) heatmap de coverage (enzima × condição)
      B) boxplot coverage global por condição
      C) boxplot nº de proteínas detectadas por condição
      D) boxplot nº de peptídeos detectados por condição

    Usa os PNGs já gerados em cov_dir. Ignora silenciosamente se algum faltar.
    Salva:
      coverage_panel_figure{view_tag}.png
    """
    # caminhos esperados
    box_cov_png     = os.path.join(cov_dir, f"coverage_boxplot_global_condition{view_tag}.svg")
    prot_count_png  = os.path.join(cov_dir, f"protein_counts_boxplot{view_tag}.svg")
    pep_count_png   = os.path.join(cov_dir, f"peptide_counts_boxplot{view_tag}.svg")

    # verifica o que existe
    paths = {
        "A": prot_count_png,
        "B": pep_count_png,
    }
    existing = {k: p for k, p in paths.items() if os.path.isfile(p)}
    if not existing:
        logging.warning(f"Nenhum dos PNGs de coverage encontrados em {cov_dir}; painel não será gerado.")
        return None

    # cria figura 2x2
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, wspace=0.25, hspace=0.25)

    panel_labels = {
        "A": "(C) Protein counts by Condition",
        "B": "(D) Peptide counts by Condition",
    }

    # helper pra desenhar cada subfig
    def add_panel(label, path, grid_spec, ax_idx):
        if label not in existing:
            return None
        img = mpimg.imread(existing[label])
        ax = fig.add_subplot(grid_spec[ax_idx])
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(panel_labels[label], loc="left", fontsize=11, fontweight="bold")
        return ax

    # layout:
    #  A | B
    #  C | D
    add_panel("A", existing.get("A"), gs, (0, 0))
    add_panel("B", existing.get("B"), gs, (0, 1))


    fig.suptitle(f"Coverage & Detection Summary {view_tag}".strip(), fontsize=13, fontweight="bold")
    panel_path = os.path.join(cov_dir, f"coverage_panel_figure{view_tag}.png")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(panel_path, dpi=300)
    plt.close(fig)
    logging.info(f"Painel de coverage salvo em: {panel_path}")
    return panel_path

def coverage_analysis_unified(paths: List[str], long: pd.DataFrame, outdir: str,
                              presence_thr: float = 0.0, protein_id_col: str = "PG.ProteinAccessions",
                              view_dir: Optional[str] = None, view_tag: str = ""):
    """
    Consolida e plota métricas de coverage e contagens, por VISÃO (GLOBAL / NO_ISO / ISOFORMS_ONLY).

    Gera:
      - coverage_detail{view_tag}.tsv
      - coverage_matrix_by_enzyme_condition{view_tag}.tsv + heatmap png
      - coverage_matrix_global_condition{view_tag}.tsv
      - protein_coverage_by_condition{view_tag}.tsv
      - boxplots de coverage (enzima×condição, global por condição)
      - peptide_coverage_detail{view_tag}.tsv
      - boxplot comparativo coverage Proteína vs Peptídeo
      - protein_counts_table{view_tag}.tsv + boxplot
      - peptide_counts_table{view_tag}.tsv + boxplot

    Retorna:
      (cov_detail_path, mat_ec_path, mat_ec_png, mat_cond_path)
    """
    view_root = view_dir or outdir
    cov_dir = os.path.join(view_root, "coverage")
    os.makedirs(cov_dir, exist_ok=True)

    # ======================================================
    # 1) Coverage em nível PROTEÍNA (PG.Coverage)
    # ======================================================
    cov_rows = [_collect_coverage_rows_from_file(p, protein_id_col=protein_id_col) for p in paths]
    cov = pd.concat(cov_rows, ignore_index=True) if cov_rows else pd.DataFrame(
        columns=["ProteinID", "Enzyme", "Condition", "CoveragePct"]
    )

    cov_detail_path = os.path.join(cov_dir, f"coverage_detail{view_tag}.tsv")

    if cov.empty:
        logging.warning("Sem coverage disponível para consolidar (nível proteína).")
        pd.DataFrame(columns=["ProteinID", "Enzyme", "Condition", "CoveragePct"]).to_csv(
            cov_detail_path, sep="\t", index=False
        )
        return None, None, None, None

    cov["ProteinID"] = cov["ProteinID"].astype(str)
    cov.to_csv(cov_detail_path, sep="\t", index=False)

    # Presença em MS2 na VISÃO (GLOBAL / NO_ISO / ISOFORMS_ONLY)
    pres = (
        long[long["Intensity"] > presence_thr]
        .groupby(["Enzyme", "Condition", "ProteinID"])
        .size()
        .reset_index(name="n_present")
    )

    # Imputar Condition quando não veio do report
    need_impute = cov["Condition"].isna()
    if need_impute.any():
        cov_na = cov[need_impute].drop(columns=["Condition"]).drop_duplicates()
        present_any = pres[["Enzyme", "Condition", "ProteinID"]].drop_duplicates()
        cov_imputed = cov_na.merge(present_any, on=["Enzyme", "ProteinID"], how="inner")
        cov = pd.concat([cov[~need_impute], cov_imputed], ignore_index=True)

    cov_pres = cov.merge(
        pres[["Enzyme", "Condition", "ProteinID"]],
        on=["Enzyme", "Condition", "ProteinID"],
        how="inner",
    )

    if cov_pres.empty:
        logging.warning("Sem valores válidos de coverage (nível proteína) — pulando plots.")
        return cov_detail_path, None, None, None

    # ======================================================
    # 2) Matriz por Enzyme × Condition + HEATMAP
    # ======================================================
    mat_ec = (
        cov_pres.groupby(["Enzyme", "Condition"])["CoveragePct"]
        .median()
        .unstack("Condition")
    )

    cols_order = [c for c in ["CTRL", "N0", "N+"] if c in mat_ec.columns]
    mat_ec = mat_ec[cols_order] if cols_order else mat_ec
    mat_ec["Total"] = mat_ec.median(axis=1, skipna=True)
    total_row = mat_ec.median(axis=0, skipna=True).to_frame().T
    total_row.index = ["Total"]
    mat_ec_full = pd.concat([mat_ec, total_row], axis=0)

    mat_ec_path = os.path.join(cov_dir, f"coverage_matrix_by_enzyme_condition{view_tag}.tsv")
    mat_ec_full.to_csv(mat_ec_path, sep="\t")

    # Heatmap
    fig, ax = plt.subplots(figsize=(8, 5 + 0.4 * len(mat_ec_full.index)))
    data = mat_ec_full.values.astype(float)
    im = ax.imshow(data, aspect="auto")
    ax.set_xticks(range(len(mat_ec_full.columns)))
    ax.set_xticklabels(mat_ec_full.columns)
    ax.set_yticks(range(len(mat_ec_full.index)))
    ax.set_yticklabels(mat_ec_full.index)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Median sequence coverage (%)")
    ax.set_title(
        f"Median sequence coverage by Enzyme × Condition (MS2-filtered) {view_tag}"
    )
    plt.tight_layout()
    mat_ec_png = os.path.join(cov_dir, f"coverage_heatmap_by_enzyme_condition{view_tag}.svg")
    plt.savefig(mat_ec_png, dpi=300)
    plt.close()

    # ======================================================
    # 3) Resumo global por condição (proteoma)
    # ======================================================
    mat_cond = (
        cov_pres.groupby(["Condition"])["CoveragePct"]
        .median()
        .to_frame(name="MedianCoveragePct")
        .T
    )
    cols_order = [c for c in ["CTRL", "N0", "N+"] if c in mat_cond.columns]
    mat_cond = mat_cond[cols_order] if cols_order else mat_cond
    mat_cond_path = os.path.join(cov_dir, f"coverage_matrix_global_condition{view_tag}.tsv")
    mat_cond.to_csv(mat_cond_path, sep="\t")

    # ======================================================
    # 4) Coverage por proteína × condição (para downstream)
    # ======================================================
    prot_cond = (
        cov_pres.groupby(["ProteinID", "Condition"])["CoveragePct"]
        .median()
        .unstack("Condition")
    )
    prot_cond = prot_cond[[c for c in ["CTRL", "N0", "N+"] if c in prot_cond.columns]]
    prot_cond_path = os.path.join(cov_dir, f"protein_coverage_by_condition{view_tag}.tsv")
    prot_cond.to_csv(prot_cond_path, sep="\t")

    # ======================================================
    # 5) BOXPLOTS de coverage (proteína)
    # ======================================================
    try:
        df_cov_plot = cov_pres.dropna(subset=["CoveragePct", "Condition"]).copy()
        order_cov = [c for c in ["CTRL", "N0", "N+"] if c in df_cov_plot["Condition"].unique()]
        palette_enzyme = {e: col_enzyme(e) for e in df_cov_plot["Enzyme"].unique()}

        # 5.1) Boxplot coverage por enzima × condição
        plt.figure(figsize=(10, 6))
        sns.boxplot(
            data=df_cov_plot,
            x="Condition",
            y="CoveragePct",
            hue="Enzyme",
            order=order_cov if order_cov else None,
            palette=palette_enzyme,
            showfliers=False,
        )
        plt.ylabel("Sequence coverage (%)")
        plt.title(f"Sequence coverage by Condition × Enzyme {view_tag}")
        plt.tight_layout()
        box_cov_ec_png = os.path.join(
            cov_dir, f"coverage_boxplot_enzyme_condition{view_tag}.svg"
        )
        plt.savefig(box_cov_ec_png, dpi=300)
        plt.close()

        # 5.2) Boxplot coverage global por condição (proteoma)
        plt.figure(figsize=(6, 6))
        sns.boxplot(
            data=df_cov_plot,
            x="Condition",
            y="CoveragePct",
            order=order_cov if order_cov else None,
            showfliers=False,
            palette=COND_COLORS,
        )
        plt.ylabel("Sequence coverage (%)")
        plt.title(f"Global sequence coverage by Condition (proteome) {view_tag}")
        plt.tight_layout()
        box_cov_cond_png = os.path.join(
            cov_dir, f"coverage_boxplot_global_condition{view_tag}.svg"
        )
        plt.savefig(box_cov_cond_png, dpi=300)
        plt.close()
    except Exception as e:
        logging.warning(f"Falha ao gerar boxplots de coverage (proteína): {e}")

    # ======================================================
    # 6) Coverage em nível PEPTÍDEO + boxplot comparativo
    # ======================================================
    try:
        pep_cov_rows = []
        for p in paths:
            enz = infer_enzyme_from_path(p) or "UNK"
            dfp = pd.read_csv(p, sep="\t", low_memory=False)

            # Coluna de coverage em nível peptídeo
            if "PEP.Coverage (Global)" in dfp.columns:
                pep_cov_col = "PEP.Coverage (Global)"
            elif "PEP.Coverage" in dfp.columns:
                pep_cov_col = "PEP.Coverage"
            else:
                continue  # sem coverage de peptídeo nesse report

            # ID do peptídeo
            pep_id_col = None
            for cand in ["PEP.StrippedSequence", "PEP.PeptideSequence"]:
                if cand in dfp.columns:
                    pep_id_col = cand
                    break
            if pep_id_col is None:
                continue

            pep_cov_vals = dfp[pep_cov_col].apply(parse_row_coverage)

            if "R.Condition" in dfp.columns:
                cond = dfp["R.Condition"].map(RUNWISE_COND_MAP)
            else:
                cond = pd.Series([None] * len(dfp))

            pep_cov_rows.append(pd.DataFrame({
                "PeptideID": dfp[pep_id_col].astype(str),
                "Enzyme": enz,
                "Condition": cond,
                "CoveragePct": pep_cov_vals
            }))

        if pep_cov_rows:
            pep_cov = pd.concat(pep_cov_rows, ignore_index=True)
            pep_cov = pep_cov.dropna(subset=["CoveragePct", "Condition"])
        else:
            pep_cov = pd.DataFrame(columns=["PeptideID", "Enzyme", "Condition", "CoveragePct"])

        pep_detail_path = os.path.join(cov_dir, f"peptide_coverage_detail{view_tag}.tsv")
        pep_cov.to_csv(pep_detail_path, sep="\t", index=False)

        # Boxplot comparativo Proteína vs Peptídeo (coverage)
        if not pep_cov.empty:
            df_prot_cov = df_cov_plot[["Condition", "CoveragePct"]].copy()
            df_prot_cov["Level"] = "Protein"

            df_pep_cov = pep_cov[["Condition", "CoveragePct"]].copy()
            df_pep_cov["Level"] = "Peptide"

            df_cmp_cov = pd.concat([df_prot_cov, df_pep_cov], ignore_index=True)
            order_cmp = [c for c in ["CTRL", "N0", "N+"] if c in df_cmp_cov["Condition"].unique()]

            plt.figure(figsize=(8, 6))
            sns.boxplot(
                data=df_cmp_cov,
                x="Condition",
                y="CoveragePct",
                hue="Level",
                order=order_cmp if order_cmp else None,
                showfliers=False,
            )
            plt.ylabel("Sequence coverage (%)")
            plt.title(f"Global coverage — Proteome vs Peptides (all enzymes) {view_tag}")
            plt.tight_layout()
            cmp_cov_png = os.path.join(
                cov_dir, f"coverage_boxplot_global_protein_vs_peptide{view_tag}.svg"
            )
            plt.savefig(cmp_cov_png, dpi=300)
            plt.close()
        else:
            logging.warning("Nenhum coverage em nível de peptídeo encontrado; pulando boxplot comparativo.")
    except Exception as e:
        logging.warning(f"Falha ao processar coverage em nível de peptídeo / boxplot comparativo: {e}")

    # ======================================================
    # 7) BOXPLOT — Quantidade de PROTEÍNAS detectadas por condição
    # ======================================================
    try:
        df_counts_prot = (
            long[long["Intensity"] > presence_thr]
            .groupby(["Condition", "SampleCol"])["ProteinID"]
            .nunique()
            .reset_index(name="ProteinCount")
        )

        if not df_counts_prot.empty:
            order_cnt = [c for c in ["CTRL", "N0", "N+"] if c in df_counts_prot["Condition"].unique()]

            plt.figure(figsize=(8, 6))
            sns.boxplot(
                data=df_counts_prot,
                x="Condition",
                y="ProteinCount",
                order=order_cnt if order_cnt else None,
                showfliers=False,
                palette=COND_COLORS,
            )
            plt.ylabel("Número de proteínas detectadas")
            plt.title(f"Proteínas detectadas por condição (todas enzimas) {view_tag}")
            plt.tight_layout()

            proteins_box_png = os.path.join(
                cov_dir, f"protein_counts_boxplot{view_tag}.svg"
            )
            plt.savefig(proteins_box_png, dpi=300)
            plt.close()

            proteins_box_tsv = os.path.join(
                cov_dir, f"protein_counts_table{view_tag}.tsv"
            )
            df_counts_prot.to_csv(proteins_box_tsv, sep="\t", index=False)
        else:
            logging.warning("Nenhuma proteína detectada para gerar boxplot de contagem.")
    except Exception as e:
        logging.warning(f"Falha ao gerar boxplot de quantidades de proteínas: {e}")

    # ======================================================
    # 8) BOXPLOT — Quantidade de PEPTÍDEOS detectados por condição
    # ======================================================
    try:
        pep_rows_counts = []
        for p in paths:
            enz = infer_enzyme_from_path(p) or "UNK"
            dfp = pd.read_csv(p, sep="\t", low_memory=False)

            has_run_pep = "PEP.MS2Quantity" in dfp.columns
            wide_pep = [c for c in dfp.columns if c.endswith(".PEP.MS2Quantity")]

            # ID de peptídeo
            pep_id_col = None
            for cand in ["PEP.StrippedSequence", "PEP.PeptideSequence"]:
                if cand in dfp.columns:
                    pep_id_col = cand
                    break
            if pep_id_col is None:
                continue

            if has_run_pep:
                intens = pd.to_numeric(dfp["PEP.MS2Quantity"], errors="coerce").fillna(0.0)
                cond_series = (
                    dfp["R.Condition"].map(RUNWISE_COND_MAP)
                    if "R.Condition" in dfp.columns
                    else pd.Series(["UNK"] * len(dfp))
                )
                sample_series = (
                    dfp["R.FileName"]
                    if "R.FileName" in dfp.columns
                    else pd.Series([f"{os.path.basename(p)}_PEP"] * len(dfp))
                )

                pep_rows_counts.append(
                    pd.DataFrame({
                        "PeptideID": dfp[pep_id_col].astype(str),
                        "Enzyme": enz,
                        "Condition": cond_series.astype(str),
                        "SampleCol": sample_series.astype(str),
                        "Intensity": intens,
                    })
                )
            elif wide_pep:
                for c in wide_pep:
                    intens = pd.to_numeric(dfp[c], errors="coerce").fillna(0.0)
                    sample_name = parse_filename_from_ms2col(c)
                    cond_token = cond_from_token(sample_name) or "UNK"
                    pep_rows_counts.append(
                        pd.DataFrame({
                            "PeptideID": dfp[pep_id_col].astype(str),
                            "Enzyme": enz,
                            "Condition": cond_token,
                            "SampleCol": sample_name,
                            "Intensity": intens,
                        })
                    )

        if pep_rows_counts:
            pep_long = pd.concat(pep_rows_counts, ignore_index=True)
            pep_counts = (
                pep_long[pep_long["Intensity"] > presence_thr]
                .groupby(["Condition", "SampleCol"])["PeptideID"]
                .nunique()
                .reset_index(name="PeptideCount")
            )
        else:
            pep_counts = pd.DataFrame(columns=["Condition", "SampleCol", "PeptideCount"])

        pep_counts_tsv = os.path.join(
            cov_dir, f"peptide_counts_table{view_tag}.tsv"
        )
        pep_counts.to_csv(pep_counts_tsv, sep="\t", index=False)

        if not pep_counts.empty:
            order_cnt_pep = [c for c in ["CTRL", "N0", "N+"] if c in pep_counts["Condition"].unique()]

            plt.figure(figsize=(8, 6))
            sns.boxplot(
                data=pep_counts,
                x="Condition",
                y="PeptideCount",
                order=order_cnt_pep if order_cnt_pep else None,
                showfliers=False,
                palette=COND_COLORS,
            )
            plt.ylabel("Número de peptídeos detectados")
            plt.title(f"Peptídeos detectados por condição (todas enzimas) {view_tag}")
            plt.tight_layout()

            peptides_box_png = os.path.join(
                cov_dir, f"peptide_counts_boxplot{view_tag}.svg"
            )
            plt.savefig(peptides_box_png, dpi=300)
            plt.close()
        else:
            logging.warning("Nenhum peptídeo detectado para gerar boxplot de contagem.")
    except Exception as e:
        logging.warning(f"Falha ao gerar boxplot de quantidades de peptídeos: {e}")

    # ======================================================
    # 9) Painel tipo figura de paper com todos os boxplots + heatmap
    # ======================================================
    try:
        build_coverage_panel(cov_dir, view_tag=view_tag)
    except Exception as e:
        logging.warning(f"Falha ao gerar painel de coverage: {e}")

    return cov_detail_path, mat_ec_path, mat_ec_png, mat_cond_path


    # Retorno original mantido
    return cov_detail_path, mat_ec_path, mat_ec_png, mat_cond_path

# -------- Heatmap FC opcional --------
def heatmap_fc(per_enzyme: pd.DataFrame, outdir: str, value_col: str = "log2FC_N+_vs_CTRL"):
    if per_enzyme.empty or value_col not in per_enzyme.columns:
        return None
    hm = (per_enzyme[["ProteinID","Enzyme",value_col]]
                    .dropna()
                    .groupby(["Enzyme","ProteinID"]).mean(numeric_only=True).reset_index())
    pivot = hm.pivot(index="ProteinID", columns="Enzyme", values=value_col)
    if pivot.empty: return None
    plt.figure(figsize=(12, max(6, len(pivot)//2)))
    sns.heatmap(pivot, cmap="vlag", center=0, cbar_kws={"label": value_col})
    plt.title(f"Heatmap of {value_col} by Protein and Enzyme")
    plt.tight_layout()
    path = os.path.join(outdir, "heatmaps", f"heatmap_{value_col.replace('+','plus')}.svg")
    plt.savefig(path, dpi=300); plt.close()
    return path

# ============================
# Construção das VISÕES (sem recursão)
# ============================
def build_views(long_ms2: pd.DataFrame,
                ann_iso: pd.DataFrame,
                protein_accessions_col: str) -> Dict[str, Dict[str, pd.DataFrame]]:
    v = {}

    # mapa de grupos mantidos e isoforms explodidos
    expl = explode_isoforms(ann_iso)  # contém isoform_kept; apenas keep_group=True
    # ids de grupos mantidos (comparação por primeira acc)
    kept_groups_first_acc = set(
        ann_iso[ann_iso["keep_group"]==True][protein_accessions_col].astype(str).str.split(";").str[0].unique()
    )

    # 1) GLOBAL_ALL
    v["GLOBAL_ALL"] = {"long": long_ms2.copy()}

    # 2) GLOBAL_NO_ISO — exclui todos os grupos mantidos
    long_no_iso = long_ms2.copy()
    mask_exclude = long_no_iso["ProteinID"].astype(str).isin(kept_groups_first_acc)
    long_no_iso = long_no_iso[~mask_exclude].copy()
    v["GLOBAL_NO_ISO"] = {"long": long_no_iso}

    # 3) ISOFORMS_ONLY — mapeia ProteinID (first acc) -> isoform_kept (explode)
    expl2 = expl[[protein_accessions_col,"isoform_kept"]].copy()
    expl2["__first__"] = expl2[protein_accessions_col].astype(str).str.split(";").str[0]
    iso_long = long_ms2.merge(expl2[["__first__","isoform_kept"]],
                              left_on="ProteinID", right_on="__first__", how="inner")
    iso_long = iso_long.drop(columns=["__first__"]).rename(columns={"isoform_kept":"ProteinID"})
    v["ISOFORMS_ONLY"] = {"long": iso_long}

    # sanitize todas
    for k in list(v.keys()):
        v[k]["long"] = sanitize_long_table(v[k]["long"])
    return v

# ============================
# MAIN
# ============================
def main():
    ap = argparse.ArgumentParser(description="Pipeline unificado (GLOBAL vs ISOFORMS) para middle-down.")
    ap.add_argument("--inputs", nargs="+", required=True, help="TSVs dos reports (um por enzima).")
    ap.add_argument("--outdir", default="outputs_WT", help="Diretório de saída.")
    ap.add_argument("--protein-id-col", default=None, help="Coluna de identificador (default: detecta).")
    ap.add_argument("--protein-accessions-col", default="PG.ProteinGroups",
                    help="Coluna com acessões separadas por ';' para mineração de isoformas.")
    ap.add_argument("--prefer", choices=["PG","PEP"], default="PG", help="Preferência MS2 (PG→fallback PEP).")
    ap.add_argument("--norm", choices=["none","median","sum"], default="median", help="Normalização por amostra (wide).")
    ap.add_argument("--agg-func", choices=["sum","median"], default="sum", help="Agregação replicatas Protein×Enzyme×Condition.")
    ap.add_argument("--pseudocount", type=float, default=1.0, help="Pseudocount para FCs.")
    ap.add_argument("--presence-thr", type=float, default=0.0, help="Limiar de presença para UpSet e coverage.")
    ap.add_argument("--custom-prefixes", default="CTRL,NMAIS,NZERO,ML,NO,LOM,TMR",
                    help="Prefixos de bancos custom para a curadoria de isoformas.")
    ap.add_argument("--top-k", type=int, default=50, help="Top K de markers globais por condição.")
    ap.add_argument("--label-k", type=int, default=15, help="Rótulos em plots rank.")
    args = ap.parse_args()

    ensure_root_dirs(args.outdir)

    # 1) Load MS2
    long_ms2 = load_reports_unified(
        args.inputs, protein_id_col=args.protein_id_col,
        prefer=args.prefer, norm=args.norm, outdir=args.outdir
    )
    if long_ms2.empty:
        raise SystemExit("Nenhuma intensidade MS2 encontrada. Verifique os reports/colunas.")

    # 2) Mineração de isoformas
    ann_iso = annotate_isoforms_table(
        args.inputs, protein_accessions_col=args.protein_accessions_col,
        custom_prefixes=[s.strip() for s in args.custom_prefixes.split(",") if s.strip()],
        outdir=args.outdir
    )

    # 3) VISÕES (sem recursão)
    views = build_views(long_ms2, ann_iso, args.protein_accessions_col)

    # 4) Para cada VISÃO: FCs, TOPs, Ranks, UpSets, Coverage
    for view_name, payload in views.items():
        view_dir = os.path.join(args.outdir, f"views/{view_name}")
        long = payload["long"].copy()  # já sanitizado

        # FCs
        per_enzyme, per_enzyme_global = summarize_by_protein(long,agg_func=args.agg_func, pseudocount=args.pseudocount)
        per_enzyme.to_csv(os.path.join(view_dir,"rank_tables","per_enzyme_log2FC.tsv"), sep="\t", index=False)
        per_enzyme_global.to_csv(os.path.join(view_dir,"rank_tables","per_enzyme_and_GLOBAL_log2FC.tsv"), sep="\t", index=False)

        # TOP markers (GLOBAL)
        export_top_markers_global(per_enzyme_global, view_dir, top_k=args.top_k)

        # Ranks
        ranks_ms2_pair(long, view_dir)
        ranks_ms2_combined(long, view_dir)
        ranks_ms2_facet(long, view_dir)

        # UpSets
        upset_by_enzyme(long, view_dir, presence_thr=args.presence_thr, view_tag=f"_{view_name}")
        upset_by_condition(long, view_dir, presence_thr=args.presence_thr, view_tag=f"_{view_name}")

        # Coverage (MS2-filtered)
        coverage_analysis_unified(
            args.inputs, long, outdir=args.outdir, presence_thr=args.presence_thr,
            protein_id_col=args.protein_id_col or "PG.ProteinAccessions",
            view_dir=view_dir, view_tag=f"_{view_name}"
        )

        # Heatmap FC (opcional: N+ vs CTRL na visão)
        heatmap_fc(per_enzyme, view_dir, value_col="log2FC_N+_vs_CTRL")

    print("✅ Concluído. Saídas em:", os.path.abspath(args.outdir))
    print("Views geradas:", ", ".join(views.keys()))
    print("Paleta fixa por enzima aplicada em todos os plots.")

if __name__ == "__main__":
    main()