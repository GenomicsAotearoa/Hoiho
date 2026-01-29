from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable, Set
import bisect, csv, os, pickle, re
from urllib.parse import unquote

try:
    import pandas as pd
except Exception:
    pd = None


@dataclass
class GFFFeature:
    seqid: str
    source: Optional[str]
    type: str
    start: int
    end: int
    score: Optional[float]
    strand: Optional[str]
    phase: Optional[str]
    attributes: Dict[str, str]

    def key_candidates(self) -> List[str]:
        raw_vals = []
        for k in ("ID","Name","gene","gene_id","locus_tag","transcript_id","protein_id","Parent","Alias","Dbxref"):
            v = self.attributes.get(k)
            if not v: continue
            parts = re.split(r"[,\s]+", v)
            for p in parts:
                if not p: continue
                raw_vals.append(p)
                if ":" in p:
                    rhs = p.split(":",1)[1]
                    if rhs: raw_vals.append(rhs)
        def sane(tok): return bool(re.fullmatch(r"[A-Za-z0-9._|\-]+", tok))
        seen, out = set(), []
        for tok in raw_vals:
            if sane(tok) and tok not in seen:
                out.append(tok); seen.add(tok)
        for tok in out:
            if "|" in tok:
                last = tok.split("|")[-1]
                if last and sane(last) and last not in seen:
                    out.append(last); seen.add(last)
        return out


class GFFIntervalIndex:
    def __init__(self, features_by_contig: Dict[str,List[GFFFeature]]):
        self.index = {}
        for contig, feats in features_by_contig.items():
            feats_sorted = sorted(feats, key=lambda f:(f.start,f.end))
            self.index[contig] = (feats_sorted, [f.start for f in feats_sorted])

    def query(self, contig, qstart, qend, types=None, how="overlap"):
        if qstart>qend: qstart,qend=qend,qstart
        feats,starts = self.index.get(contig,([],[]))
        i = bisect.bisect_right(starts, qend)
        candidates = feats[:i]
        out=[]
        for f in candidates:
            if types and f.type not in types: continue
            if how=="overlap" and not(f.end<qstart or f.start>qend): out.append(f)
            elif how=="within" and f.start>=qstart and f.end<=qend: out.append(f)
            elif how=="cover" and f.start<=qstart and f.end>=qend: out.append(f)
        return out


GFF_FEATURES_BY_CONTIG = {}
GFF_IDX = None
EGGNOG_MAP = {}
EGGNOG_ALIAS_TO_KEY = {}

def parse_gff(path, feature_types=None):
    feats_by_contig={}
    with open(path) as fin:
        for line in fin:
            if not line or line.startswith("#"): continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)!=9: continue
            seqid,source,ftype,start,end,score,strand,phase,attrs=parts
            if feature_types and ftype not in feature_types: continue
            try: start,end=int(start),int(end)
            except: continue
            d={}
            if attrs and attrs!=".":
                for kv in re.split(r";\s*",attrs):
                    if not kv: continue
                    if "=" in kv: k,v=kv.split("=",1)
                    else: k,v=kv,""
                    d[k]=unquote(v.strip('"'))
            feats_by_contig.setdefault(seqid,[]).append(
                GFFFeature(seqid, None if source=="." else source, ftype,
                           min(start,end), max(start,end),
                           None if score=="." else float(score),
                           None if strand=="." else strand,
                           None if phase=="." else phase,
                           d))
    return feats_by_contig


import csv, re, os
from typing import Dict, List, Tuple, Optional

# Column aliases we’ll accept across eggNOG-mapper variants
_EGG_ALIASES = {
    "query": ["query","query_name","qseqid","seqid","#query","#query_name","#qseqid","#seqid"],
    "preferred_name": ["preferred_name","predicted_gene_name","name","gene_name"],
    "description": ["description","annot","annotation","product"],
    "go_terms": ["gos","go_terms","go", "go-term", "go-terms", "go_annotations", "go terms", "go_terms(eggnog)"],
}

def _norm(s): 
    return s.strip().lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "")

def _detect_emapper_header(path: str) -> List[str]:
    header = None
    with open(path, "r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith("#"):
                cand = line.lstrip("#").rstrip("\n").split("\t")
                if any(_norm("query") == _norm(x) for x in cand):
                    header = cand
            else:
                break
    if header is None:
        # Fallback: first data row count – create dummy column names
        with open(path, "r", encoding="utf-8") as fin:
            for ln in fin:
                if ln.startswith("#"):
                    continue
                n = len(ln.rstrip("\n").split("\t"))
                header = [f"col{i+1}" for i in range(n)]
                break
    return [h.strip() for h in header] if header else []

def _build_alias_map(header: List[str]) -> Dict[str, str]:
    """
    Map logical keys ('query','preferred_name','description','go_terms') 
    to actual header names present in the file.
    """
    present = { _norm(h): h for h in header }
    out = {}
    for logical, candidates in _EGG_ALIASES.items():
        for c in candidates:
            key = _norm(c)
            if key in present:
                out[logical] = present[key]
                break
    return out

def parse_emapper_annotations(path):
    header=None
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                cand=line.lstrip("#").strip().split("\t")
                if any("query" in x.lower() for x in cand): header=cand
            else: break
    if not header: raise ValueError("No header found")
    norm=[h.strip() for h in header]
    lower_map={h.lower():h for h in norm}
    wanted,alias={},{}
    with open(path) as fin:
        rdr=csv.DictReader((ln for ln in fin if not ln.startswith("#")),
                           fieldnames=norm,delimiter="\t",quoting=csv.QUOTE_NONE)
        for row in rdr:
            row={k:(v.strip() if v else v) for k,v in row.items()}
            q=row.get(lower_map.get("query")) or row.get(lower_map.get("seqid"))
            if not q: continue
            desc=row.get("Description") or row.get("description")
            pref=row.get("Preferred_name") or row.get("preferred_name")
            gos=row.get("GOs") or row.get("go_terms")
            go_terms=[]
            if gos:
                go_terms=re.findall(r"GO:\d{7}",gos) or [p for p in re.split(r"[,\s|;]+",gos) if p.startswith("GO:")]
            wanted[q]={"query":q,"preferred_name":pref,"description":desc,"go_terms":sorted(set(go_terms))}
            alias.setdefault(q,q)
            if "|" in q: alias.setdefault(q.split("|")[-1],q)
    return wanted,alias

def parse_emapper_annotations_verbose(
    path: str, 
    accept_raw_go_strings: bool = True,
    strict_go_regex: bool = True
) -> Tuple[Dict[str, Dict[str, object]], Dict[str, str], Dict[str, object]]:
    """
    Returns:
      - map: query -> {query, preferred_name, description, go_terms(list[str])}
      - alias_to_key: alias -> canonical query
      - diag: diagnostics dict with coverage and column mapping info
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    header = _detect_emapper_header(path)
    if not header:
        raise ValueError("eggNOG: could not determine header; file empty or malformed?")

    alias_map = _build_alias_map(header)
    diag = {
        "path": path,
        "n_rows": 0,
        "n_with_query": 0,
        "n_with_prefname": 0,
        "n_with_description": 0,
        "n_with_go_terms": 0,
        "header": header,
        "alias_map": alias_map,  # logical -> actual column name
        "notes": [],
        "missing_keys": [k for k in ("query","preferred_name","description","go_terms") if k not in alias_map]
    }

    if "query" not in alias_map:
        diag["notes"].append("No recognizable 'query' column; nothing will map.")
    if "go_terms" not in alias_map:
        diag["notes"].append("No recognizable GO terms column; GO annotations will be empty.")

    wanted: Dict[str, Dict[str, object]] = {}
    alias_to_key: Dict[str, str] = {}

    # Stream rows
    with open(path, "r", encoding="utf-8") as fin:
        rdr = csv.DictReader(
            (ln for ln in fin if not ln.startswith("#")),
            fieldnames=header, delimiter="\t", quoting=csv.QUOTE_NONE
        )
        for row in rdr:
            diag["n_rows"] += 1

            q = row.get(alias_map.get("query", ""), "") or ""
            q = q.strip()
            if not q:
                continue
            diag["n_with_query"] += 1

            pref = row.get(alias_map.get("preferred_name", ""), "") or None
            if pref: 
                pref = pref.strip()
                if pref: diag["n_with_prefname"] += 1
            desc = row.get(alias_map.get("description", ""), "") or None
            if desc:
                desc = desc.strip()
                if desc: diag["n_with_description"] += 1

            gos_raw = row.get(alias_map.get("go_terms", ""), "") if "go_terms" in alias_map else ""
            go_terms: List[str] = []
            if gos_raw:
                # First: strict GO id extraction
                ids = re.findall(r"GO:\d{7}", gos_raw)
                if ids:
                    go_terms = sorted(set(ids))
                elif accept_raw_go_strings:
                    # Fall back to splitting and keeping things that look like GO:...
                    parts = re.split(r"[,\s|;]+", gos_raw.strip())
                    go_terms = sorted(set([p for p in parts if p.startswith("GO:")]))
                if go_terms:
                    diag["n_with_go_terms"] += 1
                elif strict_go_regex:
                    diag["notes"].append(
                        f"Row with query '{q}' had a GO column but no canonical GO IDs matched regex."
                    )

            wanted[q] = {
                "query": q,
                "preferred_name": pref if pref else None,
                "description": desc if desc else None,
                "go_terms": go_terms,
            }

            # Conservative aliases: full string and last pipe-chunk
            alias_to_key.setdefault(q, q)
            if "|" in q:
                last = q.split("|")[-1]
                if re.fullmatch(r"[A-Za-z0-9._\-]+", last or ""):
                    alias_to_key.setdefault(last, q)

    # Final warnings
    if diag["n_rows"] == 0:
        diag["notes"].append("No data rows parsed (file might be empty after comments).")
    if diag["n_with_go_terms"] == 0 and "go_terms" in alias_map:
        diag["notes"].append("Parsed zero GO terms. Check GO column content/format.")

    return wanted, alias_to_key, diag



def init_gff_index(gff_path, feature_types=("gene","mRNA")):
    global GFF_FEATURES_BY_CONTIG,GFF_IDX
    GFF_FEATURES_BY_CONTIG=parse_gff(gff_path,feature_types)
    GFF_IDX=GFFIntervalIndex(GFF_FEATURES_BY_CONTIG)
    return GFF_IDX

def init_eggnog_map(path):
    global EGGNOG_MAP,EGGNOG_ALIAS_TO_KEY
    EGGNOG_MAP,EGGNOG_ALIAS_TO_KEY=parse_emapper_annotations(path)
    return EGGNOG_MAP,EGGNOG_ALIAS_TO_KEY

def _lookup_eggnog_for_feature(f):
    for key in f.key_candidates():
        if key in EGGNOG_MAP: return EGGNOG_MAP[key]
        base=EGGNOG_ALIAS_TO_KEY.get(key)
        if base and base in EGGNOG_MAP: return EGGNOG_MAP[base]
    return None

def query_region(contig,start,end,feature_types=("gene",),how="overlap",as_dataframe=True):
    if not GFF_IDX: raise RuntimeError("Init GFF first")
    feats=GFF_IDX.query(contig,start,end,types=set(feature_types),how=how)
    rows=[]
    for f in feats:
        ann=_lookup_eggnog_for_feature(f)
        rows.append({
            "contig":f.seqid,"start":f.start,"end":f.end,"strand":f.strand,"type":f.type,
            "id":f.attributes.get("ID") or f.attributes.get("Name") or "",
            "query_key":ann["query"] if ann else None,
            "preferred_name":ann["preferred_name"] if ann else None,
            "description":ann["description"] if ann else None,
            "go_terms":";".join(ann["go_terms"]) if ann else ""
        })
    if as_dataframe and pd is not None: return pd.DataFrame(rows)
    return rows

if __name__=="__main__":
    import sys
    if len(sys.argv)<3:
        print("Usage: gff_eggnog_util.py <gff3-file> <eggnog-emapper-output>")
        sys.exit(1)
    gff_path = sys.argv[1]
    eggnog_path = sys.argv[2]
    print(f"Loading GFF from {gff_path}...")
    init_gff_index(gff_path)
    print(f"Loading eggNOG annotations from {eggnog_path}...")
    init_eggnog_map(eggnog_path)

    # Get the eggnog mapper info for g1.t1 and try for just 'g1'

    print("Done.")

