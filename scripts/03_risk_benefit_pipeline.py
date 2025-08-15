# Risk-Benefit pipeline snippet
import numpy as np
import pandas as pd
from scipy.stats import beta
import matplotlib.pyplot as plt

# ---------- USER INPUT / ADAPT ----------
# Dataframes:
# norm_expr: genes x normal_samples
# meta_norm: DataFrame with columns ['sample_id', 'tissue']
# tum_expr: genes x tumour_samples
# meta_tum: DataFrame with columns ['sample_id', 'indication']

# cutoffs: dict gene -> cutoff_value (could be global per gene or per-tissue if you adapt)
# Example: cutoffs = {'GENE1': 2.3, 'GENE2': 0.0, ...}
# If you have per-tissue cutoffs, adapt the p_normal calculation accordingly.

# Tissue risk group mapping 
tissue_to_group = {
    # High-risk examples:
    'Adipose Tissue': 'High', 'Blood': 'High', 'Blood Vessel': 'High', 'Soft tissue/Bone': 'High', 'Brain': 'High', 'Esophagus': 'High', 'Eye': 'High', 'Heart': 'High', 'Liver': 'High', 'Lung': 'High', 
    'Lymphatic tissue': 'High', 'Muscle': 'High', 'Nerve': 'High', 'Paraganglia': 'High', 'Small Intestine': 'High', 'White blood cell': 'High', 'Head and Neck region': 'High', 'Bile duct': 'High',
    # Medium-risk examples:
    'Adrenal Gland': 'Medium', 'Bladder': 'Medium', 'Colon': 'Medium', 'Kidney': 'Medium', 'Pancreas': 'Medium', 'Pituitary': 'Medium', 'Stomach': 'Medium', 'Skin': 'Medium', 'Lining of body cavities': 'Medium',
    # Low-risk examples:
    'Breast': 'Low', 'Cervix': 'Low', 'Fallopian Tube': 'Low', 'Ovary': 'Low', 'Salivary Gland': 'Low', 'Prostate': 'Low', 'Rectum': 'Low', 'Spleen': 'Low', 'Thymus': 'Low', 'Thyroid': 'Low',
    'Uterus': 'Low', 'Endometrium': 'Low', 'Vagina': 'Low'
}
group_weights = {'High': 1.0, 'Medium': 0.5, 'Low': 0.1}

# Utility lambda
lam = 3.0

# Hard reject threshold for any high-risk tissue
hard_reject_high_pct = 0.05  # 5%

# Minimum samples to trust a tissue estimate
min_samples_tissue = 8
# ----------------------------------------

def proportion_ci_wilson(k, n, alpha=0.05):
    """Wilson score interval (approx). Returns (p_hat, low, high)."""
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    from math import sqrt
    z = 1.96  # approximate for 95%
    denom = 1 + z**2 / n
    centre = p + z*z/(2*n)
    adj = z*sqrt((p*(1-p) + z*z/(4*n))/n)
    low = (centre - adj) / denom
    high = (centre + adj) / denom
    return (p, max(0, low), min(1, high))

# Precompute small helper maps for meta frames
meta_norm = meta_norm.set_index('sample_id')
meta_tum = meta_tum.set_index('sample_id')

# Map normal samples by tissue
tissue_samples = {}
for sid, row in meta_norm.iterrows():
    tissue = row['tissue']
    tissue_samples.setdefault(tissue, []).append(sid)

# Map tumour samples by indication
indication_samples = {}
for sid, row in meta_tum.iterrows():
    ind = row['indication']
    indication_samples.setdefault(ind, []).append(sid)

genes = list(cutoffs.keys())

# Output records
records = []

for gene in genes:
    cutoff = cutoffs[gene]
    # 1) Per-tissue normal % above cutoff
    tissue_info = {}
    for tissue, sids in tissue_samples.items():
        expr_vals = norm_expr.loc[gene, sids].values
        n = len(expr_vals)
        k = np.sum(expr_vals > cutoff)
        p, low, high = proportion_ci_wilson(k, n) if n>0 else (np.nan, np.nan, np.nan)
        grp = tissue_to_group.get(tissue, 'Low')
        weight = group_weights.get(grp, 0.1)
        tissue_info[tissue] = {'n': n, 'k': int(k), 'p': p, 'ci_low': low, 'ci_high': high,
                               'group': grp, 'group_weight': weight}
    # 2) RiskScore: weighted sum across tissues (group-weighted)
    # Option A: weighted sum of p across tissues
    riskscore = 0.0
    total_w = 0.0
    for t,v in tissue_info.items():
        if np.isnan(v['p']): continue
        w = v['group_weight']
        riskscore += w * v['p']
        total_w += w
    if total_w > 0:
        riskscore = riskscore / total_w  # normalize to keep scale 0..1
    else:
        riskscore = np.nan

    # Max high-risk p (for hard reject)
    max_high_p = 0.0
    for t,v in tissue_info.items():
        if v['group'] == 'High' and not np.isnan(v['p']):
            max_high_p = max(max_high_p, v['p'])

    # 3) Per-indication tumour % above cutoff -> BenefitScore per indication
    benefit = {}
    for ind, sids in indication_samples.items():
        expr_vals = tum_expr.loc[gene, sids].values
        n_t = len(expr_vals)
        k_t = np.sum(expr_vals > cutoff)
        p_t, low_t, high_t = proportion_ci_wilson(k_t, n_t) if n_t>0 else (np.nan, np.nan, np.nan)
        benefit[ind] = {'n': n_t, 'k': int(k_t), 'p': p_t, 'ci_low': low_t, 'ci_high': high_t}

    # 4) Utility per indication
    utility = {}
    for ind,v in benefit.items():
        b = v['p'] if not np.isnan(v['p']) else 0.0
        util = b / (1.0 + lam * (riskscore if not np.isnan(riskscore) else 0.0))
        utility[ind] = util

    # 5) Decision flags: hard reject or flagged for review
    hard_reject = (max_high_p > hard_reject_high_pct)
    # e.g. accept if any indication has benefit >= 0.1 and not hard_reject
    accept_any = any((v['p'] is not None and v['p'] >= 0.10) for v in benefit.values())
    decision = 'Reject' if hard_reject else ('Accept' if accept_any else 'Flag')

    records.append({
        'gene': gene,
        'cutoff': cutoff,
        'riskscore': riskscore,
        'max_high_p': max_high_p,
        'decision': decision,
        'tissue_info': tissue_info,
        'benefit': benefit,
        'utility': utility
    })

# Convert to DataFrame summary (one row per gene for ranking)
summary = []
for r in records:
    # best indication by utility
    if len(r['utility'])>0:
        best_ind = max(r['utility'].items(), key=lambda x: x[1])
        best_ind_name, best_util = best_ind
        best_benefit = r['benefit'][best_ind_name]['p']
    else:
        best_ind_name, best_util, best_benefit = (None, np.nan, np.nan)
    summary.append({
        'gene': r['gene'],
        'cutoff': r['cutoff'],
        'riskscore': r['riskscore'],
        'max_high_p': r['max_high_p'],
        'best_indication': best_ind_name,
        'best_benefit_pct': best_benefit,
        'best_utility': best_util,
        'decision': r['decision']
    })
summary_df = pd.DataFrame(summary).set_index('gene').sort_values(['best_utility'], ascending=False)

# Save or display
print(summary_df.head(30))

# ---------- PLOTS ----------
# 1) Risk vs Benefit scatter (best indication)
x = summary_df['riskscore']
y = summary_df['best_benefit_pct']
sizes = summary_df['best_benefit_pct'].fillna(0) * 200  # scale for plot point sizes
plt.figure(figsize=(6,5))
plt.scatter(x, y, s=sizes)
plt.xlabel('RiskScore (normalized)')
plt.ylabel('Best tumour % > cutoff (Benefit)')
plt.title('Risk vs Benefit (per gene)')
# Mark hard_reject threshold region (vertical line or shaded region)
plt.axhline(y=0.10, linestyle='--', label='Min benefit 10%')
plt.legend()
plt.show()

# 2) Heatmap for top-N genes: per-tissue p_normal
topn = 20
top_genes = summary_df.head(topn).index.tolist()
# Build a matrix genes x tissues
tissues = sorted(tissue_samples.keys())
mat = pd.DataFrame(index=top_genes, columns=tissues, dtype=float)
for g in top_genes:
    rec = next(r for r in records if r['gene'] == g)
    for t in tissues:
        v = rec['tissue_info'].get(t, {'p': np.nan})
        mat.loc[g, t] = v['p']
plt.figure(figsize=(max(8, len(tissues)*0.3), max(6, topn*0.2)))
plt.imshow(mat.fillna(0).values, aspect='auto', interpolation='nearest')
plt.yticks(range(len(top_genes)), top_genes)
plt.xticks(range(len(tissues)), tissues, rotation=90)
plt.title('Per-tissue percent normals > cutoff (top genes)')
plt.colorbar(label='p_normal')
plt.tight_layout()
plt.show()
