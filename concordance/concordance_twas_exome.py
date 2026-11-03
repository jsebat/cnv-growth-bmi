# Concordance of CNV effects with common-variant (TWAS) and rare-variant (exome) gene effects on height and BMI.
# JS 2026-11-02, reviewer request; goes in the supplement after the single-hit model results (Fig. 3).
#
# Inputs
#   CNV intervals, genes and class: Molly/pathway_analysis/Final/cnv_genes_annotated.csv (34 reciprocal intervals, 942 genes,
#     Dose_response_class mirror / asymmetric from Fig. 3; the class was defined on height and is used for both traits)
#   CNV effects: Sacks Supplementary Tables.xlsx, T3 (meta-analysed DEL and DUP betas, SD units, height and BMI)
#   Common variants: TWAS30 (Mancuso et al. 2017 AJHG; FUSION, >40 GTEx/CMC/METSIM/NTR/YFS panels), height.txt and bmi.txt,
#     downloaded from github.com/bogdanlab/TWAS30 -> twas30/
#   Rare variants: Kosmicki et al. 2026 (medRxiv 10.64898/2026.06.22.26355163) Supplementary Table 6, burden tests for the 207
#     exome-wide height genes -> exome/kosmicki_2026_supp_tables.xlsx; Akbari et al. 2021 (Science) Table 1, 16 exome-wide BMI
#     genes (gene, chromosome, position; transcribed below)
#
# Logic
#   A TWAS z > 0 means higher predicted expression goes with a higher trait. A duplication raises dosage of every gene in the
#   interval, a deletion lowers it, so the expected sign of the DUP effect is sign(S) and of the DEL effect is -sign(S), where
#   S = sum over genes in the interval of the gene's TWAS z (mean over tissues). Tests per trait, allele and CNV group
#   (all, mirror, asymmetric): sign concordance (one-sided binomial against 0.5) and Spearman correlation of S with the
#   observed effect (expected + for DUP, - for DEL). Sensitivity: S from transcriptome-wide significant genes only
#   (best-tissue P < 0.05 / number of gene-tissue tests), CNVs without such a gene dropped.
#   Exome: genes with an exome-wide burden association inside a CNV, paired with the CNV deletion effect (loss of function
#   ~ deletion); concordance = same sign of the gene's LoF (pLoF, else pLoF+missense) burden effect and the DEL effect.
# Nested intervals (e.g. 22q11.21 A-D and its subregions) are kept as in Fig. 3, so CNVs are not independent; P values are
# descriptive.
# Outputs: concordance_twas_by_cnv.csv, concordance_twas_tests.csv, concordance_twas_permutation.csv, concordance_exome_genes.csv
import re, numpy as np, pandas as pd
from scipy import stats

ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
HERE = ROOT + r"\CNV_ASD_growth_analysis\concordance"
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
t3 = pd.read_excel(ROOT + r"\NatureSubmission_04_13_2026\Sacks Supplementary Tables.xlsx", sheet_name="T3 CNV main effects and frequen")
norm = lambda s: re.sub(r"[^0-9a-z]", "", str(s).lower())

# ---- CNV effects (DEL, DUP) per interval
t3["allele"] = t3.locus.str.extract(r"_(DEL|DUP)$")[0]
t3["key"] = t3.locus.str.replace(r"_(DEL|DUP)$", "", regex=True).map(norm)
cnv = genes.groupby("CNV_ID").agg(cls=("Dose_response_class", "first"), n_genes=("Gene", "nunique")).reset_index()
cnv["key"] = cnv.CNV_ID.map(norm)
for a in ["DEL", "DUP"]:
    x = t3[t3.allele == a].set_index("key")
    for tr, col in [("height", "Height"), ("bmi", "BMI")]:
        cnv[f"{a}_{tr}"] = cnv.key.map(x[f"beta{col}"]); cnv[f"{a}_{tr}_p"] = cnv.key.map(x[f"p{col}"])
print("CNV intervals matched to T3:", cnv[["DEL_height", "DUP_height"]].notna().sum().to_dict(), "of", len(cnv))

# ---- TWAS gene scores
rows = []
for tr in ["height", "bmi"]:
    tw = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+")
    ntest = len(tw); thr = 0.05 / ntest
    g = tw.groupby("GENE").agg(z_mean=("TWAS.Z", "mean"), n_tissue=("TWAS.Z", "size"), p_best=("TWAS.P", "min")).reset_index()
    best = tw.loc[tw.groupby("GENE")["TWAS.P"].idxmin(), ["GENE", "TWAS.Z"]].rename(columns={"TWAS.Z": "z_best"})
    g = g.merge(best, on="GENE"); g["sig"] = g.p_best < thr
    m = genes[["CNV_ID", "Gene"]].merge(g, left_on="Gene", right_on="GENE", how="inner")
    agg = m.groupby("CNV_ID").agg(n_twas=("Gene", "nunique"), S=("z_mean", "sum"), n_sig=("sig", "sum")).reset_index()
    sig = m[m.sig].groupby("CNV_ID").agg(S_sig=("z_best", "sum")).reset_index()
    agg = agg.merge(sig, on="CNV_ID", how="left"); agg["trait"] = tr; agg["thr"] = thr
    rows.append(agg)
tw = pd.concat(rows)
by = cnv.merge(tw, on="CNV_ID", how="left")
# Mirror/asymmetric classes come from the height dose-response; for BMI there are no asymmetric CNVs (Jon 2026-11-02), so every
# CNV is mirror for BMI and BMI is analysed in all CNVs only (GROUPS_FOR below).
by.loc[by.trait == "bmi", "cls"] = "mirror"
GROUPS_FOR = {"height": ["all", "mirror", "asymmetric"], "bmi": ["all"]}
by.to_csv(HERE + r"\concordance_twas_by_cnv.csv", index=False)

# ---- tests
GROUPS = {"all": lambda d: d, "mirror": lambda d: d[d.cls == "mirror"], "asymmetric": lambda d: d[d.cls == "asymmetric"]}
out = []
for tr in ["height", "bmi"]:
    d0 = by[by.trait == tr]
    for gname, f in [(k, GROUPS[k]) for k in GROUPS_FOR[tr]]:
        for a, exp_sign in [("DEL", -1), ("DUP", +1)]:
            for score in ["S", "S_sig"]:
                d = f(d0).dropna(subset=[f"{a}_{tr}", score]); d = d[d[score] != 0]
                if len(d) < 3: out.append(dict(trait=tr, group=gname, allele=a, score=score, n_cnv=len(d))); continue
                conc = (np.sign(d[f"{a}_{tr}"]) == exp_sign * np.sign(d[score])).sum()
                bt = stats.binomtest(int(conc), len(d), 0.5, alternative="greater")
                rho, p = stats.spearmanr(d[score], d[f"{a}_{tr}"])
                p1 = p / 2 if np.sign(rho) == exp_sign else 1 - p / 2
                out.append(dict(trait=tr, group=gname, allele=a, score=score, n_cnv=len(d), n_concordant=int(conc),
                                frac_concordant=conc / len(d), p_binomial=bt.pvalue, spearman_rho=rho, expected_sign=exp_sign, p_spearman_one_sided=p1))
T = pd.DataFrame(out); T.to_csv(HERE + r"\concordance_twas_tests.csv", index=False)
pd.set_option("display.width", 220)
print(T.round(3).to_string(index=False))


# ---- gene-label permutation null for the Spearman correlation of S (keeps nested intervals, gene counts and CNV effects fixed;
#      shuffles the mean TWAS z across all genes with TWAS data, genome-wide), 10,000 permutations
rng = np.random.default_rng(20261102); NPERM = 10000
perm_rows = []
for tr in ["height", "bmi"]:
    tw0 = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+").groupby("GENE")["TWAS.Z"].mean()
    zpool = tw0.values
    gm = genes[["CNV_ID", "Gene"]][genes.Gene.isin(tw0.index)].copy()
    cnv_ids = sorted(gm.CNV_ID.unique()); idx = {c: i for i, c in enumerate(cnv_ids)}
    rowi = gm.CNV_ID.map(idx).values; zobs = gm.Gene.map(tw0).values
    eff = cnv.set_index("CNV_ID").loc[cnv_ids]
    for gname in GROUPS_FOR[tr]:
        keep = np.array([True if gname == "all" else eff.loc[c, "cls"] == gname for c in cnv_ids])
        for a, es in [("DEL", -1), ("DUP", +1)]:
            y = eff[f"{a}_{tr}"].values
            ok = keep & ~np.isnan(y)
            S_obs = np.bincount(rowi, weights=zobs, minlength=len(cnv_ids))
            r_obs = stats.spearmanr(S_obs[ok], y[ok])[0]
            yr = stats.rankdata(y[ok])
            null = np.empty(NPERM)
            for b in range(NPERM):
                zb = rng.choice(zpool, size=len(zobs), replace=False)
                Sb = np.bincount(rowi, weights=zb, minlength=len(cnv_ids))[ok]
                null[b] = np.corrcoef(stats.rankdata(Sb), yr)[0, 1]
            p_perm = (1 + np.sum(es * null >= es * r_obs)) / (NPERM + 1)
            perm_rows.append(dict(trait=tr, group=gname, allele=a, n_cnv=int(ok.sum()), spearman_rho=r_obs, p_perm_one_sided=p_perm))
P = pd.DataFrame(perm_rows); P.to_csv(HERE + r"\concordance_twas_permutation.csv", index=False)
print("\nGene-label permutation (10,000), one-sided in the expected direction:"); print(P.round(4).to_string(index=False))

# ---- exome genes inside CNVs
k6 = pd.read_excel(HERE + r"\exome\kosmicki_2026_supp_tables.xlsx", sheet_name="Table S6", header=1)
akbari = pd.DataFrame([("UHMK1", 1, 162497250, -0.06), ("GPR75", 2, 53853133, -0.34), ("ROBO1", 3, 78598912, 0.24), ("KIAA1109", 4, 122170647, 0.13),
    ("PCSK1", 5, 96393000, 0.33), ("GPR151", 5, 146514853, -0.05), ("SPARC", 5, 151663570, 0.05), ("UBR2", 6, 42564319, 0.34), ("CALCR", 7, 93426355, 0.09),
    ("PDE3B", 11, 14644075, 0.12), ("ANO4", 12, 100901785, 0.29), ("KIAA0586", 14, 58427628, -0.04), ("MC4R", 18, 60371350, 0.30), ("DPP9", 19, 4676563, -0.06),
    ("ANKRD27", 19, 32598144, 0.04), ("GIPR", 19, 45669520, -0.09)], columns=["Gene", "chr", "pos", "beta_sd"])
ex = []
for tr, gl in [("height", set(k6["Gene symbol"].dropna())), ("bmi", set(akbari.Gene))]:
    hit = genes[genes.Gene.isin(gl)][["CNV_ID", "Gene", "Dose_response_class"]]
    for _, r in hit.iterrows():
        c = cnv[cnv.CNV_ID == r.CNV_ID].iloc[0]
        if tr == "height":
            kg = k6[k6["Gene symbol"] == r.Gene]
            lof = kg[kg["Burden category"].astype(str).str.contains("pLoF", case=False)]
            pick = (lof if len(lof) else kg).sort_values("P-value").iloc[0]
            eff, cat, p = pick["Effect (SD units)"], pick["Burden category"], pick["P-value"]
        else:
            ak = akbari[akbari.Gene == r.Gene].iloc[0]; eff, cat, p = ak.beta_sd, "Akbari Table 1", np.nan
        ex.append(dict(trait=tr, CNV_ID=r.CNV_ID, cls=r.Dose_response_class, gene=r.Gene, burden_category=cat, gene_effect_sd=eff, gene_p=p,
                       DEL_effect_sd=c[f"DEL_{tr}"], DEL_p=c[f"DEL_{tr}_p"], concordant=np.sign(eff) == np.sign(c[f"DEL_{tr}"])))
E = pd.DataFrame(ex, columns=["trait", "CNV_ID", "cls", "gene", "burden_category", "gene_effect_sd", "gene_p", "DEL_effect_sd", "DEL_p", "concordant"])
E.to_csv(HERE + r"\concordance_exome_genes.csv", index=False)
print("\nExome-wide genes inside CNV intervals:"); print(E.round(4).to_string(index=False))
print("height genes tested (Kosmicki):", k6["Gene symbol"].nunique(), " BMI genes tested (Akbari):", len(akbari))
