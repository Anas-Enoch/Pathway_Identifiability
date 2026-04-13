#!/usr/bin/env python3
"""
build_ST002829_pathway_mapping.py

Reads ST002829_benchmark_ready.csv (Metabolon naming conventions),
classifies each metabolite into a KEGG-aligned pathway,
and writes ST002829_pathway_mapping.csv to results/data/.

Run from repo root:
    python3 results/scripts/build_ST002829_pathway_mapping.py
"""
import re, sys
import pandas as pd
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE
for _ in range(6):
    if (_ROOT / "results" / "data").exists():
        break
    _ROOT = _ROOT.parent

BENCH_CSV = _ROOT / "results" / "data" / "ST002829_benchmark_ready.csv"
MAP_CSV   = _ROOT / "results" / "data" / "ST002829_pathway_mapping.csv"

METABOLON_RULES = [
    ("Ether Lipid Metabolism", [
        r'1-\(1-enyl-', r'\(P-\d+:\d+', r'plasmalogen',
    ]),
    ("Glycerophospholipid Metabolism", [
        r'-GPC\b', r'-GPC\s*\(', r'-GPE\b', r'-GPE\s*\(',
        r'-GPI\b', r'-GPI\s*\(', r'-GPS\b', r'-GPS\s*\(',
        r'-GPG\b', r'-GPG\s*\(',
        r'glycerophosphocholine', r'glycerophosphoethanolamine',
        r'glycerophosphoinositol', r'glycerophosphoserine',
        r'^PC\b', r'^LPC\b', r'^PE\b', r'^LPE\b', r'^PI\b', r'^PS\b',
    ]),
    ("Glycerolipid Metabolism", [
        r'glycerol\s*\(\d+:\d+\)',
        r'1-\w+glycerol\b', r'1-\w+glycerol\s',
        r'diacylglycerol', r'\bDAG\b', r'\bTG\b', r'\bDG\b',
        r'^TG\b', r'^DG\b', r'triacylglycerol', r'triglyceride',
    ]),
    ("Sphingolipid Metabolism", [
        r'ceramide', r'\bCer\b', r'Cer[_\s]', r'HexCer',
        r'sphingomyelin', r'\bSM\b', r'SM\s+d',
        r'sphingosine', r'sphinganine', r'dihydroceramide',
        r'glucosylceramide', r'galactosylceramide', r'lactosylceramide',
        r'sphingosine-1-phosphate',
    ]),
    ("Fatty Acid Degradation", [
        r'carnitine',
    ]),
    ("Fatty Acid Metabolism", [
        r'\(a\d+:\d+\)', r'\(i\d+:\d+\)',
        r'\(12 or 13\)-methyl', r'\(14 or 15\)-methyl', r'\(16 or 17\)-methyl',
        r'oleate\b', r'palmitate\b', r'stearate\b', r'myristate\b', r'laurate\b',
        r'linoleate\b', r'linolenate\b', r'arachidonate\b', r'adrenate\b',
        r'eicosanoate\b', r'docosanoate\b', r'hexadecanoate\b',
        r'propanoate\b', r'butanoate\b', r'pentanoate\b', r'hexanoate\b',
        r'heptanoate\b', r'octanoate\b', r'nonanoate\b', r'decanoate\b',
        r'undecanoate\b', r'dodecanoate\b', r'tridecanoate\b',
        r'tetradecanoate\b', r'pentadecanoate\b', r'heptadecanoate\b',
        r'nonadecanoate\b', r'eicosenoate\b', r'palmitoleate\b',
        r'vaccenate\b', r'nervonate\b', r'lignocerate\b',
        r'behenate\b', r'arachidate\b', r'margarate\b', r'gondoate\b',
        r'mead acid', r'hydroxy.*acid',
        r'dicarboxylate\b', r'sebacate\b', r'suberate\b', r'pimelate\b',
        r'decenoate\b', r'dodecenoate\b', r'undecenoate\b',
        r'eicosapentaenoate\b', r'docosahexaenoate\b', r'docosapentaenoate\b',
        r'eicosadienoate\b', r'eicosatrienoate\b', r'docosadienoate\b',
    ]),
    ("Steroid Hormone Biosynthesis", [
        r'cholesterol\b', r'cholesteryl\b',
        r'testosterone\b', r'cortisol\b', r'cortisone\b',
        r'corticosterone\b', r'progesterone\b',
        r'estradiol\b', r'estrone\b', r'androstenedione\b',
        r'\bDHEA\b', r'pregnenolone\b', r'aldosterone\b',
        r'dihydrotestosterone\b',
    ]),
    ("Primary Bile Acid Biosynthesis", [
        r'cholate\b', r'deoxycholate\b', r'ursodeoxycholate\b',
        r'chenodeoxycholate\b', r'glycocholate\b', r'taurocholate\b',
        r'glycochenodeoxycholate\b', r'taurochenodeoxycholate\b',
        r'glycoursodeoxycholate\b', r'lithocholate\b',
        r'hyocholate\b', r'hyodeoxycholate\b', r'bile acid\b',
    ]),
    ("Tryptophan Metabolism", [
        r'kynurenine\b', r'kynurenate\b', r'kynurenic\b',
        r'anthranilic\b', r'quinolinate\b', r'picolinate\b',
        r'xanthurenate\b', r'xanthurenic\b',
        r'3-hydroxykynurenine', r'3-hydroxyanthranilate',
        r'indole', r'indoxyl', r'serotonin\b', r'tryptamine\b',
        r'5-hydroxyindole', r'5-hydroxytryptophan',
        r'skatole\b', r'melatonin\b', r'tryptophan\b',
    ]),
    ("Purine Metabolism", [
        r'\badenosine\b', r'\bguanosine\b', r'\binosine\b', r'\bxanthosine\b',
        r'hypoxanthine\b', r'\bxanthine\b', r'uric acid\b', r'allantoin\b',
        r'\bAMP\b', r'\bADP\b', r'\bATP\b', r'\bGMP\b', r'\bGDP\b', r'\bIMP\b',
        r'deoxyadenosine\b', r'deoxyguanosine\b',
        r'\badenine\b', r'\bguanine\b', r'allopurinol\b',
    ]),
    ("Modified Nucleoside Metabolism", [
        r'N1-methyladenosine', r'N6-methyladenosine', r'1-methyladenosine',
        r'7-methylguanosine', r'N2,N2-dimethylguanosine',
        r'2-methyladenosine', r'N1-methylinosine', r'pseudouridine',
        r'N4-acetylcytidine', r'5-methylcytidine', r'5-methyluridine',
        r'isopentenyladenosine', r'N2-methylguanosine', r'methylinosine',
        r'1-methylguanidine',
    ]),
    ("Pyrimidine Metabolism", [
        r'\bcytidine\b', r'\buridine\b', r'\bthymidine\b',
        r'\bcytosine\b', r'\buracil\b', r'dihydrouracil\b',
        r'orotate\b', r'orotic acid',
        r'\bUMP\b', r'\bUDP\b', r'\bCMP\b', r'\bCDP\b', r'\bCTP\b',
        r'deoxycytidine\b', r'deoxyuridine\b', r'deoxythymidine\b',
        r'N-carbamoyl-',
    ]),
    ("TCA Cycle", [
        r'\bcitrate\b', r'isocitrate\b', r'alpha-ketoglutarate\b',
        r'2-oxoglutarate\b', r'\bsuccinate\b', r'succinic acid\b',
        r'\bfumarate\b', r'fumaric acid\b', r'\bmalate\b', r'malic acid\b',
        r'oxaloacetate\b', r'citramalate\b', r'aconitate\b',
        r'2-hydroxyglutarate\b', r'3-hydroxyglutarate\b',
    ]),
    ("Glycolysis / Gluconeogenesis", [
        r'1,5-anhydroglucitol', r'1,5-AG',
        r'\bglucose\b', r'\bfructose\b', r'\bgalactose\b', r'\bmannose\b',
        r'\blactate\b', r'\bpyruvate\b', r'glucuronate\b', r'gluconate\b',
        r'\bsorbitol\b', r'\btrehalose\b',
    ]),
    ("Amino Acid Metabolism", [
        r'\balanine\b', r'\barginine\b', r'\basparagine\b', r'\baspartate\b',
        r'\bcysteine\b', r'\bglutamine\b', r'\bglutamate\b', r'\bglycine\b',
        r'\bhistidine\b', r'\bisoleucine\b', r'\bleucine\b', r'\blysine\b',
        r'\bmethionine\b', r'\bphenylalanine\b', r'\bproline\b',
        r'\bserine\b', r'\bthreonine\b', r'\btyrosine\b', r'\bvaline\b',
        r'\bcitrulline\b', r'\born[iI]thine\b', r'\bcreatine\b', r'\bcreatinine\b',
        r'\btaurine\b', r'\bbetaine\b', r'\bsarcosine\b',
        r'\bhomocysteine\b', r'\bcystathionine\b', r'ergothioneine\b',
        r'N-acetyl\w+amino', r'N-acetylaspartate\b',
        r'\bdimethylglycine\b', r'dimethylarginine\b',
        r'guanidineacetate\b', r'guanidinoacetate\b',
        r'\bpipecolic\b', r'\bhomoserine\b', r'3-methylhistidine\b',
        r'1-methylhistidine\b', r'cyclo\(pro', r'prolyl\w+',
        r'\bspermidine\b', r'\bspermine\b', r'acetylspermidine\b',
        r'putrescine\b', r'cadaverine\b',
    ]),
    ("Valine Leucine Isoleucine Degradation", [
        r'3-methyl-2-oxobutyrate', r'3-methyl-2-oxovalerate',
        r'4-methyl-2-oxopentanoate', r'4-methyl-2-oxovalerate',
        r'isovalerate\b', r'isobutyrate\b', r'2-methylbutyrate\b',
        r'3-methylglutarate\b', r'3-hydroxyisobutyrate\b',
        r'2-hydroxy-3-methylvalerate', r'2-hydroxy-3-methylbutyrate',
        r'isovalerylglycine\b', r'isobutyrylglycine\b', r'tiglylglycine\b',
        r'methylsuccinate\b', r'ethylmalonate\b',
    ]),
    ("Phenylalanine Metabolism", [
        r'phenylpyruvate\b', r'phenyllactate\b', r'phenylacetate\b',
        r'hippurate\b', r'phenylacetylglutamine\b', r'homogentisate\b',
        r'4-hydroxyphenylpyruvate', r'4-hydroxyphenyllactate',
        r'4-hydroxyphenylacetate', r'cinnamate\b', r'cinnamoyl\b',
        r'p-cresol\b', r'3,4-dihydroxyphenyl', r'dopamine\b',
    ]),
    ("One Carbon Pool by Folate", [
        r'S-adenosyl', r'\bSAM\b', r'\bSAH\b',
        r'5-methyltetrahydrofolate', r'folate\b', r'folic acid',
    ]),
    ("Nicotinate and Nicotinamide Metabolism", [
        r'nicotinamide\b', r'nicotinate\b', r'nicotinic acid',
        r'nicotinurate\b', r'1-methylnicotinamide\b',
        r'\bNAD\b', r'\bNADH\b', r'\bNMN\b',
        r'trigonelline\b', r'quinolinate\b',
        r'N-methyl-2-pyridone', r'N-methyl-4-pyridone',
    ]),
    ("Glutathione Metabolism", [
        r'glutathione\b', r'\bGSH\b', r'\bGSSG\b',
        r'cysteinylglycine\b', r'gamma-glutamyl', r'5-oxoproline\b',
        r'pyroglutamate\b',
    ]),
    ("Histidine Metabolism", [
        r'imidazole\b', r'imidazolyl\b', r'imidazoleacetate\b',
        r'imidazolelactate\b', r'urocanate\b',
        r'1-methyl-4-imidazole', r'1-methyl-5-imidazole',
        r'carnosine\b', r'anserine\b', r'homocarnosine\b',
    ]),
    ("Vitamin B6 Metabolism", [
        r'pyridoxine\b', r'pyridoxal\b', r'pyridoxamine\b',
        r'pyridoxic\b', r'4-pyridoxate\b', r'\bPLP\b',
    ]),
    ("Butanoate Metabolism", [
        r'3-hydroxybutyrate\b', r'3-hydroxybutanoate\b',
        r'acetoacetate\b', r'beta-hydroxybutyrate\b',
        r'3-hydroxybutyrylcarnitine\b',
    ]),
    ("Pentose Phosphate Pathway", [
        r'\bribose\b', r'\bribulose\b', r'\bxylulose\b',
        r'sedoheptulose\b', r'\berythrose\b', r'gluconolactone\b',
    ]),
    ("Porphyrin and Chlorophyll Metabolism", [
        r'bilirubin\b', r'biliverdin\b', r'\bheme\b', r'porphyrin\b',
    ]),
    ("Urea Cycle", [
        r'\burea\b', r'argininosuccinate\b',
    ]),
    ("Cofactor Biosynthesis", [
        r'pantothenate\b', r'pantothenic\b', r'riboflavin\b',
        r'dehydroascorbate\b', r'ascorbate\b', r'biotin\b', r'thiamine\b',
    ]),
    ("Drug Metabolism", [
        r'acetaminophen\b', r'ibuprofen\b', r'caffeine\b',
        r'theophylline\b', r'cotinine\b', r'benzoate\b',
    ]),
]

def classify(name: str) -> str | None:
    nl = name.strip().lower()
    for pathway, patterns in METABOLON_RULES:
        for pat in patterns:
            try:
                if re.search(pat, nl, re.IGNORECASE):
                    return pathway
            except re.error:
                continue
    return None


def main():
    if not BENCH_CSV.exists():
        sys.exit(f"[ERROR] Not found: {BENCH_CSV}\n"
                 "Run the benchmark-ready CSV generation first.")

    df = pd.read_csv(BENCH_CSV, nrows=0)
    mets = [c for c in df.columns if c != "condition"]
    print(f"Classifying {len(mets)} metabolites from Metabolon nomenclature...")

    rows, missed = [], []
    for m in mets:
        pw = classify(m)
        if pw:
            rows.append({"metabolite": m, "pathway": pw})
        else:
            missed.append(m)

    df_map = pd.DataFrame(rows)
    print(f"\nCoverage: {len(rows)}/{len(mets)} ({100*len(rows)/len(mets):.1f}%)")

    vc = df_map["pathway"].value_counts()
    print("\nPathway distribution:")
    print(vc.to_string())
    viable = (vc >= 3).sum()
    print(f"\nPathways with ≥3 metabolites (viable for benchmark): {viable}")

    print(f"\nUnclassified: {len(missed)}")
    if missed:
        print("First 20 unclassified:")
        for m in missed[:20]:
            print(f"  {m}")

    MAP_CSV.parent.mkdir(parents=True, exist_ok=True)
    df_map.to_csv(MAP_CSV, index=False)
    print(f"\nSaved → {MAP_CSV}")
    print("\nNext: python3 results/scripts/active/run_ST002829_benchmark.py")


if __name__ == "__main__":
    main()
