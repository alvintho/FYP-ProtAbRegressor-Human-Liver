import pandas as pd
import numpy as np
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from collections import Counter


class MRNAFeaturesCalculator:
    def __init__(self, df):
        self.df = df
        self.df_features = pd.DataFrame([])
        self.codon_table = {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "TAT": "Y",
            "TAC": "Y",
            "TAA": "*",
            "TAG": "*",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "TGT": "C",
            "TGC": "C",
            "TGA": "*",
            "TGG": "W",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }

    def __compute_features(self):
        self.df_features["gc_5utr"] = self.df["5utr"].astype(str).apply(gc_fraction)
        self.df_features["gc_cds"] = self.df["cds"].astype(str).apply(gc_fraction)
        self.df_features["gc_3utr"] = self.df["3utr"].astype(str).apply(gc_fraction)

        self.df_features["cds_mw"] = self.df["cds"].astype(str).apply(
            lambda seq: molecular_weight(seq, seq_type="DNA")
        )

        return self.df_features

    def __compute_translation_features(self):
        def calc_kozak_score(row):
            utr5 = row["5utr"]
            cds = row["cds"]
            kozak_context = utr5[-6:] + cds[:4] if len(utr5) >= 6 else utr5 + cds[:4]
            consensus = "GCCACCATGG"
            kozak_score = 0
            for i in range(min(len(kozak_context), len(consensus))):
                if kozak_context[i].upper() == consensus[i]:
                    kozak_score += 1
            return kozak_score

        self.df_features["kozak_score"] = self.df.apply(calc_kozak_score, axis=1)

        return self.df_features

    def __compute_codon_usage_features(self):
        def get_codons(seq):
            seq = Seq(seq.upper())
            if len(seq) % 3 != 0:
                seq = seq[: -(len(seq) % 3)]
            return [str(seq[i : i + 3]) for i in range(0, len(seq), 3)]

        def calculate_cai(codons):
            if not codons:
                return 0

            codon_counts = Counter(codons)
            total_codons = sum(codon_counts.values())

            if total_codons == 0:
                return 0

            aa_codon_counts = {}
            for codon, count in codon_counts.items():
                if codon in self.codon_table:
                    aa = self.codon_table[codon]
                    if aa != "*":  # Exclude stop codons
                        if aa not in aa_codon_counts:
                            aa_codon_counts[aa] = {}
                        aa_codon_counts[aa][codon] = count

            weights = {}
            for aa in aa_codon_counts:
                max_count = max(aa_codon_counts[aa].values())
                for codon in aa_codon_counts[aa]:
                    weights[codon] = (
                        aa_codon_counts[aa][codon] / max_count if max_count > 0 else 0
                    )

            valid_codons = [
                c for c in codons if c in weights and self.codon_table[c] != "*"
            ]

            if not valid_codons:
                return 0

            weight_values = [weights.get(codon, 0) for codon in valid_codons]
            weight_values = [w for w in weight_values if w > 0]

            if not weight_values:
                return 0

            return np.exp(np.mean(np.log(weight_values)))

        def calculate_rare_codon_freq(codons):
            if not codons:
                return pd.Series({"rare_codon_freq": 0})

            # Calculate rare codon frequency
            codon_counts = Counter(codons)
            total_codons = sum(codon_counts.values())
            if total_codons > 0:
                codon_freq = {c: n / total_codons for c, n in codon_counts.items()}
                if codon_freq:
                    threshold = np.percentile(list(codon_freq.values()), 10)
                    rare_codon_freq = sum(
                        freq for freq in codon_freq.values() if freq <= threshold
                    )
                else:
                    rare_codon_freq = 0
            else:
                rare_codon_freq = 0

            return pd.Series({"rare_codon_freq": rare_codon_freq})

        codons_list = self.df["cds"].apply(get_codons)

        self.df_features["cai"] = codons_list.apply(calculate_cai)

        rare_codon_freq = codons_list.apply(calculate_rare_codon_freq)
        self.df_features = pd.concat([self.df_features, rare_codon_freq], axis=1)

        return self.df_features

    def compute_mrna_features(self):
        self.__compute_features()
        self.__compute_translation_features()
        self.__compute_codon_usage_features()
        return self.df_features


class ProteinFeaturesCalculator:
    def __init__(self, df):
        self.df = df
        self.df_features = pd.DataFrame([])

        self.aa_groups = {
            "hydrophobic": set("AVILMFWY"),
            "hydrophilic": set("RKDENQHS"),
            "charged": set("RKDEH"),
            "aromatic": set("FYW"),
            "polar_uncharged": set("STNQ"),
        }

    def __compute_physicochemical_properties(self):

        def calculate_properties(seq):
            try:
                seq = seq.replace("*", "")
                if not seq or not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in seq):
                    return pd.Series(
                        {
                            "protein_mw": 0,
                            "protein_pi": 0,
                            "protein_aromaticity": 0,
                            "protein_instability": 0,
                            "protein_gravy": 0,
                            "protein_flexibility": 0,
                            "protein_helix_frac": 0,
                            "protein_turn_frac": 0,
                            "protein_sheet_frac": 0,
                            "protein_charge_ph7": 0,
                            "protein_aliphatic_index": 0,
                        }
                    )

                analysis = ProteinAnalysis(seq)

                ss_fractions = (
                    analysis.secondary_structure_fraction()
                )

              
                aliphatic_index = (
                    (
                        1.8 * (seq.count("A") / len(seq))
                        + 4.2 * (seq.count("V") / len(seq))
                        + 5.9 * (seq.count("I") / len(seq))
                        + 3.8 * (seq.count("L") / len(seq))
                    )
                    * 100
                    if seq
                    else 0
                )

                return pd.Series(
                    {
                        "protein_mw": analysis.molecular_weight(),
                        "protein_pi": analysis.isoelectric_point(),
                        "protein_aromaticity": analysis.aromaticity(),
                        "protein_instability": analysis.instability_index(),
                        "protein_gravy": analysis.gravy(),
                        "protein_flexibility": (
                            np.mean(analysis.flexibility())
                            if analysis.flexibility()
                            else 0
                        ),
                        "protein_helix_frac": ss_fractions[0],
                        "protein_turn_frac": ss_fractions[1],
                        "protein_sheet_frac": ss_fractions[2],
                        "protein_charge_ph7": analysis.charge_at_pH(7.0),
                        "protein_aliphatic_index": aliphatic_index,
                    }
                )
            except:
                return pd.Series(
                    {
                        "protein_mw": 0,
                        "protein_pi": 0,
                        "protein_aromaticity": 0,
                        "protein_instability": 0,
                        "protein_gravy": 0,
                        "protein_flexibility": 0,
                        "protein_helix_frac": 0,
                        "protein_turn_frac": 0,
                        "protein_sheet_frac": 0,
                        "protein_charge_ph7": 0,
                        "protein_aliphatic_index": 0,
                    }
                )

        physchem_features = self.df["protein_sequence"].apply(calculate_properties)
        self.df_features = pd.concat([self.df_features, physchem_features], axis=1)
        return self.df_features

    def __compute_amino_acid_composition(self):
        def calculate_composition(seq):
            try:
                seq = seq.replace("*", "")
                if not seq or not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in seq):
                    aa_features = {
                        f"protein_aa_freq_{aa}": 0 for aa in "ACDEFGHIKLMNPQRSTVWY"
                    }
                    aa_features.update(
                        {
                            "protein_hydrophobic_frac": 0,
                            "protein_hydrophilic_frac": 0,
                            "protein_charged_frac": 0,
                            "protein_aromatic_frac": 0,
                            "protein_polar_uncharged_frac": 0,
                        }
                    )
                    return pd.Series(aa_features)
                aa_counts = Counter(seq)
                total_aa = len(seq)
                aa_freq = {
                    f"protein_aa_freq_{aa}": aa_counts.get(aa, 0) / total_aa
                    for aa in "ACDEFGHIKLMNPQRSTVWY"
                }

                # Calculate amino acid group frequencies
                hydrophobic_count = sum(
                    aa_counts.get(aa, 0) for aa in self.aa_groups["hydrophobic"]
                )
                hydrophilic_count = sum(
                    aa_counts.get(aa, 0) for aa in self.aa_groups["hydrophilic"]
                )
                charged_count = sum(
                    aa_counts.get(aa, 0) for aa in self.aa_groups["charged"]
                )
                aromatic_count = sum(
                    aa_counts.get(aa, 0) for aa in self.aa_groups["aromatic"]
                )
                polar_uncharged_count = sum(
                    aa_counts.get(aa, 0) for aa in self.aa_groups["polar_uncharged"]
                )

                aa_features = aa_freq
                aa_features.update(
                    {
                        "protein_hydrophobic_frac": hydrophobic_count / total_aa,
                        "protein_hydrophilic_frac": hydrophilic_count / total_aa,
                        "protein_charged_frac": charged_count / total_aa,
                        "protein_aromatic_frac": aromatic_count / total_aa,
                        "protein_polar_uncharged_frac": polar_uncharged_count
                        / total_aa,
                    }
                )
                return pd.Series(aa_features)
            except:
                aa_features = {
                    f"protein_aa_freq_{aa}": 0 for aa in "ACDEFGHIKLMNPQRSTVWY"
                }
                aa_features.update(
                    {
                        "protein_hydrophobic_frac": 0,
                        "protein_hydrophilic_frac": 0,
                        "protein_charged_frac": 0,
                        "protein_aromatic_frac": 0,
                        "protein_polar_uncharged_frac": 0,
                    }
                )
                return pd.Series(aa_features)

        aa_composition = self.df["protein_sequence"].apply(calculate_composition)
        self.df_features = pd.concat([self.df_features, aa_composition], axis=1)
        return self.df_features

    def compute_protein_features(self):
        self.__compute_physicochemical_properties()
        self.__compute_amino_acid_composition()
        return self.df_features
