import numpy as np
import math
import warnings
from Bio.PDB import PDBParser
from typing import Dict, Optional, Any, List, Set, Union
from pathlib import Path
warnings.filterwarnings("ignore")

"""
Calculates the iPSAE metric, linked to Dunbrack's preprint publication https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2, posted in December 2025.
Reused by during Nipah competition (with slight adaptations): https://github.com/adaptyvbio/nipah_ipsae_pipeline/blob/main/ipsae.py.

The Dunbrack's script is adapted here to the current evaluation pipeline workflow.

"""

import numpy as np
import math
import warnings
from Bio.PDB import PDBParser
from typing import Dict, Optional, Any, List, Set, Union
from pathlib import Path

# Suppress Biopython PDB warnings
warnings.filterwarnings("ignore")


# taken from Evaluation Engine

class IPSAECalculator:
    """
    Calculates iPSAE, pDockQ, and LIS metrics following the Dunbrack Lab methodology (v4).
    Ref: https://github.com/DunbrackLab/IPSAE (version of January 2026: https://github.com/DunbrackLab/IPSAE/commit/6174cf9e71cb1bd660cc805856a18c4871a6dec3)
    """
    # 8 angstrom as per evaluation engine "interactions.contact_cutoff" distance = 8 (see configs)
    def __init__(self, pae_cutoff: float = 8.0, dist_cutoff: float = 8.0, pdockq_cutoff: float = 8.0):
        self.pae_cutoff = pae_cutoff
        self.dist_cutoff = dist_cutoff
        self.pdockq_cutoff = pdockq_cutoff
        self.parser = PDBParser(QUIET=True)

    def _calc_d0(self, L: float, pair_type: str = 'protein') -> float:
        """Calculates d0 normalization factor based on protein length (copy pasted from calc_d0 function of Dunbrack lab).
        The more the protein/complex is large, the more the error is big and should be normalized to avoid penalizing too much the large proteins.

        Args:
            L: float, total length of protein chain(s).
            pair_type: str, type of considered pair (default: protein)
        Returns:
            d0 value (float): the calculated d0 normalization factor

        """
        L = float(L)
        min_value = 1.0
        if pair_type == 'nucleic_acid':
            min_value = 2.0

        # TM-score formula proposed by Zhang and Skolnick (https://pubmed.ncbi.nlm.nih.gov/15849316/)
        if L > 27:
            d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
        else:
            d0 = 1.0

        return max(min_value, d0)

    def _calc_d0_array(self, L: np.ndarray, pair_type: str = 'protein') -> np.ndarray:
        """Vectorized d0 calculation."""
        L = np.array(L, dtype=float)
        L = np.maximum(26, L)  # Fixed based on script comment 01.03.2026
        min_value = 1.0
        if pair_type == 'nucleic_acid':
            min_value = 2.0

        return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)

    def _ptm_func(self, x: np.ndarray, d0: float) -> np.ndarray:
        """Vectorized PTM function: 1 / (1 + (x/d0)^2)."""
        return 1.0 / (1 + (x / d0) ** 2.0)

    def compute(
            self,
            pdb_path: str,
            pae_matrix: np.ndarray,
            plddt_vector: np.ndarray,
            binder_chain: str = "B",
            target_chain: str = "A"
    ) -> Dict[str, float]:
        try:
            # 1. Parse Structure
            structure = self.parser.get_structure("complex", pdb_path)
            model = structure[0]

            ca_coords = []
            chains_list = []
            plddts = []

            current_idx = 0
            for chain in model:
                for res in chain:
                    if 'CA' in res and res.id[0] == ' ':
                        chains_list.append(chain.id)
                        ca_coords.append(res['CA'].get_coord())
                        if current_idx < len(plddt_vector):
                            plddts.append(plddt_vector[current_idx])
                        current_idx += 1

            numres = len(chains_list)
            chains = np.array(chains_list)
            coordinates = np.array(ca_coords)
            plddt_final = np.array(plddts)

            # Dimension synchronization
            if pae_matrix.shape[0] != numres:
                pae_matrix = pae_matrix[:numres, :numres]
            if plddt_final.shape[0] != numres:
                plddt_final = plddt_final[:numres]

            # Matrix of distances
            diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
            distances = np.sqrt(np.nansum(diff ** 2, axis=2))

            # Chain indices
            c1, c2 = target_chain, binder_chain
            c1_indices = np.where(chains == c1)[0]
            c2_indices = np.where(chains == c2)[0]

            if len(c1_indices) == 0 or len(c2_indices) == 0:
                return {"max_ipsae": 0.0, "min_ipsae": 0.0, "pdockq": 0.0, "pdockq2": 0.0, "lis": 0.0}

            # --- pDockQ / pDockQ2 ---
            pdockq, pdockq2 = 0.0, 0.0

            # Mask contacts between C1 and C2
            # Matrix of booleans (numres, numres)
            contact_mask = (distances <= self.pdockq_cutoff)
            # Isolation of C1 (rows) vs C2 (columns)
            interface_mask = contact_mask[np.ix_(c1_indices, c2_indices)]

            npairs = np.sum(interface_mask)

            if npairs > 0:
                # Identification of interface residues (at least one contact between them)
                c1_if_idx = c1_indices[np.any(interface_mask, axis=1)]
                c2_if_idx = c2_indices[np.any(interface_mask, axis=0)]

                # Mean pLDDT of interface residues
                interface_plddts = np.concatenate([plddt_final[c1_if_idx], plddt_final[c2_if_idx]])
                mean_plddt = np.mean(interface_plddts)

                # pDockQ
                x_q1 = mean_plddt * math.log10(npairs)
                pdockq = 0.724 / (1 + math.exp(-0.052 * (x_q1 - 152.611))) + 0.018

                # pDockQ2
                # Mean PTM on contact pairs
                pae_contacts = pae_matrix[np.ix_(c1_indices, c2_indices)][interface_mask]
                mean_ptm = np.mean(self._ptm_func(pae_contacts, 10.0))
                x_q2 = mean_plddt * mean_ptm
                pdockq2 = 1.31 / (1 + math.exp(-0.075 * (x_q2 - 84.733))) + 0.005

            # --- LIS ---
            sub_pae = pae_matrix[np.ix_(c1_indices, c2_indices)]
            valid_pae = sub_pae[sub_pae <= 12.0]
            lis_score = np.mean((12.0 - valid_pae) / 12.0) if valid_pae.size > 0 else 0.0

            # --- ipSAE (directional) ---
            def calc_directional_scores(idxA, idxB):
                # PAE A (rows) vs B (columns)
                sub_pae_ab = pae_matrix[np.ix_(idxA, idxB)]
                valid_mask = sub_pae_ab < self.pae_cutoff

                n0chn = len(idxA) + len(idxB)
                d0chn = self._calc_d0(n0chn)

                # n0dom : residues involved in the interface
                n0dom = np.sum(np.any(valid_mask, axis=1)) + np.sum(np.any(valid_mask, axis=0))
                d0dom = self._calc_d0(n0dom)

                # n0res per residue from A
                n0res_list = np.sum(valid_mask, axis=1)
                d0res_list = self._calc_d0_array(n0res_list)

                scores_chn, scores_dom, scores_res = [], [], []

                for i_local in range(len(idxA)):
                    if n0res_list[i_local] > 0:
                        pae_vals = sub_pae_ab[i_local, valid_mask[i_local]]

                        scores_chn.append(np.mean(self._ptm_func(pae_vals, d0chn)))
                        scores_dom.append(np.mean(self._ptm_func(pae_vals, d0dom)))
                        scores_res.append(np.mean(self._ptm_func(pae_vals, d0res_list[i_local])))
                    else:
                        scores_chn.append(0.0);
                        scores_dom.append(0.0);
                        scores_res.append(0.0)

                return {
                    "chn": np.max(scores_chn) if scores_chn else 0.0,
                    "dom": np.max(scores_dom) if scores_dom else 0.0,
                    "res": np.max(scores_res) if scores_res else 0.0
                }

            res_AB = calc_directional_scores(c1_indices, c2_indices)
            res_BA = calc_directional_scores(c2_indices, c1_indices)

            return {
                # Local: normalized by d0res (this is the min ipsae value we can use for min)
                # Interface: normalized by d0dom
                # Global: normalized by d0chn
                "ipsae_local_max": float(max(res_AB["res"], res_BA["res"])),
                "ipsae_global_max": float(max(res_AB["chn"], res_BA["chn"])),
                "ipsae_interface_max": float(max(res_AB["dom"], res_BA["dom"])),
                "iptm_global_max": float(max(res_AB["chn"], res_BA["chn"])),
                "ipsae_local_min": float(min(res_AB["res"], res_BA["res"])),
                "ipsae_global_min": float(min(res_AB["chn"], res_BA["chn"])),
                "ipsae_interface_min": float(min(res_AB["dom"], res_BA["dom"])),
                "iptm_global_min": float(min(res_AB["chn"], res_BA["chn"])),
                "pdockq": float(pdockq),
                "pdockq2": float(pdockq2),
                "lis": float(lis_score),
                "npairs_pdockq": int(npairs),
                "mean_plddt_interface": float(mean_plddt) if 'mean_plddt' in locals() else 0.0
            }

        except Exception as e:
            print(f"  [ERROR] Failed to compute IPSAE metrics: {e}")
            return {"ipsae": 0.0, "pdockq": 0.0, "pdockq2": 0.0, "lis": 0.0}

    def compute_from_boltz(
            self,
            pdb_path: str,
            pae_npz_path: str,
            plddt_npz_path: str,
            binder_chain: str = "H",
            target_chain: str = "T"
    ) -> Dict[str, float]:

        # 1. Keep model 0 from complex
        structure = self.parser.get_structure("complex", pdb_path)
        model = structure[0]

        num_res_complex = 0
        for chain in model:
            # Filter on standards residues (alpha carbon)
            for residue in chain:
                if 'CA' in residue and residue.id[0] == ' ':
                    num_res_complex += 1

        print(f"  [DEBUG] PDB residues (standard CA): {num_res_complex}")

        # 2. Load PAE matrix and PLDDT vector, truncate with a warning if there is a mismatch of residues number
        # should not happen but can indicate something went wrong with boltz output usage (mismatch between pdb and matrix for example)
        pae = load_boltz_npz_filtered(pae_npz_path, num_res_complex, 'pae')
        plddt = load_boltz_npz_filtered(plddt_npz_path, num_res_complex, 'plddt')

        print(f"  [DEBUG] NPZ Shapes after filter: PAE {pae.shape}, pLDDT {plddt.shape}")

        # 3. Run metrics calculation
        return self.compute(
            pdb_path=pdb_path,
            pae_matrix=pae,
            plddt_vector=plddt,
            binder_chain=binder_chain,
            target_chain=target_chain
        )

    def get_metrics_description(self) -> Dict[str, str]:
        """Returns a brief description of the calculated metrics, if needed for printing/info."""
        return {
            "ipsae": "ipSAE with d0 based on n0res (per-residue local interface size).",
            "ipsae_d0chn": "ipSAE with d0 based on total residues in the chain pair.",
            "ipsae_d0dom": "ipSAE with d0 based on residues actually in the interface (PAE < cutoff).",
            "pdockq": "Standard pDockQ (Bryant et al. 2022).",
            "pdockq2": "pDockQ2 using PAE-based PTM (Zhu et al. 2023).",
            "lis": "Logit-Isomorphic Score (Kim et al. 2024)."
        }

