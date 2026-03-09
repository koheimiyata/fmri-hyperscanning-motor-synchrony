# InterSync: fMRI Hyperscanning Study of Interpersonal Motor Synchrony

This repository contains analysis code and behavioral data for an fMRI hyperscanning study of interpersonal motor synchrony.
The corresponding fMRI dataset is available on OpenNeuro.

---

## Dependencies

- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) or later (analyses were conducted with SPM25)
- [fMRIPrep](https://fmriprep.org/) v21+
- MATLAB R2020b or later
- `first_eigenvariate.m` (included in `Functions/`)
- `r_to_zfisher.m` (included in `Functions/`)
- `predict_behavior.m` (included in `Functions/`)
- `fdr_bh.m` (included in `Functions/`)

---

## Directory Structure

The scripts assume the following directory structure:

```
(root)/
├── BIDS/
│   └── InterSync/
│       └── derivatives/
│           └── fMRIPrep/          # fMRIPrep output
└── Analysis/
    ├── Group/
    │   ├── InterBrain/
    │   │   └── CPM/               # 4D NIfTI files for CPM analysis
    │   └── OneT/                  # Second-level SPM results
    ├── Mat_file/                  # Intermediate .mat files
    │   └── CPM/
    ├── Mfile/                     # Analysis scripts (this repository)
    │   └── Functions/             # Helper functions
    ├── ROI_temp/                  # ROI mask files (e.g., Shen268)
    └── SPM_temp/                  # SPM batch template .mat files
```

Each participant pair has its own folder under `Analysis/`:

```
Analysis/
└── P02/
    ├── subA/
    │   ├── R1st/                  # 3D NIfTI volumes (run 1)
    │   ├── R2nd/                  # 3D NIfTI volumes (run 2)
    │   ├── T1/                    # Anatomical images
    │   ├── Results_CPM/           # First-level SPM results (CPM)
    │   ├── Results_GLM/           # First-level SPM results (GLM)
    │   └── batch/                 # Saved SPM batch files
    └── subB/
        └── ...
```

---

## Script Overview

Scripts are prefixed with a letter indicating the analysis stage and should be run in numerical order within each stage.

### P — Preprocessing

| Script | Description |
|--------|-------------|
| `P1_fMRIPrep2SPM.m` | Convert fMRIPrep output to SPM format; split 4D NIfTI into 3D volumes |
| `P2_Discard_initialFiles.m` | Move the first 3 volumes of each run to a `discarded/` subfolder (T1 equilibration) |
| `P3_SaveConfoundsAsMat.m` | Extract confound regressors from fMRIPrep TSV files; save as `.txt` and `.mat` |
| `P4_Smooth.m` | Apply spatial smoothing (5 mm FWHM Gaussian kernel) using SPM |

### C — Connectome-based Predictive Modeling (CPM)

| Script | Description |
|--------|-------------|
| `C1_FLevel_CPM.m` | Build and run SPM first-level GLM for CPM (nuisance regression only, no task regressors) |
| `C2_ResidT2F.m` | Concatenate SPM residuals (Res_*.nii) into 4D NIfTI files per run |
| `C3_ConMatrix_Shen268.m` | Compute intra- and inter-brain ROI-to-ROI connectivity matrices using Shen 268 atlas |
| `C4_Permutation_CPM.m` | Run CPM with leave-one-out cross-validation to predict behavioral scores |
| `C5_SummarizeMatrix.m` | Summarize CPM positive-mask edges at the canonical network level |

### G — GLM Analysis

| Script | Description |
|--------|-------------|
| `G1_FLevel_GLM.m` | Build and run SPM first-level GLM with task and nuisance regressors (FAST autocorrelation) |
| `G2_CManager.m` | Run SPM Contrast Manager for each subject |
| `G3_SLevel_GLM.m` | Run SPM second-level one-sample t-test for each contrast |

---

## Execution Order

```
[Preprocessing]
P1_fMRIPrep2SPM
  → P2_Discard_initialFiles
  → P3_SaveConfoundsAsMat
  → P4_Smooth

[CPM analysis]
C1_FLevel_CPM
  → C2_ResidT2F
  → C3_ConMatrix_Shen268
  → C4_Permutation_CPM
  → C5_SummarizeMatrix

[GLM analysis]
G1_FLevel_GLM
  → G2_CManager
  → G3_SLevel_GLM
```

---

## Notes

- Pair IDs P01 and P19 are excluded from the analysis (see manuscript for details).
- The prefix `R` in run names (`R1st`, `R2nd`) denotes the circle-drawing task runs, distinguishing them from other tasks in this project.
- All scripts use relative paths based on the `Analysis/Mfile/` working directory. Before running, set the MATLAB working directory to `Analysis/Mfile/`.
- Helper functions (`Functions/`) are added to the MATLAB path at the start of scripts that require them via `addpath(fullfile(pwd, 'Functions'))`.


---

## Citation

If you use this code, please cite:

> [Authors]. ([Year]). [Title]. [Journal]. [DOI]

This code is partially based on:

> Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang J, Chun MM, Papademetris X & Constable RT. (2015). Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity. *Nature Neuroscience*, 18, 1664–1671.

---

## License

This code is released under the [MIT License](LICENSE).
Portions based on Finn & Shen (2015) are released under the GNU GPL v2 (see individual script headers).
