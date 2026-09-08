# Fresh 41-gene DINP–CRC single-cell localization

- Source object: `16023185-de21-4c0d-a9c8-73abdd52d142.h5ad`; shape `370,115 cells × 38,361 features`.
- Source title: All Cells
- Gene input: `query_41_genes.csv` (41 genes).
- Source coverage: **40/41** genes; missing: `MIR675`.
- Score genes: **40**; the score is the mean gene-wise z score over the observed score genes.
- Eligible scope: **223,160 cells**, **62 donors**, primary CRC/normal cells with mapped compartments.

## Compartment localization

```
compartment                group  n_cells  mean_program_score  median_program_score  sd_cell_program_score
endothelial colon adenocarcinoma     4912            0.156664              0.123634               0.316100
endothelial               normal     2608            0.043470              0.039336               0.199298
 epithelial colon adenocarcinoma   113216           -0.048571             -0.047940               0.181898
 epithelial               normal    55079           -0.058347             -0.077831               0.193787
 fibroblast colon adenocarcinoma     3031            0.388460              0.383741               0.283272
 fibroblast               normal     2200            0.333672              0.292705               0.300090
    myeloid colon adenocarcinoma    39167            0.153187              0.140935               0.261732
    myeloid               normal     2947           -0.027658             -0.020513               0.222504
```

## Donor-aware paired contrasts

```
                                            contrast      family                          target  n_paired_donors  mean_delta_tumor_minus_normal  median_delta_tumor_minus_normal  sd_delta  mean_delta_95ci_low  mean_delta_95ci_high  paired_t_statistic   paired_t_p  paired_wilcoxon_statistic  paired_wilcoxon_p direction  paired_t_BH_FDR  paired_wilcoxon_BH_FDR
                          epithelial tumor vs normal compartment                      epithelial               36                      -0.018559                        -0.030086  0.074259            -0.043685              0.006566           -1.499556 1.426965e-01                      243.0       1.615267e-01      down     1.902620e-01            1.884478e-01
                             myeloid tumor vs normal compartment                         myeloid               36                       0.100197                         0.105463  0.117569             0.060417              0.139977            5.113442 1.138456e-05                       75.0       1.328049e-05        up     3.035883e-05            4.648170e-05
                          fibroblast tumor vs normal compartment                      fibroblast                0                            NaN                              NaN       NaN                  NaN                   NaN                 NaN          NaN                        NaN                NaN      flat              NaN                     NaN
                         endothelial tumor vs normal compartment                     endothelial                0                            NaN                              NaN       NaN                  NaN                   NaN                 NaN          NaN                        NaN                NaN      flat              NaN                     NaN
source tumor-labeled epithelial vs normal epithelial     subtype source_tumor_labeled_epithelial               36                       0.003166                         0.015405  0.077839            -0.023171              0.029503            0.244053 8.086150e-01                      299.0       6.028147e-01        up     8.086150e-01            6.028147e-01
         other tumor epithelial vs normal epithelial     subtype          other_tumor_epithelial               36                      -0.040284                        -0.037142  0.093262            -0.071840             -0.008729           -2.591690 1.383839e-02                      187.0       2.096025e-02      down     2.214142e-02            2.934435e-02
                          macrophage tumor vs normal     subtype                      macrophage               35                       0.170227                         0.151274  0.114969             0.130734              0.209720            8.759599 3.090659e-10                        6.0       8.149073e-10        up     2.472527e-09            5.704351e-09
                            monocyte tumor vs normal     subtype                        monocyte               33                       0.115411                         0.105897  0.127264             0.070285              0.160537            5.209531 1.080206e-05                       59.0       1.999666e-05        up     3.035883e-05            4.665887e-05
                           dendritic tumor vs normal     subtype                       dendritic               34                       0.075712                         0.077788  0.122392             0.033007              0.118416            3.607035 1.010742e-03                      111.0       9.691622e-04        up     2.021485e-03            1.696034e-03
                         granulocyte tumor vs normal     subtype                     granulocyte                3                       0.073911                         0.079461  0.282457            -0.627752              0.775574            0.453229 6.948086e-01                        NaN                NaN        up     7.940670e-01                     NaN
```

## Interpretation boundary

This is a localization analysis of the fresh 41-gene program in a versioned CRC single-cell reference. It does not establish that DINP exposure causes the program, that the program mediates the epidemiologic association, or that source-labeled tumor epithelial cells are definitively malignant without independent CNV/malignancy validation. The absent MIR675 feature is not treated as a negative expression measurement.

## Reproducibility

The complete input/output hashes and source metadata are in `singlecell_41_manifest.json`; cell-level expression was read from the source H5AD without re-querying Census or reprocessing the whole object.
