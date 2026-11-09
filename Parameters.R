# ============================================================
# Parameter file  •  Pre-lab + C-peptide module
# ============================================================

# ---------- 1.  Biomarker log-likelihood ratios  ----------
# Provide LR(+) and LR(–) for each auto-antibody, then convert
# to natural-log space (base-e) for use in additive models.
#
# GAD: from the PUBLISHED meta-analysis (19 studies) — Front Endocrinol 2026,
#   DOI 10.3389/fendo.2026.1771950; pooled sens 0.531, spec 0.926
#   => LR+ 0.531/(1-0.926) = 7.18 ; LR- (1-0.531)/0.926 = 0.51
#   These MUST match Table 3 of the manuscript. Superseded values, do not reuse:
#     9.15 / 0.484 = the 2025 DTM poster (21 studies)
#     8.9  / 0.47  = an intermediate pre-publication draft (19 studies)

lr_table <- data.frame(
  biomarker = c("GAD", "IA2", "ZnT8", "IAA"),
  LR_pos    = c(  7.18,  ## <<< EDIT HERE  (LR+ for GAD)
                  21.4,  ## <<< EDIT HERE  (LR+ for IA-2)
                  11.7,  ## <<< EDIT HERE  (LR+ for ZnT8)
                  5  ## <<< EDIT HERE  (LR+ for IAA judgement-based)
  ),
  LR_neg    = c( 0.51,  ## <<< EDIT HERE  (LR– for GAD)
                 0.626, ## <<< EDIT HERE  (LR– for IA-2)
                 0.759, ## <<< EDIT HERE  (LR– for ZnT8)
                 0.92  ## <<< EDIT HERE  (LR– for IAA, judgement-based)
  )
)

lr_table$logLR_pos <- log(lr_table$LR_pos)
lr_table$logLR_neg <- log(lr_table$LR_neg)

# ---------- 2.  Pooled log-likelihood ratios (demographic/clinical)  ----------
# Log likelihood ratios for male sex and PARENTAL history of T2D, pooled across
# the regions that measured each term (inverse-variance, common effect).
#
# Revised 2026-07-28 for the BMJ R1 revision. Previously -0.03 / -0.553, taken
# from a brms random-effects fit (Model/Model.R) that also pooled the Scania
# ANDIS parental-history coefficient.  ANDIS recorded FAMILY history of diabetes
# in ANY relative, not parental history -- a different construct.  It was the
# sole source of between-region heterogeneity (Q = 110.8, p < 1e-4 with it;
# Q = 0.06, p = 0.81 without), and carried 33% of the inverse-variance weight,
# biasing the pooled estimate toward the null.  Exeter never reported either
# coefficient, so only Mohan Clinics and Ukraine T1D Exomes contribute to fh.
#   fh:       Mohan -0.850 (0.046), Ukraine -0.894 (0.180)  -> -0.853 (0.045)
#   male_sex: ANDIS +0.046 (0.063), Mohan -0.108 (0.049),
#             Ukraine -0.172 (0.190)                        -> -0.055 (0.038)
# See Model/R1_repool_parental_history.R in the manuscript folder.
pooled_coef <- list(
  male_sex = -0.055, ## <<< EDIT HERE (log LR for male sex)      Table 2 prints -0.05
  fh = -0.853        ## <<< EDIT HERE (log LR for parental T2D history)  prints -0.85
)

# ---------- 3.  C-peptide gamma distribution parameters  ----------
# MODIFY THESE PARAMETERS TO CHANGE C-PEPTIDE LIKELIHOOD RATIOS (fasting only)
cpeptide_params <- list(
  stimulated = list(
    t1d = list(alpha = NA, lambda = NA),
    t2d = list(alpha = NA, lambda = NA)
  ),
  fasting = list(
    # Gamma parameters derived from the South Asia fasting C-peptide means/SDs:
    #   alpha = (mean/SD)^2 ; lambda = mean/SD^2
    #   T1D  mean 0.327, SD 0.197 -> 2.7553, 8.4259
    #   T2D  mean 1.059, SD 0.468 -> 5.1204, 4.8351
    # These MUST match Table 4 of the manuscript. Superseded: t2d alpha = 5.13.
    t1d = list(alpha = 2.76, lambda = 8.43),
    t2d = list(alpha = 5.12, lambda = 4.84)
  )
)

# ---------- 4.  Bundle for downstream use ----------
params <- list(
  lr = lr_table,
  pooled = pooled_coef,
  cpeptide = cpeptide_params
)

