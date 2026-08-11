## Initial cluster annotation settings ---------------------------------------

initial_cluster_annotation_list <- list(
  Acinar = c(0, 13, 24),
  Ductal_like_acinar = c(33, 1, 12, 17, 30),
  Normal_ductal = 23,
  PanIN = c(10, 27),
  PDAC = 20,
  Islet = 14,
  CAF = c(9, 6, 7),
  Mural = 11,
  Endothelial = c(34, 3),
  Lymph_T = c(2, 16),
  Lymph_B = 15,
  Plasma = c(18, 37),
  Myeloid = c(5, 4, 31),
  Mast = 26,
  Nerve = 29,
  Unclassified = c(8, 19, 21, 22, 25, 28, 32, 35, 36)
)

initial_cluster_annotation <- unlist(
  lapply(names(initial_cluster_annotation_list), function(ct) {
    setNames(
      rep(ct, length(initial_cluster_annotation_list[[ct]])),
      initial_cluster_annotation_list[[ct]]
    )
  })
)


## CAF subcluster settings ---------------------------------------------------

caf_subcluster_order <- c(4, 8, 0, 1, 2, 5, 3, 7, 6)


## Spatial niche annotation settings -----------------------------------------

niche_order <- paste0(
  "niche_",
  c(2, 3, 8, 9, 5, 6, 1, 4, 7, 10)
)

niche_annotation <- c(
  "niche_1"  = "N7. Plasma cell niche",
  "niche_2"  = "N1. Normal acinar niche",
  "niche_3"  = "N2. Atrophic acinar niche",
  "niche_4"  = "N8. Lymphoid niche",
  "niche_5"  = "N5. Islet niche",
  "niche_6"  = "N6. Vascular niche",
  "niche_7"  = "N9. Lymphoid–myeloid niche",
  "niche_8"  = "N3. Normal–PanIN niche",
  "niche_9"  = "N4. Cancer niche",
  "niche_10" = "N10. Neuronal niche"
)

normal_acinar_niche <- "niche_2"