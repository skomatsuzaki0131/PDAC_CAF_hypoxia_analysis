## Initial cluster annotation settings ---------------------------------------

List.Annot <- list(
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

Vec.Annot <- unlist(
  lapply(names(List.Annot), function(ct) {
    setNames(
      rep(ct, length(List.Annot[[ct]])),
      List.Annot[[ct]]
    )
  })
)