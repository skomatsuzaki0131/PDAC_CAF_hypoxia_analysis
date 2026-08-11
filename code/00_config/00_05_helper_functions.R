## Helper functions ----------------------------------------------------------

suppressPackageStartupMessages({
  library(stringr)
  library(stringi)
  library(magrittr)
})

format_gene_set_label <- function(x) {
  x %>%
    stringr::str_remove("_ssGSEA$") %>%
    str_remove("^(HALLMARK_|KEGG_)") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower() %>%
    str_to_sentence() %>%
    str_replace_all("\\bTnfa\\b", "TNF-α") %>%
    str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
    str_replace_all("\\bTgf beta\\b", "TGF-β") %>%
    str_replace_all("\\bMycaf\\b", "myCAF") %>%
    str_replace_all("\\bKras\\b", "KRAS") %>%
    str_replace_all("\\bP53\\b", "p53") %>%
    str_replace_all("\\bIcaf\\b", "iCAF") %>%
    str_replace_all("\\bDna\\b", "DNA") %>%
    str_replace_all("\\bUv\\b", "UV") %>%
    str_replace_all("\\bMyc\\b", "MYC") %>%
    str_replace_all("\\bWnt beta catenin\\b", "Wnt/β-catenin") %>%
    str_replace_all("\\bdn\\b", "down") %>%
    str_replace_all("\\bnfkb\\b", "NF-κB") %>%
    str_replace_all("\\bIl2 stat5\\b", "IL-2/STAT5") %>%
    str_replace_all("\\bG2m\\b", "G2/M") %>%
    str_replace_all("\\bE2f\\b", "E2F") %>%
    str_replace_all("\\bMtorc1\\b", "mTORC1") %>%
    str_replace_all("\\bPi3k akt mtor\\b", "PI3K-AKT-mTOR") %>%
    str_replace_all("targets v", "targets V") %>%
    str_replace_all("orig\\b", " hypoxia") %>%
    str_replace_all("\\biCAF signature\\b", "iCAF") %>%
    str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    str_replace_all("\\bInterferon alpha\\b", "IFN-α") %>%
    str_replace_all("\\bInterferon gamma\\b", "IFN-γ") %>%
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}

format_gene_set_label_for_table <- function(x) {
  x %>%
    stringr::str_remove("_ssGSEA$") %>%
    stringi::stri_trans_general("Any-ASCII") %>%
    stringr::str_replace_all("[–—−]", "-") %>%
    stringr::str_remove("^(HALLMARK_|KEGG_)") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence() %>%
    stringr::str_replace_all("\\bTnfa\\b", "TNF-alpha") %>%
    stringr::str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
    stringr::str_replace_all("\\bTgf beta\\b", "TGF-beta") %>%
    stringr::str_replace_all("\\bMycaf\\b", "myCAF") %>%
    stringr::str_replace_all("\\bKras\\b", "KRAS") %>%
    stringr::str_replace_all("\\bP53\\b", "p53") %>%
    stringr::str_replace_all("\\bIcaf\\b", "iCAF") %>%
    stringr::str_replace_all("\\bDna\\b", "DNA") %>%
    stringr::str_replace_all("\\bUv\\b", "UV") %>%
    stringr::str_replace_all("\\bMyc\\b", "MYC") %>%
    stringr::str_replace_all("\\bWnt beta catenin\\b", "Wnt/beta-catenin") %>%
    stringr::str_replace_all("\\bdn\\b", "down") %>%
    stringr::str_replace_all("\\bnfkb\\b", "NF-kappaB") %>%
    stringr::str_replace_all("\\bIl2 stat5\\b", "IL-2/STAT5") %>%
    stringr::str_replace_all("\\bG2m\\b", "G2/M") %>%
    stringr::str_replace_all("\\bE2f\\b", "E2F") %>%
    stringr::str_replace_all("\\bMtorc1\\b", "mTORC1") %>%
    stringr::str_replace_all("\\bPi3k akt mtor\\b", "PI3K-AKT-mTOR") %>%
    stringr::str_replace_all("targets v", "targets V") %>%
    stringr::str_replace_all("orig\\b", " hypoxia") %>%
    stringr::str_replace_all("\\biCAF signature\\b", "iCAF") %>%
    stringr::str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    stringr::str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    stringr::str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    stringr::str_replace_all("Epithelial mesenchymal transition", "EMT")
}