## Helper functions ----------------------------------------------------------

clean_gs_label_ForFig <- function(x) {
  x %>%
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

clean_gs_label_ForTable <- function(x) {
  x %>%
    stringi::stri_trans_general("Any-ASCII") %>%
    str_replace_all("[–—−]", "-") %>%
    str_remove("^(HALLMARK_|KEGG_)") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower() %>%
    str_to_sentence() %>%
    str_replace_all("\\bTnfa\\b", "TNF-alpha") %>%
    str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
    str_replace_all("\\bTgf beta\\b", "TGF-beta") %>%
    str_replace_all("\\bMycaf\\b", "myCAF") %>%
    str_replace_all("\\bKras\\b", "KRAS") %>%
    str_replace_all("\\bP53\\b", "p53") %>%
    str_replace_all("\\bIcaf\\b", "iCAF") %>%
    str_replace_all("\\bDna\\b", "DNA") %>%
    str_replace_all("\\bUv\\b", "UV") %>%
    str_replace_all("\\bMyc\\b", "MYC") %>%
    str_replace_all("\\bWnt beta catenin\\b", "Wnt/beta-catenin") %>%
    str_replace_all("\\bdn\\b", "down") %>%
    str_replace_all("\\bnfkb\\b", "NF-kappaB") %>%
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
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}