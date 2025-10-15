# Define mapping of class to TCGA
CLASS_TO_TCGA <- list(
  "Hyper_Colon"       = "TCGA_COAD",
  "Hyper_Lung_LUAD"   = "TCGA-LUAD",
  "Hyper_Lung_LUSC"   = "TCGA-LUSC",
  "Hyper_Breast"      = "TCGA-BRCA",
  "Hyper_Prostate"    = "TCGA-PRAD",
  "Hyper_Melanoma"    = "TCGA-SKCM",
  "Hyper_Mesotelioma" = "TCGA-MESO",
  "Hyper_Urothelial"  = "TCGA-BLCA"
)


# Define mapping of class to frequency and file path
CLASS_PARAMS_WGS <- list(
  Colon          = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Colorectal_Carcinoma_NCIT_C2955.tsv")),
  Lung           = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Lung_Carcinoma_NCIT_C4878.tsv")),
  Prostate       = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Prostate_Carcinoma_NCIT_C4863.tsv")),
  Urothelial     = list(freq = 22.5,file = paste0(PATH_INITIAL, "data/progenetix/Urothelial_Carcinoma_NCIT_C4030.tsv")),
  EwS            = list(freq = 10,  file = paste0(PATH_INITIAL, "data/progenetix/Ewing_Sarcoma_NCIT_C4817.tsv")),
  Mesotelioma    = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Malignant_Mesothelioma_NCIT_C4456.tsv")),
  Melanoma       = list(freq = 20,  file = paste0(PATH_INITIAL, "data/progenetix/Melanoma_NCIT_C3224.tsv")),
  Breast         = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Breast_Carcinoma_NCIT_C4872.tsv")),
  MM             = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Multiple_Myeloma_NCIT_C3242.tsv")),
  Pancreatic     = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Pancreatic_Adenocarcinoma_NCIT_C8294.tsv"))
)