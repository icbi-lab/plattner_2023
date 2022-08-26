import pandas as pd

# Cell type definitions

# map cell phenotypes as named in the phenotype data, to idx ordered as used
# in manuscript heatmaps
cellTypeMap = {
    "Tcells_undefined": 0,
    "Prol_Tcells": 1,
    "CD8+_Tcells": 2,
    "Intraepithelial_CD8+_Tcells": 3,
    "Intraepithelial_CD4+_PD1+_Tcells": 4,
    "CD4+_Tcells": 5,
    "CD8-_Tcells": 6,
    "Tregs": 7,
    "ICOS+_Tregs": 8,
    "CD57+_cells": 42,
    "CD45RO_undefined": 9,
    "CD11c+_cells": 10,
    "CD38+_cells": 43,
    "VISTA+_CD31+_CD38+_cells": 11,
    "ILC": 12,
    "CD163+_HLA-DR+_macrophages": 13,
    "HLA-DR+_macrophages": 14,
    "CD11c+_HLA-DR+_macrophages": 15,
    "Macrophages_undefined": 16,
    "HLA-DR+_monocytes": 17,
    "Monocytes": 18,
    "Granulocytes": 19,
    "Fibroblasts": 20,
    "HLA-DR+_fibroblasts": 21,
    "D2-40_fibroblasts": 22,
    "D2-40+_TGFb+_fibroblasts": 23,
    "TGFb+_tum": 24,
    "Bcat+_FOXP3+_PDL1+_tum": 25,
    "Bcat+_CD15+_tum": 26,
    "Bcat+_Prol_tum": 27,
    "Bcat+_HLA-DR+_Prol_tum": 28,
    "Bcat+_HLA-DR+_tum": 29,
    "Bcat+_TGFb+_tum": 30,
    "Bcat+_tum": 31,
    "Bcat+_HLA-DR+_TGFb+_tum": 32,
    "Bcat+_Prol_TGFb+_tum": 33,
    "Apoptotic_Bcat+_Prol_tum": 34,
    "Apoptotic_Bcat+_tum": 35,
    "Apoptotic_Bcat+_HLA-DR+_tum": 36,
    "Vessels": 37,
    "D2-40+_vessels": 38,
    "TGFb+_vessels": 39,
    "VISTA+_vessels": 40,
    "Unknown": 41,
}

# map cell phenotypes from phenotype data to names used in analysis and manuscript
cellTypeNameMap = {
    "Tcells_undefined": "Tcells undefined",
    "Prol_Tcells": "proliferating Tcells",
    "CD8+_Tcells": "CD8+ Tcells",
    "Intraepithelial_CD8+_Tcells": "Intraepithelial CD8+ Tcells",
    "Intraepithelial_CD4+_PD1+_Tcells": "Intraepithelial CD4+PD1+ Tcells",
    "CD4+_Tcells": "CD4+ Tcells",
    "CD8-_Tcells": "CD8- Tcells",
    "Tregs": "Tregs",
    "ICOS+_Tregs": "ICOS+ Tregs",
    "CD163+_HLA-DR+_macrophages": "CD163+HLA-DR+ macrophages",
    "HLA-DR+_macrophages": "HLA-DR+ macrophages",
    "CD11c+_HLA-DR+_macrophages": "CD11c+HLA-DR+ macrophages",
    "Macrophages_undefined": "Macrophages undefined",
    "HLA-DR+_monocytes": "HLA-DR+ monocytes",
    "Monocytes": "Monocytes",
    "Granulocytes": "Granulocytes",
    "Fibroblasts": "Fibroblasts",
    "HLA-DR+_fibroblasts": "HLA-DR+ fibroblasts",
    "D2-40_fibroblasts": "D2-40 fibroblasts",
    "D2-40+_TGFb+_fibroblasts": "D2-40+TGFb+ fibroblasts",
    "Bcat+_FOXP3+_PDL1+_tum": "Bcat+FOXP3+PDL1+ tumor cells",
    "TGFb+_tum": "TGFb+ tumor cells",
    "Bcat+_CD15+_tum": "Bcat+CD15+ tumor cells",
    "Bcat+_Prol_tum": "Bcat+ proliferating tumor cells",
    "Bcat+_HLA-DR+_Prol_tum": "Bcat+HLA-DR+ proliferating tumor cells",
    "Bcat+_HLA-DR+_tum": "Bcat+HLA-DR+ tumor cells",
    "Bcat+_TGFb+_tum": "Bcat+TGFb+ tumor cells",
    "Bcat+_tum": "Bcat+ tumor cells",
    "Bcat+_HLA-DR+_TGFb+_tum": "Bcat+HLA-DR+TGFb+ tumor cells",
    "Bcat+_Prol_TGFb+_tum": "proliferating Bcat+TGFb+ tumor cells",
    "Apoptotic_Bcat+_Prol_tum": "Apoptotic Bcat+ proliferating tumor cells",
    "Apoptotic_Bcat+_tum": "Apoptotic Bcat+ tumor cells",
    "Apoptotic_Bcat+_HLA-DR+_tum": "Apoptotic Bcat+HLA-DR+ tumor cells",
    "CD57+_cells": "CD57+ cells",
    "CD45RO_undefined": "CD45RO undefined",
    "CD11c+_cells": "CD11c+ cells",
    "CD38+_cells": "CD38+ cells",
    "VISTA+_CD31+_CD38+_cells": "VISTA+CD31+CD38+ cells",
    "ILC": "ILC",
    "Vessels": "Vessels",
    "D2-40+_vessels": "D2-40+ vessels",
    "TGFb+_vessels": "TGFb+ vessels",
    "VISTA+_vessels": "VISTA+ vessels",
    "Unknown": "Unknown",
}

# make sure we have the correct sorting, so apply same sorting as in cellTypeMap
cellTypeNameMap = {n: cellTypeNameMap[n] for n in cellTypeMap}

# assign a unique color to each cell type (colors as used in manuscript)
cellTypeColorMap = {
    "Tcells_undefined": "#53F328",
    "Prol_Tcells": "#62D243",
    "CD8+_Tcells": "#37FB03",
    "Intraepithelial_CD8+_Tcells": "#299C0A",
    "Intraepithelial_CD4+_PD1+_Tcells": "#71BD5D",
    "CD4+_Tcells": "#1D590D",
    "CD8-_Tcells": "#0F590D",
    "Tregs": "#478A44",
    "ICOS+_Tregs": "#043A02",
    "CD57+_cells": "#94FF33",
    "CD45RO_undefined": "#6AAB2F",
    "CD11c+_cells": "#5AA615",
    "CD38+_cells": "#4BD82E",
    "VISTA+_CD31+_CD38+_cells": "#1F820A",
    "ILC": "#04F90E",
    "CD163+_HLA-DR+_macrophages": "#0ACC55",
    "HLA-DR+_macrophages": "#13A14A",
    "CD11c+_HLA-DR+_macrophages": "#53D786",
    "Macrophages_undefined": "#53AF77",
    "HLA-DR+_monocytes": "#46FBC3",
    "Monocytes": "#73D5F1",
    "Granulocytes": "#2896E5",
    "Fibroblasts": "#9873F1",
    "HLA-DR+_fibroblasts": "#754CDA",
    "D2-40_fibroblasts": "#5E3DAF",
    "D2-40+_TGFb+_fibroblasts": "#062697",
    "TGFb+_tum": "#B61818",
    "Bcat+_FOXP3+_PDL1+_tum": "#FB1313",
    "Bcat+_CD15+_tum": "#BD2207",
    "Bcat+_Prol_tum": "#FF5004",
    "Bcat+_HLA-DR+_Prol_tum": "#F3842C",
    "Bcat+_HLA-DR+_tum": "#F99E0B",
    "Bcat+_TGFb+_tum": "#BB4B0B",
    "Bcat+_tum": "#F77B13",
    "Bcat+_HLA-DR+_TGFb+_tum": "#F7A413",
    "Bcat+_Prol_TGFb+_tum": "#F7C613",
    "Apoptotic_Bcat+_Prol_tum": "#C6A839",
    "Apoptotic_Bcat+_tum": "#C67939",
    "Apoptotic_Bcat+_HLA-DR+_tum": "#F39E0B",
    "Vessels": "#C6398C",
    "D2-40+_vessels": "#F70795",
    "TGFb+_vessels": "#9F0560",
    "VISTA+_vessels": "#B831F7",
    "Unknown": "#C2C2C2",
}

# In []:
cellTypeClassMap = {
    "Tcells_undefined": "lymphocytes",
    "Prol_Tcells": "lymphocytes",
    "CD8+_Tcells": "lymphocytes",
    "Intraepithelial_CD8+_Tcells": "lymphocytes",
    "Intraepithelial_CD4+_PD1+_Tcells": "lymphocytes",
    "CD4+_Tcells": "lymphocytes",
    "CD8-_Tcells": "lymphocytes",
    "Tregs": "lymphocytes",
    "ICOS+_Tregs": "lymphocytes",
    "CD57+_cells": "lymphocytes",
    "CD45RO_undefined": "lymphocytes",
    "CD11c+_cells": "macrophages",
    "CD38+_cells": "macrophages",
    "VISTA+_CD31+_CD38+_cells": "macrophages",
    "ILC": "macrophages",
    "CD163+_HLA-DR+_macrophages": "macrophages",
    "HLA-DR+_macrophages": "macrophages",
    "CD11c+_HLA-DR+_macrophages": "macrophages",
    "Macrophages_undefined": "macrophages",
    "HLA-DR+_monocytes": "macrophages",
    "Monocytes": "macrophages",
    "Granulocytes": "fibroblasts",
    "Fibroblasts": "fibroblasts",
    "HLA-DR+_fibroblasts": "fibroblasts",
    "D2-40_fibroblasts": "fibroblasts",
    "D2-40+_TGFb+_fibroblasts": "fibroblasts",
    "TGFb+_tum": "tumorcells",
    "Bcat+_FOXP3+_PDL1+_tum": "tumorcells",
    "Bcat+_CD15+_tum": "tumorcells",
    "Bcat+_Prol_tum": "tumorcells",
    "Bcat+_HLA-DR+_Prol_tum": "tumorcells",
    "Bcat+_HLA-DR+_tum": "tumorcells",
    "Bcat+_TGFb+_tum": "tumorcells",
    "Bcat+_tum": "tumorcells",
    "Bcat+_HLA-DR+_TGFb+_tum": "tumorcells",
    "Bcat+_Prol_TGFb+_tum": "tumorcells",
    "Apoptotic_Bcat+_Prol_tum": "tumorcells",
    "Apoptotic_Bcat+_tum": "tumorcells",
    "Apoptotic_Bcat+_HLA-DR+_tum": "tumorcells",
    "Vessels": "endothelial",
    "D2-40+_vessels": "endothelial",
    "TGFb+_vessels": "endothelial",
    "VISTA+_vessels": "endothelial",
    "Unknown": "endothelial",
}

cellTypeClassColors = {
    "lymphocytes": "#53F328",
    "macrophages": "#53AF77",
    "fibroblasts": "#9873F1",
    "tumorcells": "#F99E0B",
    "endothelial": "#9F0560",
}

# map cell types to the cell type call color
cellTypeClassColorMap = {n: cellTypeClassColors[cellTypeClassMap[n]] for n in cellTypeClassMap}

# In[4]:

# map cell type colors to cell type names used in manuscript
cellTypeNameColorMap = {cellTypeNameMap[n]: cellTypeColorMap[n] for n in cellTypeMap}


# In[5]:

# map cell type names as used in manuscript to cell type numbers (01, 02, 03 ...) as used in manuscript
cellTypeNameColorMapIdx = {
    cellTypeNameMap[n]: str(cellTypeMap[n] + 1).zfill(2) for n in cellTypeMap
}


# In[6]:

# map cell type numbers (01, 02, 03 ...) as used in manuscript to cell type colors used in manuscript
cellTypeIdxColorMap = {cellTypeNameColorMapIdx[x]: cellTypeNameColorMap[x] for x in cellTypeNameColorMap}


# In[7]:

# make pandas series
cellTypeIdxColorMap_s = pd.Series(cellTypeIdxColorMap)

# In []:

# map cell type numbers (01, 02, 03 ...) as used in manuscript to cell type colors used in manuscript
cellTypeClassIdxColorMap = {
    str(list(cellTypeColorMap.keys()).index(n) + 1).zfill(2): cellTypeClassColors[cellTypeClassMap[n]]
    for n in cellTypeClassMap
}
cellTypeDummyIdxColorMap = {str(list(cellTypeColorMap.keys()).index(n) + 1).zfill(2): "white" for n in cellTypeClassMap}

cellTypeClassIdxColorMap_s = pd.Series(cellTypeClassIdxColorMap)
cellTypeDummyIdxColorMap_s = pd.Series(cellTypeDummyIdxColorMap)

# In[8]:

# make pandas series
cellTypeColorMap_s = pd.Series(cellTypeNameColorMap)


# In[9]:

# swap key/value of custom_cellTypeMapOrder (custom order number is key and name is value now)
inv_cellTypeMapOrder = {cellTypeMap[x]: x for x in cellTypeMap}

# get custom sorted cell type numbers as used in manuscript
cellTypeMapOrder_idx = [str(i + 1).zfill(2) for i in list(dict(sorted(inv_cellTypeMapOrder.items())).keys())]

# In[10]:
# swap key/value of cellTypeNameMap (cellTypeNameMap name is key and data celltype name is value now)
inv_cellTypeNameMap = {cellTypeNameMap[x]: x for x in cellTypeNameMap}
