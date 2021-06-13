**U**nified **H**ierarchical Cell **A**nnotation **F**ramework (uHAF)
========

# Introduction

We designed the structured vocabulary in macroscope level and the microscope level. The macroscope level annotations describe the anatomical information including the relationship between the human body, system, organ, the region and the subregion. The microscope level annotations describe the histological type (the “tissue_type”) and the cellular identities (the “cell_type”) that arise from the molecular heterogeneities of the cells. The macroscope annotation terms form composition ("part_of") relationships while the terms in the microscope level possess inheritance relationship ("is_a"). 


### Macroscope ontology

The macroscope level mainly containing anatomical information including system, organ, region, subregion.
- Human body
  - System
    - Organ
      - Suborgan 
> For example:
- Human body
  - Cardiovascular system
    - Heart
      - Atrium
        - Left atrium
        - Right atrium
-------


### Microscope ontology

The microscope level annotations describe the histological type (the “tissue_type”) and the cellular identities (the “cell_type”) that arise from the molecular heterogeneities of the cells.
- Cell
  - Tissue_type
    - Cell_type
> For example:
- Cell
  - Connective tissue
    - Hematocyte
      - Leukocyte
        - Lymphocyte
           - T cell
             - CD4 T cell
             - CD8 T cell
-------


# File table 
| File name | Description |
| --------- | ----------- |
| "uHAF macro ontology.owl" | Macroscopic ontology |
| "uHAF micro ontology.owl" | Microscopic ontology |
| "uHAF macro-micro map.csv" | The connection between macroscopic and microscopic ontology, only containing the observed cell type in hECA v1.0 |
| "uHAF marker reference.xlsx" | The marker reference list for UHAF annotation, each sheet is one organ |
| "uHAF-CL tree.csv" | The uHAF-CL tree only contains the observed cell type in hECA v1.0 |


# Visulization

Use [Protege](https://protege.stanford.edu/) to visulize the tree structure of UHAF.


# Reference

1. *柏树令，应大军 (2015). 系统解剖学, 3 edn (人民卫生出版社).*
2. *李继承，曾园山 (2018). 组织学与胚胎学, 9 edn (人民卫生出版社).*
3. *唐军民，张雷 (2013). 组织学与胚胎学, 3 edn (北京大学医学出版社).*
4. *Mescher, A.L. (2016). Junqueira’s Basic Histology Text and Atlas, 13 edn (McGraw-Hill).*

