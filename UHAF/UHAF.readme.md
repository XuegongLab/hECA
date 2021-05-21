**U**nified **H**ierarchical Cell **A**nnotation **F**ramework (UHAF)
========

# Introduction

We designed the structured vocabulary in macroscope level and the microscope level. The macroscope level annotations describe the anatomical information including the relationship between the human body, system, organ, the region and the subregion. The microscope level annotations describe the histological type (the “tissue_type”) and the cellular identities (the “cell_type”) that arise from the molecular heterogeneities of the cells. The macroscope annotation terms form composition ("part_of") relationships while the terms in the microscope level possess inheritance relationship ("is_a"). 


### Macroscope ontology

The macroscope level mainly containing anatomical information including system, organ, region, subregion:
- Human body
  - System
    - Organ
      - Region
        - Subregion 
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


# Visulization

Use [Protege](https://protege.stanford.edu/) to visulize the tree structure of UHAF.


# Reference

1. *Junqueira's Basic Histology text & atlas (13th EDITON).* 
2. *唐君民, 张雷等. 组织学与胚胎学(第3版)[M]. 北京: 北京大学医学出版社, 2013.* 
3. *柏树令, 应大军等. 系统解剖学(第3版)[M]. 北京: 人民卫生出版社, 2015.* 
4. *《组织学与胚胎学》人卫第9版教材--高清彩色*

