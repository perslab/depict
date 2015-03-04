#! /bin/bash

/home/tools/plink/plink_v1-90_stable_beta_3f_2-Mar/plink --bfile /home/data/1000G/data/phase1/bed_CEU_GBR_TSI_unrelated/CEU_GBR_TSI_unrelated.phase1_release_v3.20101123.snps_indels_svs.genotypes --extract /tmp/ldl_teslovich_nature2010.rsID --make-bed --out CEU_GBR_TSI_unrelated.phase1_release_v3.20101123.snps_indels_svs.genotypes_ldl_teslovich_nature2010
