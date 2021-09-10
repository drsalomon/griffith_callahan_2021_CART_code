source("/windir/c/Users/redas/desktop/jupyter_directory/salomon_lab_folder/car_t_project/r_scripts/ssGSEA2.0.R")
gct_file <- "/windir/c/Users/redas/desktop/jupyter_directory/salomon_lab_folder/car_t_project/pepdep_data/tio2_data/ptmsea_output/raw/bcell/qvalue/qvalue_data.gct"
output_loc <- "/windir/c/Users/redas/desktop/jupyter_directory/salomon_lab_folder/car_t_project/pepdep_data/tio2_data/ptmsea_output/raw/bcell/qvalue/output"
sig_db <- "/windir/c/Users/redas/desktop/jupyter_directory/salomon_lab_folder/car_t_project/database/ptm.sig.db.all.flanking.human.v1.9.0.gmt"
signat.all <- unlist(lapply(sig_db, readLines))
signat.all <- strsplit(signat.all, '	')
names(signat.all) <- sapply(signat.all, function(x)x[1])
signat.all <- lapply(signat.all, function(x) x[-c(1,2)])
names(gct_file) <- paste(  sub('\\.gct$', '', sub('.*/','', gct_file)), 'ssGSEA', sep='_' )
input.ds <- gct_file
i <- 1
gsea.res <-
  ssGSEA2(
    "/windir/c/Users/redas/desktop/jupyter_directory/salomon_lab_folder/car_t_project/pepdep_data/tio2_data/ptmsea_output/raw/bcell/qvalue/qvalue_data.gct",
    gene.set.databases=sig_db,
    sample.norm.type="none",
    weight=1,
    statistic="area.under.RES",
    output.score.type="NES",
    nperm=10000,
    min.overlap=3,
    correl.type="z.score",
    output.prefix=output_loc,
    par=T,
    spare.cores=1,
    export.signat.gct=F,
    extended.output=T,
    param.file=T,
    global.fdr=TRUE,
  )
