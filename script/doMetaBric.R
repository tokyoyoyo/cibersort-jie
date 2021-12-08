source("/home/Zyc0/zjj/Scrip/CibersortTool/DATACIBER.R")

file_need_pre <- "/home/Zyc0/zjj/Data/UnArchive/brca_metabric/data_expression_median.txt"

file_preProcessed <- "/home/Zyc0/zjj/Data/PreProcessed/metabric_mRNA_preProcessed.txt"

matrix_file <- "/home/Zyc0/zjj/Data/CellMatrix/metabric_cell_matrix.txt"

pre_process(file_need_pre, file_preProcessed)

do_cibersort(file_preProcessed, matrix_file)