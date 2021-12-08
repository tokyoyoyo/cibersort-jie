#导入依赖的包

library("tibble")





#预处理函数，功能如下
#       1、基因去重，多个基因重名时保留表达量大的基因，
#          表达量越小，越有可能是噪声，保留平均表达量最大（行平均值最大）的一行
#
#       2、去除不存在gene_id的基因
#
#       3、去除不再需要的gene_id列
#
#'      @param file_path 待处理文件的路径
#'      @param result_data_path 计算结果写入路径
#
pre_process <- function(file_path, result_data_path) {

  cat("\n\n\n=========================",
      paste(strptime(date(),'%a %b %d %H:%M:%S %Y') ),
      "  数据预处理程序启动=======================\n\n",
      file = "log.txt", append = TRUE)

  mrna_data <- read.table(file_path, header = T, sep = "\t", check.names = F)
  #读取下载的数据，不考虑行名

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "     数据读取完成，开始处理\n",
      file = "log.txt", append = TRUE)


  mrna_data <- mrna_data[!is.na(mrna_data[, "Entrez_Gene_Id"]), ]
  mrna_data <- mrna_data[mrna_data[, "Entrez_Gene_Id"] != "", ]
  #去除gene_id为Na或空字符串的数据，只保留确实存在gene_id的数据

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  已删除gene_id字段不合法基因数据\n",
      file = "log.txt", append = TRUE)


  colnames(mrna_data)[1] <- "Gene symbol"
  #cibersort处理的数据的基因名所在列需要命名为"Gene symbol"

  mrna_data <- mrna_data[, -2]
  #删掉不再需要的基因id列

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  已删除gene_id列\n", file = "log.txt", append = TRUE)


  index <- order(rowMeans(mrna_data[, -1]), decreasing = T)
  #计算行平均值，按降序排列，获得降序的索引列表

  mrna_data <- mrna_data[index, ]
  #按降序的行序列order调整数据集的基因顺序

  keep <- !duplicated(mrna_data$"Gene symbol")
  #选择保留的行，对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
  #理由如下
  #表达量越大越有可能是我们期盼的那个基因，因为在生物学研究中，
  #首先关注到的往往是表达量大的基因

  mrna_data <- mrna_data[keep, ]
  #得到最后处理之后的表达谱矩阵

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  根据基因平均表达量去重\n\n",
      "开始保存预处理结果\n",file = "log.txt", append = TRUE)

  
  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "数据中含有 ", sum(is.na(mrna_data))," 个NA值。\n",
      "  去除含有NA值的基因\n\n"
      ,file = "log.txt", append = TRUE)
  
  mrna_data <-na.omit(mrna_data)
  #去除含有空值的行
  
  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      " 去除含有NA值的基因\n\n",
      file = "log.txt", append = TRUE)

  write.table(mrna_data,
                file = result_data_path,
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
  #将格式化的数据保存为cibersort能够接受的数据

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  ==========数据已完全写入", result_data_path, "==========\n",
      file = "log.txt", append = TRUE)


  return(result_data_path)
  #返回经过处理的数据写入的文件路径
}



#CIBERSORT函数，功能如下
#       1、计算样本免疫细胞矩阵
#
#       2、将结果矩阵写入传入的路径
#
#'       @Param file_path 待处理数据的路径
#'       @param result_data_path 指定文件写入路径
do_cibersort <- function(file_path, result_data_path) {

  source("../depdence/CIBERSORT4.R")
  #导入CIBERSORT脚本,在此通过相对路径导入，如果导入失败就试下用绝对路径

  cat("\n\n\n=========================",
        paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
        "  CIBERSORT程序启动=======================\n\n",
        file = "log.txt", append = TRUE)


  lm22 <- "/home/Zyc0/zjj/Scrip/CibersortTool/LM22.txt"
  lm22 <- read.table(lm22,
                      header = T,
                      sep = "\t",
                      row.names = 1,
                      check.names = F)
  #读取lm22免疫细胞基因表达特征矩阵

cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
"  读取lm22免疫细胞基因表达特征矩\n", file = "log.txt", append = TRUE)


  data_for_cibersort <- read.table(file_path,
                                  header = T,
                                  sep = "\t",
                                  row.names = 1,
                                  check.names = F)
  #读取待处理的mRNA表达矩阵

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  读取待处理的mRNA表达矩阵： ",
      file_path, "\n\n", file = "log.txt", append = TRUE)
  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  基因交集数为： ",
      length(intersect(rownames(lm22), rownames(data_for_cibersort))), "\n\n",
      file = "log.txt", append = TRUE)


cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  ==========开始反卷积mRNA数据计算免疫细胞矩阵==========\n",
      file = "log.txt", append = TRUE)

  cell_matrix <- CIBERSORT(lm22,
                  data_for_cibersort,
                  perm = 1000,
                  QN = TRUE,
                  absolute = FALSE)
  #CIBERSORT运算

cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  ==========计算完成,开始写入数据==========\n",
      file = "log.txt", append = TRUE)

  write.table(cell_matrix,
              file = result_data_path,
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE,
              quote = FALSE)
  #将计算出的免疫细胞数据保存下来

  cat(paste(strptime(date(),'%a %b %d %H:%M:%S %Y')),
      "  ==========数据已完全写入",
      result_data_path, "==========\n",
      file = "log.txt", append = TRUE)
}