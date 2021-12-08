# cibersort的简单封装

### 文件结构
ciberSortEnc(cibersort封装)
|
|----depdence
|     |----CIBERSORT3.R
|     |----CIBERSORT4.R
|     |----LM22.txt
|
|----backup
|     |----CIBERSORT3_copy.R
|     |----CIBERSORT4_copy.R
|     |----LM22_copy.txt
|     |----LM22_source.xls
|
|----script
|     |----doMetaBric.R
|
|----src
|     |----DATACIBER.R
|
|----README.md

### 使用简介
 ##### depdence
 CIBERSORT3.R 与 CIBERSORT4.R 是cibersort的两个版本，4版本比较新。LM22是22种免疫细胞的特征基因矩阵，也可以自制特征基因矩阵进行处理。

 ##### backup
文件夹中的文件是depdence中的备份，您有需要时可以直接在cibersort源码中添加注释或日志等语句，如果搞出什么问题，难以恢复就可以从backup中复制一份到depdence替换原本的代码。其中 LM22_source.xls 是我cibersort文章链接中下载的，如果您想自行创建LM22矩阵的txt文件，可以复制表格部分粘贴到txt文件中。

 ##### script
doMetaBric.R 是以metabric数据集为例的一个脚本文件，可以在PC手动运行或在服务器以下面的方式进行后台运行，运算结果会被持久化到文件中，并且可以在log.txt文件中查看日志（强烈建议您根据自己的需求添加或更改日志）
```bash

nohup Rscript <脚本路径> &

# nohup Rscript doMetaBric.R &
#可以在脚本所在文件夹下运行上述语句调用，关于nohup工具可以参考

https://www.runoob.com/linux/linux-comm-nohup.html

```

 ##### src
对cibersort的使用进行了封装。
DATACIBER.R中只有两个方法，

pre_process是预处理函数，功能如下
       1、基因去重，多个基因重名时保留表达量大的基因，
            表达量越小，越有可能是噪声，保留平均表达量最大（行平均值最大）的一行

       2、去除不存在gene_id的基因

       3、去除不再需要的gene_id列

    该方法有两个参数
      @param file_path 待处理文件的路径
      @param result_data_path 计算结果写入路径


do_cibersort是调用cibersort进行反卷积运算的函数，功能如下
       1、计算样本免疫细胞矩阵

       2、将结果矩阵写入传入的路径

    该方法有两个参数
       @Param file_path 待处理数据的路径
       @param result_data_path 指定文件写入路径

预处理方法针对不同数据集可能需要您手动进行一些调整，我在测试时使用的数据集市metabric数据集，下载地址为
```bash

https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz

```
通过数据集中的说明文件您可以大概了解该数据集的存放形式，由于说明文件过于简短，非常容易忽略。
  