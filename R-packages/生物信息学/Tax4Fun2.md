# 基本原理
16S rRNA基因测序已经成为研究微生物群落及其对各种生态系统中不断变化的环境条件的反应的有力技术。Wemheuer等开发了Tax4Fun2包，用于从16S rRNA基因序列中预测原核生物群落的功能谱和功能基因冗余。

Tax4Fun2以R包的形式发布，当前最新的版本是1.1.5（作者已删库跑路）。

Tax4Fun2提供了两个参考数据库（Ref99NR和Ref100NR），同Ref100NR相比，Ref99NR体积更小，对硬件的要求更低，预测速度更快。为了获取默认的参考数据库，可以使用`buildReferenceData`函数自动下载和构建参考数据。

还可以通过自己提供基因组数据来自定义构建参考数据集，核糖体RNA序列（16S rRNA或18S rRNA）通过BLAST搜索使用SILVA SSURef 132数据库（以uclust）进行鉴定。功能配置文件由BLASTp与Diamond对KEGG KO数据库生成。使用prodigal 2.6.3在功能注释之前预测蛋白质序列。目前功能注释仅适用于原核生物。提取的rRNA序列和功能图谱随后使用`addUserDataByClustering`或`addUserData`函数构建参考数据集。`addUserDataByClustering`函数需要使用vsearch中的功能，当参考数据较少时应使用`addUserData`函数。

![Tax4Fun2工作流程](/R-packages/生物信息学/plots/Tax4Fun2工作流程.webp)

# 使用方法
## 1 安装 Tax4Fun2 包，构建默认参考数据并安装所有依赖项
### 1.1 设置工作目录
```R
setwd("D:/16sdata/Tax4Fun2_test")
getwd()
```

### 1.2 将Tax4Fun2安装包下载在工作目录下，安装Tax4Fun2
原下载连接已失效（作者删库跑路？），可以在我的GitHub中[下载](https://github.com/ZihaoShu/Tax4Fun2/raw/main/Tax4Fun2_1.1.5.tar.gz)。

所有的官方参考数据均可在我的GitHub中[下载](https://github.com/ZihaoShu/Tax4Fun2)。

```R
install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)
library(Tax4Fun2)
```
### 1.3 构建默认数据库
```R
buildReferenceData(path_to_working_directory = ".", 
                   use_force = FALSE, 
                   install_suggested_packages = TRUE)
```
### 1.4 下载测试数据
```R
getExampleData(path_to_working_directory = ".")
```
### 1.5安装依赖
进入NCBI官方网站下载最新版本的[blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-win64.exe)，截至2021/11/4已更新到v2.12.0，
将二进制文件放置在由buildReferenceData() ['Tax4Fun2_ReferenceData_v2'] 生成的文件夹。

## 2 生成自己的参考数据集
### 2.1 提取单个基因组中的SSU序列（16S rRNA 和 18S rRNA）
```R
extractSSU(genome_file = "OneProkaryoticGenome.fasta", 
           file_extension = "fasta", 
           path_to_reference_data = "Tax4Fun2_ReferenceData_v2")
```
### 2.2 为单个原核基因组分配功能
```R
assignFunction(genome_file = "OneProkaryoticGenome.fasta", 
               file_extension = "fasta", 
               path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
               num_of_threads = 1, fast = TRUE)
```

### 2.3 生成参考数据
```
generateUserData(path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                 path_to_user_data = ".", 
                 name_of_user_data = "User_Ref0", 
                 SSU_file_extension = "_16SrRNA.ffn", 
                 KEGG_file_extension = "_funPro.txt")
```

## 3 进行功能预测
### 3.1 仅使用默认参考数据进行功能预测
#### 3.1.2 运行RefBlast工具参考Ref99NR进行序列同源比对
```R
runRefBlast(path_to_otus = "KELP_otus.fasta", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
            path_to_temp_folder = "Kelp_Ref99NR", 
            database_mode = "Ref99NR", 
            use_force = T, num_threads = 6)
```
### 3.1.3 功能注释
```R
makeFunctionalPrediction(path_to_otu_table = "KELP_otu_table.txt", 
                         path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                         path_to_temp_folder = "Kelp_Ref99NR", 
                         database_mode = "Ref99NR", 
                         normalize_by_copy_number = TRUE, 
                         min_identity_to_reference = 0.97, 
                         normalize_pathways = FALSE)
```
### 3.2 使用默认数据库和用户生成的数据库（非集群）进行功能预测
#### 3.2.1 生成用户数据库
```R
generateUserData(path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                 path_to_user_data = "KELP_UserData", 
                 name_of_user_data = "KELP1", 
                 SSU_file_extension = ".ffn", 
                 KEGG_file_extension = ".txt")
```
#### 3.2.2 运行RefBlast工具
```R
runRefBlast(path_to_otus = "KELP_otus.fasta", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
            path_to_temp_folder = "Kelp_Ref99NR_withUser1", 
            database_mode = "Ref99NR", 
            use_force = T, num_threads = 6, include_user_data = T, 
            path_to_user_data = "KELP_UserData", 
            name_of_user_data = "KELP1")
```
#### 3.2.3 功能注释
```R
makeFunctionalPrediction(path_to_otu_table = "KELP_otu_table.txt", 
                         path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                         path_to_temp_folder = "Kelp_Ref99NR_withUser1", 
                         database_mode = "Ref99NR", 
                         normalize_by_copy_number = T, 
                         min_identity_to_reference = 0.97, 
                         normalize_pathways = F, include_user_data = T, 
                         path_to_user_data = "KELP_UserData", 
                         name_of_user_data = "KELP1")
```

## 4 计算（多）功能冗余指数（实验性）
计算 KEGG 函数的系统发育分布（高 FRI -> 高冗余，低 FRI -> 函数冗余较少，可能会随着社区变化而丢失）
### 4.1 运行RefBlast工具
```R
runRefBlast(path_to_otus = "Water_otus.fna", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
            path_to_temp_folder = "Water_Ref99NR", 
            database_mode = "Ref99NR", 
            use_force = T, num_threads = 6)
```
### 4.2 计算FRIs
```R
calculateFunctionalRedundancy(path_to_otu_table = "Water_otu_table.txt", 
                              path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                              path_to_temp_folder = "Water_Ref99NR", 
                              database_mode = "Ref99NR", 
                              min_identity_to_reference = 0.97)
```

---
# 参考资料
[Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences](https://doi.org/10.1186/s40793-020-00358-7)