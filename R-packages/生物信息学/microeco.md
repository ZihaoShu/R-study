# 基本原理
高通量测序技术在微生物生态学中的应用产生了大量的测序数据。微生物群落生态学的中的测序数据可以任意分为生物信息学分析和统计分析。生物信息学分析时一项典型的计算密集型工作。基于扩增子测序数据获得必要的操作分类单元（OTU）表、分类分配和系统发育树后，下游数据分析趋于多样化、复杂和耗时。

R语言中除了经典的生态包vegan外，还有一些包用来执行微生物群落数据的复杂分析，如phyloseq、microbiome等。

microeco基于R6开发，提供了封装的面向对象编程范式。为分析方法的每个部分创建类，使包框架模块化、清晰和简洁。整合了一些常用的方法，使包支持广泛的微生物群落分析，如LEfSe、RDA、共现网络分析、功能预测和空模型分析等。

![microeco工作流程](/R-packages/生物信息学/plots/microeco工作流程.jpeg)

# 使用方法
## microtable
是为创建储存输入文件对象和后续数据预处理而设计的基本类，数据预处理的内置函数包括
- `tidy_dataset()` - 修剪对象中的所有基本文件并使其信息一致
- `tidy_taxonomy()` - 使分类表具有统一的格式
- `rarefy_samples()` - 当样本间的序列数差异很大，可以使用该函数使得每个样本的序列号相等，以减少测序深度对$\alpha$-和$\beta$-多样性计算的影响。
- `merge_samples()` - 根据特定的一组合并样本
- `merge_taxa()` - 根据特定的分类等级合并样本
- `cal_abund()` - 可以自动计算每个分类等级的分类单元丰度并返回所包含所有表的列表（taxa_abund）。
- `cal_alphadiv()` - 可用于计算alpha多样性，包括Chao1、Shannon-Wiener指数、Simpson指数和系统发育alpha多样性。
- `cal_betadiv()` - 可用于获取beta多样性的距离矩阵，包括Bray-Curtis、Jaccard和UniFrac。

## trans_aund和trans_venn
对样本数据进行可视化
- `trans_abund` - 可以用来实现条形图、箱线图、热图和饼图的绘制。
- `trans_venn` - 用于韦恩图的绘制，可以对超过五个组进行分析。
- `trans_venn_com()` - 用来分析每个部分的分类群组成，该函数可以将每个分数中的OTU组成转换成群落组成表，并返回一个新的微表对象，用于快速绘制每个分数的分类组成。

## trans_alpha
- `trans_alpha` - 将microtable对象中的alpha多样性表转换为其他格式，用于统计和绘图。
- `cal_diff()` - 可以通过Kruskal-Wallis秩和检验或ANOVA（方差分析）对所有 alpha 多样性指标进行多重比较来执行组间差异检验。
- `plot_alpha()` - 用于显示alpha多样性数据，并通过多重比较或配对比较直接添加显着性。

## trans_beta
- `trans_beta` - 目前实现了几种常用的无约束排序方法，例如主成分分析、主坐标分析 (PCoA) 和非度量多维缩放。
- `cal_group_distance()` - 可以分别用于组距离变换。
- `plot_group_distance()` - 可用于绘图。
- `plot_clustering()` - 实现聚类分析和绘图。
- `cal_manova()` - 用于置换多变量方差分析（perMANOVA）的整体比较、配对组比较和多因素比较。

## trans_diff
指示分类单元的识别对于解释不同群体之间群落结构差异的生物学机制具有重要意义。随着微生物群落生态学中测序数据复杂性的增加，评估组间差异显着的分类群是一项挑战。
- `trans_diff` - 目前实现了三种众所周知的方法：LEfSe、随机森林和metastat。

## trans_envs
评估环境因素对微生物群落结构的影响对于推断支配群落组装的潜在机制至关重要。
- `trans_env` - 数据转化。
- `cal_rda()` - 实现冗余分析(RDA)和基于距离的RDA(db-RDA)。
- `trans_rda()` - 加强转化。
- `plot_rda()` - 加强绘图方法。
- `cal_mantel()` - 计算所有环境变量和所有距离矩阵。

相关热图可以通过以下两个步骤来完成：
- `cal_cor()` - 计算分类群和环境因素之间的相关性。
- `plot_corr()` - 绘图。

## trans_null
- `trans_nullmodel` - 提供了一个封装，包括系统发育信号的计算，包括计算phylogenetic signal, beta mean pairwise phylogenetic distance (betaMPD), beta mean nearest taxon distance (betaMNTD), beta nearest taxon index (betaNTI), beta net relatedness index (betaNRI) 和 Bray-Curtis-based Raup-Crick (RCbray).。
- `cal_process()` - RCbray 和 betaNTI（或 betaNRI）之间的组合可用于推断在特定假设下主导群落组装的每个生态过程的强度。可用该函数来解析每个推断进程的百分比来实现。

## trans_network
微生物生态学中的共现模式分析是一个热门话题，通常通过应用网络分析来进行分析。
- `trans_network` - 提供了三种网络构建方法：相关网络、SPIIEC-EASI网络和基于概率图模型 (PGM)网络。

## trans_func
功能分析是微生物群落数据分析中一个有吸引力的部分，主要来自其对生物学问题的可解释性。
- `trans_func` - 创建一个trans_func项目
- `cal_spe_func()` - 可以将OTU的分类分配与原核生物的FAPROTAX数据库的分类信息或与原核生物的分类信息进行匹配。FUNGuild 数据库用于真菌。
- `cal_spe_func_perc()` - 计算与群落或网络模块中的每个性状相关的分类群（丰度未加权）或个体（丰度加权）的百分比。

![trans_func工作流程](/R-packages/生物信息学/plots/trans_func工作流程.jpeg)

---
# 参考资料
Liu C, Cui Y, Li X, et al. [microeco: an R package for data mining in microbial community ecology](https://doi.org/10.1093/femsec/fiaa255)[J]. FEMS microbiology ecology, 2021, 97(2): fiaa255.

[Tutorial for R microeco package](https://chiliubio.github.io/microeco_tutorial/)
