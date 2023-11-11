# ---- 批量下载GWAS catalog数据 ----
setwd("~") # 设置下载文件位置
options(timeout=300)  # 设置超时为300秒
files <- paste0("33462482-", "GCST", 90011301:90011730, "-EFO_0007874",".h.tsv.gz") # 设置下载文件名
urls <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", # 自行设置文件下载地址
               "GCST90011001-GCST90012000/",
               "GCST",
               90011301:90011730,
               "/",
               "harmonised/",
               files)
purrr::map(1:430, ~download.file(destfile = files[.x], url = urls[.x], mode="wb"))


# ---- 加载包 ----
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments")
# devtools::install_github("explodecomputer/plinkbinr")
# devtools::install_github("rondolab/MR-PRESSO")
# devtools::install_github("WSpiller/RadialMR")
library(devtools)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(MendelianRandomization)
library(gridExtra)
library(grid)
library(lattice)
library(LDlinkR)
library(ggpubr)
library(simex)
library(MRPRESSO)
library(TwoSampleMR)
library(data.table)
library(R.utils)
library(ieugwasr)
library(plinkbinr)
library(dplyr)
library(data.table)
library(doParallel)
library(RadialMR)


# ---- 读取本地clump所需文件路径 ----
get_plink_exe() # 获取plink.exe本地位置，并填入下行。
plink_path = "C:/Users/TCMzh/AppData/Local/R/win-library/4.3/plinkbinr/bin/plink_Windows.exe"
bfile_path = "E:/OneDrive/mendelR/Dependent_files/EUR" # 欧洲人群基因组本地数据
maffile_path = "E:/OneDrive/mendelR/Dependent_files/" # fileFrequency.frq文件路径
setwd("E:/OneDrive/mendelR/") # 设置本地路径
getwd()


# ---- 定义get_eaf_from_1000G函数功能：提取SNP的EAF值 ----
## 来自github("HaobinZhou/Get_MR")（该包已停止维护）的get_eaf_from_1000G函数，用于从1000G的MAF文件中提取EAF并将其与输入数据匹配。
get_eaf_from_1000G<-function(dat,path,type="exposure"){
  corrected_eaf_expo<-function(data_MAF){
    effect=data_MAF$effect_allele.exposure
    other=data_MAF$other_allele.exposure
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    cor1<-which(data_MAF$effect_allele.exposure !=data_MAF$A1)
    data_MAF$eaf.exposure=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.exposure[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.exposure ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.exposure[error]=NA
    data_MAF$type[error]="error"
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    return(data_MAF)
  }
  corrected_eaf_out<-function(data_MAF){
    effect=data_MAF$effect_allele.outcome
    other=data_MAF$other_allele.outcome
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    cor1<-which(data_MAF$effect_allele.outcome !=data_MAF$A1)
    data_MAF$eaf.outcome=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.outcome[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.outcome ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.outcome[error]=NA
    data_MAF$type[error]="error"
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    return(data_MAF)
  }
  if(type=="exposure" && (is.na(dat$eaf.exposure[1])==T || is.null(dat$eaf.exposure)==T)){
    r<-nrow(dat)
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    dat<-corrected_eaf_expo(dat)
    cor1<-dat$cor1
    harna<-dat$harna
    error<-dat$error
    dat<-dat$data_MAF
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    print(paste0("一共有",length(harna),"个SNP在1000G中找不到，占比",length(harna)/r*100,"%"))
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    return(dat)
  }
  if(type=="outcome" && (is.na(dat$eaf.outcome[1])==T || is.null(dat$eaf.outcome)==T)){
    r<-nrow(dat)
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    dat<-corrected_eaf_out(dat)
    cor1<-dat$cor1
    harna<-dat$harna
    error<-dat$error
    dat<-dat$data_MAF
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    print(paste0("一共有",length(harna),"个SNP在1000G找不到，占比",length(harna)/r*100,"%"))
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    return(dat)
  }
  else{return(dat)}
}


# ---- 定义get_f函数：估计工具变量的强度----
## 来自github("HaobinZhou/Get_MR")（该包已停止维护）的get_f函数，用于估计工具变量的强度，获取R2（这里的R2与Clumping过程中的r2不同）和F统计量。
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("数据不包含eaf，无法计算F统计量")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)}
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}


# ---- 通过样本量计算功率 ----
n <- (16587+155860)  # 表示结局（疾病）16587 cases, 155860 controls
ratio <- 16587/155860

sig <- 0.05

Betas <- seq(from=0, to=0.5, by=0.0005)
Powers <- as.data.frame(Betas)
Powers$ORs <- exp(Betas)

Powers$Var0.5 <- (pnorm(sqrt(n*0.005*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var1 <- (pnorm(sqrt(n*0.01*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var2.5 <- (pnorm(sqrt(n*0.025*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var5 <- (pnorm(sqrt(n*0.05*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var10 <- (pnorm(sqrt(n*0.1*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100

PowerPlot <- ggplot(Powers, aes(ORs)) +
  geom_line(aes(y = Var0.5, colour = "9")) +  
  geom_line(aes(y = Var1, colour = "7")) +
  geom_line(aes(y = Var2.5, colour = "5")) + 
  geom_line(aes(y = Var5, colour = "3")) +
  geom_line(aes(y = Var10, colour = "1")) +
  xlab("Odds ratio per unit increase in risk factor") +
  scale_y_continuous("Power (%)", limits=c(0, 100), breaks=c(20,40,60,80,100)) +
  theme(axis.title.x = element_text(face="bold", size=20), axis.text.x  = element_text(vjust=0.5, size=16)) +
  theme(axis.title.y = element_text(face="bold", size=20), axis.text.y  = element_text(vjust=0.5, size=16)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_colour_discrete("% Variance", labels= c("10.0", "5.0", "2.5", "1.0", "0.5")) +
  geom_hline(yintercept=70) +
  ggtitle("Power for analysing all GAME-ON cases & controls") +
  theme(plot.title=element_text(lineheight=5, size= rel(1.2), face="bold")) +
  theme_classic()

dest.plot <- "bacterial_pneumonia_power_2SMR.png"

png(dest.plot, width = 10*500, height = 5*500, res=500)
PowerPlot
dev.off()


# ---- 读取本地暴露数据并进行本地clump ----
exp_MiBioGen <- fread("Taxa_abundance/MiBioGen.allHits.p1e4.txt",header = T) # 读取本地数据

View(exp_MiBioGen[1:10,])
sum(is.na(exp_MiBioGen$rsID))

exp_MiBioGen_p <- subset(exp_MiBioGen, P.weightedSumZ<1e-5) # p值过滤

exp_MiBioGen_format <- format_data( # 本地数据修改为TwoSampleMR格式
  exp_MiBioGen_p,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "bac", 
  snp_col = "rsID",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "eff.allele",
  other_allele_col = "ref.allele",
  pval_col = "P.weightedSumZ",
  samplesize_col = "N",
  min_pval = 1e-200,
  z_col = "Z.weightedSumZ",
  chr_col = "chr",
  pos_col = "bp",
  log_pval = FALSE
)
exp_MiBioGen_eaf <- get_eaf_from_1000G(exp_MiBioGen_format, maffile_path, type = "exposure")  # 补充eaf值

exp_MiBioGen_clump <- ld_clump( # 本地clump
  dplyr::tibble(
    rsid=exp_MiBioGen_eaf$SNP, 
    pval=exp_MiBioGen_eaf$pval.exposure,
    id=exp_MiBioGen_eaf$id.exposure
  ),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1,
  pop = "EUR",
  plink_bin = plink_path,
  bfile = bfile_path
)

ivs_MiBioGen <- subset(exp_MiBioGen_eaf,SNP %in% exp_MiBioGen_clump$rsid) # clump后工具变量个snps数据表


# ---- 计算F统计量 ----
ivs_MiBioGen_f <- get_f(ivs_MiBioGen, F_value = 10) # F>10表示工具变量足够强
summary(ivs_MiBioGen_f$F)
summary(ivs_MiBioGen_f$R2)
write.csv(ivs_MiBioGen_f, "ivs_MiBioGen.csv") # 剔除弱工具变量后的工具变量表


# ---- 读取结局并与暴露进行harmonization ----
out_bacterial_pneumonia_orig <- fread("diseases/f~~~.gz",  # 读取结局，以finngen数据库数据为例
                                      header = T)
View(out_bacterial_pneumonia_orig[1:10,])
out_bacterial_pneumonia_format <- format_data( # 本地数据修改为TwoSampleMR格式
  out_bacterial_pneumonia_orig,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "rsids", 
  beta_col = "beta", 
  se_col = "sebeta", 
  eaf_col = "af_alt", 
  effect_allele_col = "alt", 
  other_allele_col = "ref", 
  pval_col = "pval",
  ncase_col = "af_alt_cases", 
  ncontrol_col = "af_alt_controls", 
  gene_col = "nearest_genes", 
  min_pval = 1e-200,
  chr_col = "#chrom", 
  pos_col = "pos", 
  log_pval = FALSE
)
out_bacterial_pneumonia_format$outcome <- "~" # 定义outcome列为某结局名称

harmonise_data <- TwoSampleMR::harmonise_data(
  exposure_dat = ivs_MiBioGen_f, 
  outcome_dat = out_bacterial_pneumonia_format
)

table(harmonise_data_final$mr_keep) # 查看harmonise后数据mr_keep值的数量，值为ture即为参与mr分析的snp

harmonise_data <- harmonise_data[harmonise_data$mr_keep == TRUE,] # 保留mr_keep为true的SNP

snp_lesthan_3 <- data.frame(table(harmonise_data$id.exposure)) # 去除小于3个SNP的暴露变量
snp_lesthan_3 <- subset(snp_lesthan_3, Freq < 3)
harmonise_data_final <- subset(harmonise_data, 
                               !harmonise_data$id.exposure %in% snp_lesthan_3$Var1)
write.csv(harmonise_data_final, "harmonise_data_final.csv") # 导出
rm(list = ls()) # 删除当前工作环境中所有对象


# ---- MR分析 ----
harmonise_data_final <- fread("harmonise_data_final.csv", header = T)

dat_list <- split(harmonise_data_final,harmonise_data_final$id.exposure) # 转为list,进行mr并行

mr_Speed <- function(dat) {
  res=TwoSampleMR::mr(
    dat,
    method_list =c("mr_egger_regression",
                   "mr_weighted_median",
                   "mr_ivw",
                   "mr_simple_mode",
                   "mr_weighted_mode"
                   # "mr_ivw_mre",
                   # "mr_ivw_fe"
                   )) 
  return(res)
}

cl <- makeCluster(4)  # 设置CPU数、线程数(多线程、多指标，测试：2475个snp经4核心跑了14分钟) 
registerDoParallel(cl) # 登记

clusterEvalQ(cl,library(TwoSampleMR)) # 添加并行计算中用到的包

clusterExport(cl = cl ,varlist = c("dat_list","mr_Speed"))  # 要到的数据与function填入
start1=Sys.time()
res=parLapply(cl = cl,
              X = dat_list, 
              fun = mr_Speed)
end1=Sys.time()
print(paste0("mr()耗时",end1-start1))

res <- do.call(rbind,res)  # 合并为data.frame
write.csv(res, "res.csv")  # 导出mr()结果
# # 主要结果，进行FDR矫正
# res$p_FDR <- p.adjust(res$pval,method = "bonferroni")
# res <- dplyr::arrange(res, p_FDR)

# mr_report(harmonise_data, output_path = "mr",author="Tong Zhou", study = paste("Gut microbiota", "Bacterial pneumoniae",sep=""))
or_results <- generate_odds_ratios(res)
results <- cbind.data.frame(
  outcome = or_results$outcome,
  exposure = or_results$exposure,
  nsnp = or_results$nsnp,
  method = or_results$method,
  b = or_results$b,
  se = or_results$se,
  pval = or_results$pval,
  or = or_results$or,
  or_lci95 = or_results$or_lci95,
  or_uci95 = or_results$or_uci95
)
write.csv(results, "or_res.csv") # 导出带有odds_ratios的mr结果
rm(list = ls()) # 删除当前工作环境中所有对象


# ---- 可视化准备：MR结果读取 ----
harmonise_data <- fread("harmonise_data_final.csv", header = TRUE)
res_data <- fread("res.csv", header = TRUE)

res_single <- mr_singlesnp(harmonise_data)
res_ivw_p5e2 <- res_data %>% # 过滤出IVW_p5E-2的显著IVs结果
  filter(method == 'Inverse variance weighted' & pval < 0.05) %>%
  select(id.exposure, exposure, id.outcome, outcome) %>%
  as.data.frame()

write.csv(res_single, "res_single_snp.csv")
write.csv(res_ivw_p5e2, "res_ivw_p5e2.csv")


# ---- 每个IVs_p5e2的因果效应散点图 ----
mr_scatter <- mr_scatter_plot(res_data, harmonise_data)
# 设置文件夹路径
scatter_filepath <- "res_scatter"
if (!dir.exists(scatter_filepath)) {
  dir.create(scatter_filepath)
}
for(i in 1:nrow(res_ivw_p5e2)) {
  tryCatch({
    # 提取特定的id.exposure和id.outcome组合
    id_exp <- res_ivw_p5e2$id.exposure[i]
    id_out <- res_ivw_p5e2$id.outcome[i]
    plot_key <- paste0(id_exp, ".", id_out)
    # 检查是否存在对应的图表
    if (!is.null(mr_scatter[[plot_key]])) {
      # 文件名使用exposure和outcome的值
      file_name <- paste0(scatter_filepath, "/", res_ivw_p5e2$exposure[i], "_", res_ivw_p5e2$outcome[i], "_scatter.pdf")
      # 保存PDF
      pdf(file_name, width = 15, height = 20)
      # 使用双括号动态引用列表中的图表
      print(ggarrange(mr_scatter[[plot_key]], ncol = 2, nrow = 2, widths = c(2, 2), heights = c(1, 1)))
      dev.off()
    } else {
      message("No plot for ", plot_key)
    }
  }, error = function(e) {
    message("An error occurred: ", e$message)
    dev.off() # make sure to close the device in case of error
  })
}


# ---- 每个IVs_p5e2的SNP效应forest pot ----
mr_forest <- mr_forest_plot(res_single)
mr_forest
forest_filepath <- "res_forest"
if (!dir.exists(forest_filepath)) {
  dir.create(forest_filepath)
}
for(i in 1:nrow(res_ivw_p5e2)) {
  tryCatch({
    id_exp <- res_ivw_p5e2$id.exposure[i]
    id_out <- res_ivw_p5e2$id.outcome[i]
    plot_key <- paste0(id_exp, ".", id_out)
    if (!is.null(mr_forest[[plot_key]])) {
      file_name <- paste0(forest_filepath, "/", res_ivw_p5e2$exposure[i], "_", res_ivw_p5e2$outcome[i], "_scatter.pdf")
      pdf(file_name, width = 15, height = 20)
      print(ggarrange(mr_forest[[plot_key]], ncol = 2, nrow = 2, widths = c(2, 2), heights = c(1, 1)))
      dev.off()
    } else {
      message("No plot for ", plot_key)
    }
  }, error = function(e) {
    message("An error occurred: ", e$message)
    dev.off() 
  })
}


# ---- 每个IVs_p5e2的SNP效应funnel pot ----
mr_funnel <- mr_funnel_plot(res_single)
funnel_filepath <- "res_funnel"
if (!dir.exists(funnel_filepath)) {
  dir.create(funnel_filepath)
}
for(i in 1:nrow(res_ivw_p5e2)) {
  tryCatch({
    id_exp <- res_ivw_p5e2$id.exposure[i]
    id_out <- res_ivw_p5e2$id.outcome[i]
    plot_key <- paste0(id_exp, ".", id_out)
    if (!is.null(mr_funnel[[plot_key]])) {
      file_name <- paste0(funnel_filepath, "/", res_ivw_p5e2$exposure[i], "_", res_ivw_p5e2$outcome[i], "_scatter.pdf")
      pdf(file_name, width = 15, height = 20)
      print(ggarrange(mr_funnel[[plot_key]], ncol = 2, nrow = 2, widths = c(2, 2), heights = c(1, 1)))
      dev.off()
    } else {
      message("No plot for ", plot_key)
    }
  }, error = function(e) {
    message("An error occurred: ", e$message)
    dev.off() 
  })
}


# ---- 每个IVs_p5e2的SNP留一法forest pot ----
res_loo <- mr_leaveoneout(harmonise_data)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

loo_filepath <- "res_loo"
if (!dir.exists(loo_filepath)) {
  dir.create(loo_filepath)
}
for(i in 1:nrow(res_ivw_p5e2)) {
  tryCatch({
    id_exp <- res_ivw_p5e2$id.exposure[i]
    id_out <- res_ivw_p5e2$id.outcome[i]
    plot_key <- paste0(id_exp, ".", id_out)
    if (!is.null(mr_loo[[plot_key]])) {
      file_name <- paste0(loo_filepath, "/", res_ivw_p5e2$exposure[i], "_", res_ivw_p5e2$outcome[i], "_scatter.pdf")
      pdf(file_name, width = 15, height = 20)
      print(ggarrange(mr_loo[[plot_key]], ncol = 2, nrow = 2, widths = c(2, 2), heights = c(1, 1)))
      dev.off()
    } else {
      message("No plot for ", plot_key)
    }
  }, error = function(e) {
    message("An error occurred: ", e$message)
    dev.off() 
  })
}


# ---- 异质性分析和多效性分析 ----
Het<-mr_heterogeneity(harmonise_data)
Plt<-mr_pleiotropy_test(harmonise_data)

write.csv(Het, "res_heterogeneity.csv") # 异质性：Cochrance's Q值
write.csv(Plt, "res_pleiotropy.csv") # 多效性：Egger回归的截距

res_single2 <- res_single[grep("^rs", res_single$SNP), ] # 以下为求出每一个IVs的I2循环

res_isq <- data.frame()

unique_combinations <- unique(res_single2[c("id.exposure", "id.outcome", "exposure", "outcome")])

for (i in 1:nrow(unique_combinations)) {

  current_combination <- unique_combinations[i, ]

  subset_res_single2 <- subset(res_single2, 
                               id.exposure == current_combination$id.exposure & 
                                 id.outcome == current_combination$id.outcome & 
                                 exposure == current_combination$exposure & 
                                 outcome == current_combination$outcome)
  
  subset_harmonise_data <- subset(harmonise_data, 
                                  id.exposure == current_combination$id.exposure & 
                                    id.outcome == current_combination$id.outcome & 
                                    exposure == current_combination$exposure & 
                                    outcome == current_combination$outcome)

  res_meta <- metafor::rma(yi = subset_res_single2$b, 
                           sei = subset_res_single2$se,
                           weights = 1/subset_harmonise_data$se.outcome^2,
                           data = subset_res_single2,
                           method = "FE")

  I2_value <- as.numeric(res_meta$I2)
  H2_value <- as.numeric(res_meta$H2)
  QE_value <- as.numeric(res_meta$QE)
  QEp_value <- as.numeric(res_meta$QEp)

  temp_df <- data.frame(
    id.exposure = current_combination$id.exposure,
    id.outcome = current_combination$id.outcome,
    exposure = current_combination$exposure,
    outcome = current_combination$outcome,
    I2 = I2_value,
    H2 = H2_value,
    QE = QE_value,
    QEp = QEp_value
  )

  res_isq <- rbind(res_isq, temp_df)
}

write.csv(res_isq, "res_isq.csv") # 异质性：I2、H2、IVW法的Cochrance's Q值、Q_pval值


# ---- 径向MR分析检测异常值 ----
harmonise_single <- harmonise_data[harmonise_data$SNP%in%res_single$SNP,]
raddat <- format_radial(harmonise_single$beta.exposure, 
                        harmonise_single$beta.outcome, 
                        harmonise_single$se.exposure, 
                        harmonise_single$se.outcome, 
                        harmonise_single$SNP)
# ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
# dim(ivwrad$outliers)[1] 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]  # 查看ivw异常值数量
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] # 查看egger异常值数量
#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE) # 绘制径向图（snp数量过多无法出图，可跳过）
ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "res_outliers.csv", row.names=F, quote=F) # 导出异常值表格


# ---- 去除异常值后再次MR分析 ----
harmonise_filer <- harmonise_data[!harmonise_data$SNP %in% ivwrad$outliers$SNP,] # 清理异常值snp

dat_list <- split(harmonise_filer,harmonise_filer$id.exposure) # mr分析
# install.packages("doParallel")
library(doParallel)
mr_Speed <- function(dat) {
  res=TwoSampleMR::mr(
    dat,
    method_list =c("mr_egger_regression",
                   "mr_weighted_median",
                   "mr_ivw",
                   "mr_simple_mode",
                   "mr_weighted_mode"
                   # "mr_ivw_mre",
                   # "mr_ivw_fe"
    )) 
  return(res)
}
cl <- makeCluster(4)  # 设置CPU数\线程数 
registerDoParallel(cl) # 登记
clusterEvalQ(cl,library(TwoSampleMR)) # 添加并行计算中用到的包
clusterExport(cl = cl ,varlist = c("dat_list","mr_Speed")) # 要到的数据与function填入
start1=Sys.time()
res=parLapply(cl = cl,
              X = dat_list,
              fun = mr_Speed)
end1=Sys.time()
print(paste0("mr()耗时",end1-start1))
res <- do.call(rbind,res) # 合并为data.frame
write.csv(res, "filter_res.csv")

or_results <- generate_odds_ratios(res) # 提取od值
results <- cbind.data.frame(
  outcome = or_results$outcome,
  exposure = or_results$exposure,
  nsnp = or_results$nsnp,
  method = or_results$method,
  b = or_results$b,
  se = or_results$se,
  pval = or_results$pval,
  or = or_results$or,
  or_lci95 = or_results$or_lci95,
  or_uci95 = or_results$or_uci95
)
write.csv(results, "filter_or_res.csv")


# ---- MR-PRESSO ----
harmonise_single = as.data.frame(harmonise_single)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", 
                              BetaExposure = "beta.exposure", 
                              SdOutcome = "se.outcome", 
                              SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE, 
                              DISTORTIONtest = TRUE, 
                              data = harmonise_single, 
                              NbDistribution = 3000,   # ND值1000到10000，根据SNP数量和电脑配置自行调节
                              SignifThreshold = 0.05)
write.table(mr_presso, "bmi_mr_presso.txt")

