####################################################################

#install package dependencies
install.packages(c('quantreg', 'mvtnorm', 'expm', 'bigsnpr', 'xgboost'))

# install the package
# setwd to the root of the folder path of "QMVtest_1.0.zip" and run the following command
install.packages("QMVtest_1.0.zip",repos = NULL)

# or setwd to the root of the folder path of "QMVtest_1.0.tar.gz" and run the following command
install.packages("QMVtest_1.0.tar.gz",repos = NULL,
                 type="source")

# load the package
library(QMVtest)

# check the help files
help(package='QMVtest')

################examples for QMV_test()#############

#Phedata: phenotype(Y) and covariates(Sex,age and BMI) data for 4000 individuals
data(Phedata)
#Genotype: a data for 4000 individuals and 31 SNPs
data(Genotype)

#location tests: QXcat and QZmax
loc_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=10,
                    method='location')
#result
head(loc_res)
# SNP        QXcat        QZmax
# rs001_1 1.000000e+00 8.594667e-01
# rs001_2 3.922845e-02 4.584644e-02
# rs001_3 7.148378e-07 2.427947e-05
# rs001_4 2.929903e-01 3.395235e-01
# rs001_5 8.756749e-09 1.970099e-06
# rs001_6 2.841031e-06 8.609484e-07

#scale test: wM3VNA3.3
scale_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=10,
                      method='scale')

#result
head(scale_res)
# SNP    wM3VNA3.3
# rs001_1 5.240164e-02
# rs001_2 1.395192e-02
# rs001_3 8.769083e-06
# rs001_4 5.063405e-01
# rs001_5 1.978776e-02
# rs001_6 6.181299e-03

#joint tests: QMVXcat and QMVZmax
joint_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=10,
                      method='joint')

#result
head(joint_res)
# SNP      QMVXcat      QMVZmax
# rs001_1 2.069245e-01 1.846653e-01
# rs001_2 4.657895e-03 5.343983e-03
# rs001_3 1.679668e-10 4.954418e-09
# rs001_4 4.314340e-01 4.746143e-01
# rs001_5 4.067860e-09 7.040526e-07
# rs001_6 3.311627e-07 1.067094e-07

#All of the above association tests
all_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=10,
                    method='all')

#result
head(all_res)
# SNP        QXcat        QZmax    wM3VNA3.3      QMVXcat      QMVZmax
# rs001_1 1.000000e+00 8.594667e-01 5.240164e-02 2.069245e-01 1.846653e-01
# rs001_2 3.922845e-02 4.584644e-02 1.395192e-02 4.657895e-03 5.343983e-03
# rs001_3 7.148378e-07 2.427947e-05 8.769083e-06 1.679668e-10 4.954418e-09
# rs001_4 2.929903e-01 3.395235e-01 5.063405e-01 4.314340e-01 4.746143e-01
# rs001_5 8.756749e-09 1.970099e-06 1.978776e-02 4.067860e-09 7.040526e-07
# rs001_6 2.841031e-06 8.609484e-07 6.181299e-03 3.311627e-07 1.067094e-07

################examples for QMV_test_ped()#############

path <- system.file("extdata", package = "QMVtest")

# Reading the pedfile and storing the data in temporary directory
ped_res <- QMV_test_ped(paste(path,"example.raw",sep = '/'),trait_missing=NA,
                        Genotype_missing=NA,
                        ped_raw.header=FALSE,
                        Covariate_path=paste(path,"Covariate.txt",sep = '/'),
                        Covariate_missing = NA,
                        Covariate.header=TRUE,
                        map_path=paste(path,"example.map",sep = '/'),
                        map_header=FALSE,
                        missing_cutoff=0.15,
                        MAF_Cutoff=NULL,
                        MGC_Cutoff=30,
                        method='location')

#result
head(ped_res)
# marker.ID chromosome genetic.dist physical.pos        QXcat        QZmax
#   rs001_1          X            0            1 9.901308e-01 8.528473e-01
#   rs001_2          X            0            2 5.241943e-02 1.565153e-01
#   rs001_3          X            0            3 4.810147e-01 3.653642e-01
#   rs001_4          X            0            4 1.000000e+00 9.994096e-01
#   rs001_5          X            0            5 2.348808e-03 9.227721e-04
#   rs001_6          X            0            6 3.358941e-07 2.991962e-06

################examples for QMV_test_bed()#############

path <- system.file("extdata", package = "QMVtest")

# Reading the bedfile and storing the data in temporary directory
bed_res <- QMV_test_bed(paste(path,"example-missing.bed",sep = '/'),
                        fastImpute=FALSE,
                        Covariate_path=paste(path,"Covariate.txt",sep = '/'),
                        Covariate_missing = NA,
                        Covariate.header=TRUE,
                        missing_cutoff=0.15,
                        MAF_Cutoff=NULL,
                        MGC_Cutoff=30,
                        method='location')

#result
head(bed_res)
# marker.ID chromosome genetic.dist physical.pos allele1 allele2        QXcat        QZmax
#   rs001_1         23            0            1       2       1 1.000000e+00 0.8753876284
#   rs001_2         23            0            2       2       1 5.046415e-02 0.1523528667
#   rs001_3         23            0            3       2       1 5.376889e-01 0.4021731516
#   rs001_4         23            0            4       2       1 1.000000e+00 0.9815029317
#   rs001_5         23            0            5       2       1 2.687029e-03 0.0010289593
#   rs001_6         23            0            6       2       1 3.854863e-07 0.0000036112
