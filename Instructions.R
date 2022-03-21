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

#Phedata: phenotype (Y) and covariates (Sex, age and BMI) data for 4000 unrelated individuals
data(Phedata)
#Genotype: a data for 4000 unrelated individuals and 31 SNPs
data(Genotype)

#the location tests (i.e., QXcat and QZmax)
loc_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=30,
                    method='location')
#result
head(loc_res)
# SNP     QXcat        QZmax
# rs001_1 1.000000e+00 8.588311e-01
# rs001_3 7.150318e-07 2.428027e-05
# rs001_4 2.884434e-01 3.344079e-01
# rs001_5 8.681130e-09 1.949980e-06
# rs001_6 2.476245e-06 7.409794e-07
# rs001_7 5.027080e-01 5.699959e-01

#the scale test (i.e., wM3VNA3.3)
scale_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=30,
                      method='scale')

#result
head(scale_res)
# SNP     wM3VNA3.3
# rs001_1 5.240164e-02
# rs001_3 8.769083e-06
# rs001_4 5.063405e-01
# rs001_5 1.978776e-02
# rs001_6 6.181299e-03
# rs001_7 6.849780e-01

#the joint tests (i.e., QMVXcat and QMVZmax)
joint_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=30,
                      method='joint')

#result
head(joint_res)
# SNP     QMVXcat      QMVZmax
# rs001_1 2.069245e-01 1.845621e-01
# rs001_3 1.680106e-10 4.954574e-09
# rs001_4 4.270229e-01 4.700339e-01
# rs001_5 4.034222e-09 6.972588e-07
# rs001_6 2.907451e-07 9.252731e-08
# rs001_7 7.114539e-01 7.576364e-01

#All of the above association tests
all_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=30,
                    method='all')

#result
head(all_res)
# SNP      QXcat        QZmax        wM3VNA3.3   QMVXcat      QMVZmax
# rs001_1 1.000000e+00 8.588311e-01 5.240164e-02 2.069245e-01 1.845621e-01
# rs001_3 7.150318e-07 2.428027e-05 8.769083e-06 1.680106e-10 4.954574e-09
# rs001_4 2.884434e-01 3.344079e-01 5.063405e-01 4.270229e-01 4.700339e-01
# rs001_5 8.681130e-09 1.949980e-06 1.978776e-02 4.034222e-09 6.972588e-07
# rs001_6 2.476245e-06 7.409794e-07 6.181299e-03 2.907451e-07 9.252731e-08
# rs001_7 5.027080e-01 5.699959e-01 6.849780e-01 7.114539e-01 7.576364e-01

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
# marker.ID  chromosome  genetic.dist  physical.pos    QXcat        QZmax
# rs001_1    X           0             1               9.535619e-01 8.193516e-01
# rs001_2    X           0             2               5.076249e-02 1.538114e-01
# rs001_3    X           0             3               4.802256e-01 3.652537e-01
# rs001_4    X           0             4               1.000000e+00 9.993240e-01
# rs001_5    X           0             5               1.582337e-03 6.297673e-04
# rs001_6    X           0             6               3.015064e-07 2.761678e-06

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
# marker.ID   chromosome  genetic.dist  physical.pos  allele1 allele2    QXcat         QZmax
# rs001_1     23          0             1             2       1          9.923633e-01  8.408148e-01
# rs001_2     23          0             2             2       1          4.883641e-02  1.496640e-01
# rs001_3     23          0             3             2       1          5.365879e-01  4.016883e-01
# rs001_4     23          0             4             2       1          1.000000e+00  9.850508e-01
# rs001_5     23          0             5             2       1          1.827692e-03  7.086208e-04
# rs001_6     23          0             6             2       1          3.439666e-07  3.310289e-06
