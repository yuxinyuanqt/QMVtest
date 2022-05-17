####################################################################

#install package dependencies
install.packages(c('quantreg', 'mvtnorm', 'expm', 'bigsnpr', 'xgboost'))

#install the package
#setwd to the root of the folder path of "QMVtest_1.0.zip" and run the following command
install.packages("QMVtest_1.0.zip",repos = NULL)

#or setwd to the root of the folder path of "QMVtest_1.0.tar.gz" and run the following command
install.packages("QMVtest_1.0.tar.gz",repos = NULL,
                 type="source")

#load the package
library(QMVtest)

#check the help files
help(package='QMVtest')

################examples for QMV_test()#############

##A function to obtain the p-values and the test statistics of the location tests, the scale test, the location-scale tests or all

#Phedata: a dataset for 4000 unrelated individuals and four variables, including one phenotype and three covariates (i.e., Sex,age and BMI).
data(Phedata)
#Genotype: a dataset for 4000 unrelated individuals and 31 SNPs.
data(Genotype)

####the location tests (i.e., QXcat and QZmax)
loc_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=30,
                    method='location')

##return the first six columns of the result
# $Tstat
# SNP     QXcat     QZmax
# rs001_1  2.830649 0.2450410
# rs001_3 35.553406 4.2797979
# rs001_4  6.846574 1.0307625
# rs001_5 44.816709 4.8158328
# rs001_6 32.923876 5.0072004
# rs001_7  5.370449 0.6344447
# $pvalue
# SNP        QXcat        QZmax
# rs001_1 1.000000e+00 8.588311e-01
# rs001_3 7.150318e-07 2.428027e-05
# rs001_4 2.884434e-01 3.344079e-01
# rs001_5 8.681130e-09 1.949980e-06
# rs001_6 2.476245e-06 7.409794e-07
# rs001_7 5.027080e-01 5.699959e-01


####the scale test (i.e., wM3VNA3.3)
scale_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=30,
                      method='scale')

##return the first six columns of the result
# $Tstat
# SNP     wM3VNA3.3
# rs001_1 2.57215320
# rs001_3 8.75227708
# rs001_4 0.77759966
# rs001_5 3.29056092
# rs001_6 4.13339700
# rs001_7 0.49612353
#
# $pvalue
# SNP     wM3VNA3.3
# rs001_1 5.240164e-02
# rs001_3 8.769083e-06
# rs001_4 5.063405e-01
# rs001_5 1.978776e-02
# rs001_6 6.181299e-03
# rs001_7 6.849780e-01

####the joint tests (i.e., QMVXcat and QMVZmax)
joint_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=30,
                      method='joint')

##return the first six columns of the result
# $Tstat
# SNP     QMVXcat    QMVZmax
# rs001_1  5.8976347  6.2020007
# rs001_3 51.5904344 44.5402499
# rs001_4  3.8476051  3.5518796
# rs001_5 44.9696121 34.1407671
# rs001_6 35.9899885 38.4030396
# rs001_7  2.1322285  1.8809893
#
# $pvalue
# SNP      QMVXcat      QMVZmax
# rs001_1 2.069245e-01 1.845621e-01
# rs001_3 1.680106e-10 4.954574e-09
# rs001_4 4.270229e-01 4.700339e-01
# rs001_5 4.034222e-09 6.972588e-07
# rs001_6 2.907451e-07 9.252731e-08
# rs001_7 7.114539e-01 7.576364e-01

##all of the above association tests
all_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=30,
                    method='all')

##return the first six columns of the result
# $Tstat
# SNP     QXcat     QZmax     wM3VNA3.3  QMVXcat    QMVZmax
# rs001_1  2.830649 0.2450410 2.57215320  5.8976347  6.2020007
# rs001_3 35.553406 4.2797979 8.75227708 51.5904344 44.5402499
# rs001_4  6.846574 1.0307625 0.77759966  3.8476051  3.5518796
# rs001_5 44.816709 4.8158328 3.29056092 44.9696121 34.1407671
# rs001_6 32.923876 5.0072004 4.13339700 35.9899885 38.4030396
# rs001_7  5.370449 0.6344447 0.49612353  2.1322285  1.8809893
#
# $pvalue
# SNP        QXcat        QZmax    wM3VNA3.3      QMVXcat      QMVZmax
# rs001_1 1.000000e+00 8.588311e-01 5.240164e-02 2.069245e-01 1.845621e-01
# rs001_3 7.150318e-07 2.428027e-05 8.769083e-06 1.680106e-10 4.954574e-09
# rs001_4 2.884434e-01 3.344079e-01 5.063405e-01 4.270229e-01 4.700339e-01
# rs001_5 8.681130e-09 1.949980e-06 1.978776e-02 4.034222e-09 6.972588e-07
# rs001_6 2.476245e-06 7.409794e-07 6.181299e-03 2.907451e-07 9.252731e-08
# rs001_7 5.027080e-01 5.699959e-01 6.849780e-01 7.114539e-01 7.576364e-01

################examples for QMV_test_ped()#############

##Read PLINK raw files into R, and compute the p-values and the test statistics of the location tests, the scale test, the location-scale tests or all

path <- system.file("extdata", package = "QMVtest")

#Reading the pedfile and storing the data in temporary directory
#the joint tests (i.e., QMVXcat and QMVZmax)
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
                        method='joint')

##return the first six columns of the result
# $Tstat
# marker.ID chromosome genetic.dist physical.pos   QMVXcat    QMVZmax
# rs001_1          X            0            1    1.3989611  1.7023433
# rs001_2          X            0            2    6.1107536  3.8936150
# rs001_3          X            0            3    5.9060861  6.4534139
# rs001_4          X            0            4    0.4464889  0.4478413
# rs001_5          X            0            5   13.0847191 14.9273343
# rs001_6          X            0            6   37.3398049 32.9102002
#
# $pvalue
# marker.ID chromosome genetic.dist physical.pos      QMVXcat      QMVZmax
# rs001_1          X            0            1     8.443756e-01 7.902919e-01
# rs001_2          X            0            2     1.910284e-01 4.205949e-01
# rs001_3          X            0            3     2.062725e-01 1.677495e-01
# rs001_4          X            0            4     9.784968e-01 9.783759e-01
# rs001_5          X            0            5     1.086917e-02 4.854328e-03
# rs001_6          X            0            6     1.533086e-07 1.246129e-06

################examples for QMV_test_bed()#############

##Read PLINK bed files into R, and compute the p-values and the test statistics of the location tests, the scale test, the location-scale tests or all

path <- system.file("extdata", package = "QMVtest")

#Reading the bedfile and storing the data in temporary directory
#the joint tests (i.e., QMVXcat and QMVZmax)
#no imputation for missing genotypes
bed_res <- QMV_test_bed(paste(path,"example-missing.bed",sep = '/'),
                        fastImpute=FALSE,
                        Covariate_path=paste(path,"Covariate.txt",sep = '/'),
                        Covariate_missing = NA,
                        Covariate.header=TRUE,
                        missing_cutoff=0.15,
                        MAF_Cutoff=NULL,
                        MGC_Cutoff=30,
                        method='joint')

##return the first six columns of the result
# $Tstat
# marker.ID chromosome genetic.dist physical.pos allele1 allele2    QMVXcat    QMVZmax
# rs001_1         23            0            1       2       1     1.2383211  1.5697568
# rs001_2         23            0            2       2       1     6.1956130  3.9557804
# rs001_3         23            0            3       2       1     5.2496780  5.8287861
# rs001_4         23            0            4       2       1     0.4105257  0.4406500
# rs001_5         23            0            5       2       1    12.7042006 14.5991775
# rs001_6         23            0            6       2       1    37.2744346 32.7459421
# $pvalue
# marker.ID chromosome genetic.dist physical.pos allele1 allele2      QMVXcat      QMVZmax
# rs001_1         23            0            1       2       1      8.717499e-01 8.142172e-01
# rs001_2         23            0            2       2       1      1.850083e-01 4.120234e-01
# rs001_3         23            0            3       2       1      2.626247e-01 2.123047e-01
# rs001_4         23            0            4       2       1      9.816061e-01 9.790155e-01
# rs001_5         23            0            5       2       1      1.281532e-02 5.609000e-03
# rs001_6         23            0            6       2       1      1.581391e-07 1.346428e-06
