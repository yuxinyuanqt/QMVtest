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
                    MGC_Cutoff=20,
                    method='location')

##return the first six columns of the result
# $Tstat
# SNP     QXcat     QZmax
# rs001_1 14.858562 2.7868118
# rs001_2 18.329556 3.3692829
# rs001_3  4.144115 0.5351290
# rs001_4  3.897154 0.2378622
# rs001_5  9.035181 1.4266495
# rs001_6 56.639895 6.4432892
#
# $pvalue
# SNP        QXcat        QZmax
# rs001_1 1.000748e-02 6.426004e-03
# rs001_2 2.127724e-03 9.364946e-04
# rs001_3 7.737065e-01 6.394159e-01
# rs001_4 8.402068e-01 8.644602e-01
# rs001_5 1.204522e-01 1.731820e-01
# rs001_6 2.944397e-11 1.666498e-10


####the scale test (i.e., wM3VNA3.3)
scale_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      method='scale')

##return the first six columns of the result
# $Tstat
# SNP     wM3VNA3.3
# rs001_1 2.74847284
# rs001_2 1.63325427
# rs001_3 0.01996321
# rs001_4 0.21477120
# rs001_5 3.55160246
# rs001_6 7.52380846
#
# $pvalue
# SNP     wM3VNA3.3
# rs001_1 4.134161e-02
# rs001_2 1.794681e-01
# rs001_3 9.961706e-01
# rs001_4 8.862103e-01
# rs001_5 1.382821e-02
# rs001_6 5.106703e-05

####the joint tests (i.e., QMVXcat and QMVZmax)
joint_res <- QMV_test(Genotype,Phedata$Y,
                      Phedata$Sex,
                      Covariate=Phedata[,c(-1,-2)],
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      method='joint')

##return the first six columns of the result
# $Tstat
# SNP     QMVXcat    QMVZmax
# rs001_1 15.5806163 16.466576
# rs001_2 15.7409207 17.382249
# rs001_3  0.5207989  0.902074
# rs001_4  0.5898163  0.532902
# rs001_5 12.7950931 12.068913
# rs001_6 68.2618071 64.794996
#
# $pvalue
# SNP      QMVXcat      QMVZmax
# rs001_1 3.636774e-03 2.452928e-03
# rs001_2 3.387261e-03 1.628796e-03
# rs001_3 9.714449e-01 9.242631e-01
# rs001_4 9.641886e-01 9.702199e-01
# rs001_5 1.232164e-02 1.684610e-02
# rs001_6 5.282339e-14 2.842231e-13

##all of the above association tests
all_res <- QMV_test(Genotype,Phedata$Y,
                    Phedata$Sex,
                    Covariate=Phedata[,c(-1,-2)],
                    missing_cutoff=0.15,
                    MAF_Cutoff=NULL,
                    MGC_Cutoff=20,
                    method='all')

##return the first six columns of the result
# $Tstat
# SNP     QXcat     QZmax     wM3VNA3.3  QMVXcat    QMVZmax
# rs001_1 14.858562 2.7868118 2.74847284 15.5806163 16.466576
# rs001_2 18.329556 3.3692829 1.63325427 15.7409207 17.382249
# rs001_3  4.144115 0.5351290 0.01996321  0.5207989  0.902074
# rs001_4  3.897154 0.2378622 0.21477120  0.5898163  0.532902
# rs001_5  9.035181 1.4266495 3.55160246 12.7950931 12.068913
# rs001_6 56.639895 6.4432892 7.52380846 68.2618071 64.794996
#
# $pvalue
# SNP        QXcat        QZmax    wM3VNA3.3      QMVXcat      QMVZmax
# rs001_1 1.000748e-02 6.426004e-03 4.134161e-02 3.636774e-03 2.452928e-03
# rs001_2 2.127724e-03 9.364946e-04 1.794681e-01 3.387261e-03 1.628796e-03
# rs001_3 7.737065e-01 6.394159e-01 9.961706e-01 9.714449e-01 9.242631e-01
# rs001_4 8.402068e-01 8.644602e-01 8.862103e-01 9.641886e-01 9.702199e-01
# rs001_5 1.204522e-01 1.731820e-01 1.382821e-02 1.232164e-02 1.684610e-02
# rs001_6 2.944397e-11 1.666498e-10 5.106703e-05 5.282339e-14 2.842231e-13

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
                        MGC_Cutoff=20,
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
                        MGC_Cutoff=20,
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
