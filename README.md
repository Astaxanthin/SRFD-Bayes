# Tumor fractions deciphered from circulating cell-free DNA methylation for cancer early diagnosis

Official MATLAB code for the paper "Tumor fractions deciphered from circulating cell-free DNA methylation for cancer early diagnosis".

## Prerequisites

MATLAB R2018b

## Experimental settings

The experiments consist of two parts: deconvolution and diagnosis, each is performed on both simulation dataset and real dataset.

In the deconvolution step, the semi-reference-deconvolution (SRFD) is first implemented on cfDNA methylation data to obtain a reference database, which is then utilized to deconvolve the test samples to decipher their fraction vectors. To estimate the tumor fractions of real samples, the reference database learned from simulation dataset is directly utilized for the deconvolution of real cfDNA methylation data from cancer patients.

In the diagnosis step, the diagnostic prior is first obtained from the machine learning based classifiers, and then the conditional probability distribution is computed from the tumor components in the fraction vectors of each test sample. The prior and the conditional probability distribution are combined to make the final Bayesian diagnostic decision.



## Usage

```
  run main.m
```

## File instruction


The 'data' directory contains two type of datasets: simulation dataset and real dataset.

+ *simulation_dataset*

  The simulation data contains four datasets with the CNV event probabilities of 0, 10%, 30% and 50%. Each dataset consists of 2400 training (400 for each category) and 2400 test (400 for each category) samples.

  *train_data* is a $K_s \times N$ matrix. $N=2400$ suggests the number of samples, each with $K_s=350$ dimensional methylation levels, i.e. $\beta$ value. Similarly, *test_data* is also a $K_s \times N$ matrix. *train_theta* and *test_theta* denote the simulated tumor fraction of training and test samples.

+ *real_dataset*
  
  The real data mainly contains three datasets: chip-based methylation data from GSE122126, GSE108462 and GSE129374, Xu et al. data and Chen et al. data. The chip-based data is formulated in file *validation_real_data.mat* .


## References
If you find this work or code useful, please cite this study. The citation will be updated soon. If you have any questions about this code, please contact zhouxiao17@mails.tsinghua.edu.cn

