#!/usr/bin env python3

'''
Timothy Isonio
July/August 2016
'''

import sys, math
import numpy as np
import matplotlib.pyplot as plt

def y(maf, h, n):
  y_hat = np.random.normal(0, h**2/n) * (np.random.binomial(2, maf) - 2*maf) / np.sqrt(2*maf*(1-maf))
  e = np.random.normal(0, 1-h**2)

  return(y_hat, y_hat + e)

# number of individuals
m = 1000

# trait heritability
h = .4

# number of SNPs
n = 1000

# each SNP has allele frequency~uniform(0.05, 0.95)
a = 0.05

# read in from a file
# with open(sys.argv[1], 'r') as in_tsv:
#   mafs = [float(line.strip().split('\t')[-1]) for line in in_tsv.readlines()]  # MAF is always last column
#   mafs = [maf for maf in mafs if maf < 1]  # ignore rows where the last column > 1



# mafs = [np.random.uniform(a, 1-a) for i in range(n)]
# beta_stddev = math.sqrt(h**2/n)
# betas = [np.random.normal(0, beta_stddev) for i in range(n)]

mafs = np.random.uniform(a, 1-a, size=n)
beta_stddev = math.sqrt(h**2/n)
betas = np.random.normal(0, beta_stddev, size=n)

# mafs = np.fromiter( (np.random.uniform(a, 1-a) for i in range(n)), dtype=float)
# betas = np.fromiter( (np.random.normal(0, beta_stddev) for i in range(n)), dtype=float)

y_hats = []
y_trues = []
#betas = []
gs = []

for ind_index in range(m):
  # y_vectorize = np.vectorize(y)
  # y_hats, y_trues = y_vectorize(mafs, h, n)

  # y_hat = np.sum(y_hats)
  # y_true = np.sum(y_trues)
  
  # print(y_hat, y_true)
  # sys.exit()

  y_hat = 0
  y_true = 0
  for snp_index in range(n):

    # generate betas for each person
    #beta = np.random.normal(0, math.sqrt(h**2/n))
    # betas.append(beta)
    beta = betas[snp_index] # grab this SNP's beta
    p = mafs[snp_index]  # grab this SNP's allele frequency
    
    # generate genotype for each person at each SNP
    g = np.random.binomial(2, p)
    gs.append(g)
    
    # standarize genotype
    x = (g - 2*p) / math.sqrt(2*p*(1-p)) 
        
    y_hat += beta*x
    y_true += beta*x
    
  e = np.random.normal(0, math.sqrt(1-h**2))
  y_true += e
  y_hats.append(y_hat)
  y_trues.append(y_true)

y_trues_minus_y_hats = [a - b for a, b in zip(y_trues, y_hats)]


print("mean, var y_hats: ", np.mean(y_hats), np.var(y_hats))
print("mean, var y_trues: ", np.mean(y_trues), np.var(y_trues))
print("mean, var y_trues_minus_y_hats: ", np.mean(y_trues_minus_y_hats), np.var(y_trues_minus_y_hats))

# plt.hist(y_hats)
# plt.show()

# plt.hist(y_trues)
# plt.show()

# plt.hist(y_trues_minus_y_hats)
# plt.show()



# for entry in y_results:
#   print(entry)
