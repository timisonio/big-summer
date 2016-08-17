#!/usr/bin env python3

'''
Timothy Isonio
July/August 2016
'''
import sys, math
import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt

# number of individuals
m = 1000

# number of SNPs
n = 1000

# each SNP has allele frequency~uniform(0.05, 0.95)
a = 0.05

# trait heritability
h2 = 0.9

m_sample = 1000 # number of people to sample

def generate_beta_hats(x, y):
  # x is a 2D numpy array (holding genotypes)
  # y is a vector (holding phenotypes)

  x_height, x_width = x.shape

  clf = linear_model.LinearRegression(n_jobs=3) # use multicore
  clf.fit(x, y)
  return clf.coef_

def generate_genotypes(mafs, m):
  # generate genotypes from MAFs using binomial, standardize, for each SNP for each person
  return np.array(
    [
      [
        (np.random.binomial(2, p) - 2*p) / math.sqrt(2*p*(1-p)) for p in mafs
      ] for ind_index in range(m)
    ]
  )

def generate_phenotypes(betas, genotypes, n, m, h2):
  # generate phenotypes using model y=sum(b*x)+e
  return [
    np.sum(
      [
        betas[snp_index]*genotypes[ind_index][snp_index] for snp_index in range(n)
      ]
    ) + np.random.normal(0, math.sqrt(1-h2))
    for ind_index in range(m)
  ]

def generate_yhats(phenotypes, genotypes, n, m, h2):
  beta_hats = generate_beta_hats(genotypes, phenotypes)
  return [
    np.sum(
      [
        beta_hats[snp_index]*genotypes[ind_index][snp_index] for snp_index in range(n)
      ]
    ) + np.random.normal(0, math.sqrt(1-h2))
    for ind_index in range(m)
  ]

def generate_yhats_with_betas(betas, genotypes, n, m, h2):
    return [
    np.sum(
      [
        betas[snp_index]*genotypes[ind_index][snp_index] for snp_index in range(n)
      ]
    ) + np.random.normal(0, math.sqrt(1-h2))
    for ind_index in range(m)
  ]

mafs = np.random.uniform(a, 1-a, size=n)
beta_stddev = math.sqrt(h2/n)
betas = np.random.normal(0, beta_stddev, size=n)

genotypes = generate_genotypes(mafs, m)
phenotypes = generate_phenotypes(betas, genotypes, n, m, h2)


# generate a new pool of individuals not in the larger pool
new_genotypes = generate_genotypes(mafs, m_sample)
new_phenotypes = generate_phenotypes(betas, new_genotypes, n, m_sample, h2)
#new_yhats = generate_yhats(new_phenotypes, new_genotypes, n, m_sample, h2)
new_yhats = generate_yhats_with_betas(betas, new_genotypes, n, m_sample, h2)

new_mse = np.mean([(a - b)**2 for a, b in zip(new_yhats, new_phenotypes)])
new_mpe = 100 / n * np.sum((a - b)/a for a, b in zip(new_phenotypes, new_yhats))

# print("E[y true]:", np.mean(phenotypes))
# print("Var[y true]:", np.var(phenotypes))
# print("New MSE, New MPE: ", new_mse, new_mpe)



print(new_mse, new_mpe, sep=',', end=',')


# pull out some people from the larger pool
inside_genotypes = genotypes[0:m_sample]
inside_phenotypes = phenotypes[0:m_sample]
#inside_yhats = generate_yhats(inside_phenotypes, inside_genotypes, n, m_sample, h2)
inside_yhats = generate_yhats_with_betas(betas, inside_genotypes, n, m_sample, h2)

inside_mse = np.mean([(a - b)**2 for a, b in zip(inside_yhats, inside_phenotypes)])
inside_mpe = 100 / n * np.sum((a - b)/a for a, b in zip(inside_phenotypes, inside_yhats))

# print("Inside MSE, Inside MPE: ", inside_mse, inside_mpe)
print(inside_mse, inside_mpe, sep=',')



# beta_hats = generate_beta_hats(genotypes, phenotypes)
# print(betas[0:5], beta_hats[0:5])

# beta_mse = [(a - b)**2 for a, b in zip(beta_hats, betas)]

# beta_mpe = 100 * np.sum([(a - b)/a for a, b in zip(betas, beta_hats)]) / n
# print(beta_mpe)

# print(np.mean(beta_mse))
