import scanpy as sc
import random
from models import*
'''
define three methods
m1 = 50 PCs taken from all_genes
m2 = 50 PCs taken from HVGs
m3 = no PCs (2500 HVGs only)
'''
seed = random.randint(1000, 9999)
print("seed is "+ str(seed))
'''
run & train models for Astrocytes
'''
#m1
PC_all_genes = sc.read('/Users/lily/PycharmProjects/scRNA/PC_all_genes.h5ad')
astro_lr_1 = logistic_regression(PC_all_genes, seed)
astro_rf_1 = random_forest(PC_all_genes, seed)
astro_dnn_1 = DeepNeuralNet(PC_all_genes,seed)

#m2
PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/PC_HVGs.h5ad')
astro_lr_2 = logistic_regression(PC_HVGs, seed)
astro_rf_2 = random_forest(PC_HVGs, seed)
astro_dnn_2 = DeepNeuralNet(PC_HVGs,seed)

#m3
no_PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/no_PC_HVGs.h5ad')
astro_lr_3 = logistic_regression(no_PC_HVGs, seed)
astro_rf_3 = random_forest(no_PC_HVGs, seed)
astro_dnn_3 = DeepNeuralNet(no_PC_HVGs,seed)

#print('check:', PC_HVGs == no_PC_HVGs)
'''
run & train models for DAN
'''
#m1
DAN_PC_all_genes = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_PC_all_genes.h5ad')
DAN_lr_1 = logistic_regression(DAN_PC_all_genes, seed)
DAN_rf_1 = random_forest(DAN_PC_all_genes, seed)
DAN_dnn_1 = DeepNeuralNet(DAN_PC_all_genes, seed)

#m2
DAN_PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_PC_HVGs.h5ad')
DAN_lr_2 = logistic_regression(DAN_PC_HVGs, seed)
DAN_rf_2 = random_forest(DAN_PC_HVGs, seed)
DAN_dnn_2 = DeepNeuralNet(DAN_PC_HVGs,seed)

#m3
DAN_no_PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_no_PC_HVGs.h5ad')
DAN_lr_3 = logistic_regression(DAN_no_PC_HVGs, seed)
DAN_rf_3 = random_forest(DAN_no_PC_HVGs, seed)
DAN_dnn_3 = DeepNeuralNet(DAN_no_PC_HVGs,seed)

####result disply###
print('Here is the display of results')
print('Astrocytes: ')
print('m1:')
print('logistic regression acc: ',astro_lr_1,'\n',
      'random forest acc: ', astro_rf_1,'\n',
      'deep neural net train, valid, test acc: ', astro_dnn_1
      )

print('m2:')
print('logistic regression acc: ',astro_lr_2,'\n',
      'random forest acc: ', astro_rf_2,'\n',
      'deep neural net train, valid, test acc: ', astro_dnn_2
      )

print('m3:')
print('logistic regression acc: ',astro_lr_3,'\n',
      'random forest acc: ', astro_rf_3,'\n',
      'deep neural net train, valid, test acc: ', astro_dnn_3
      )

print('DA neurons: ')
print('m1:')
print('logistic regression acc: ',DAN_lr_1,'\n',
      'random forest acc: ', DAN_rf_1,'\n',
      'deep neural net train, valid, test acc: ',  DAN_dnn_1
      )
print('m2:')
print('logistic regression acc: ',DAN_lr_2, '\n',
      'random forest acc: ', DAN_rf_2, '\n',
      'deep neural net train, valid, test acc: ',  DAN_dnn_2
      )
print('m3:')
print('logistic regression acc: ',DAN_lr_3, '\n',
      'random forest acc: ', DAN_rf_3, '\n',
      'deep neural net train, valid, test acc: ',  DAN_dnn_3
      )