from preprocess_indel_files import preprocess_indel_files
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

def compute_summary_statistics(name_genes_unique, name_genes_grna_unique, name_indel_type_unique, indel_count_matrix):
  # Compute TSNE on indels that commonly occur across all genes and grna
  number_of_files_per_indel = []
  for i in range(indel_count_matrix.shape[0]):
    number_of_files_per_indel.append(np.count_nonzero(indel_count_matrix[i]))
  print "The twenty most commonly occurring indels are:"
  for i in range(20):
    print name_indel_type_unique[np.argsort(number_of_files_per_indel)[::-1][i]]
  #
  indel_count_matrix_small = np.copy(indel_count_matrix)
  indel_count_matrix_small = indel_count_matrix_small[np.argsort(number_of_files_per_indel)[::-1][0:200]]
  X = TSNE(n_components=2, random_state=0).fit_transform(np.transpose(indel_count_matrix_small))
  plt.scatter(X[:, 0], X[:, 1])
  plt.savefig('all-genes.pdf')
  plt.clf()

  # Plot heat map of cosine distances
  heat_map_inner_prod = np.matmul(np.transpose(indel_count_matrix_small), indel_count_matrix_small)
  # The next 2 lines change the plot from inner product to cosine distance
  #  For just the inner product, comment these two lines
  col_wise_norms = np.expand_dims(np.linalg.norm(indel_count_matrix_small, axis = 0), axis = 1)
  heat_map_inner_prod = np.divide(heat_map_inner_prod, np.matmul(col_wise_norms, np.transpose(col_wise_norms)))
  fig, axis = plt.subplots()
  heat_map = axis.pcolor(heat_map_inner_prod, cmap = plt.cm.Blues)
  axis.set_yticks(np.arange(heat_map_inner_prod.shape[0])+0.5, minor=False)
  axis.set_xticks(np.arange(heat_map_inner_prod.shape[1])+0.5, minor=False)
  axis.invert_yaxis()
  axis.set_yticklabels(name_indel_type_unique, minor=False)
  axis.set_xticklabels(name_indel_type_unique, minor=False)
  plt.xticks(fontsize=2.5, rotation=90)
  plt.yticks(fontsize=2.5)
  plt.colorbar(heat_map)
  plt.savefig('inner_product_heat_map.pdf')
  plt.clf()

  ##
  # Threshold indels for each column
  threshold = 0.05
  indel_count_matrix_threshold = indel_count_matrix/np.reshape(np.sum(indel_count_matrix, axis = 0), (1, -1))
  indel_count_matrix_threshold[indel_count_matrix_threshold <= threshold] = 0.0
  indel_count_matrix_threshold[indel_count_matrix_threshold == np.inf] = 0.0
  indel_count_matrix_threshold[indel_count_matrix_threshold == -np.inf] = 0.0
  #
  # For each gene-grna pair, count the number of indels above threshold
  count_insertions_gene_grna = np.zeros(len(name_genes_grna_unique), dtype = int)
  count_deletions_gene_grna = np.zeros(len(name_genes_grna_unique), dtype = int)
  for i in range(len(name_genes_grna_unique)):
    for j in range(indel_count_matrix_threshold.shape[0]):
      if np.sum(indel_count_matrix_threshold[j][3*i:3*i+3]) > 0:
        if name_indel_type_unique[j].find('I') != -1:
          count_insertions_gene_grna[i] += 1
        if name_indel_type_unique[j].find('D') != -1:
          count_deletions_gene_grna[i] += 1

  # Save the output
  indel_family_count_gene_grna = open('indel_family_count_gene_grna.txt', 'w')
  indel_family_count_gene_grna.write("gene_grna_name,insertion_count,deletion_count\n")
  for i in range(len(name_genes_grna_unique)):
    indel_family_count_gene_grna.write("%s,%d,%d\n" % (name_genes_grna_unique[i], count_insertions_gene_grna[i], count_deletions_gene_grna[i]))

  return count_insertions_gene_grna, count_deletions_gene_grna