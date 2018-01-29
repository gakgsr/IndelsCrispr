from preprocess_indel_files import preprocess_indel_files
from compute_summary_statistic import compute_summary_statistics
import numpy as np
from sklearn import linear_model
from sklearn import metrics

def convert_to_int(nucleotide):
  nucleotide_array = ['A', 'C', 'G', 'T']
  return nucleotide_array.index(nucleotide) + 1

def load_gene_sequence(sequence_file_name, name_genes_grna_unique):
  # Create numpy matrix of size len(name_genes_grna_unique) * 23, to store the sequence
  sequence_pam_per_gene_grna = np.zeros((len(name_genes_grna_unique), 23), dtype = int)
  # Obtain the grna and PAM sequence corresponding to name_genes_grna_unique
  with open(sequence_file_name) as f:
    for line in f:
      line = line.replace('"', '')
      line = line.replace(' ', '')
      line = line.replace('\n', '')
      l = line.split(',')
      if l[1] + '-' + l[0] in name_genes_grna_unique:
        index_in_name_genes_grna_unique = name_genes_grna_unique.index(l[1] + '-' + l[0])
        for i in range(20):
          sequence_pam_per_gene_grna[index_in_name_genes_grna_unique, i] = convert_to_int(l[2][i])
        for i in range(3):
          sequence_pam_per_gene_grna[index_in_name_genes_grna_unique, 20 + i] = convert_to_int(l[3][i])
  print "Sum of elements in input to regression is %d" % np.sum(sequence_pam_per_gene_grna)
  return sequence_pam_per_gene_grna

def perform_logistic_regression(sequence_pam_per_gene_grna, count_insertions_gene_grna, count_deletions_gene_grna):
  log_reg = linear_model.LogisticRegression(C=1e5)
  threshold_insertions = 1
  count_insertions_gene_grna_binary = np.copy(count_insertions_gene_grna)
  count_insertions_gene_grna_binary[count_insertions_gene_grna >= threshold_insertions] = 1
  count_insertions_gene_grna_binary[count_insertions_gene_grna < threshold_insertions] = 0
  print "Number of positive testing samples in insertions is %f" % np.sum(count_insertions_gene_grna_binary[50:])
  print "Number of positive training samples in insertions is %f" % np.sum(count_insertions_gene_grna_binary[:50])
  log_reg.fit(sequence_pam_per_gene_grna[:50], count_insertions_gene_grna_binary[:50])
  log_reg_pred = log_reg.predict(sequence_pam_per_gene_grna[50:])
  print "Accuracy Score for Insertions: %f" % metrics.accuracy_score(count_insertions_gene_grna_binary[50:], log_reg_pred)
  threshold_deletions = 3
  count_deletions_gene_grna_binary = np.copy(count_deletions_gene_grna)
  count_deletions_gene_grna_binary[count_deletions_gene_grna >= threshold_deletions] = 1
  count_deletions_gene_grna_binary[count_deletions_gene_grna < threshold_deletions] = 0
  print "Number of positive testing samples in deletions is %f" % np.sum(count_deletions_gene_grna_binary[50:])
  print "Number of positive training samples in deletions is %f" % np.sum(count_deletions_gene_grna_binary[:50])
  log_reg.fit(sequence_pam_per_gene_grna[:50], count_deletions_gene_grna_binary[:50])
  log_reg_pred = log_reg.predict(sequence_pam_per_gene_grna[50:])
  print "Accuracy Score for Deletions: %f" % metrics.accuracy_score(count_deletions_gene_grna_binary[50:], log_reg_pred)

data_folder = "../IndelsData/"
name_genes_unique, name_genes_grna_unique, name_indel_type_unique, indel_count_matrix = preprocess_indel_files(data_folder)
count_insertions_gene_grna, count_deletions_gene_grna = compute_summary_statistics(name_genes_unique, name_genes_grna_unique, name_indel_type_unique, indel_count_matrix)
sequence_file_name = "sequence_pam_gene_grna.csv"
sequence_pam_per_gene_grna = load_gene_sequence(sequence_file_name, name_genes_grna_unique)
perform_logistic_regression(sequence_pam_per_gene_grna, count_insertions_gene_grna, count_deletions_gene_grna)