import numpy as np
import glob

def preprocess_indel_files(data_folder):
  # First process the files to glean the names of the genes and the different indels
  name_genes = []
  name_genes_grna = []
  name_indel_type = []
  for each_file in glob.glob(data_folder + "counts-*.txt"):
    with open(each_file) as f:
      i = 0
      process_file = False
      add_file = False
      for line in f:
        line = line.replace('"', '')
        line = line.replace('\n', '')
        l = line.split(',')
        if i == 0:
          if len(l) == 4:
            process_file = True
            curr_gene_name = each_file[len(data_folder) + 7:-4].split('-')[0]
            curr_gene_grna_name = curr_gene_name + '-' + l[0].split('-')[-2]
        if i > 0 and process_file:
          indel_type = ''
          # Some positions are of the form: "-23:-21D,-19:-15D", which get split by the process when we call split()
          # We try to account for such things in this space
          for j in range(0, len(l) - 4):
            indel_type += l[j]
          # We only consider I or D
          if line.find('I') != -1 or line.find('D') != -1:
            name_indel_type.append(indel_type)
            if not add_file:
              name_genes.append(curr_gene_name)
              name_genes_grna.append(curr_gene_grna_name)
              add_file = True
        i += 1

  # Take the unique values, in sorted order
  name_genes_unique = list(set(name_genes))
  name_genes_grna_unique = list(set(name_genes_grna))
  name_indel_type_unique = list(set(name_indel_type))
  name_genes_unique.sort()
  name_genes_grna_unique.sort()
  name_indel_type_unique.sort()

  ##
  # Then process the files again to get the actual counts from only the desired files, and from the desired rows and columns
  indel_count_matrix = np.zeros((len(name_indel_type_unique), 3*len(name_genes_grna_unique)))
  for each_file in glob.glob(data_folder + "counts-*.txt"):
    with open(each_file) as f:
      i = 0
      process_file = False
      for line in f:
        line = line.replace('"', '')
        line = line.replace('\n', '')
        l = line.split(',')
        if i == 0:
          if len(l) == 4 and each_file[len(data_folder) + 7:-4].split('-')[0] + '-' + l[0].split('-')[-2] in name_genes_grna_unique:
            process_file = True
            curr_gene_name = each_file[len(data_folder) + 7:-4].split('-')[0]
            col_index = name_genes_grna_unique.index(curr_gene_name + '-' + l[0].split('-')[-2])
        if i > 0 and process_file:
          indel_type = ''
          # Some positions are of the form: "-23:-21D,-19:-15D", which get split by the process when we call split()
          # We try to account for such things in this space
          for j in range(0, len(l) - 4):
            indel_type += l[j]
          # We ignore SNV, others, and no variants
          if line.find('I') != -1 or line.find('D') != -1:
            row_index = name_indel_type_unique.index(indel_type)
            for j in range(3):
              if l[j + len(l) - 4] != 'NA':
                indel_count_matrix[row_index, 3*col_index + j] = float(l[j + len(l) - 4])
        i += 1

  ##
  # Process the proportions file to get the proportions data
  indel_prop_matrix = np.zeros((len(name_indel_type_unique), 3*len(name_genes_grna_unique)))
  for each_file in glob.glob(data_folder + "proportions-*.txt"):
    with open(each_file) as f:
      i = 0
      process_file = False
      for line in f:
        line = line.replace('"', '')
        line = line.replace('\n', '')
        l = line.split(',')
        if i == 0:
          if len(l) == 4 and each_file[len(data_folder) + 12:-4].split('-')[0] + '-' + l[0].split('-')[-2] in name_genes_grna_unique:
            process_file = True
            curr_gene_name = each_file[len(data_folder) + 12:-4].split('-')[0]
            col_index = name_genes_grna_unique.index(curr_gene_name + '-' + l[0].split('-')[-2])
        if i > 0 and process_file:
          indel_type = ''
          # Some positions are of the form: "-23:-21D,-19:-15D", which get split by the process when we call split()
          # We try to account for such things in this space
          for j in range(0, len(l) - 4):
            indel_type += l[j]
          # We ignore SNV, others, and no variants
          if line.find('I') != -1 or line.find('D') != -1:
            row_index = name_indel_type_unique.index(indel_type)
            for j in range(3):
              if l[j + len(l) - 4] != 'NA':
                indel_prop_matrix[row_index, 3*col_index + j] = float(l[j + len(l) - 4])
        i += 1

  ##
  # Save the indel counts, indel type, and gene-grna name information
  indel_type_file = open('indel_type.txt', 'w')
  for indel_type in name_indel_type_unique:
    indel_type_file.write("%s\n" % indel_type)
  genes_grna_file = open('genes_grna.txt', 'w')
  for genes_grna in name_genes_grna_unique:
    genes_grna_file.write("%s\n" % genes_grna)
  np.savetxt("indel_count_matrix.txt", indel_count_matrix)
  np.savetxt("indel_prop_matrix.txt", indel_prop_matrix)

  return name_genes_unique, name_genes_grna_unique, name_indel_type_unique, indel_count_matrix, indel_prop_matrix