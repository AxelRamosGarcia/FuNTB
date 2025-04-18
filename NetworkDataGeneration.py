import re
from itertools import permutations

def update_positions(A,B):
  '''
  The purpose of this function is to return a dictionary that contains every gene variation information gene:{position: frequency}, calculated along with every individual within the sample

  Args:
    'A' (list): Current dicitonary of position's variations
    'B' (list): New dictionary of positions to update A

  Returns:
    dict: Comprised by position's variations and its frequency, updating A list keys and values using B list.
  '''
  return {x: A.get(x, 0) + B.get(x, 0) for x in set(A).union(B)}

def generate_groups_permutations(groups):
  """
    Generates all possible permutations of pairs within a given list of groups.

    Parameters:
        groups (list): A list of groups.

    Returns:
        list: A list of tuples representing all possible permutations of pairs within the groups.

    Example:
        groups = ['A', 'B', 'C']
        permutations = generate_groups_permutations(groups)
        print(permutations)
        # Output: [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]
  """
  return list(permutations([i for i in range(len(groups))], 2))

def GenesData(samples, samples_variations):
  """
    Retrieves gene data from samples_variations for the given samples.

    Parameters:
        samples (list): A list of sample names.
        samples_variations (dict): A dictionary containing sample variations data.

    Returns:
        dict: A dictionary containing gene data aggregated from samples_variations.

    Example:
        samples = ['Sample1', 'Sample2', 'Sample3']
        variations = {
            'Sample1': {
                'GeneA': {'variation positions': [1, 2, 3]},
                'GeneB': {'variation positions': [4, 5]},
            },
            'Sample2': {
                'GeneA': {'variation positions': [2, 3]},
                'GeneC': {'variation positions': [6]},
            },
            'Sample3': {
                'GeneB': {'variation positions': [4, 5]},
                'GeneD': {'variation positions': [7, 8]},
            },
        }
        gene_data = GenesData(samples, variations)
        print(gene_data)
        # Output: {'GeneA': [1, 2, 3], 'GeneB': [4, 5], 'GeneC': [6], 'GeneD': [7, 8]}

  """
  genes_data = {}
  for sample in samples:
    genes = list(samples_variations[sample].keys())
    for gene in genes:
      if gene.startswith("ppe") or gene.startswith("pp_") or gene.startswith("pe_") or gene.startswith("PPE") or gene.startswith("PP_") or gene.startswith("PE_"):
        pass
      else:
        if gene not in genes_data:
          genes_data[gene] = samples_variations[sample][gene]['variation_positions']
        else:
          genes_data[gene] = update_positions(genes_data[gene],samples_variations[sample][gene]['variation_positions'])
  return genes_data

def relevant_positions(experimentalPositions,controlPositions):
  """
    Filters out positions from experimentalPositions that are not present in controlPositions.

    Parameters:
        experimentalPositions (dict): A dictionary of experimental positions.
        controlPositions (set): A set of control positions.

    Returns:
        dict: A dictionary containing relevant positions from experimentalPositions.

    Example:
        experimental = {'A': 1, 'B': 2, 'C': 3, 'D': 4}
        control = {'B', 'D'}
        relevant = relevant_positions(experimental, control)
        print(relevant)
        # Output: {'A': 1, 'C': 3}

  """
  rpos = {}
  for position in experimentalPositions.keys():
    if position not in controlPositions:
      rpos[position] = experimentalPositions[position]
  return rpos

def uncommon_positions(experimental_group,control_group):
  """
    Returns a dictionary containing the positions of genes in the experimental group that are uncommon or absent
    in the control group. If a gene is present in both groups, the positions of that gene are compared and stored.

    Args:
        experimental_group (dict): Dictionary representing the experimental group.
            Keys represent genes, and values contain relevant information including positions.
        control_group (dict): Dictionary representing the control group.
            Keys represent genes, and values contain relevant information including positions.

    Returns:
        dict: A dictionary containing the positions of genes in the experimental group that are uncommon or absent
        in the control group. If a gene is present in both groups, the positions of that gene are compared and stored.
        The keys of the dictionary are genes, and the values are either the positions or the result of the comparison.
  """
  gene_pos = {}
  genes_of_interes = experimental_group.keys()
  for gene in genes_of_interes:
    if gene in control_group.keys():
      gene_pos[gene] = relevant_positions(experimental_group[gene],control_group[gene])
    else:
      gene_pos[gene] = experimental_group[gene]
  return gene_pos

def relevant_genes(uncommon_genes):
  """
    Filters out genes from uncommon_genes that have no variations.

    Parameters:
        uncommon_genes (dict): A dictionary of uncommon genes.

    Returns:
        dict: A dictionary containing relevant genes from uncommon_genes along with additional information.

    Example:
        uncommon = {'GeneA': {'Pos1': 2, 'Pos2': 1}, 'GeneB': {}, 'GeneC': {'Pos3': 4}}
        relevant = relevant_genes(uncommon)
        print(relevant)
        # Output: {'GeneA': [2, {'Pos1': 2, 'Pos2': 1}, 3], 'GeneC': [1, {'Pos3': 4}, 4]}

  """
  relevant_genes = {}
  for gene in uncommon_genes.keys():
    if len(uncommon_genes[gene]) > 0:
      relevant_genes[gene] = uncommon_genes[gene]
  return relevant_genes

def remove_pp_and_pe(gene_nodes):
  """
  Removes gene nodes that start with "ppe", "pp_", "pe_", "PPE", "PP_", or "PE_" from a list.

  Parameters:
    gene_nodes (list): A list of gene nodes.

  Returns:
    list: A list containing the gene nodes that do not start with the specified prefixes.

  Example:
    >>> gene_list = ['ppeGene', 'pp_Gene', 'pe_Gene', 'PPEGene', 'PP_Gene', 'PE_Gene', 'OtherGene']
    >>> remove_pp_and_pe(gene_list)
    ['OtherGene']
  """
  nodes2remove = []
  for gene in gene_nodes:
    if gene.startswith("ppe") or gene.startswith("pp_") or gene.startswith("pe_") or gene.startswith("PPE") or gene.startswith("PP_") or gene.startswith("PE_"):
      nodes2remove.append(gene)
  return nodes2remove