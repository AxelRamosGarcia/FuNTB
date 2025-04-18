import numpy as np
import networkx as nx
import random
from scipy import stats
from NetworkDataGeneration import remove_pp_and_pe, update_positions

def GetNetwork(groups_names, relevant_genes_data):
  """
  Constructs a network graph using the provided group names and relevant genes data.

  Parameters:
    group_names (list): A list of group names.
    relevant_genes_data (list): A list of dictionaries containing relevant genes data for each group.

  Returns:
    networkx.Graph: A network graph constructed from the group names and relevant genes data.

  Example:
    >>> groups = ['Group A', 'Group B', 'Group C']
    >>> genes_data = [
    ...   {'Gene1': 0.5, 'Gene2': 0.8, 'Gene3': 0.3},
    ...   {'Gene2': 0.6, 'Gene3': 0.2, 'Gene4': 0.9},
    ...   {'Gene1': 0.4, 'Gene3': 0.7, 'Gene4': 0.5}
    ... ]
    >>> network = GetNetwork(groups, genes_data)
    >>> network.nodes()
    ['Group A', 'Group B', 'Group C', 'Gene1', 'Gene2', 'Gene3', 'Gene4']
    >>> network.edges()
    [('Group A', 'Gene2'), ('Group A', 'Gene3'), ('Group B', 'Gene3'), ('Group B', 'Gene4'), ('Group C', 'Gene3'), ('Group C', 'Gene4')]

  Note:
    - The function uses the 'networkx' library to construct the network graph.
    - The 'relevant_genes_data' list contains dictionaries where the keys are gene names and the values are relevant data associated with each gene.
    - The function connects the group names with genes that are common across multiple groups based on the 'relevant_genes_data'.
    - The function removes nodes without edges and nodes starting with "ppe", "pp_", "pe_", "PPE", "PP_", or "PE_" using the 'remove_pp_and_pe' function.

  """
  Variation_Graph = nx.Graph()
  Variation_Graph.add_nodes_from(groups_names)

  gropus_genes = []

  for i in range(len(groups_names)):
    lower_bound = (len(groups_names)-1)*(i)
    upper_bound = (i+1)*(len(groups_names)-1)
    genes = set(relevant_genes_data[lower_bound].keys())
    for j in range(lower_bound + 1, upper_bound): # upper_bound-1
      genes = set(genes.intersection(relevant_genes_data[j].keys()))
    gropus_genes.append(genes)

    edges = [(groups_names[i], gene) for gene in genes]
    Variation_Graph.add_edges_from(edges)
    Variation_Graph.add_nodes_from(genes)

  # Find nodes without edges
  nodes_without_edges = [node for node in Variation_Graph.nodes() if Variation_Graph.degree[node] == 0]

  # Remove nodes without edges
  Variation_Graph.remove_nodes_from(nodes_without_edges)

  # Find nodes without edges
  pp_ppe_nodes = remove_pp_and_pe(Variation_Graph.nodes())

  # Remove nodes without edges
  Variation_Graph.remove_nodes_from(pp_ppe_nodes)

  return Variation_Graph, gropus_genes

def GetGroupsGenesInformation(groups_names, Groups_genes_network, relevant_genes_data):
  """
    Retrieve information about genes in different groups.

    Args:
        groups_names (list): A list of group names.
        Groups_genes_network (list): A list of lists representing the gene network for each group.
        relevant_genes_data (list): A list of dictionaries containing relevant gene information.

    Returns:
        dict: A dictionary mapping group names to gene information.

  """
  group_genes_information = {}
  for i, group_name in enumerate(groups_names):
    group_genes_information[group_name] = {}
    for gene in Groups_genes_network[i]:
      lower_bound = (len(groups_names)-1)*(i) # Start-reading position
      upper_bound = (i+1)*(len(groups_names)-1) # End-reading position
      group_genes_information[group_name][gene] = {}
      for j in range(lower_bound, upper_bound):
        group_genes_information[group_name][gene] = update_positions(group_genes_information[group_name][gene], relevant_genes_data[j][gene])
  return group_genes_information

def FitnessScore(Group_genes_information, Variations_network):
    """
    Calculates the fitness score based on the group gene information and updates the variation network.

    Parameters:
        Group_genes_information (dict): A dictionary containing information about group genes.
        Variations_network (networkx.Graph): The variation network to be updated.

    Returns:
        None
    """

    # Initialize an empty dictionary to store fitness score for each gene
    Fitness_score_genes_info = {}

    # Calculate fitness score for each gene in each group
    for group, genes in Group_genes_information.items():
        for gene in genes:
            if gene not in Fitness_score_genes_info:
                Fitness_score_genes_info[gene] = Group_genes_information[group][gene]
            else:
                Fitness_score_genes_info[gene] = update_positions(Fitness_score_genes_info[gene], Group_genes_information[group][gene])
    
    # Update the variation network with the calculated fitness scores
    for gene, variations in Fitness_score_genes_info.items():
        Variations_network.nodes[gene]['size'] = sum(variations.values()) / len(variations)

def GetPSVComponents(nodes, groups_names):
  """
    Calculates the z-scores of node sizes for nodes that are not in groups_names.

    Parameters:
        nodes (networkx.Graph): The graph object containing the nodes.

    Returns:
        fs (list): A list of node sizes for nodes not in groups_names.
        z_scores (numpy.ndarray): An array of z-scores calculated from fs.
  """
  fs = []
  zsc = []
  new_nodes = []

  for node in nodes:
    if node not in groups_names:
      new_nodes.append(node)
      fs.append(nodes[node]['size'])
    else:
       pass
  z_scores = stats.zscore(fs)
  return fs, z_scores, new_nodes

def normalize(values_to_normalize, method):
  """
    Normalizes a given array of values using different normalization methods.

    Parameters:
        values_to_normalize (list or numpy.ndarray): The array of values to be normalized.
        method (str): The normalization method to be used. Supported options are 'min_max', 'zsc', 'log', and any other value for 'robust' scaling.

    Returns:
        numpy.ndarray: The normalized values based on the specified method.

    Raises:
        None

    Example:
        values = [1, 2, 3, 4, 5]
        normalized_values = normalize(values, 'min_max')
        print(normalized_values)
        # Output: [0.0, 0.25, 0.5, 0.75, 1.0]

  """
  values_to_normalize = np.array(values_to_normalize)
  # Min-Max Scaling
  if method == 'min_max':
    normalized_ps_values = (values_to_normalize - np.min(values_to_normalize)) / (np.max(values_to_normalize) - np.min(values_to_normalize))
  # Z-score Standardization
  elif method == 'zsc':
    normalized_ps_values = (values_to_normalize - np.mean(values_to_normalize)) / np.std(values_to_normalize)
  # Log Transformation
  elif method == 'log':
    normalized_ps_values = np.log(values_to_normalize + 1)
  # Robust Scaling
  else:
    normalized_ps_values = (values_to_normalize - np.median(values_to_normalize)) / (np.percentile(values_to_normalize, 75) - np.percentile(values_to_normalize, 25))

  return normalized_ps_values

def CalculatePSVScore(Varitations_network_nodes, fitess_Scores, z_scores, alpha_factor, groups_names, nodes):
  """
    Calculates the PSV (Protein-Specific Variation) score for each gene in the variation network.

    Parameters:
        Variations_network_nodes (list): A list of nodes in the variation network.
        fitness_Scores (list): A list of fitness scores for the genes.
        z_scores (list): A list of z-scores for the genes.
        alpha_factor (float): A weight factor for combining fitness scores and z-scores.
        groups_names (list): A list of group names to exclude from calculation.

    Returns:
        None
  """
  #nodes  = [node for node in list(Varitations_network_nodes) if node not in groups_names]
  pvs_scores = [(alpha_factor*fitess_Scores[i]) + ((1-alpha_factor)*z_scores[i]) for i in range(len(fitess_Scores))]
  pvs_scores = normalize(pvs_scores, 'log')
  z_scores = normalize(z_scores, 'log')

  Pareto_Front_nodes = {}

  for i, gene_node in enumerate(nodes):
    Varitations_network_nodes[gene_node]['size'] = pvs_scores[i]
    Pareto_Front_nodes[gene_node] = {}
    Pareto_Front_nodes[gene_node]['pvs'] = pvs_scores[i]
    Pareto_Front_nodes[gene_node]['sigma'] = z_scores[i]

  return Pareto_Front_nodes

# Alteration density score
def AlterationDensityScore(number_of_variations, number_of_positions):
   return round(number_of_variations/number_of_positions, 4)

# Dominant Altered gene score
def DominantAlteredGeneScore(max_gene_variations, number_of_variations):
   return round(max_gene_variations/number_of_variations, 4)

# Cluter Diversity Alteration Score
def ClusterDiversityAlterationScore(number_of_variations, total_groups_genes):
   return round(number_of_variations/total_groups_genes, 4)

# Comprehensive Alteration Impact Score (CAIS)
def CAIS(alteration_density_score, dominant_altered_gene_score, custer_diversity_score, alpha, beta, gamma):
  alpha = float(alpha)
  beta = float(beta)
  gamma = float(gamma)

  #assert alpha + beta + gamma == 1,"Coefficients addition is different to one"
  cais = alpha*(alteration_density_score) + beta*(dominant_altered_gene_score) + gamma*(custer_diversity_score)
  return round(cais, 4)


def GetTotalGroupGenes(gene, Variations_network):
   # Gtting the gene neighbours
   neighbours = set(list(Variations_network.neighbors(gene)))

   # Counting total genes neighbours
   gene_neighbours = set()

   for neighbour in neighbours:
      new_neighbours = set(list(Variations_network.neighbors(neighbour)))
      gene_neighbours = gene_neighbours.union(new_neighbours)

   return len(gene_neighbours)

# Get CAIS per gene
def GetCAIS(Group_genes_information, Variations_network, alpha, beta, gamma):
    # Initialize an empty dictionary to store fitness score for each gene
    genes_updated_info = {}

    # Genes positions update for score calculation
    for group, genes in Group_genes_information.items():
        for gene in genes:
            if gene not in genes_updated_info:
                genes_updated_info[gene] = Group_genes_information[group][gene]
            else:
                genes_updated_info[gene] = update_positions(genes_updated_info[gene], Group_genes_information[group][gene])
    
    # Scores calculation
    genes_scores = {}
    for gene, variations in genes_updated_info.items():
       genes_scores[gene] = {}
       # Alteration Density Score
       genes_scores[gene]['ads'] = AlterationDensityScore(sum(variations.values()), len(variations)) # number_of_variations, number_of_positions
       
       # Dominant Altered gene score
       genes_scores[gene]['dags'] = DominantAlteredGeneScore(max(variations.values()), sum(variations.values())) # max_gene_variations, number_of_variations

       # Cluter Diversity Alteration Score
       total_groups_genes = GetTotalGroupGenes(gene, Variations_network)
       genes_scores[gene]['cdas'] = ClusterDiversityAlterationScore(sum(variations.values()), total_groups_genes) # number_of_variations, total_groups_genes

       # Comprehensive Alteration Impact Score (CAIS)
       genes_scores[gene]['cais'] = CAIS(genes_scores[gene]['ads'], genes_scores[gene]['dags'], genes_scores[gene]['cdas'], alpha, beta, gamma)

    # Scores normalization
    for gene, scores in genes_scores.items():
       pass
    
    # Update the variation network with the calculated fitness scores
    '''
    for gene, variations in genes_updated_info.items():
      Variations_network.nodes[gene]['size'] = sum(variations.values()) / len(variations)
    '''
    return genes_scores

# Groups nodes colors
get_colors = lambda n: list(map(lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF),range(n)))