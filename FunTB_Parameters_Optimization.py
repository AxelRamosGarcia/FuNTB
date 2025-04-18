# α, 𝛽, and 𝛾 Optimization parameters

def calculate_fitness_score(pareto_fronts):
    """
    Calculate the fitness score based on the number of genes of interest found in each Pareto front.
    The score is weighted to prioritize genes found in earlier Pareto fronts.

    :param pareto_fronts: A dataframe containing genes and their associated Pareto front information. 
                          Each dataframe should represent genes in one Pareto front (PF1 to PF5).
    :return: The fitness score used for optimizing α, β, and γ parameters.

    Fitness score is calculated based on the formula:
    Fitness Score = 0.5 * PF1GFs + 0.3 * PF2GFs + 0.2 * PF3GFs + 0.07 * PF4GFs + 0.03 * PF5GFs
    
    where:
    - PF1GFs is the proportion of genes of interest found in Pareto Front 1,
    - PF2GFs is the proportion of genes of interest found in Pareto Front 2, 
    - PF3GFs is the proportion of genes of interest found in Pareto Front 3, and so on.

    The genes of interest (GoI) are predefined based on domain knowledge related to antibiotic resistance.
    These include: 
    ["inhA", "katG", "fabG1", "rpoA", "rpoB", "rpoC", "rrl", "atpE", "ethA", "pncA", "fbiA", "fbiB", 
     "fgd1", "ddn", "gibB", "rsmG", "rrs", "gyrB", "rpsL", "gidB", "embA", "embB", "embC"].

    """
    # Genes of interest

    GoI = set(["inhA", "katG", "fabG1", "rpoA", "rpoB", "rpoC", "rrl", "atpE", "ethA", "pncA", "fbiA","fbiB", "fgd1", "ddn", "gibB", "rsmG", "rrs",
               "gyrB", "rpsL", "gidB", "embA", "embB", "embC"])

    # Coefficients for each Pareto front
    weights = [0.5, 0.3, 0.2, 0.07, 0.03]

    # Total number of genes found across all Pareto fronts and fitness score
    total_genes = 0
    genes_found = []
    genes_of_interest_list = []

    for pareto_front in pareto_fronts:
        genes_of_interest =  pareto_front[pareto_front["Gene"].isin(GoI)]["Gene"].tolist()
        genes_of_interest_list.append(genes_of_interest)
        genes_found.append(len(genes_of_interest))
        total_genes += len(genes_of_interest)

    fitness_score = sum([coefficient*genes for coefficient,genes in zip(weights,genes_found)])/total_genes

    return fitness_score, genes_found, genes_of_interest_list