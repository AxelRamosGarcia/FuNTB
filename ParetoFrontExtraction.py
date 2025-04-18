# Find from the first to the n-th Pareto front

# Importing libraries

from paretoset import paretoset


# Get Desired Pareto Fronts

def GetParetoFronts(genes_scores_df, num_pareto_sets):

    # Initialize a list to store all Pareto fronts
    pareto_fronts = []

    # Keep track of the current dataset to apply Pareto set on
    current_data = genes_scores_df.copy()

    # Loop to get multiple Pareto fronts
    for i in range(num_pareto_sets):
        # Get the Pareto set for the current data
        mask = paretoset(current_data[['ADS', 'DAGS', 'CDAS', 'CAIS']], sense=["min", "max", "max", "max"])
        
        # Add the Pareto front to the list
        pareto_front = current_data[mask]
        pareto_fronts.append(pareto_front)
        
        # Remove the current Pareto set from the data
        current_data = current_data[~mask]
        
        # Stop if there is no more data
        if current_data.empty:
            break
    
    return pareto_fronts

def ResistanceGeneAssociation(pareto_fronts):

    GoI = set(["inhA", "katG", "fabG1", "rpoA", "rpoB", "rpoC", "rrl", "atpE", "ethA", "pncA", "fbiA",
           "fbiB", "fgd1", "ddn", "gibB", "rsmG", "rrs", "gyrB", "rpsL", "gidB", "embA", "embB",
           "embC"])
    
    resistance_dictionary = {"inhA": "isoniazid", "katG": "isoniazid", "fabG1": "isoniazid",
                         "rpoA": "rifampicin", "rpoB": "rifampicin", "rpoC": "rifampicin",
                         "gyra": "ofloxacin", "rrl": "aminoglycosides", "atpE": "bedaquiline",
                         "ethA": "ethionamide", "pncA": "pyrazinamide", "fbiA": "Delamid/Pretomanid",
                         "fbiB": "Delamid/Pretomanid", "fgd1": "Delamid/Pretomanid", "ddn": "delamid/pretomanid", 
                         "gibB": "streptomycin", "rsmG": "streptomycin", "rrs": "streptomycin",
                         "gyrB": "ofloxacin", "rpsL": "streptomycin", "gidB": "streptomycin",
                         "embA": "ethambutol", "embB": "ethambutol", "embC": "ethambutol"}
    
    DrugResistanceAssociationDictionary = {}

    for pareto_front_number, pareto_front in enumerate(pareto_fronts):
        current_genes = set(pareto_front["Gene"])
        resistance_and_sensitivity = GoI.intersection(current_genes)
        GoI = GoI.difference(resistance_and_sensitivity)

        # Drug resistance association

        for gene in resistance_and_sensitivity:
            DrugResistanceAssociationDictionary[gene] = resistance_dictionary[gene]
    
    return DrugResistanceAssociationDictionary