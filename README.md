# Protein evolution

This folder contains code developed for our recent paper available at https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1007010

    Shou-Wen Wang, Anne-Florence Bitbol, Ned S. Wingreen, Revealing evolutionary constraints on proteins through sequence analysis, PLOS Comp. Bio. (2019)
    
    

## Structure

### Binary_Model

The single-mutation effect generated from our elastic network model of PDZ is provided in the Data folder.  An example is provided to illustrate how to run the code.

### PottsModel_with_21_states_plus_real_data_analysis

You can play around with the 21-state Potts model, or or identify the single-mutation effects (or sector sites) for a particular protein family.  The PDZ homolog sequences as well as the experimentally measured single-mutation effects from Ranganathan's group is provided in the Data folder. You can also use your own data.  An example is provided to illustrate how to run the code.
     
For analyzing the real protein data, run “dealing_Protein_data” for the ICOD method; and run “dealing_Protein_data_SCA” for the SCA method.  Options are detailed in the beginning of these files.  The default data to use is the PDZ data contained in the folder “./Data/PDZal.mat”

### Normal_mode_analysis_for_PDZ_elastic_network_model

The folder gives our code for analyzing the open-to-close conformational change of an elastic network model for PDZ domain, and compute the single-point mutation effect associated with this conformational change.  
