%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: PathwayAnalysis
% Author: Johanna Rapp
% Date: 20-Aug-2026
% Description: Combines step01_Add_Reactants_HitList.m and
% step02_Add_Pathways_toHitList.m. See detialed instructions in the
% functions.
% Required Input files: 
% hitList.xlsx,
% TableS10_Genes_Reactants_iML1515.xlsx,
% TableS11_Ecocyc_Pathways.xlsx
% TableS12_EcocycPathways_Reactants_iML1515.xlsx
% Table S10-12 from Rapp et al. 2026, Cell Systems
% (https://doi.org/10.1016/j.cels.2025.101518)
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PathwayAnalysis
hitList = readtable("hitList.xlsx");

% Add reactants from iML1515 model
[hitList] = step01_Add_Reactants_HitList(hitList);

% Add Pathway from Ecocyc database
[hitList] = step02_Add_Pathways_toHitList(hitList);

writetable(hitList, 'HitList_withPathways.xlsx')

%% Explanations of the Variable Names in Table "HitList_withPathways.xlsx"
% Pathways:	Subpathways from Ecocyc where gene is associated with. A gene 
% can be part of different pathways (pathways are separated with "/"). 
% "Not in PL" ("Not in Pathway List") - Gene is not included in Ecocyc subpathways.

% Pathways Abbreviation: Subpathways abbreviation from Ecocyc where gene is 
% associated with. A gene can be part of different pathways (pathways are 
% separated with "/").

% Position:	Position of the gene within the pathway (according to Table S11).

% Upstream Accumulation: 
% 0 - accumulating metabolite is not part of reactions upstream in the pathway(s). 
% 1- accumulating metabolite is part of a reaction upstream in the pathway(s), 
% but it is not a reactant of the CRISPRi targeted gene. 
% 2 - accumulating metabolite is part of reactions upstream in the pathway(s), 
% but it is also a reactant of the CRISPRi targeted gene.

% Downstream Accumulation:	
% 0 - accumulating metabolite is not part of reactions downstream in the pathway(s). 
% 1- accumulating metabolite is part of a reaction downstream in the pathway(s), 
% but it is not a reactant of the CRISPRi targeted gene. 
% 2 - accumulating metabolite is part of reactions downstream in the pathway(s), 
% but it is also a reactant of the CRISPRi targeted gene.
