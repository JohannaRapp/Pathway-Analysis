%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step01_Add_Reactants_HitList
% Author: Johanna Rapp
% Date: 20-Aug-2026
% Description: Uses differential metabolite-strain pairs from FI-MS screen.
% Analyses if accumulating metabolite is a reactant of the corresponding
% enzyme  targetd by CRISPRi knockdown.
% Uses gene–protein–reaction associations from genome scale model iML1515
% (Supplementary Table S10, Rapp et al. 2026, Cell Systems). 
% If FI-MS hit consists of isobaric metabolites, it extracts the isobar
% whoch is a reactant  the corresponding enzyme  targetd by CRISPRi knockdown.
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hitList] = step01_Add_Reactants_HitList(hitList)
reacts = readtable("TableS10_Genes_Reactants_iML1515.xlsx");

for i = 1:size(hitList,1)
    strain_temp = string(hitList.Gene{i});
    MonoMass_temp = hitList.MonoMass(i);
    idx_reactants = find(reacts.gene==strain_temp);
    reacts_temp = reacts(idx_reactants,:);
    idx_match = find(abs(double(reacts_temp.MonoisotopicMass)-MonoMass_temp)<0.001);

    if ~isempty(idx_match)
        hitList.Reactant(i) = 1;
        hitList.ReactantAbb{i} = reacts_temp.Abbreviation{idx_match};
    else
        hitList.Reactant(i) = 0;
        hitList.ReactantAbb{i} = 'NaN';
    end

end
end