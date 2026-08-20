%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step02_Add_Pathways_toHitList
% Author: Johanna Rapp
% Date: 20-Aug-2026
% Description: Uses differential metabolite-strain pairs from FI-MS screen.
% Pathways of Ecocyc are added to CRISPRi targets. If an enzyme targeted by
% CRISPRi via gene knockdown is involved in different pathways all pathways
% are shown.
% Uses pathway list from Ecocyc summarized in Supplementary Table S12, Rapp
% et al. 2026, Cell Systems). Enzymes in the pathwy list were assigned with
% the corresponding reactants (substrates and products) of the iML1515
% model.
% Analyses if a metabolite accumulation is a reactant response, a response
% inside or outside the targeted pathway.
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hitList] = step02_Add_Pathways_toHitList(hitList)
pathways_reactant = readtable("TableS12_EcocycPathways_Reactants_iML1515.xlsx");
pathways = readtable("TableS11_Ecocyc_Pathways.xlsx");

for i = 1:size(hitList,1)
    clearvars -except hitList pathways_reactant pathways i
    idx = find(string(hitList.Gene{i})==pathways.GeneName);

    if ~isempty(idx)
        % extract unique pathways
        pathways_temp = pathways(idx,:);
        hitList.Pathways {i} = join(pathways_temp.Pathway,'/ ');
        hitList.PathwaysAbb {i} = join(pathways_temp.Abb,'/ ');
        hitList.Position {i} = join(string(pathways_temp.Position), '/ '); % position of gene in the pathway
        
        for a = 1:size(pathways_temp,1) % loop over all pathways where enzyme is part of
           
            idx_reactants = find(pathways_reactant.Pathway==string(pathways_temp.Pathway(a)));
           
            idx_upstream_reactants = find(pathways_reactant.Position(idx_reactants)<pathways_temp.Position(a));
            upstream = pathways_reactant(idx_reactants(idx_upstream_reactants),:);
  
            idx_downstream_reactants = find(pathways_reactant.Position(idx_reactants)>pathways_temp.Position(a));
            downstream = pathways_reactant(idx_reactants(idx_downstream_reactants),:);
            
            mass_temp = hitList.MonoMass(i);

            % check if mass of hitList entry matches upstream metabolites
            if ~isempty(idx_upstream_reactants)
                idx_match = find(abs(double(upstream.mass)-mass_temp)<0.001);
    
                if ~isempty(idx_match)
                    upstream_acc(a,1) = 1;
                    position(a,1) = join(string(upstream.Position(idx_match)),'-'); % position of gene in pathway where hit is an reactant
                    upstream_gene(a,1) = join(string(upstream.GeneName(idx_match)),'-'); % gene where hit is an reactant
                else
                    upstream_acc(a,1) = 0;
                    position(a,1) = "no up acc";
                    upstream_gene(a,1) = "no up acc";
                end
            else
                upstream_gene(a,1) = "no up gene";
                position(a,1) = "no up gene";
                upstream_acc(a,1) = 0;
            end

            % check if mass of hitList entry matches downstream metabolites
            if ~isempty(idx_downstream_reactants)
                idx_match_down = find(abs(double(downstream.mass)-mass_temp)<0.001);
    
                if ~isempty(idx_match_down)
                    downstream_acc(a,1) = 1;
                    position_down(a,1) = join(string(downstream.Position(idx_match_down)),'-'); % position of gene in pathway where hit is an reactant
                    downstream_gene(a,1) = join(string(downstream.GeneName(idx_match_down)),'-'); % gene where hit is an reactant
                else
                    downstream_acc(a,1) = 0;
                    position_down(a,1) = "no down acc";
                    downstream_gene(a,1) = "no down acc";
                end
           else
                downstream_gene(a,1) = "no down gene";
                position_down(a,1) = "no down gene";
                downstream_acc(a,1) = 0;
            end
            clearvars idx_reactants idx_upstream_reactants upstream idx_downstream_reactants...
                downstream mass_temp idx_match idx_match_down
        end % end pathway loop
        
        hitList.UpstreamAcc(i) = join(string(upstream_acc),'/ ');
        hitList.UpstreamPos(i) = join(string(position),'/ ');
        hitList.UpstreamGene(i) = join(string(upstream_gene),'/');
        hitList.DownstreamAcc(i) = join(string(downstream_acc),'/ ');
        hitList.DownstreamPos(i) = join(string(position_down),'/ ');
        hitList.DownstreamGene(i) = join(string(downstream_gene),'/');
        
    else % if gene from hitList is not found in pathway list
        hitList.Pathways {i} = "Not in PL";
        hitList.PathwaysAbb {i} = "Not in PL";
        hitList.UpstreamAcc(i) = "Not in PL";
        hitList.UpstreamPos(i) = "Not in PL";
        hitList.UpstreamGene(i) = "Not in PL";
        hitList.DownstreamAcc(i) = "Not in PL";
        hitList.DownstreamPos(i) = "Not in PL";
        hitList.DownstreamGene(i) = "Not in PL";
    end
end % end hitList loop

clearvars -except hitList
%%
for i = 1:size(hitList,1)
    temp_up = strsplit(hitList.UpstreamAcc(i),'/ ');
    unique_temp_up = unique(temp_up);
    zero_up = find(unique_temp_up=="0");
    one_up = find(unique_temp_up=="1");
    if ~isempty(zero_up) & ~isempty(one_up)
        hitList.UpstreamAccAll(i,1) = "1";
    elseif ~isempty(zero_up)
        hitList.UpstreamAccAll(i,1) = "0";
    elseif ~isempty(one_up)
        hitList.UpstreamAccAll(i,1) = "1";
    elseif isempty(zero_up) & isempty(one_up)
        hitList.UpstreamAccAll(i,1) = "Not in PL";
    end

    temp_down = strsplit(hitList.DownstreamAcc(i),'/ ');
    unique_temp_down = unique(temp_down);
    zero_down = find(unique_temp_down=="0");
    one_down = find(unique_temp_down=="1");
    if ~isempty(zero_down) & ~isempty(one_down)
        hitList.DownstreamAccAll(i,1) = "1";
    elseif ~isempty(zero_down)
        hitList.DownstreamAccAll(i,1) = "0";
    elseif ~isempty(one_down)
        hitList.DownstreamAccAll(i,1) = "1";
    elseif isempty(zero_down) & isempty(one_down)
        hitList.DownstreamAccAll(i,1) = "Not in PL";
    end

end
end
