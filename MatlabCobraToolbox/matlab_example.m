%% use E. coli iJO1366 model
if exist('iJO1366.mat', 'file')
    fprintf(' iJO1366.mat found.\n');
else
    % download if not in the search path (need curl)
    [status_curl, result_curl] = system('curl --max-time 15 -s -k -L --head http://bigg.ucsd.edu/static/models/iJO1366.mat');
    % check if the URL exists
    if status_curl == 0 && ~isempty(strfind(result_curl, ' 200'))
        status_curlDownload = system('curl --max-time 60 -O -L http://bigg.ucsd.edu/static/models/iJO1366.mat');
        if status_curlDownload == 0
            fprintf(' + Downloaded:      http://bigg.ucsd.edu/static/models/iJO1366.mat\n');
        end
    else
        fprintf(' > The URL http://bigg.ucsd.edu/static/models/iJO1366.mat cannot be reached.\n');
    end
end

model = load('iJO1366.mat');
model = model.iJO1366;
% add biomass metabolites
model = addMetabolite(model, 'biomass_core[e]');
model.S(findMetIDs(model, 'biomass_core[e]'), findRxnIDs(model, 'BIOMASS_Ec_iJO1366_core_53p95M')) = 1;
model = addMetabolite(model, 'biomass_wt[e]');
model.S(findMetIDs(model, 'biomass_wt[e]'), findRxnIDs(model, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = 1;
% add export reactions
model = addReaction(model, {'EX_biomass_core(e)', 'Core biomass export'}, 'biomass_core[e] ->');
model = addReaction(model, {'EX_biomass_wt(e)', 'WT biomass export'}, 'biomass_wt[e] ->');

%% compute chemical formulae
% get all extracellular metabolites
metEx = any(model.S(:, sum(model.S ~= 0, 1) <= 1), 2);
% original chemical formulae
metFormulae = model.metFormulas;
% use all extracellular metabolites to infer all other metabolites
% (the best practice is to use as many known metabolites as one can to
%  reveal all hidden inconsistencies, but all extracellular metabolites + 
%  sink/demand metabolites is sufficient to give defined results for all metabolites)
model.metFormulas(~metEx) = {''};
[model2, metCompute, ele, metEle, rxnBal, S_fill, solInfo, N, LP] = computeMetFormulae(model, [], model.rxns(sum(model.S ~= 0, 1) > 1));
% compare with the original results
yn = compareMetFormulas(metFormulae, model2.metFormulas);
fprintf('%d / %d metabolites have the same chemical formulae computed.\n', sum(yn), numel(model.mets));
fprintf('Mets with formulae different from the original model:\n');
format short
disp([{'mets', 'Original formulae', 'Computed formulae'}; model.mets(~yn), metFormulae(~yn), model2.metFormulas(~yn)])
% **** It can be seen that the differences all originate from the specific elements used to name
%      each conserved moieties in the network. 

%% Match the computed formulae with the original chemical formulae
% Try to match the generic elements defined by computeMetFormulae to those defined in the orginal model

% Different ordering of the columns for conserved moieties when using EFM
% tool or the matlab null.m. Change manually if you choose a particular
% method rather than decided by the function

CMchoice = 3;  % if solved by the default matlab null.m
if exist('CalculateFluxModes.m', 'file')
    CMchoice = 2;  % if solved using EFMtool
end
fprintf('\nTry matching the original formulae ... \n')
match = {'Conserve_a_', 'R', 'R';
    'Conserve_aa', 'X', 'R';
    'Conserve_ab', 'X', 'R';
    'Conserve_ac', 'R', 'R';
    'Conserve_ad', 'R', 'R';
    'Conserve_ae', 'C22H31N7O17P3S', 'R';
    'Conserve_af', 'XH', 'R';
    'Conserve_ag', 'R', 'R';
    'Conserve_ah', 'C3H9N', 'R';
    'Conserve_ai', 'R', 'X';
    'Conserve_aj', 'C6H6N4', 'X';
    'Conserve_ak', 'R', 'X';
    'Conserve_al', 'COX', 'X';
    'Conserve_a', 'X', 'X';
    'Conserve_b', 'X', 'HX';
    'Conserve_c', 'X', 'X';
    'Conserve_d', 'R', 'HSR';
    'Conserve_e', 'R', 'O2S6R';
    'Conserve_f', 'R', 'COX';
    'Conserve_g', 'X', 'C22H31N7O17P3S';
    'Conserve_h', 'O2S6R', 'O2S6R';
    'Conserve_i', 'R', 'HSR';
    'Conserve_j', 'R', 'CN';
    'Conserve_k', 'R', 'C6H6N4';
    'Conserve_l', 'R', 'C3H9N';
    'Conserve_m', 'HSR', 'HOR';
    'Conserve_n', 'R', 'X';
    'Conserve_o', 'R', 'R';
    'Conserve_p', 'X', 'R';
    'Conserve_q', 'HSR', 'C10H17O10PR2';
    'Conserve_r', 'RHO', 'R';
    'Conserve_s', 'O2S6R', 'R';
    'Conserve_t', 'R', 'R';
    'Conserve_u', 'R', 'R';
    'Conserve_v', 'R', 'R';
    'Conserve_w', 'R', 'R';
    'Conserve_x', 'R', 'R';
    'Conserve_y', 'CN', 'R';
    'Conserve_z', 'C10H17O10PR2', 'R'};
metFormulae2 = model2.metFormulas;
for j = 1:size(match, 1)
    metFormulae2 = strrep(metFormulae2, match{j, 1}, match{j, CMchoice});
end
% standardize both sets of formulae
[~, ele1, metEle1] = checkEleBalance(metFormulae);
[~, ele2, metEle2] = checkEleBalance(metFormulae2);
metFormulae1 = convertMatrixFormulas(ele1, metEle1);
metFormulae2 = convertMatrixFormulas(ele2, metEle2);
% compare
yn2 = compareMetFormulas(metFormulae1, metFormulae2);
fprintf('%d / %d metabolites have the same chemical formulae computed.\n', sum(yn2), numel(model.mets));
fprintf('Mets with formulae different from the original model:\n');
disp([{'mets', 'Original formulae', 'Computed formulae'}; model.mets(~yn2), metFormulae1(~yn2), metFormulae2(~yn2)])

%% Look at the inconsistencies found in the results
fprintf('Minimum inconsistencies found when balancing each element:\n')
for j = 1:numel(solInfo.sol)
    fprintf('%s: %.4f\n', strjoin(solInfo.ele(solInfo.eleConnect(:, j)), ', '), solInfo.sol(j).minIncon.obj)
end
% **** should see zero inconsistency for all elements

%% Check molecular weights of the biomass metabolites
bm = findMetIDs(model,{'biomass_core[e]'; 'biomass_wt[e]'});
fprintf('Molecular weight (g/mmol):\nbiomass_core[e]\t%.12f\nbiomass_wt[e]\t%.12f\n', ...
    MW(model2.metFormulas{bm(1)}) / 1000, MW(model2.metFormulas{bm(2)}) / 1000)
% **** You shall see that the molecular weight of the biomass produced in silico in the E. coli
%      model is pretty accurate.

%% Compute the range for the molecular weight of a metabolite

% Create an inconsistency in the model.

% Take for example the metabolite man6pglyc_c which has only two connecting reactions 
% so that the algorithm would deem that it is equally likely for the inconsistency
% to lie in either of the reactions, but not only one of them.
r = findRxnIDs(model, 'MANPGH');
h = findMetIDs(model, 'h_c');
model.S(h, r) = 1;
[metMWrange,metForm,ele,metEle,rxnBal,infeasibility,inconUB,sol] = computeMetMWrange(model, [], 'man6pglyc_c', model.rxns(sum(model.S ~= 0, 1) > 1));
model.S(h, r) = 0;
% Look at the inconsistencies found in the results
fprintf('Minimum inconsistencies found when balancing each element:\n')
for j = 1:numel(sol.minIncon)
    fprintf('%s: %.4f\n', ele{j}, sol.minIncon(j).obj)
end
fprintf('man6pglyc_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrange(1), metForm{1}, strjoin(model.rxns(any(abs(rxnBal.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrange(2), metForm{2}, strjoin(model.rxns(any(abs(rxnBal.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
% Note that the the difference between metMWrange(1) and metMWrange(2) is
% 1.0XXX. 1 due to the difference in the hydrogen atom, .0XXX due to the
% tolerance allowed when optimizing for each element.
% And the unbalanced reaction in the min case is different from that in the max case.

% But if we take a metabolite, e.g. glc__D_c (glucose) or fru_c (fructose) with many connecting
% reactions, the reaction with the one single inconsistency will be pinned down
r = findRxnIDs(model, 'XYLI2');
h = findMetIDs(model, 'h_c');
model.S(h, r) = 1;
[metMWrangeGlc,metFormGlc,ele,metEle,rxnBalGlc, ~, ~, solGlc] = computeMetMWrange(model, [], 'glc__D_c', model.rxns(sum(model.S ~= 0, 1) > 1));
[metMWrangeFru,metFormFru,~,~,rxnBalFru, ~, ~, solFru] = computeMetMWrange(model, [], 'fru_c', model.rxns(sum(model.S ~= 0, 1) > 1));
model.S(h, r) = 0;
fprintf('Minimum inconsistencies found:\nelement\tglc__D_c\tfru_c\n')
for j = 1:numel(sol.minIncon)
    fprintf('%s:\t%.4f\t%.4f\n', ele{j}, solGlc.minIncon(j).obj, solFru.minIncon(j).obj)
end
fprintf('\nglc__D_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeGlc(1), metFormGlc{1}, strjoin(model.rxns(any(abs(rxnBalGlc.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeGlc(2), metFormGlc{2}, strjoin(model.rxns(any(abs(rxnBalGlc.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('\nfru_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeFru(1), metFormFru{1}, strjoin(model.rxns(any(abs(rxnBalFru.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeFru(2), metFormFru{2}, strjoin(model.rxns(any(abs(rxnBalFru.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
% Note that here the differences between max and min MW are smaller and << 1, 
% only due to the tolerance allowed. And the unbalanced reaction is unique.

% One can use the above method to check the MW of the biomass produced in
% silico regardless of the existing inconsistency