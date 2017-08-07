function [yn,element,metEle] = compareMetFormulas(metFormulas_1, metFormulas_2)
% Given two cell arrays of the same length of chemical formulas 'form_1' and
% 'form_2', determine if they are the same.
%
% USAGE:
%    [yn,ele,metEle] = compareMetFormulas(metFormulas_1, metFormulas_2)
%
% INPUTS:
%    metFormulas_1:   a cell array of chemical formulae 
%    metFormulas_2:   a second cell array of chemical formulae for comparison, same length as metFormulas_1
% 
% OUTPUTS:
%    yn:        logical vector, true if equal, false if unequal
%    ele:       cell array of elements in the formulas
%    metEle:    numel(form1) x numel(ele) x 2 matrix of metFormulas_1 and metFormulas_2
%               metEle(i,j,k) is the stoichiometry of element{j} in form_i

if iscell(metFormulas_1) && iscell(metFormulas_2) && numel(metFormulas_1) ~= numel(metFormulas_2)
    error('The number of chemical formulas in each cell is different!')
end
[~,ele1,metEle1] = checkEleBalance(metFormulas_1);
[~,ele2,metEle2] = checkEleBalance(metFormulas_2);
if iscell(metFormulas_1) && ischar(metFormulas_2)
    metEle2 = repmat(metEle2,numel(metFormulas_1),1);
elseif ischar(metFormulas_1) && iscell(metFormulas_2)
    metEle1 = repmat(metEle1,numel(metFormulas_2),1);
end
nForm = size(metEle1,1);
    


element = union(ele1, ele2);
[~,id] = ismember({'C';'H';'N';'O';'P';'S'}, element);
id2 = setdiff(1:numel(element), id);
element = element([id(id~=0); id2(:)]);
id = strcmp(element,'Charge');
if any(id)
    element = [element(~id); element(id)];
end

metEle = zeros(nForm,numel(element), 2);
[~,idEl1] = ismember(ele1,element);
[~,idEl2] = ismember(ele2,element);
metEle(:,idEl1,1) = metEle1;
metEle(:,idEl2,2) = metEle2;

yn = all(abs(metEle(:,:,1) - metEle(:,:,2)) < 1e-6, 2);
