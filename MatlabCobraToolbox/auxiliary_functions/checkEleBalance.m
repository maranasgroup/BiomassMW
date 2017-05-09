function [EleBal, element, metEle] = checkEleBalance(model, element, metEle,selfCall)
%Check the elemental balance of the reactions in the COBRA mode 'model'.
%
%[EleBal, element, metEle] = checkEleBalance(model)
%'metFormulas' must be a field in 'model'. Support parentheses/brackets.
%Each element must start with a capital letter followed by lowercase letters or '_', followed by a number. 
%For charges, just include them in the formulas as 'Charge-1' or 'Charge2' etc.
%  'EleBal' is the elemental balance of each element in 'element' for each
%    reaction, (nE x nR matrix)
%  'element' is nE x 1 cell array of elements detected in model.metFormulas
%  'metEle' is the elemental composition of metabolites (nM x nE matrix)
%     metEle(m,e) is the stoichiometry of element e in metabolite m
%
%[EleBal, element, metEle] = checkEleBalance(metFormulas)
%The first argument can also be 'metFormulas', nM x 1 cell array of strings and 'S' nM x nR matrix 
% In this case, only 'element' and 'metEle' is returned 'EleBal' would be empty.
%
%[EleBal, element, metEle] = checkEleBalance(model/metFormulas, element, metEle)
% If you already have 'element' and 'metEle' from previous runs, this
% includes the newly calculated formulas for metabolites in 'model' into 
% 'element' and 'metEle'
%
%e.g., [~, element,metEle] = checkEleBalance({'H2O'; '[H2O]2(CuSO4)'}) would return:
%  element = {'H';'O';'Cu';'S'}
%  metEle = [2, 1, 0, 0; 4, 6, 1, 1]
%
%Siu Hung Joshua Chan May 2017

if nargin < 4
    selfCall = false;
end
%for recalling the original formula at the top level if there are parentheses
%in the formula leading to iterative calling
persistent formTopLv
if ~selfCall
    formTopLv = '';
end
nE_max = 150;
if ~isstruct(model)
    form = model;
    if ~iscell(form)
        form = {form};
    end
else
    form = model.metFormulas;
end
if ~exist('element', 'var')
    element = {};
elseif isempty(element)
    element = {};
elseif exist('metEle', 'var') && ~isempty(metEle) && numel(element) ~= size(metEle, 2)
    error('number of elements in metEle and element are not the same.')
elseif numel(unique(element)) < numel(element)
    error('Repeated elements in the input ''element'' array.')
else
    element = element(:);
end

calc = true(numel(form),1);
if ~exist('metEle', 'var') || isempty(metEle)
    metEle = zeros(numel(form), max([nE_max,numel(element)]));
elseif size(metEle,1) ~= numel(form)
    error('Incorrect size of the input element vector for metabolites')
else
    %calculate only those entries in the input metEle with NaN
    calc(~any(isnan(metEle), 2)) = false;
end
form = strrep(form,'[','(');
form = strrep(form,']',')');
form = strrep(form,'{','(');
form = strrep(form,'}',')');
nE = numel(element);
for j = 1:numel(form)
    form{j} = strtrim(form{j});
    if ~isempty(form{j}) && calc(j)
        %get all outer parentheses
        parenthesis = [];
        stP = [];
        stPpos = [];
        lv = 0;
        k = 1;
        while k <= length(form{j})
            if strcmp(form{j}(k),'(')
                if lv == 0
                    pStart = k;
                end
                lv = lv + 1;
            elseif strcmp(form{j}(k),')')
                if lv == 1
                    parenthesis = [parenthesis; [pStart k]];
                    %closed parenthesis, get the following stoichiometry if any
                    stPre = regexp(form{j}(k+1:end), '^(\+|\-)?\d*\.?\d*', 'match');
                    if isempty(stPre)
                        stP = [stP; 1];
                        stPpos = [stPpos; k k];
                    elseif isnan(str2double(stPre{1}))
                        f = [form{j}(parenthesis(end,1):parenthesis(end,2)), stPre{1}];
                        if isempty(formTopLv)
                            s2 = sprintf('the input formula ''%s''',form{j});
                        else
                            s2 = sprintf('''(%s)'' in the input formula ''%s''',form{j},formTopLv);
                        end
                        error(['#%d: Invalid stoichiometry in chemical formula. ''%s''',...
                            ' from %s'],j,f,s2);
                    else
                        %stoichiometry
                        stP = [stP; str2double(stPre{1})];
                        %position of the stoichiometry in the text
                        stPpos = [stPpos; k+1, k+length(stPre{1})];
                        k = k + length(stPre{1});
                    end
                end
                lv = lv - 1;
            end
            k = k + 1;
        end
        if isempty(parenthesis)
        %if isempty(strfind(form{j}, '(')) %&& ~strcmp(form{j}, 'none')
        %No parenthesis, parse the formula
            re = regexp(form{j}, '([A-Z][a-z_]*)((?:\+|\-)?\d*\.?\d*)', 'tokens');
            s = strjoin(cellfun(@(x) strjoin(x,''),re,'UniformOutput',false),'');
            errorFlag = 0;
            if ~strcmp(s,form{j})
                errorFlag = 1;
                
            else
                f = cellfun(@(x) ~isempty(x{2}) & isnan(str2double(x{2})),re);
                if any(f)
                    f = strjoin(cellfun(@(x) strjoin(x,''),re(f),'UniformOutput',false),''', ''');
                    errorFlag = 2;
                end
            end
            if errorFlag > 0
                if isempty(formTopLv)
                    s2 = sprintf('the input formula ''%s''',form{j});
                else
                    s2 = sprintf('''(%s)'' in the input formula ''%s''',form{j},formTopLv);
                end
                if errorFlag == 1
                    error(['#%d: Invalid chemical formula. Only ''%s'' can be recognized' ...
                        ' from %s.\n'...
                        'Each element should start with a capital letter followed by lower case letters '...
                        'or ''_'' with indefinite length and followed by a number.\n'],...
                        j, s, s2)
                elseif errorFlag == 2
                    error(['#%d: Invalid stoichiometry in chemical formula. ''%s''',...
                        ' from %s'],j,f,s2);
                end
            end
            elementJ = repmat({''},numel(re), 1);
            nEj = 0;
            stoichJ = zeros(numel(re),1);
            for k = 1:numel(re)
                [ynK,idK] = ismember(re{k}(1),elementJ(1:nEj));
                if ynK
                    k2 = idK;
                else
                    nEj = nEj + 1;
                    elementJ{nEj} = re{k}{1};
                    k2 = nEj;
                end
                if isempty(re{k}{2})
                    stoichJ(k2) = stoichJ(k2) + 1;
                else
                    stoichJ(k2) = stoichJ(k2) + str2double(re{k}{2});
                end
            end
            elementJ = elementJ(1:nEj);
            stoichJ = stoichJ(1:nEj);
            [ynE, ~] = ismember(elementJ, element);
            nE = nE + sum(~ynE);
            element = [element; elementJ(~ynE)];
            [~, idE] = ismember(elementJ, element);
            metEle(j, :) = 0;
            metEle(j, idE) = stoichJ;
        else
            %parentheses found. iteratively get the formula inside
            %parentheses
            rest = true(length(form{j}),1);
            if isempty(formTopLv)
                formTopLv = form{j};
            end
            for k = 1:size(parenthesis,1)
                [~, element,metEleJ] = checkEleBalance(form{j}(...
                    (parenthesis(k,1)+1):(parenthesis(k,2)-1)),element,[],true);
                metEle(j,1:numel(element)) = metEle(j,1:numel(element)) + metEleJ * stP(k);
                rest(parenthesis(k,1):stPpos(k,2)) = false;
            end
            if any(rest)
                [~, element,metEleJ] = checkEleBalance(form{j}(rest),element,[],true);
                metEle(j,1:numel(element)) = metEle(j,1:numel(element)) + metEleJ;
            end
            formTopLv = '';
        end
    elseif calc(j)
        metEle(j,:) = NaN;
    end
end
metEle = metEle(:, 1:numel(element));
if isstruct(model)
    EleBal = metEle' * model.S;
else
    EleBal = [];
end
end