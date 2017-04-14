function [model,metCompute,S_fill,rxnBal,ele,metEle,N,LP] = computeFormulasFillMets(model,metKnown,rxns,metFill,findCM,nameCM,param)
% Compute the chemical formulas of the unknown metabolites 
% using a set of metabolites with known formulae and a set of reactions.
% To include charge balance in the computation, simply in all formulas, add
% e.g. 'Charge2' for charge +2 or 'Charge-1' for charge -1. 'C' must be
% capitalized and 'harge' must be in lower case.
% The minimum conflict is found by allowing filling up by metFill (e.g. proton) in the reaction stoichiometry
%
% [model,metCompute,S_fill,rxnBal,ele,metEle,N] = computeFormulasFillMets(model,metKnown,rxns,metFill,findCM);
% Input:
%   model:              COBRA model
%   metKnown:           known metabolites (character array or IDs)
%   rxns:               the set of reactions for inferring formulae (character array or IDs)
%   metFill:            the chemical formulas for compounds for freely filling the
%                       imbalance, e.g. {'HCharge1', 'H2O'}
%   findCM = 'efmtool': find conserved moieties using EFMtool (most comprehensive, 
%                       but computational cost may be high.)  
%                       (default, if EFMtool not in path, switch to rational basis)
%   findCM = true  :    find conserved moieties and add into formulas using the rational basis of model.S'
%   findCM = N     :    directly supply the rational basis for N(S') or 
%                       the set of extreme rays for {n | n * S = 0, n >= 0}
%   findCM = false :    not to find conserved moieties and return minimal formulas
% Output:
%   model:          COBRA model with updated formulas
%   metCompute:     M_unknown x 2 array of cell with [mets | computed formulas]
%   S_fill:         mF x n matrix representing the adjustment of the S-matrix 
%                   if you supply mF mets as 'metFill'
%   rxnBal:         N by E matrix for the elemental balance of rxns, the ij-th
%                   entry is the imbalance of the i-th rxn regarding the j-th element in 'ele'
%   ele             elements corresponding to the row of rxnBal, as well as
%                   to the coloumn of metEle
%   metEle:         chemical formulas in a M by E matrix.
%   N:              the extreme ray or rational null space matrix

%[model,metCompute,S_fill,rxnBal,ele,metEle,N] = computeFormulasFillMets(model,metKnown,rxns,metFill,findCM,nameCM);
% nameCM = 0 the program assigns default names for conserved moieties
%            (Conserve_a, Conserve_b, ...)
% nameCM = 1 to name true conserved moieties interactively (exclude dead end mets). 
% nameCM = 2 to name all (including dead end)
%
% [model,metCompute,rxnBal,ele,metEle,N,LP] = computeFormulas(...)
% Return also the CPLEX optimization object 'LP'

%check mets with formulas that need to be transformed
%form0 = model.metFormulas;
if ~isfield(model,'metFormulas')
    error('model does not have the field ''metFormulas.''')
end
if nargin < 7 || isempty(param) 
    param = struct();
end
if nargin < 6 || isempty(nameCM)
    nameCM = 0;
end
if nargin < 5 || isempty(findCM)
    findCM = 'efmtool';
end
if nargin < 4 || isempty(metFill)
    %defaulted proton for filling imbalance
    metFill = {'HCharge1'};
elseif ischar(metFill)
    metFill = {metFill};
end
if nargin < 3 || isempty(rxns)
    rxns = find(sum(model.S~=0,1) > 1 & (model.lb ~= 0 | model.ub ~= 0)');
end
if ischar(rxns)
    rxns = {rxns};
end
if iscell(rxns)
    rxnC = findRxnIDs(model,rxns);
else
    rxnC = rxns;
end
if any(rxnC == 0)
    if iscell(rxns)
        error('%s in rxns is not in the model.', rxns{find(rxnC==0,1)});
    else
        error('rxns index must be positive integer.')
    end
end
if nargin < 2 || isempty(metKnown)
    metKnown = model.mets(~cellfun(@isempty,model.metFormulas));
end
if ischar(metKnown)
    metKnown = {metKnown};
end
if iscell(metKnown)
    metK = findMetIDs(model,metKnown);
else
    metK = metKnown;
end
if any(metK == 0) 
    if iscell(metKnown)
        error('%s in metKnown is not in the model.', metKnown{find(metK==0,1)});
    else
        error('metKnown index must be positive integer.')
    end
end
metKform = cellfun(@isempty,model.metFormulas(metK));
if any(metKform)
    warning('Some mets in metKnown do not have formulas in the model. Ignore them.');
end
%All formulas must be in the form of e.g. Abc2Bcd1. Elements are
%represented by one capital letter followed by lower case letter or
%underscore, followed by a number for the stoichiometry. No brackets or
%other symbols allowed.
[metK,metKform] = deal(metK(~metKform), model.metFormulas(metK(~metKform)));
re = regexp(metKform,'[A-Z][a-z_]*(\-?\d+\.?\d*)?','match');
re = cellfun(@(x) strjoin(x,''),re,'UniformOutput',false);
goodForm = strcmp(strtrim(re), strtrim(metKform));
if ~all(goodForm)
    goodForm = find(~goodForm,1);
    error('%s has an invalid formula %s\n',metKnown{goodForm},metKform{goodForm});
end
% [ynF,idF] = ismember(metFill,metK);
% if ~all(ynF)
%     error('''%s'' in ''metFill'' is not in ''metKnown''.',strjoin(model.mets(metFill(~ynF)),''', '''));
% end
%% minimum inconsistency
%formulas for known metabolites
[~,eleK,metEleK] = checkEleBalance(metKform);
%formulas for filling metabolites
[~,eleK,metEleF] = checkEleBalance(metFill, eleK);
if numel(eleK) > size(metEleK,2)
    metEleK = [metEleK, zeros(size(metEleK,1), numel(eleK) - size(metEleK,2))];
end
eleCh = strcmp(eleK,'Charge');
m = size(model.S,1);
nE = numel(eleK);
mK = numel(metK);
mU = m - mK;
mF = numel(metFill);
metU = setdiff((1:m)',metK);
nR = numel(rxnC);
nameM = strcat(repmat(model.mets(metU), nE,1),'_', reshape(repmat(eleK(:)',mU,1),nE*mU,1));
nameR = strcat(repmat(model.rxns(rxnC), nE,1),'_', reshape(repmat(eleK(:)',nR,1),nE*nR,1));
if mF > 0
    nameF = strcat(repmat(metFill(:), nR,1),'_', reshape(repmat(model.rxns(rxnC)',mF,1),mF*nR,1));
end
LP = Cplex();
%n_ik, stoichoimetry for element k in met i
LP.addCols(zeros(mU*nE,1),[],zeros(mU*nE,1),inf(mU*nE,1),[],char(nameM));
if any(eleCh)
    %if charge is in the formula, allow negative charges
    eleChK = find(eleCh);
    LP.Model.lb((mU*(eleChK-1) + 1):(mU*eleChK)) = -inf;
end
%x^+_jk, positive inconsistency of element k in reaction j
LP.addCols(1000*ones(nR*nE,1),[],zeros(nR*nE,1),inf(nR*nE,1),[],char(strcat(nameR,'_pos')));
%x^-_jk, negative inconsistency of element k in reaction j
LP.addCols(1000*ones(nR*nE,1),[],zeros(nR*nE,1),inf(nR*nE,1),[],char(strcat(nameR,'_neg')));
if mF > 0
    %z^+_ij, mets for filling inconsistency
    LP.addCols(0.1*ones(nR*mF,1),[],zeros(nR*mF,1),inf(nR*mF,1),[],char(strcat(nameF,'_pos')));
    %z^-_ij, mets for filling inconsistency
    LP.addCols(0.1*ones(nR*mF,1),[],zeros(nR*mF,1),inf(nR*mF,1),[],char(strcat(nameF,'_neg')));
end
%RHS
b = -model.S(metK,rxnC)' * metEleK;
b = b(:);
%constraint matrix
[row,col,entry] = find(model.S(metU,rxnC)');
row = repmat(row(:),nE,1)+reshape(repmat(0:nR:nR*(nE-1),numel(row),1),numel(row)*nE,1);
col = repmat(col(:),nE,1)+reshape(repmat(0:mU:mU*(nE-1),numel(col),1),numel(col)*nE,1);
entry = repmat(entry,nE,1);
row2 = reshape(repmat(1:nR*nE, mF, 1),nR*nE*mF,1);
col2 = repmat((1:mF*nR)',nE,1);
entry2 = reshape(repmat(metEleF,nR,1),nR*nE*mF,1);
A = [sparse(row,col,entry,nR*nE,mU*nE), ... % sum(i, S_ij * n_ie)
    sparse(1:nR*nE,1:nR*nE,ones(nR*nE,1),nR*nE,nR*nE), ... % + x^+_je
    sparse(1:nR*nE,1:nR*nE,-ones(nR*nE,1),nR*nE,nR*nE), ... % - x^-_je
    sparse(row2,col2,entry2), -sparse(row2,col2,entry2)]; % + z^+_ij * n_ie - z^-_ij * n_ie
% S_unknown' n_unknown = - S_known' * n_known
LP.addRows(b, A, b, char(strcat(nameR,'_bal')));
%Handle parameters
cobraParam = struct();
cobraParamList = {'feasTol','optTol'};
CplexParamList = {'feasibility','optimality'};
[cobraParam.feasTol, cobraParam.optTol] = getCobraSolverParams('LP',cobraParamList, param);
for j = 1:numel(CplexParamList)
    if isfield(param,CplexParamList{j})
        cobraParam.(cobraParamList{j}) = param.(CplexParamList{j});
    elseif isfield(param,'simplex') && isfield(param.simplex,CplexParamList{j})
        cobraParam.(cobraParamList{j}) = param.simplex.(CplexParamList{j});
    elseif isfield(param,'simplex') && isfield(param.simplex,'tolerances') && isfield(param.simplex.tolerances,CplexParamList{j})
        cobraParam.(cobraParamList{j}) = param.simplex.tolerances.(CplexParamList{j});
    else
        param.simplex.tolerances.(CplexParamList{j}) = cobraParam.(cobraParamList{j});
    end
end
LP = setCplexParam(LP,param);
feasTol = LP.Param.simplex.tolerances.feasibility.Cur;
%Solve
LP.solve;
%manually check feasibility
feas = checkSolFeas(LP);

if ~(feas <= feasTol)
    if isfield(LP.Solution,'x')
        metEle = zeros(numel(model.mets),numel(eleK));
        metEle(metK,:) = metEleK;
        metEle(metU,:) = reshape(LP.Solution.x(1:mU*nE),mU,nE);
        rxnBal = metEle' * model.S;
        S_fill = sparse(repmat((1:mF)',nR,1),reshape(repmat(rxnC(:)',mF,1),mF*nR,1), ...
            LP.Solution.x(((mU+nR*2)*nE+1):((mU+nR*2)*nE+nR*mF)) ...
            - LP.Solution.x(((mU+nR*2)*nE+nR*mF+1):((mU+nR*2)*nE+nR*mF*2)), ...
            mF, size(model.S,2));
        % S_fill(metFill,rxnC) = reshape(LP.Solution.x(((mU+nR*2)*nE+1):((mU+nR*2)*nE+nR*mF*2)),mF,nR);
    else
        metEle = [];
        rxnBal = [];
        S_fill = [];
    end
    ele = eleK;
    %terminate if infeasible (should not happen) 
    fprint('Infeasible during optimization for minimal inconsistency.');
    return
end
x = LP.Solution.x;

%% minimal formulas
%bound the inconsistency with the above solution
ind = (mU*nE+1):((mU+2*nR)*nE);
LP.Model.ub(ind(x(ind) <= 0)) = 0;
LP.Model.ub(ind(x(ind) > 0)) = x(ind(x(ind) > 0)) * (1 + 1e-12);
LP.Model.lb(ind(x(ind) > feasTol)) = x(ind(x(ind) > feasTol)) * (1 - 1e-12);
ind = ((mU+2*nR)*nE+1):((mU+2*nR)*nE+mF*nR*2);
LP.Model.ub(ind(x(ind) <= 0)) = 0;
LP.Model.ub(ind(x(ind) > 0)) = x(ind(x(ind) > 0)) * (1 + 1e-12);
LP.Model.lb(ind(x(ind) > feasTol)) = x(ind(x(ind) > feasTol)) * (1 - 1e-12);
%Change objective
LP.Model.obj(:) = 0;
LP.Model.obj(1:mU*nE) = 1;
if any(eleCh)
    %decompose charge variables into +ve and -ve part if exist
    orderCh = find(eleCh);
    LP.Model.obj((mU*(orderCh - 1) + 1):(mU*orderCh)) = 0;
    LP.addCols(ones(mU,1),[],zeros(mU,1),inf(mU,1),[],char(strcat(model.mets(metU),'_charge+')));
    LP.addCols(ones(mU,1),[],zeros(mU,1),inf(mU,1),[],char(strcat(model.mets(metU),'_charge-')));
    LP.addRows(zeros(mU,1),...
        [sparse(mU,mU*(orderCh - 1)), sparse(1:mU,1:mU,ones(mU,1),mU,mU), ...
        sparse(mU, mU*(nE-orderCh) + nR*2*nE), sparse(mU, mF*nR*2), ...
        sparse(1:mU,1:mU,-ones(mU,1),mU,mU),sparse(1:mU,1:mU,ones(mU,1),mU,mU)],...
        zeros(mU,1),char(strcat(model.mets(metU),'_charge_decomp')));
end
%solve
LP.solve;
%check feasibility
feas = checkSolFeas(LP);
metEle = zeros(numel(model.mets),numel(eleK));
metEle(metK,:) = metEleK;
if feas <= feasTol
    metEle(metU,:) = reshape(LP.Solution.x(1:mU*nE),mU,nE);
    S_fill = sparse(repmat((1:mF)',nR,1),reshape(repmat(rxnC(:)',mF,1),mF*nR,1), ...
            LP.Solution.x(((mU+nR*2)*nE+1):((mU+nR*2)*nE+nR*mF)) ...
            - LP.Solution.x(((mU+nR*2)*nE+nR*mF+1):((mU+nR*2)*nE+nR*mF*2)), ...
            mF, size(model.S,2));
else
    %use the previous solution if currectly infeasible
    metEle(metU,:) = reshape(x(1:mU*nE),mU,nE);
    fprint('Infeasible during optimization for minimal formulas.');
end


%% find conserved moieties
if nargin < 4 || isempty(findCM)
    findCM = true;
end
cont = true;
if ischar(findCM) && strcmpi(findCM,'efmtool')
    %use EFMtool
    pathEFM = which('CalculateFluxModes.m');
    if isempty(pathEFM)
        warning('EFMtool not in Matlab path. Use rational basis.');
    else
        dirEFM = strsplit(pathEFM,filesep);
        dirEFM = strjoin(dirEFM(1:end-1),filesep);
        dirCur = pwd;
        cd(dirEFM);
        %will very probably fail due to lack of memory if there are many
        %dead end metabolites, may add code to remove deadend mets first
        %         [~,removedMets] = removeDeadEnds(model);
        %         metDead = findMetIDs(model,removedMets);
        N = CalculateFluxModes(full(model.S'),zeros(numel(model.mets),1));
        N = N.efms;
        cd(dirCur);
        cont = false;
    end
    findCM = true;
end
if cont
    if size(findCM,1) == numel(model.mets)
        %input is the null space matrix / set of extreme rays
        N = findCM;
        findCM = true;
    elseif numel(findCM) == 1 && findCM
        N = null(full(model.S'),'r');
    else
        N = [];
    end
end
if findCM
    fprintf('Find conserved moieties...\n');
    
    %clear close-to-zero values
    N(abs(N) < 1e-8) = 0;
    N = sparse(N);
    %true generic conserved moieties, positive and not involving known mets
    Ncm = N(:,~any(N < 0, 1) & ~any(N(metK,:),1));
    %add them into formulas
    metEle = [metEle, Ncm];
    ele = [eleK(:); cell(size(Ncm,2),1)];
    j2 = 1;
    for j = 1:size(Ncm,2)
        while any(strcmp(ele(1:nE),['Conserve_' num2alpha(j2)]))
            j2 = j2 + 1;
        end
        ele{nE+j} = ['Conserve_' num2alpha(j2)];
        j2 = j2 + 1;
    end
else
    ele = eleK(:);
end
rxnBal = metEle' * model.S;
model.metFormulas = convertMatrixFormulas(ele,metEle,10);
if nameCM > 0 && findCM ~= 0
    %manually name conserved moieties
    ele0 = ele;
    nDefault = 0;
    nCM = size(Ncm,2);
    eleDel = false(nE + nCM, 1);
    if nameCM == 1
        %get dead end metatbolites
        [~,removedMets] = removeDeadEnds(model);
        metDead = findMetIDs(model,removedMets);
    end
    for j = 1:nCM
        fprintf('\n');
        writeCell2Text([model.mets(Ncm(:,j)~=0),model.metFormulas(Ncm(:,j)~=0),...
            model.metNames(Ncm(:,j)~=0)]);
        fprintf('\n');
        if nameCM == 1 && any(Ncm(metDead,j),1)
            %use the defaulted for dead end mets
            nDefault = nDefault + 1;
            ele{nE+j} = ele0{nE + nDefault};
        else
            cont = false;
            while true
                s = input(['Enter the formula for the conserved moiety (e.g. OHRab_cd):\n',...
                    '(hit return to use default name ''Conserve_xxx'')\n'],'s');
                if isempty(s)
                    %use the defaulted
                    nDefault = nDefault + 1;
                    ele{nE+j} = ele0{nE + nDefault};
                    break
                end
                re = regexp(s,'[A-Z][a-z_]*(\-?\d+\.?\d*)?','match');
                if strcmp(strjoin(re,''),s)
                    %manual input formula, continue to checking
                    cont = true;
                    break
                end
            end
            if cont
                %get the matrix for the input formula
                nEnew = numel(ele) - nE - nCM;
                [~, eleJ, metEleJ] = checkEleBalance(s,ele([1:nE, (nE+nCM+1):end]));
                metEle(:,[1:nE, (nE+nCM+1):end]) = metEle(:,[1:nE, (nE+nCM+1):end])...
                        + metEle(:,nE+j) * metEleJ(1,1:(nE+nEnew));
                if numel(eleJ) > nE + nEnew
                    %there are new elements
                    ele = [ele(:); eleJ((numel(ele)-nCM+1):end)];
                    metEle = [metEle, metEle(:,nE+j) * metEleJ(1,(nE+nEnew+1):end)];
                end
%                 [ynJ,idJ] = ismember(eleJ,ele(1:nE));
%                 if any(ynJ)
%                     metEle(:,idJ(ynJ)) = metEle(:,idJ(ynJ)) + metEle(:,nE+j) * metEleJ(1,ynJ);
%                 end
%                 eleJ = eleJ(~ynJ);
%                 metEleJ = metEleJ(1,~ynJ);
%                 if numel(ele) > nE + nCM
%                     [ynJ,idJ] = ismember(eleJ,ele(nE+nCM+1:end));
%                     if any(ynJ)
%                         metEle(:,nE+nCM+idJ(ynJ)) = metEle(:,nE+nCM+idJ(ynJ)) ...
%                             + metEle(:,nE+j) * metEleJ(1,ynJ);
%                     end
%                     eleJ = eleJ(~ynJ);
%                     metEleJ = metEleJ(1,~ynJ);
%                 end
%                 if ~isempty(eleJ)
%                     ele = [ele(:);eleJ(:)];
%                     metEle = [metEle, metEleJ];
%                 end
                eleDel(nE + j) = true;
            end
        end
    end
    %del defaulted but replaced columns
    if any(eleDel)
        eleDel = find(eleDel);
        ele(eleDel) = [];
        metEle(:,eleDel) = [];
    end
    %1:nE              : real elements
    %nE+1:nE+nDefault  : default generic elements (Conserve_xxx)
    %nE+nDeafult+1:end : generic element by user's input
    %Change if names of default generic elements are mixed up with user
    %input
    j0 = 0;
    for j = 1:nDefault
        j0 = j0 + 1;
        nameJ = ['Conserve_' num2alpha(j0)];
        while any(strcmp(ele([1:nE, nE+nDefault+1:end]),nameJ))
            j0 = j0 + 1;
            nameJ = ['Conserve_' num2alpha(j0)];
        end
        ele{nE+j} = nameJ;
    end
end
model.metFormulas = convertMatrixFormulas(ele,metEle,10);
metCompute = [model.mets(metU) model.metFormulas(metU)];
end

function s = num2alpha(index,charSet)
%s = num2alpha(j,charSet)
%Given a nonzero integer j and a character set charSet, convert j into
%a string formed from the characters in charSet having order j.
%'charSet' defaulted to be '_abcdefghijklmnopqrstuvwxyz' in which '_' acts
%like 0 and 'z' acts like the largest digit 9 in decimal expression
%e.g. num2slpha(0) is '_' , num2slpha(1) is 'a', num2slpha(27^2) is 'a__'
%
if nargin < 2
    charSet = ['_' char(97:122)];
end
if numel(index) > 1
    s = cell(numel(index),1);
    for j = 1:numel(index)
        s{j} = num2alpha(index(j),charSet);
    end
    return
end
N = length(charSet);
s = '';
k = floor(index/N);
while k > 0
    s = [s charSet(index - k*N + 1)];
    index = k;
    k = floor(index/N);
end
s = [charSet(mod(index,N) + 1) s];
end