function [metMWrange,metForm,metFeas,rxnBal,ele,metEle,LP,sol] = computeMetMWrange(model,metKnown,metInterest,rxns,percent,varargin)
% Compute the minimum and maximum molecular weight (MW) of a metabolie
% given a set of formulas for the known metabolites using a set of reactions. 
% Done by first minimizing the mass-and-charge imbalance. Then fixing the 
% imbalance at the minimum or higher level, minimize and maximize the molecular
% weight of the met of interest.
%
%[metMWrange,metForm,metFeas,rxnBal,ele,metEle,LP] = computeMetMWrange(model,metKnown,metInterest,rxns,percent,param)
% Input: 
%    model         COBRA model
%    metKnown      A known set of metabolites with formulas (character array or IDs)
%    metInterest   The metabolite of interest (character or ID)
%    rxns          The set of reactions used to compute the formula of metInterest 
%                  (character array or IDs)
%    (optional)
%    percent       The percent of inconsistency allowed when calculating the range
%                  for the molecular weight. The constraints added are
%                  sum(inconsistency_e) <= min(sum(inconsistency_e)) * (1 + percent)
%                  for each element e. Default percent = 0
%    parameters    Parameters for Cobra Toolbox (see help solveCobraLP), 
%                  parameter name followed by value, or parameter structre
%
% Output:
%    metMWrange    [min MW, max MW] (1 x 2 vector)
%    metForm       The corresponding empirical formulas (1 x 2 cell array)
%    metFeas       Infeasibility from the corresponding optimization (1 x 2 vector)
%    rxnBal        imbalance of each reaction (#element x #rxns matrix)
%    ele           elements corresponding to the row of rxnBal, as well as
%                  to the coloumn of metEle
%    metEle        chemical formulas in a M by E matrix
%    LP            Cplex LP structure

%check mets with formulas that need to be transformed
%form0 = model.metFormulas;

if nargin < 5 || isempty(percent)
    percent = 0;
end
if nargin < 4 || isempty(rxns)
    rxnC = find(sum(model.S~=0,1)>1 & (model.lb ~=0 | model.ub ~= 0)')';
elseif iscell(rxns) || ischar(rxns)
    rxnC = findRxnIDs(model,rxns);
else
    rxnC = rxns;
end
if any(rxnC == 0)
    error('%s in rxns is not in the model.', rxns{find(rxnC==0,1)});
end
if nargin < 2 || isempty(metKnown)
    metKnown = model.mets(~cellfun(@isempty,model.metFormulas));
end
if iscell(metKnown) || ischar(metKnown)
    metK = findMetIDs(model,metKnown);
else
    metK = metKnown;
end
if any(metK == 0)
    error('%s in metKnown is not in the model.', metKnown{find(metK==0,1)});
end
if iscell(metInterest) || ischar(metInterest)
    metInterest = findMetIDs(model,metInterest);
end
metInterest0 = metInterest;
if metInterest == 0
    error('The biomass met ID is incorrect.');
elseif ismember(metInterest,metK)
    metK(metK == metInterest) = [];
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
goodForm = strcmp(re, strtrim(metKform));
if ~all(goodForm)
    goodForm = find(~goodForm,1);
    error('%s has an invalid formula %s\n',metKnown{goodForm},metKform{goodForm});
end
feasTol = getCobraSolverParams('LP',{'feasTol'});
%% minimum inconsistency
[~,ele,metEleK] = checkEleBalance(metKform);
eleCh = strcmp(ele,'Charge');
m = size(model.S,1);
nE = numel(ele);
mK = numel(metK);
mU = m - mK;
metU = setdiff((1:m)',metK);
metInterest = find(metU == metInterest);
nR = numel(rxnC);
nameM = strcat(repmat(model.mets(metU), nE,1),'_', reshape(repmat(ele(:)',mU,1),nE*mU,1));
nameR = strcat(repmat(model.rxns(rxnC), nE,1),'_', reshape(repmat(ele(:)',nR,1),nE*nR,1));

ind = struct(); %index structure for variables and constraints
nVar = 0;%number of variables
nCon = 0;%number constraints
%n_ik, stoichoimetry for element k in met i
for j = 1:nE
    ind.var.(['m' ele{j}]) = [nVar + 1, nVar + mU];
    nVar = nVar + mU;
end
%x^+_jk, positive inconsistency of element k in reaction j
for j = 1:nE
    ind.var.(['xPos' ele{j}]) = [nVar + 1, nVar + nR];
    nVar = nVar + nR;
end
%x^-_jk, negative inconsistency of element k in reaction j
for j = 1:nE
    ind.var.(['xNeg' ele{j}]) = [nVar + 1, nVar + nR];
    nVar = nVar + nR;
end
%index for reaction balance constraints
for j = 1:nE
    ind.con.(['rBal' ele{j}]) = [nCon + 1, nCon + nR];
    nCon = nCon + nR;
end
%constraint matrix: S_unknown' n_unknown
[row,col,entry] = find(model.S(metU,rxnC)');
row = repmat(row(:),nE,1)+reshape(repmat(0:nR:nR*(nE-1),numel(row),1),numel(row)*nE,1);
col = repmat(col(:),nE,1)+reshape(repmat(0:mU:mU*(nE-1),numel(col),1),numel(col)*nE,1);
entry = repmat(entry,nE,1);
LP = struct();
LP.A = [sparse(row,col,entry,nR*nE,mU*nE), sparse(1:nR*nE,1:nR*nE,ones(nR*nE,1),nR*nE,nR*nE), ...
    sparse(1:nR*nE,1:nR*nE,-ones(nR*nE,1),nR*nE,nR*nE)];
%RHS: - S_known' * n_known
LP.b = -model.S(metK,rxnC)' * metEleK;
LP.b = LP.b(:);
[LP.c, LP.lb] = deal(zeros(nVar,1));
LP.c(ind.var.(['xPos' ele{1}])(1) : ind.var.(['xNeg' ele{end}])(2)) = 1;
if any(eleCh)
    %if charge is in the formula, allow negative charges
    LP.lb(ind.var.mCharge(1):ind.var.mCharge(2)) = -inf;
end
LP.ub = inf(nVar,1);
LP.osense = 1; %minimize
LP.csense = char('E' * ones(1, nR*nE)); %all equality constraints
if nargin < 6
    sol = solveCobraLP(LP);
else
    sol = solveCobraLP(LP, varargin);
end
%manually check feasibility
feas = checkSolFeas(LP, sol);
if ~(feas <= feasTol)
    %terminate if infeasible (should not happen)
    if isfield(sol,'full') && numel(sol.full) == nVar
        metEle = zeros(numel(model.mets),numel(ele));
        metEle(metK,:) = metEleK;
        metEle(metU,:) = reshape(sol.full(1:mU*nE),mU,nE);
        rxnBal = metEle' * model.S;
    else
        metEle = [];
        rxnBal = [];
    end
    fprint('Infeasible during optimization for minimal inconsistency.');
    return
end
x = sol.full;
%reuse the basis if it exists.
if isfield(sol, 'basis')
    LP.basis = sol.basis;
end
%% metabolite MW range under minimum inconsistency
%give bounds for the met of interest
LP.ub(metInterest:mU:(metInterest+mU*(nE-1))) = 1e7;
% LP.Model.ub(metInterest:mU:(metInterest+mU*(nE-1))) = 1e7;
if any(eleCh)
    LP.lb(metInterest+mU*(find(eleCh)-1)) = -1e7;
    LP.Model.lb(metInterest+mU*(find(eleCh)-1)) = -1e7;
end
% Constraint element-wise, sum(inconsist_e) <= min_inconsist_e * (1+percent)
row = [];
col = [];
entry = [];
b2 = zeros(nE,1);
for k = 1:nE
    ind.con.(['xMax' ele{k}]) = nCon + k;
    indK = (mU*nE+nR*(k-1)+1):(mU*nE+nR*k);
    row = [row; k*ones(nR*2, 1)];
    col = [col; (ind.var.(['xPos' ele{k}])(1):ind.var.(['xPos' ele{k}])(2))';...
        (ind.var.(['xNeg' ele{k}])(1):ind.var.(['xNeg' ele{k}])(2))'];
    entry = [entry; ones(nR*2, 1)];
    %total inconsistency for the element (small tolerance allowed)
    b2(k) = sum(x([(ind.var.(['xPos' ele{k}])(1):ind.var.(['xPos' ele{k}])(2))';...
        (ind.var.(['xNeg' ele{k}])(1):ind.var.(['xNeg' ele{k}])(2))'])) * (1+1e-7+abs(percent))+1e-5;
    %     LP.addRows(0,...
    %         sparse(ones(nR*2,1),[indK, (indK + nR*nE)], ones(nR*2,1),1,size(LP.Model.A,2)),...
    %         sum(x([indK, (indK + nR*nE)]))*(1+1e-7+abs(percent))+1e-5, strcat('min_incon_',ele{k}));
end
nCon = nCon + nE;
LP.A = [LP.A; sparse(row, col, entry, nE, nVar)];
LP.b = [LP.b; b2];
LP.csense = [LP.csense char('L' * ones(1, nE))];
%special care for gurobi cbasis as the number of rows has changed
if isfield(LP, 'basis') && isfield(LP.basis, 'cbasis')
    LP.basis.cbasis(end + 1 : end + nE) = 0;
end
%molecular weight of each element
c = MW(ele);
%exclude generic groups, which have value NaN in c
c2 = c;
c2(isnan(c)) = 0;
%change objective to the molecular weight of the met of interest
LP.c(:) = 0;
LP.c(metInterest:mU:(metInterest+mU*(nE-1))) = c2;
LP.osense = 1;
if nargin < 6
    sol(2) = solveCobraLP(LP);
else
    sol(2) = solveCobraLP(LP, varargin);
end
feas = checkSolFeas(LP, sol(2));
%declare output variables
metFeas(1) = feas;
metMWrange = zeros(1,2);
metEle = NaN(numel(model.mets),numel(ele),2);

[metEle(metK,:,1),metEle(metK,:,2)] = deal(metEleK);
if isfield(sol(2), 'full') && numel(sol(2).full) == nVar
    metEle(metU,:,1) = reshape(sol(2).full(1:mU*nE),mU,nE);
    metMWrange(1) = sol(2).obj;
end
if ~(feas <= feasTol)
    %may happen if scaling is allowed and meanwhile the met of interest can
    %have unbounded molecular weight
    fprintf('Infeasible during minimization for molecular weight.\n');
end
LP.osense = -1;
% If infeasible, there are numerical problems. Try using the basis from
%(default not reuse, seems slower)
for j = 1:2
    if nargin < 6
        sol(3) = solveCobraLP(LP);
    else
        sol(3) = solveCobraLP(LP, varargin);
    end
    if j == 1
        if isempty(sol(3).full) || sol(3).stat ~= 1
            if isfield(sol(2), 'basis') && ~isempty(sol(2).basis)
                LP.basis = sol(2).basis;
            end
        else
            %good solution
            break
        end
    end
end
feas = checkSolFeas(LP, sol(3));
metFeas(2) = feas;
if isfield(sol(3), 'full') && numel(sol(3).full) == nVar
    metEle(metU,:,2) = reshape(sol(3).full(1:mU*nE),mU,nE);
    metMWrange(2) = sol(3).obj;
end
rxnBal(:,:,1) = metEle(:,:,1)' * model.S;
rxnBal(:,:,2) = metEle(:,:,2)' * model.S;
if ~(feas <= feasTol)
    %may happen if scaling is allowed and meanwhile the met of interest can
    %have unbounded molecular weight
    fprintf('Infeasible during maximization for molecular weight.\n');
end
formBM0 = [metEle(metInterest0,:,1);metEle(metInterest0,:,2)];
metForm = convertMatrixFormulas(ele,formBM0);
metForm = metForm(:)';

eleG = isnan(c) & ~eleCh(:);
if any(abs(formBM0(:,eleG)) > 1e-6,1)
    fprintf('Biomass contains some generic groups.\n');
end

end
