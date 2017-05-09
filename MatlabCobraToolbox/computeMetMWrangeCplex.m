function [metMWrange,metForm,metFeas,rxnBal,ele,metEle,LP] = computeMetMWrangeCplex(model,metKnown,metInterest,rxns,percent,param)
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
%    param         parameter structure for Cplex
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
if nargin < 6
    param = struct();
end
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
%underscore, followed by a number for the stoichiometry. 
%Brackets/parentheses are also supported.
[metK,metKform] = deal(metK(~metKform), model.metFormulas(metK(~metKform)));
%Now handled by checkEleBalance
% re = regexp(metKform,'[A-Z][a-z_]*(\-?\d+\.?\d*)?','match');
% re = cellfun(@(x) strjoin(x,''),re,'UniformOutput',false);
% goodForm = strcmp(re, strtrim(metKform));
% if ~all(goodForm)
%     goodForm = find(~goodForm,1);
%     error('%s has an invalid formula %s\n',metKnown{goodForm},metKform{goodForm});
% end

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
LP = Cplex();
%n_ik, stoichoimetry for element k in met i
LP.addCols(zeros(mU*nE,1),[],zeros(mU*nE,1),inf(mU*nE,1),[],char(nameM));
if any(eleCh)
    %if charge is in the formula, allow negative charges
    eleChK = find(eleCh);
    LP.Model.lb((mU*(eleChK-1) + 1):(mU*eleChK)) = -inf;
end
%x^+_jk, positive inconsistency of element k in reaction j
LP.addCols(ones(nR*nE,1),[],zeros(nR*nE,1),inf(nR*nE,1),[],char(strcat(nameR,'_pos')));
%x^-_jk, negative inconsistency of element k in reaction j
LP.addCols(ones(nR*nE,1),[],zeros(nR*nE,1),inf(nR*nE,1),[],char(strcat(nameR,'_neg')));
%RHS
b = -model.S(metK,rxnC)' * metEleK;
b = b(:);
%constraint matrix
[row,col,entry] = find(model.S(metU,rxnC)');
row = repmat(row(:),nE,1)+reshape(repmat(0:nR:nR*(nE-1),numel(row),1),numel(row)*nE,1);
col = repmat(col(:),nE,1)+reshape(repmat(0:mU:mU*(nE-1),numel(col),1),numel(col)*nE,1);
entry = repmat(entry,nE,1);
A = [sparse(row,col,entry,nR*nE,mU*nE), sparse(1:nR*nE,1:nR*nE,ones(nR*nE,1),nR*nE,nR*nE), ...
    sparse(1:nR*nE,1:nR*nE,-ones(nR*nE,1),nR*nE,nR*nE)];
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
    %terminate if infeasible (should not happen)
    if isfield(LP.Solution,'x')
        metEle = zeros(numel(model.mets),numel(ele));
        metEle(metK,:) = metEleK;
        metEle(metU,:) = reshape(LP.Solution.x(1:mU*nE),mU,nE);
        rxnBal = metEle' * model.S;
    else
        metEle = [];
        rxnBal = [];
    end
    fprint('Infeasible during optimization for minimal inconsistency.');
    return
end
x = LP.Solution.x;

%% metabolite MW range under minimum inconsistency
%give bounds for the met of interest
LP.Model.ub(metInterest:mU:(metInterest+mU*(nE-1))) = 1e7;
if any(eleCh)
    LP.Model.lb(metInterest+mU*(find(eleCh)-1)) = -1e7;
end
% Constraint element-wise, sum(inconsist_e) <= min_inconsist_e * (1+percent)
for k = 1:nE
    ind = (mU*nE+nR*(k-1)+1):(mU*nE+nR*k);
    LP.addRows(0,...
        sparse(ones(nR*2,1),[ind, (ind + nR*nE)], ones(nR*2,1),1,size(LP.Model.A,2)),...
        sum(x([ind, (ind + nR*nE)]))*(1+1e-7+abs(percent))+1e-5, strcat('min_incon_',ele{k}));
end
%molecular weight of each element
c = MW(ele);
c2 = c;
c2(isnan(c)) = 0;
%change objective to the molecular weight of the met of interest
LP.Model.obj(:) = 0;
LP.Model.obj(metInterest:mU:(metInterest+mU*(nE-1))) = c2;
LP.Model.sense = 'minimize';
LP.solve;
feas = checkSolFeas(LP);
%declare output variables
metFeas(1) = feas;
metMWrange = zeros(1,2);
metEle = NaN(numel(model.mets),numel(ele),2);

[metEle(metK,:,1),metEle(metK,:,2)] = deal(metEleK);
if isfield(LP.Solution,'x')
    metEle(metU,:,1) = reshape(LP.Solution.x(1:mU*nE),mU,nE);
    metMWrange(1) = LP.Model.obj(:)' * LP.Solution.x;
end
if ~(feas <= feasTol)
    %may happen if scaling is allowed and meanwhile the met of interest can
    %have unbounded molecular weight
    fprintf('Infeasible during minimization for biomass weight.\n');
end
LP.Model.sense = 'maximize';
LP.solve;
feas = checkSolFeas(LP);
metFeas(2) = feas;
if isfield(LP.Solution,'x')
    metEle(metU,:,2) = reshape(LP.Solution.x(1:mU*nE),mU,nE);
    metMWrange(2) = LP.Model.obj(:)' * LP.Solution.x;
end
rxnBal(:,:,1) = metEle(:,:,1)' * model.S;
rxnBal(:,:,2) = metEle(:,:,2)' * model.S;
if ~(feas <= feasTol)
    %may happen if scaling is allowed and meanwhile the met of interest can
    %have unbounded molecular weight
    fprintf('Infeasible during maximization for biomass weight.\n');
end
formBM0 = [metEle(metInterest0,:,1);metEle(metInterest0,:,2)];
metForm = convertMatrixFormulas(ele,formBM0);
metForm = metForm(:)';

eleG = isnan(c) & ~eleCh(:);
if any(abs(formBM0(:,eleG)) > 1e-6,1)
    fprintf('Biomass contains some generic groups.\n');
end

end

function [mw, elRes, stRes] = MW(form)
%Return the molecular weight for the cell array or string of formula 'form'
%The formulas must not have '(' and ')'.
%Must be an elemental symbol followed by a number.
%
%Output:
%  mw:      molecular weight, in g/mol
%  elRes:   Unrecognized elements in form. MW assumes zero weight for them.
%  stRes:   the elemental matrix for the columns corresponding to elRes
%
%Siu Hung Joshua Chan Nov 2016

[~,element, metEle] = checkEleBalance(form);
H = 1.00794; He = 4.002602; Li = 6.941; Be = 9.012182; B = 10.811; C = 12.011; N = 14.00674;
O = 15.9994; F = 18.9984032; Ne = 20.1797; Na = 22.989768; Mg = 24.305; Al = 26.981539; Si = 28.0855;
P = 30.973762; S = 32.066; Cl = 35.4527; Ar = 39.948; K = 39.0983; Ca = 40.078; Sc = 44.95591;
Ti = 47.88; V = 50.9415; Cr = 51.9961; Mn = 54.93805; Fe = 55.847; Co = 58.9332; Ni = 58.69;
Cu = 63.546; Zn = 65.39; Ga = 69.723; Ge = 72.61; As = 74.92159; Se = 78.96; Br = 79.904; Kr = 83.8;
Rb = 85.4678; Sr = 87.62; Y = 88.90585; Zr = 91.224; Nb = 92.90638; Mo = 95.94; Tc = 98.9063;
Ru = 101.07; Rh = 102.9055; Pd = 106.42; Ag = 107.8682; Cd = 112.411; In = 114.82; Sn = 118.71;
Sb = 121.75; Te = 127.6; I = 126.90447; Xe = 131.29; Cs = 132.90543; Ba = 137.327; La = 138.9055;
Ce = 140.115; Pr = 140.90765; Nd = 144.24; Pm = 146.9151; Sm = 150.36; Eu = 151.965; Gd = 157.25;
Tb = 158.92534; Dy = 162.5; Ho = 164.93032; Er = 167.26; Tm = 168.93421; Yb = 173.04; Lu = 174.967;
Hf = 178.49; Ta = 180.9479; W = 183.85; Re = 186.207; Os = 190.2; Ir = 192.22; Pt = 195.08;
Au = 196.96654; Hg = 200.59; Tl = 204.3833; Pb = 207.2; Bi = 208.98037; Po = 208.9824; At = 209.9871;
Rn = 222.0176; Ac = 223.0197; Th = 226.0254; Pa = 227.0278; U = 232.0381; Np = 231.0359; Pu = 238.0289;
Am = 237.0482; Cm = 244.0642; Bk = 243.0614; Cf = 247.0703; Es = 247.0703; Fm = 251.0796; Md = 252.0829;
No = 257.0951; Lr = 258.0986; Rf = 259.1009; Db = 260.1053; Sg = 261.1087; Bh = 262.1138; Hs = 263.1182;
Mt = 262.1229;
% Unamed = NaN;
residue = false(numel(element),1);
for j = 1:numel(element)
    if ~exist(element{j}, 'var')
        residue(j) = true;
    end
end
if all(residue)
    mw = NaN(size(metEle,1),1);
else
    v = ['[ ' strjoin(element(~residue), ' ') ' ]'];
    eval(['mw = metEle(:,~residue) * ' v ''';']);
    mw(mw == 0) = NaN;
end
    elRes = element(residue);
    stRes = metEle(:, residue);
end