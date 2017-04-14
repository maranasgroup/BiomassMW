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