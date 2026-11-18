% =========================================================================
% LISTING C.4: solve_Eigen_FGM.m (Phu luc C.2, trang 133-135)
% Bai toan gia tri rieng: Dao dong tu do (omega_bar) & Mat on dinh (P_bar)
% =========================================================================
function [omega_bar, P_bar, modes_vib, modes_buck] = solve_Eigen_FGM( ...
    L, b, h, Ec, Em, rhoc, rhom, nu, n, bc, Ne, numModes)

if nargin < 12 || isempty(numModes)
    numModes = 4;
end

% 1. TICH PHAN DO CUNG (A11, B11, D11, A55) VA QUAN TINH KHOI LUONG (I0, I1, I2)
Gc = Ec / (2 * (1 + nu));
Gm = Em / (2 * (1 + nu));
Ks = 5 * (1 + nu) / (6 + 5 * nu);

dE = Ec - Em; dG = Gc - Gm; dRho = rhoc - rhom;
A11 = b * h * (Em + dE / (n + 1));
B11 = b * (h^2) * dE * (1 / (n + 2) - 1 / (2 * (n + 1)));
D11 = b * (h^3) * (Em / 12 + dE * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))));
Cs  = Ks * b * h * (Gm + dG / (n + 1));

I0 = b * h * (rhom + dRho / (n + 1));
I1 = b * (h^2) * dRho * (1 / (n + 2) - 1 / (2 * (n + 1)));
I2 = b * (h^3) * (rhom / 12 + dRho * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))));

% 2. THIET LAP MA TRAN PHAN TU 6x6: [ke], [me], [kge]
numNodes = Ne + 1;
numDof   = 3 * numNodes;
le       = L / Ne;
invLe    = 1 / le;

% Ma tran do cung phan tu ke (SRI cho cat)
ke_mb = zeros(6, 6);
gp = [-1/sqrt(3), 1/sqrt(3)];
for i = 1:2
    B_mb = [-invLe, 0, 0, invLe, 0, 0;
             0,     0, -invLe, 0, 0, invLe];
    ke_mb = ke_mb + (B_mb' * [A11, B11; B11, D11] * B_mb) * (le / 2);
end
Bs  = [0, -invLe, -0.5, 0, invLe, -0.5];
ke  = ke_mb + (Bs' * Cs * Bs) * le;

% Ma tran khoi luong nhat quan me (Consistent Mass Matrix)
cM = le / 6;
me = zeros(6, 6);
me([1,4], [1,4]) = I0 * cM * [2, 1; 1, 2];
me([2,5], [2,5]) = I0 * cM * [2, 1; 1, 2];
me([3,6], [3,6]) = I2 * cM * [2, 1; 1, 2];
me([1,4], [3,6]) = -I1 * cM * [2, 1; 1, 2];
me([3,6], [1,4]) = me([1,4], [3,6])';

% Ma tran do cung hinh hoc kge = Bg' * Bg * le
Bg  = [0, -invLe, 0, 0, invLe, 0];
kge = (Bg' * Bg) * le;

% 3. GHEP NOI TOAN CUC
K  = zeros(numDof, numDof);
M  = zeros(numDof, numDof);
Kg = zeros(numDof, numDof);

for e = 1:Ne
    edof = [3*e-2, 3*e-1, 3*e, 3*e+1, 3*e+2, 3*e+3];
    K(edof, edof)  = K(edof, edof)  + ke;
    M(edof, edof)  = M(edof, edof)  + me;
    Kg(edof, edof) = Kg(edof, edof) + kge;
end

% 4. GIAI BAI TOAN DAO DONG TU DO: eig(K_ff, M_ff)
switch upper(bc)
    case 'S-S'
        fixedVib = [1, 2, 3*Ne + 2];
    case 'C-C'
        fixedVib = [1, 2, 3, 3*Ne + 1, 3*Ne + 2, 3*Ne + 3];
    case 'C-F'
        fixedVib = [1, 2, 3];
end
freeVib = setdiff(1:numDof, fixedVib);

[V_vib, D_vib] = eig(K(freeVib, freeVib), M(freeVib, freeVib));
[lam_vib, idx_vib] = sort(real(diag(D_vib)));
lam_vib = lam_vib(lam_vib > 1e-6);
omega = sqrt(lam_vib(1:min(numModes, length(lam_vib))));

% Chuan hoa tan so rieng khong thu nguyen omega_bar
omega_bar = omega * (L^2 / h) * sqrt(rhom / Em);
modes_vib = V_vib(:, idx_vib(1:length(omega)));

% 5. GIAI BAI TOAN MAT ON DINH TINH (BUCKLING): eig(K_ff, Kg_ff)
%    Luu y: Voi bien ngam C-C trong buckling, giai phong bac tu do doc truc u(L)
switch upper(bc)
    case 'S-S'
        fixedBuck = [1, 2, 3*Ne + 2];
    case 'C-C'
        % Giai phong DOF 3*Ne+1 = u(L) de thanh chiu nen tu do doc truc
        fixedBuck = [1, 2, 3, 3*Ne + 2, 3*Ne + 3];
    case 'C-F'
        fixedBuck = [1, 2, 3];
end
freeBuck = setdiff(1:numDof, fixedBuck);

[V_buck, D_buck] = eig(K(freeBuck, freeBuck), Kg(freeBuck, freeBuck));
[P_cr_all, idx_buck] = sort(real(diag(D_buck)));
valid_idx = find(P_cr_all > 1e-3 & isfinite(P_cr_all));
P_cr = P_cr_all(valid_idx(1:min(numModes, length(valid_idx))));

% Chuan hoa tai nen toi han khong thu nguyen P_bar = P_cr * 12*L^2 / (Em*b*h^3)
P_bar = P_cr * (12 * L^2) / (Em * b * h^3);
modes_buck = V_buck(:, idx_buck(valid_idx(1:length(P_cr))));
end
