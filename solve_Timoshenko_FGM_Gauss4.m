% =========================================================================
% MUC C.1: solve_Timoshenko_FGM_Gauss4.m (Phu luc C.1, trang 127)
% Phien ban thuan Gauss 4 diem qua chieu day (kiem chung Sai so khi n > 5)
% =========================================================================
function [x, w, theta, u, w_max] = solve_Timoshenko_FGM_Gauss4( ...
    L, b, h, Ec, Em, nu, n, q0, P, bc, Ne, loadType)

Gc = Ec / (2 * (1 + nu));
Gm = Em / (2 * (1 + nu));
Ks = 5 * (1 + nu) / (6 + 5 * nu);

% CO DINH 4 DIEM GAUSS-LEGENDRE QUA CHIEU DAY CHO MOI n
xi_g = [-0.861136311594053, -0.339981043584856, ...
         0.339981043584856,  0.861136311594053];
w_g  = [ 0.347854845137454,  0.652145154862546, ...
         0.652145154862546,  0.347854845137454];

A11 = 0; B11 = 0; D11 = 0; A55 = 0;
for k = 1:4
    z   = 0.5 * h * xi_g(k);
    Vc  = (z / h + 0.5)^n;
    Ez  = (Ec - Em) * Vc + Em;
    Gz  = (Gc - Gm) * Vc + Gm;
    jac = 0.5 * h * w_g(k);
    A11 = A11 + b * Ez * jac;
    B11 = B11 + b * Ez * z * jac;
    D11 = D11 + b * Ez * (z^2) * jac;
    A55 = A55 + b * Gz * jac;
end
Cs = Ks * A55;

numNodes = Ne + 1; numDof = 3 * numNodes; le = L / Ne; invLe = 1 / le;
K_global = zeros(numDof, numDof);
F_global = zeros(numDof, 1);

% Uon-Mang 2 Gauss + Cat SRI 1 Gauss
gp = [-1 / sqrt(3), 1 / sqrt(3)];
Ke_mb = zeros(6, 6);
for i = 1:2
    B_mb = [-invLe, 0, 0, invLe, 0, 0;
             0,     0, -invLe, 0, 0, invLe];
    Ke_mb = Ke_mb + (B_mb' * [A11, B11; B11, D11] * B_mb) * (le / 2);
end
Bs = [0, -invLe, -0.5, 0, invLe, -0.5];
Ke = Ke_mb + (Bs' * Cs * Bs) * le;

for e = 1:Ne
    edof = [3*e-2, 3*e-1, 3*e, 3*e+1, 3*e+2, 3*e+3];
    K_global(edof, edof) = K_global(edof, edof) + Ke;
    if strcmpi(loadType, 'Uniform')
        F_global(edof) = F_global(edof) + [0; q0*le/2; 0; 0; q0*le/2; 0];
    end
end

if strcmpi(loadType, 'Point')
    if strcmpi(bc, 'C-F')
        F_global(3*Ne + 2) = F_global(3*Ne + 2) + P;
    else
        midNode = floor(numNodes / 2) + 1;
        F_global(3*(midNode-1) + 2) = F_global(3*(midNode-1) + 2) + P;
    end
end

switch upper(bc)
    case 'S-S', prescribedDof = [1, 2, 3*Ne + 2];
    case 'C-C', prescribedDof = [1, 2, 3, 3*Ne + 1, 3*Ne + 2, 3*Ne + 3];
    case 'C-F', prescribedDof = [1, 2, 3];
end
activeDof = setdiff(1:numDof, prescribedDof);

Q = zeros(numDof, 1);
Q(activeDof) = K_global(activeDof, activeDof) \ F_global(activeDof);

x = linspace(0, L, numNodes)';
u = Q(1:3:end); w = Q(2:3:end); theta = Q(3:3:end);
w_max = max(abs(w));
end
