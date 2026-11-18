% =========================================================================
% LISTING C.3: analytical_FGM.m (Phu luc C.1, trang 132-133)
% Nghiem giai tich Navier cho dam FGM tua don (S-S) theo EBT va TBT
% =========================================================================
function [w_mid_TBT, w_mid_EBT, w_bar_TBT, w_bar_EBT, D_star, z_na] = ...
    analytical_FGM(L, b, h, Ec, Em, nu, n, q0, P, loadType, numTerms)

if nargin < 11 || isempty(numTerms)
    numTerms = 200;
end

% 1. DO CUNG MAT CAT GIAI TICH DONG
Gc = Ec / (2 * (1 + nu));
Gm = Em / (2 * (1 + nu));
dE = Ec - Em;
dG = Gc - Gm;

A11 = b * h * (Em + dE / (n + 1));
B11 = b * (h^2) * dE * (1 / (n + 2) - 1 / (2 * (n + 1)));
D11 = b * (h^3) * (Em / 12 + dE * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))));
A55 = b * h * (Gm + dG / (n + 1));

D_star = D11 - (B11^2) / A11;       % Do cung uon hieu dung
z_na   = B11 / A11;                 % Vi tri truc trung hoa vat ly
Ks     = 5 * (1 + nu) / (6 + 5 * nu); % He so cat Reddy
Cs     = Ks * A55;

% 2. TINH TONG CHUOI FOURIER NAVIER TAI TAM NHIP x = L/2
w_mid_EBT = 0;
w_mid_TBT = 0;

for m = 1:2:(2*numTerms - 1)
    alpha_m = m * pi / L;
    if strcmpi(loadType, 'Uniform')
        qm = 4 * q0 / (m * pi);
        q_ref = q0;
    else
        % Tai tap trung P dat tai x = L/2
        qm = (2 * P / L) * sin(alpha_m * L / 2);
        q_ref = P / L;
    end
    
    w_ebt_m = qm / ((alpha_m^4) * D_star);
    w_tbt_m = w_ebt_m + qm / ((alpha_m^2) * Cs);
    
    s_mid = sin(alpha_m * L / 2);
    w_mid_EBT = w_mid_EBT + w_ebt_m * s_mid;
    w_mid_TBT = w_mid_TBT + w_tbt_m * s_mid;
end

% 3. CHUAN HOA DO VONG KHONG THU NGUYEN (w_bar * 10^2)
%    Theo quy uoc Reddy: w_bar = 100 * Em * b * h^3 / (q0 * L^4) * w_mid
I_ref = b * (h^3) / 12;
if strcmpi(loadType, 'Uniform')
    scaleFactor = 100 * Em * b * (h^3) / (q0 * L^4);
else
    scaleFactor = 100 * Em * I_ref / (P * L^3);
end

w_bar_EBT = scaleFactor * w_mid_EBT;
w_bar_TBT = scaleFactor * w_mid_TBT;
end
