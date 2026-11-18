% =========================================================================
% MUC C.1: main_benchmark_bending.m (Phu luc C.1, trang 127)
% Chuong trinh dieu phoi kiem chung Uon tinh, Khoa cat & Tich phan Lai
% =========================================================================
clear; clc; close all;

% 1. THONG SO CHUAN REDDY (E1 = 14.4 GPa, E2 = 1.44 GPa, nu = 0.38, n = 5)
L = 1.0; b = 1.0; L_h = 20; h = L / L_h;
Ec = 14.4e9; Em = 1.44e9; nu = 0.38; q0 = 1.0; P = 1.0;
scale = 100 * Em * b * (h^3) / (q0 * L^4);

[~, ~, w_bar_TBT_ex, w_bar_EBT_ex] = analytical_FGM(L, b, h, Ec, Em, nu, 5, q0, P, 'Uniform');
fprintf('=== 1. HOI TU LUOI TIMOSHENKO & EULER (BANG 3.2, n=5, L/h=20) ===\n');
fprintf('Nghiem giai tich Navier: EBT = %.4f | TBT = %.4f\n', w_bar_EBT_ex, w_bar_TBT_ex);

N_list = [2, 5, 10, 20, 30, 40, 50, 100];
for N = N_list
    [~, ~, ~, ~, w_tbt] = solve_Timoshenko_FGM(L, b, h, Ec, Em, nu, 5, q0, P, 'S-S', N, 'Uniform', 'SRI');
    [~, ~, ~, ~, w_ebt] = solve_Euler_FGM(L, b, h, Ec, Em, 5, q0, P, 'S-S', N, 'Uniform');
    wb_tbt = scale * w_tbt; wb_ebt = scale * w_ebt;
    err_tbt = abs(wb_tbt - w_bar_TBT_ex) / w_bar_TBT_ex * 100;
    err_ebt = abs(wb_ebt - w_bar_EBT_ex) / w_bar_EBT_ex * 100;
    fprintf('N = %3d | EBT = %.4f (%.2f%%) | TBT = %.4f (%.2f%%)\n', ...
        N, wb_ebt, err_ebt, wb_tbt, err_tbt);
end

fprintf('\n=== 2. KHAO SAT KHOA CAT FI VS SRI (BANG D.5, N=20, n=5) ===\n');
for lh = [4, 10, 20, 50, 100]
    h_k = L / lh;
    sc_k = 100 * Em * b * (h_k^3) / (q0 * L^4);
    [~, ~, ~, ~, w_fi]  = solve_Timoshenko_FGM(L, b, h_k, Ec, Em, nu, 5, q0, P, 'S-S', 20, 'Uniform', 'FI');
    [~, ~, ~, ~, w_sri] = solve_Timoshenko_FGM(L, b, h_k, Ec, Em, nu, 5, q0, P, 'S-S', 20, 'Uniform', 'SRI');
    fprintf('L/h = %3d | FI = %.4f | SRI = %.4f\n', lh, sc_k*w_fi, sc_k*w_sri);
end

fprintf('\n=== 3. GAUSS 4 DIEM VS TICH PHAN LAI (BANG D.6, N=100, L/h=20) ===\n');
for n_val = [1, 2, 5, 6, 8, 10, 50]
    [~, ~, w_ex] = analytical_FGM(L, b, h, Ec, Em, nu, n_val, q0, P, 'Uniform');
    [~, ~, ~, ~, w_g4]  = solve_Timoshenko_FGM_Gauss4(L, b, h, Ec, Em, nu, n_val, q0, P, 'S-S', 100, 'Uniform');
    [~, ~, ~, ~, w_hyb] = solve_Timoshenko_FGM(L, b, h, Ec, Em, nu, n_val, q0, P, 'S-S', 100, 'Uniform', 'SRI');
    err_g4  = abs(scale*w_g4  - w_ex) / w_ex * 100;
    err_hyb = abs(scale*w_hyb - w_ex) / w_ex * 100;
    fprintf('n = %2d | Gauss4 err = %6.2f%% | Hybrid err = %5.2f%%\n', n_val, err_g4, err_hyb);
end
