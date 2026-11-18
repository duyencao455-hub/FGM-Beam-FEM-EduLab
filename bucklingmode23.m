%% ========================================================================
%% BUCKLING MODE 2 VÀ MODE 3 DẦM FGM (Al2O3/Al, nu = 0.23)
%% ĐẦY ĐỦ 2 NHÁNH: TIMOSHENKO (TBT) VÀ EULER-BERNOULLI (EBT)
%% ========================================================================
clear; clc; close all;

% 1. CẤU HÌNH THÔNG SỐ VẬT LIỆU VÀ HÌNH HỌC CHUẨN
Ec = 380e9;   % Gốm Al2O3 (Pa)
Em = 70e9;    % Kim loại Al (Pa)
nu = 0.23;    % Hệ số Poisson
Ks = 5/6;     % Hệ số hiệu chỉnh cắt
b = 1.0;      % Chiều rộng (m)
h = 0.1;      % Chiều dày (m)
Ne = 50;      % Mật độ lưới 50 phần tử

p_list = [0, 1, 2, 5, 10];
L_over_h_list = [5, 10];

% 2. VÒNG LẶP KHẢO SÁT
for L_h = L_over_h_list
    L = L_h * h;
    fprintf('\n========================================================================================================\n');
    fprintf('  BẢNG KIỂM CHỨNG BUCKLING MODE 2 & MODE 3 (TỶ SỐ L/h = %d)\n', L_h);
    fprintf('========================================================================================================\n');
    
    for BC_cell = {'C-C', 'S-S'}
        BC = BC_cell{1};
        
        for mode_num = [2, 3]
            fprintf('\n>>> ĐIỀU KIỆN BIÊN: %s | KHẢO SÁT: MODE %d <<<\n', BC, mode_num);
            fprintf('--------------------------------------------------------------------------------------------------------\n');
            fprintf('| %-4s | %-12s | %-12s | %-13s | %-12s | %-12s | %-13s |\n', ...
                'p', 'Exact (TBT)', 'FEM (TBT)', 'Sai số TBT(%)', 'Exact (EBT)', 'FEM (EBT)', 'Sai số EBT(%)');
            fprintf('--------------------------------------------------------------------------------------------------------\n');
            
            for p = p_list
                % A. ĐẶC TRƯNG ĐỘ CỨNG MẶT CẮT FGM
                A = b * h * (Em + (Ec - Em)/(p + 1));
                B = b * h^2 * (Ec - Em) * p / (2 * (p + 1) * (p + 2));
                D = b * h^3 * (Em/12 + (Ec - Em) * (1/(p + 3) - 1/(p + 2) + 1/(4*(p + 1))));
                D_star = D - (B^2)/A;
                A55 = A / (2 * (1 + nu));
                cs = 1 / (Ks * A55);
                
                % B. NGHIỆM GIẢI TÍCH (EXACT) CHO MODE 2 VÀ MODE 3
                if strcmp(BC, 'S-S')
                    beta_sq_list = [(2*pi)^2, (3*pi)^2]; % Mode 2, Mode 3
                else % C-C
                    beta_sq_list = [80.762899, 157.913670]; % Mode 2, Mode 3
                end
                beta_sq = beta_sq_list(mode_num - 1);
                
                PE_exact = (beta_sq * D_star) / (L^2);
                PT_exact = PE_exact / (1 + cs * PE_exact);
                
                Pbar_EBT_exact = PE_exact * (12 * L^2) / (Em * b * h^3);
                Pbar_TBT_exact = PT_exact * (12 * L^2) / (Em * b * h^3);
                
                % C. GIẢI PHẦN TỬ HỮU HẠN (FEM TBT & FEM EBT)
                Pbar_TBT_fem = solve_FEM_TBT(Ec, Em, nu, p, b, h, L, Ne, BC, Ks, mode_num);
                Pbar_EBT_fem = solve_FEM_EBT(D_star, Em, b, h, L, Ne, BC, mode_num);
                
                % D. TÍNH SAI SỐ
                err_TBT = abs(Pbar_TBT_fem - Pbar_TBT_exact) / Pbar_TBT_exact * 100;
                err_EBT = abs(Pbar_EBT_fem - Pbar_EBT_exact) / Pbar_EBT_exact * 100;
                
                fprintf('| %-4d | %-12.3f | %-12.3f | %-13.2f | %-12.3f | %-12.3f | %-13.2f |\n', ...
                    p, Pbar_TBT_exact, Pbar_TBT_fem, err_TBT, Pbar_EBT_exact, Pbar_EBT_fem, err_EBT);
            end
            fprintf('--------------------------------------------------------------------------------------------------------\n');
        end
    end
end

%% ==================== BỘ GIẢI FEM TIMOSHENKO (TBT + SRI) ====================
function P_bar = solve_FEM_TBT(Ec, Em, nu, p, b, h, L, Ne, BC, Ks, mode_target)
    A = b * h * (Em + (Ec - Em)/(p + 1));
    B = b * h^2 * (Ec - Em) * p / (2 * (p + 1) * (p + 2));
    D = b * h^3 * (Em/12 + (Ec - Em) * (1/(p + 3) - 1/(p + 2) + 1/(4*(p + 1))));
    A55 = A / (2 * (1 + nu));
    
    Le = L / Ne;
    num_dofs = (Ne + 1) * 3;
    K_glob = zeros(num_dofs, num_dofs);
    Kg_glob = zeros(num_dofs, num_dofs);
    
    Bm = [-1/Le, 0, 0, 1/Le, 0, 0];
    Bb = [0, 0, 1/Le, 0, 0, -1/Le];
    Bs_0 = [0, -1/Le, -0.5, 0, 1/Le, -0.5]; % SRI tại tâm phần tử
    Bg = [0, -1/Le, 0, 0, 1/Le, 0];
    
    Ke = (A*(Bm'*Bm) + D*(Bb'*Bb) + B*(Bm'*Bb + Bb'*Bm) + (Ks*A55)*(Bs_0'*Bs_0)) * Le;
    Kge = (Bg' * Bg) * Le;
    
    for e = 1:Ne
        dofs = (e-1)*3 + (1:6);
        K_glob(dofs, dofs) = K_glob(dofs, dofs) + Ke;
        Kg_glob(dofs, dofs) = Kg_glob(dofs, dofs) + Kge;
    end
    
    if strcmp(BC, 'S-S')
        fixed_dofs = [1, 2, num_dofs-1];
    else % C-C
        fixed_dofs = [1, 2, 3, num_dofs-1, num_dofs];
    end
    active_dofs = setdiff(1:num_dofs, fixed_dofs);
    
    [~, D_eig] = eig(K_glob(active_dofs, active_dofs), Kg_glob(active_dofs, active_dofs));
    eig_vals = diag(D_eig);
    pos_idx = find(eig_vals > 1e-3);
    sorted_vals = sort(eig_vals(pos_idx), 'ascend');
    
    P_cr = sorted_vals(mode_target);
    P_bar = P_cr * (12 * L^2) / (Em * b * h^3);
end

%% ==================== BỘ GIẢI FEM EULER-BERNOULLI (HERMITE C1) ====================
function P_bar = solve_FEM_EBT(D_star, Em, b, h, L, Ne, BC, mode_target)
    Le = L / Ne;
    num_dofs = (Ne + 1) * 2; % Mỗi nút 2 bậc tự do: [w, phi]
    K_glob = zeros(num_dofs, num_dofs);
    Kg_glob = zeros(num_dofs, num_dofs);
    
    % Ma trận độ cứng uốn Hermite chuẩn
    Ke = (D_star / Le^3) * [ ...
        12,     6*Le,   -12,     6*Le;
        6*Le,   4*Le^2, -6*Le,   2*Le^2;
       -12,    -6*Le,    12,    -6*Le;
        6*Le,   2*Le^2, -6*Le,   4*Le^2];
    
    % Ma trận độ cứng hình học Hermite chuẩn (Phương trình 3.94 trong luận văn)
    Kge = (1 / (30*Le)) * [ ...
        36,     3*Le,   -36,     3*Le;
        3*Le,   4*Le^2, -3*Le,  -Le^2;
       -36,    -3*Le,    36,    -3*Le;
        3*Le,  -Le^2,   -3*Le,   4*Le^2];
    
    for e = 1:Ne
        dofs = (e-1)*2 + (1:4);
        K_glob(dofs, dofs) = K_glob(dofs, dofs) + Ke;
        Kg_glob(dofs, dofs) = Kg_glob(dofs, dofs) + Kge;
    end
    
    if strcmp(BC, 'S-S')
        fixed_dofs = [1, num_dofs-1]; % w1 = 0, w_end = 0 (góc xoay tự do)
    else % C-C
        fixed_dofs = [1, 2, num_dofs-1, num_dofs]; % w1=0, phi1=0, w_end=0, phi_end=0
    end
    active_dofs = setdiff(1:num_dofs, fixed_dofs);
    
    [~, D_eig] = eig(K_glob(active_dofs, active_dofs), Kg_glob(active_dofs, active_dofs));
    eig_vals = diag(D_eig);
    pos_idx = find(eig_vals > 1e-3);
    sorted_vals = sort(eig_vals(pos_idx), 'ascend');
    
    P_cr = sorted_vals(mode_target);
    P_bar = P_cr * (12 * L^2) / (Em * b * h^3);
end