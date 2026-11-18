% =========================================================================
% BÀI TOÁN 4: KHẢO SÁT CHỈ SỐ VẬT LIỆU n (BẢNG 1.9 & 1.10)
% =========================================================================
clc; clear; close all;

% =========================================================================
% 1. THÔNG SỐ ĐẦU VÀO CƠ BẢN
% =========================================================================
h = 88e-6;              % Chiều cao (m)
L = 20 * h;             % Chiều dài dầm (m)
b = 2 * h;              % Bề rộng (m)
Ec = 14.4e9;            % Gốm (Pa)
Em = 1.44e9;            % Kim loại (Pa)
nu = 0.38;              % Hệ số Poisson
N_elem = 100;           % Số phần tử
load_val = 1;           % Giá trị tải (1 N/m hoặc 1 N)
BC = 2;                 % Điều kiện biên: 2 = Tựa đơn (S-S)

% BẠN CÓ THỂ THAY ĐỔI MẢNG NÀY TÙY Ý (Ví dụ: [0; 0.5; 2; 5; 100])
n_array = [0; 1; 5; 6; 7; 10]; 

% Khởi tạo mảng lưu trữ kết quả TÍNH TOÁN THỰC TẾ (Giải tích & PTHH)
w_Ref_BET_q = zeros(length(n_array), 1);
w_Ref_TBT_q = zeros(length(n_array), 1);
w_Ref_BET_P = zeros(length(n_array), 1);
w_Ref_TBT_P = zeros(length(n_array), 1);

w_FEM_Gauss4_q = zeros(length(n_array), 1); 
w_FEM_BET_q    = zeros(length(n_array), 1);
w_FEM_TBT_q    = zeros(length(n_array), 1);
w_FEM_BET_P    = zeros(length(n_array), 1);
w_FEM_TBT_P    = zeros(length(n_array), 1);

% =========================================================================
% 2. VÒNG LẶP TÍNH TOÁN (GIẢI TÍCH CHÍNH XÁC & PHẦN TỬ HỮU HẠN)
% =========================================================================
disp('Đang chạy giải thuật Toán học và Phần tử hữu hạn, vui lòng đợi...');
for i = 1:length(n_array)
    n = n_array(i);
    
    % --- A. TÍNH NGHIỆM GIẢI TÍCH CHÍNH XÁC (Dùng làm Tham chiếu) ---
    [w_Ref_BET_q(i), w_Ref_TBT_q(i), w_Ref_BET_P(i), w_Ref_TBT_P(i)] = ...
        analytical_FGM(L, h, b, Ec, Em, nu, n, load_val);

    % --- B. CHẠY MÔ HÌNH THUẦN GAUSS 4 ĐIỂM (Để báo cáo sai số) ---
    [~, w_FEM_Gauss4_q(i), ~] = solve_Timoshenko_FGM_Gauss4(L, h, b, Ec, Em, nu, n, N_elem, load_val, 1, BC);

    % --- C. CHẠY MÔ HÌNH HYBRID TẢI PHÂN BỐ (load_type = 1) ---
    [~, w_FEM_BET_q(i), ~] = solve_Euler_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 1, BC);
    [~, w_FEM_TBT_q(i), ~] = solve_Timoshenko_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 1, BC);
    
    % --- D. CHẠY MÔ HÌNH HYBRID TẢI TẬP TRUNG (load_type = 2) ---
    [~, w_FEM_BET_P(i), ~] = solve_Euler_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 2, BC);
    [~, w_FEM_TBT_P(i), ~] = solve_Timoshenko_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 2, BC);
end
disp('Hoàn tất tính toán! Đang xuất dữ liệu...');

% Hàm vô danh định dạng số tiếng Việt (Dấu phẩy thập phân)
fmt = @(num, dec) strrep(sprintf(sprintf('%%.%df', dec), num), '.', ',');

% =========================================================================
% 3. XUẤT BẢNG 1.9 & 1.10 RA COMMAND WINDOW
% =========================================================================
fprintf('\nBảng 1.9: Sai số của sơ đồ Tích phân Gauss 4 điểm theo chỉ số n.\n');
fprintf('----------------------------------------------------------------------------------\n');
fprintf(' %-9s | %-20s | %-22s | %-20s \n', 'Chỉ số n', 'Gauss 4 điểm (FEM)', 'Giải tích (Ferreira)', 'Sai số tương đối (%)');
fprintf('----------------------------------------------------------------------------------\n');
for i = 1:length(n_array)
    err_G4 = abs(w_FEM_Gauss4_q(i) - w_Ref_TBT_q(i)) / w_Ref_TBT_q(i) * 100;
    fprintf(' %-9s | %-20s | %-22s | %-20s \n', ...
        num2str(n_array(i)), fmt(w_FEM_Gauss4_q(i),4), fmt(w_Ref_TBT_q(i),4), fmt(err_G4,2));
end
fprintf('----------------------------------------------------------------------------------\n\n');

fprintf('Bảng 1.10: Kiểm chứng Thuật toán Hybrid - Độ võng chuẩn hóa w.\n');
fprintf('------------------------------------------------------------------------------------------\n');
fprintf(' %-4s | %-15s | %-19s | %-13s | %-18s \n', 'n', 'Gauss 4 điểm', 'Hybrid (đề xuất)', 'Tham chiếu', 'Sai số Hybrid (%)');
fprintf('------------------------------------------------------------------------------------------\n');
for i = 1:length(n_array)
    err_hyb = abs(w_FEM_TBT_q(i) - w_Ref_TBT_q(i)) / w_Ref_TBT_q(i) * 100; 
    fprintf(' %-4s | %-15s | %-19s | %-13s | %-18s \n', ...
        num2str(n_array(i)), fmt(w_FEM_Gauss4_q(i),4), fmt(w_FEM_TBT_q(i),4), fmt(w_Ref_TBT_q(i),4), fmt(err_hyb,2));
end
fprintf('------------------------------------------------------------------------------------------\n\n');

%% ========================================================================
% HÀM CON 3: GIẢI DẦM TIMOSHENKO (THUẦN GAUSS 4 ĐIỂM)
% ========================================================================
function [w_max, w_star, U] = solve_Timoshenko_FGM_Gauss4(L, h, b, E1, E2, nu, n, ne, q0, load_type, BC)
    nn = ne + 1; ks = 5 * (1 + nu) / (6 + 5 * nu); % He so Cowper (Reddy 2017); 
    A = 0; B = 0; D = 0; S_val = 0;
    xi_z = [-0.861136, -0.339981, 0.339981, 0.861136]; wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];
    for i = 1:4
        z = (h/2) * xi_z(i); J = h/2; Vc = (z/h + 0.5)^n; Ez = (E1 - E2) * Vc + E2; Gz = Ez / (2 * (1 + nu));
        A = A + wi_z(i) * (Ez) * J; B = B + wi_z(i) * (z * Ez) * J; D = D + wi_z(i) * (z^2 * Ez) * J; S_val = S_val + wi_z(i) * (ks * Gz) * J;
    end
    Cb = [A*b, B*b; B*b, D*b]; Cs = S_val * b;
    GDof = 3 * nn; K = zeros(GDof); F = zeros(GDof, 1); le = L / ne;
    for e = 1:ne
        idx_e = (3*e-2):(3*e+3); ke = zeros(6);
        for k = [-0.57735, 0.57735]
            dN_dx = [-0.5, 0.5] * (2/le); Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 dN_dx(1) 0 0 dN_dx(2)];
            ke = ke + Bb' * Cb * Bb * (le/2);
        end
        Bs = [0 -1/le -0.5 0 1/le -0.5]; ke = ke + Bs' * Cs * Bs * le; 
        fe = zeros(6,1); if load_type == 1, fe(2) = q0 * le / 2; fe(5) = q0 * le / 2; end
        K(idx_e, idx_e) = K(idx_e, idx_e) + ke; F(idx_e) = F(idx_e) + fe;
    end
    if load_type == 2, center_node = floor(ne / 2) + 1; F(3 * center_node - 1) = q0; end
    if BC == 1, fixedNodes = [1, 2, 3, GDof-2, GDof-1, GDof]; elseif BC == 2, fixedNodes = [1, 2, 3*nn-1]; elseif BC == 3, fixedNodes = [1, 2, 3]; end
    activeDof = setdiff(1:GDof, fixedNodes); U = zeros(GDof, 1); U(activeDof) = K(activeDof, activeDof) \ F(activeDof);
    w_max = max(abs(U(2:3:end))); I0 = (b * h^3) / 12;
    if load_type == 1, w_star = (w_max * E2 * I0) / (abs(q0) * L^4) * 100; 
    else,              w_star = (w_max * E2 * I0) / (abs(q0) * L^3) * 100; end
end

%% ========================================================================
% HÀM CON 4: TÍNH NGHIỆM GIẢI TÍCH (ANALYTICAL EXACT SOLUTION)
% Dựa trên chuỗi Navier cho dầm Tựa Đơn (Simply Supported)
% ========================================================================
function [w_BET_q, w_TBT_q, w_BET_P, w_TBT_P] = analytical_FGM(L, h, b, E1, E2, nu, n, q0)
    % Tính toán chính xác các hệ số ma trận vật liệu (A, B, D)
    E_diff = E1 - E2;
    A = b * h * (E_diff / (n + 1) + E2);
    B = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
    D = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
    
    % Độ cứng cắt
    ks = 5 * (1 + nu) / (6 + 5 * nu); % He so Cowper (Reddy 2017);
    S_val = ks * A / (2 * (1 + nu));
    
    % Độ cứng uốn tương đương (Khử hiệu ứng ghép màng-uốn B^2/A)
    EI_eq = D - (B^2 / A);
    GA_eq = S_val;
    
    % Tính độ võng cực đại (w_max) theo lý thuyết đàn hồi kinh điển
    % 1. Tải phân bố đều (q)
    w_max_BET_q = (5 * abs(q0) * L^4) / (384 * EI_eq);
    w_max_TBT_q = w_max_BET_q + (abs(q0) * L^2) / (8 * GA_eq);
    
    % 2. Tải tập trung giữa nhịp (P)
    w_max_BET_P = (abs(q0) * L^3) / (48 * EI_eq);
    w_max_TBT_P = w_max_BET_P + (abs(q0) * L) / (4 * GA_eq);
    
    % Chuẩn hóa kết quả
    I0 = (b * h^3) / 12;
    w_BET_q = (w_max_BET_q * E2 * I0) / (abs(q0) * L^4) * 100;
    w_TBT_q = (w_max_TBT_q * E2 * I0) / (abs(q0) * L^4) * 100;
    
    w_BET_P = (w_max_BET_P * E2 * I0) / (abs(q0) * L^3) * 100;
    w_TBT_P = (w_max_TBT_P * E2 * I0) / (abs(q0) * L^3) * 100;
end
