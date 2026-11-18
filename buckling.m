% =========================================================================
% CHƯƠNG TRÌNH TỰ ĐỘNG TÍNH TOÁN, SO SÁNH VÀ XUẤT BẢNG SAI SỐ KÈM ĐỒ THỊ
% Bài toán: Lực tới hạn Buckling dầm FGM (Timoshenko vs EBT vs HOBT)
% Reference: Luan C. Trinh et al. (2016)
% =========================================================================
clear; clc; close all;
%% 1. KHAI BÁO THÔNG SỐ ĐẦU VÀO CHUNG
E_c = 380e9;      % Gốm (Pa)
E_m = 70e9;       % Kim loại (Pa)
nu  = 0.23;       % Hệ số Poisson
L = 1.0;          % Chiều dài dầm (m)
ratio = 5;        % Tỷ lệ L/h = 5
h = L / ratio;    % Chiều dày (m)
b = 1.0;          % Chiều rộng (m)
% Các chỉ số vật liệu cần khảo sát cho vòng lặp chính
p_vals = [0, 1, 2, 5, 10]; 
num_cases = length(p_vals);
% KẾT QUẢ THAM KHẢO TỪ BÀI BÁO (Trinh 2016 - FOBT)
% Giá trị P_bar cho dầm L/h=5
Ref_CC = [154.350, 80.498, 62.614, 50.384, 44.267]; % Ngàm-Ngàm
Ref_PP = [48.835 , 24.687, 19.245 , 16.024 , 14.427]; % Khớp-Khớp
% Tích phân Gauss cho vật liệu
xi_z = [-0.861136, -0.339981, 0.339981, 0.861136];
wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];

%% ========================================================================
%% 2. KHẢO SÁT SỰ HỘI TỤ CỦA LƯỚI (MESH CONVERGENCE STUDY)
%% ========================================================================
fprintf('====================================================================================\n');
fprintf('           PHẦN 1: KHẢO SÁT HỘI TỤ MẬT ĐỘ LƯỚI (ne)\n');
fprintf('====================================================================================\n');
p_test = 1; % Chọn p=1 làm mẫu khảo sát
ne_array = [2, 4, 10, 20, 30, 50, 100]; % Các mật độ lưới cần thử
P_conv_results = zeros(length(ne_array), 4); % Lưu kết quả cho 4 cột (TBT C-C, EBT C-C, TBT P-P, EBT P-P)

% Tính toán trước đặc trưng vật liệu cho trường hợp p_test = 1
A_c = 0; B_c = 0; D_c = 0; As_c = 0; ks = 5/6;
for i = 1:4
    z = (h/2) * xi_z(i); 
    wt = (h/2) * wi_z(i); 
    Vc = (z/h + 0.5)^p_test;              
    Ez = (E_c - E_m) * Vc + E_m;     
    Gz = Ez / (2 * (1 + nu));        
    A_c  = A_c  + Ez * wt * b;
    B_c  = B_c  + Ez * z * wt * b;
    D_c  = D_c  + Ez * z^2 * wt * b;
    As_c = As_c + ks * Gz * wt * b;
end
Cb_c = [A_c, B_c; B_c, D_c]; Cs_c = As_c;

fprintf('| Số phần tử (ne) |  TBT (C-C)  |  EBT (C-C)  |  TBT (P-P)  |  EBT (P-P)  |\n');
fprintf('------------------------------------------------------------------------------------\n');

for r = 1:length(ne_array)
    ne_t = ne_array(r);
    nn_t = ne_t + 1;      
    le_t = L / ne_t;      
    GDof_t = 3 * nn_t;
    
    for BC = 1:2
        if BC == 1
            fixed_t = [1, 2, 3, GDof_t-1, GDof_t]; % C-C
        else
            fixed_t = [1, 2, GDof_t-1];          % P-P
        end
        active_t = setdiff(1:GDof_t, fixed_t);
        
        K_TBT  = zeros(GDof_t, GDof_t); Kg_TBT = zeros(GDof_t, GDof_t);
        K_EBT  = zeros(GDof_t, GDof_t); Kg_EBT = zeros(GDof_t, GDof_t);
        
        xi_L = [-0.7745966692, 0, 0.7745966692]; wi_L = [5/9, 8/9, 5/9];
        
        for e = 1:ne_t
            dof_idx = (3*e-2):(3*e+3); 
            
            % --- TBT Element Stiffness ---
            ke_tbt = zeros(6, 6);
            for xi_x = [-1/sqrt(3), 1/sqrt(3)]
                dN_dx = [-0.5, 0.5] * (2/le_t);
                Bb_tbt = [dN_dx(1), 0, 0, dN_dx(2), 0, 0; 0, 0, dN_dx(1), 0, 0, dN_dx(2)];
                ke_tbt = ke_tbt + Bb_tbt' * Cb_c * Bb_tbt * (le_t/2); 
            end
            Bs_mat = [0, -1/le_t, -0.5, 0, 1/le_t, -0.5];       
            ke_tbt = ke_tbt + Bs_mat' * Cs_c * Bs_mat * le_t;
            Bg_tbt = [0, -1/le_t, 0, 0, 1/le_t, 0];
            keg_tbt = (Bg_tbt' * Bg_tbt) * le_t; 
            K_TBT(dof_idx, dof_idx)  = K_TBT(dof_idx, dof_idx) + ke_tbt; 
            Kg_TBT(dof_idx, dof_idx) = Kg_TBT(dof_idx, dof_idx) + keg_tbt;
            
            % --- EBT Element Stiffness ---
            ke_ebt = zeros(6,6); detJ = le_t / 2; 
            for i = 1:3
                xi = xi_L(i); wt = wi_L(i);
                Bm = zeros(1, 6); Bm(1, [1, 4]) = [-1/le_t, 1/le_t];    
                Bb = zeros(1, 6);
                Bb(1, 2) = - (6 / le_t^2) * xi; Bb(1, 3) = (1 / le_t) * (1 - 3*xi);       
                Bb(1, 5) = (6 / le_t^2) * xi;   Bb(1, 6) = - (1 / le_t) * (1 + 3*xi);       
                ke_ebt = ke_ebt + (A_c*(Bm'*Bm) + B_c*(Bm'*Bb + Bb'*Bm) + D_c*(Bb'*Bb)) * detJ * wt;
            end
            keg_ebt = zeros(6,6);
            kg_w = (1/(30*le_t)) * [ 36, 3*le_t, -36, 3*le_t; 3*le_t, 4*le_t^2, -3*le_t, -le_t^2;
                                    -36, -3*le_t, 36, -3*le_t; 3*le_t, -le_t^2, -3*le_t, 4*le_t^2 ];
            keg_ebt([2,3,5,6], [2,3,5,6]) = kg_w;
            K_EBT(dof_idx, dof_idx)  = K_EBT(dof_idx, dof_idx) + ke_ebt;
            Kg_EBT(dof_idx, dof_idx) = Kg_EBT(dof_idx, dof_idx) + keg_ebt;
        end
        
        [~, D_eig_TBT] = eig(K_TBT(active_t, active_t), Kg_TBT(active_t, active_t));
        eig_TBT = diag(D_eig_TBT);
        P_cr_TBT = min(eig_TBT(eig_TBT > 1e-5 & isreal(eig_TBT)));
        P_conv_results(r, (BC-1)*2 + 1) = (P_cr_TBT * 12 * L^2) / (E_m * b * h^3);
        
        [~, D_eig_EBT] = eig(K_EBT(active_t, active_t), Kg_EBT(active_t, active_t));
        eig_EBT = diag(D_eig_EBT);
        P_cr_EBT_FEM = min(eig_EBT(eig_EBT > 1e-5 & isreal(eig_EBT)));
        P_conv_results(r, (BC-1)*2 + 2) = (P_cr_EBT_FEM * 12 * L^2) / (E_m * b * h^3);
    end
    fprintf('|        %-8d | %-11.3f | %-11.3f | %-11.3f | %-11.3f |\n', ...
        ne_t, P_conv_results(r,1), P_conv_results(r,2), P_conv_results(r,3), P_conv_results(r,4));
end
fprintf('====================================================================================\n');
fprintf('=> KẾT LUẬN: Nghiệm hội tụ ổn định từ ne = 50. Chọn ne = 50 cho các tính toán sau.\n\n');

% Vẽ đồ thị hội tụ
figure('Name', 'Khao sat hoi tu luoi', 'NumberTitle', 'off', 'Position', [100, 200, 1100, 450]);
subplot(1, 2, 1);
plot(ne_array, P_conv_results(:,1), '-bs', ne_array, P_conv_results(:,2), '-r^', 'LineWidth', 1.5, 'MarkerSize', 6); grid on;
xlabel('Số lượng phần tử (ne)', 'FontWeight', 'bold'); ylabel('Lực tới hạn P_{bar}');
title('Hội tụ ĐK biên Ngàm - Ngàm (C-C)'); legend('TBT FEM', 'EBT FEM', 'Location', 'northeast');

subplot(1, 2, 2);
plot(ne_array, P_conv_results(:,3), '-bs', ne_array, P_conv_results(:,4), '-r^', 'LineWidth', 1.5, 'MarkerSize', 6); grid on;
xlabel('Số lượng phần tử (ne)', 'FontWeight', 'bold'); ylabel('Lực tới hạn P_{bar}');
title('Hội tụ ĐK biên Khớp - Khớp (P-P)'); legend('TBT FEM', 'EBT FEM', 'Location', 'northeast');

%% ========================================================================
%% 3. VÒNG LẶP CHÍNH XUẤT BẢNG SO SÁNH (SỬ DỤNG LƯỚI ĐÃ HỘI TỤ ne = 50)
%% ========================================================================
ne = 50;          % Chốt số phần tử đã khảo sát hội tụ
nn = ne + 1;      
le = L / ne;      
GDof = 3 * nn; 
fprintf('====================================================================================================================\n');
fprintf('           PHẦN 2: BẢNG SO SÁNH LỰC TỚI HẠN BUCKLING KHÔNG THỨ NGUYÊN (P_bar)\n');
fprintf('       Mô hình: Timoshenko (TBT), Euler-Bernoulli (EBT) vs HOBT (Trinh 2016)\n');
fprintf('====================================================================================================================\n');
% Khởi tạo Figure 2 cho đồ thị kết quả
figure('Name', 'Bieu do P_bar theo p', 'NumberTitle', 'off', 'Position', [150, 150, 1100, 500]);
for BC_type = 1:2
    if BC_type == 1
        fprintf('\n[1] ĐIỀU KIỆN BIÊN: NGÀM - NGÀM (C-C)\n');
        Ref_Data = Ref_CC;
        fixed = [1, 2, 3, GDof-1, GDof]; % Khóa u,w,th trái | w,th phải
    else
        fprintf('\n[2] ĐIỀU KIỆN BIÊN: KHỚP - KHỚP (P-P)\n');
        Ref_Data = Ref_PP;
        fixed = [1, 2, GDof-1];          % Khóa u,w trái | w phải
    end
    
    fprintf('--------------------------------------------------------------------------------------------------------------------\n');
    fprintf('| Chỉ số p |      FOBT        |   FEM (TBT)  | Sai số TBT(%%) | Giải Tích (EBT) |   FEM (EBT)  | Sai số EBT(%%) |\n');
    fprintf('--------------------------------------------------------------------------------------------------------------------\n');
    
    active = setdiff(1:GDof, fixed);
    
    P_FEM_TBT_plot = zeros(1, num_cases);
    P_FEM_EBT_plot = zeros(1, num_cases);
    P_Exact_EBT_plot = zeros(1, num_cases);
    
    for idx = 1:num_cases
        p = p_vals(idx);
        
        % --- TÍNH TÍCH PHÂN VẬT LIỆU ---
        A = 0; B = 0; D = 0; As = 0; ks = 5/6;
        for i = 1:4
            z = (h/2) * xi_z(i); 
            wt = (h/2) * wi_z(i); 
            Vc = (z/h + 0.5)^p;              
            Ez = (E_c - E_m) * Vc + E_m;     
            Gz = Ez / (2 * (1 + nu));        
            A  = A  + Ez * wt * b;
            B  = B  + Ez * z * wt * b;
            D  = D  + Ez * z^2 * wt * b;
            As = As + ks * Gz * wt * b;
        end
        Cb = [A, B; B, D]; 
        Cs = As;           
        
        % --- TẠO MA TRẬN ĐỘ CỨNG K VÀ Kg CHO TBT VÀ EBT ---
        K_TBT  = zeros(GDof, GDof); Kg_TBT = zeros(GDof, GDof);
        K_EBT  = zeros(GDof, GDof); Kg_EBT = zeros(GDof, GDof);
        
        xi_L = [-0.7745966692, 0, 0.7745966692]; wi_L = [5/9, 8/9, 5/9];
        
        for e = 1:ne
            dof_idx = (3*e-2):(3*e+3); 
            
            %% === PHẦN TIMOSHENKO (TBT) ===
            ke_tbt = zeros(6, 6);
            for xi_x = [-1/sqrt(3), 1/sqrt(3)]
                dN_dx = [-0.5, 0.5] * (2/le);
                Bb_tbt = [dN_dx(1), 0, 0, dN_dx(2), 0, 0; 0, 0, dN_dx(1), 0, 0, dN_dx(2)];
                ke_tbt = ke_tbt + Bb_tbt' * Cb * Bb_tbt * (le/2); 
            end
            Bs_mat = [0, -1/le, -0.5, 0, 1/le, -0.5];       
            ke_tbt = ke_tbt + Bs_mat' * Cs * Bs_mat * le; 
            Bg_tbt = [0, -1/le, 0, 0, 1/le, 0];
            keg_tbt = (Bg_tbt' * Bg_tbt) * le; 
            K_TBT(dof_idx, dof_idx)  = K_TBT(dof_idx, dof_idx) + ke_tbt; 
            Kg_TBT(dof_idx, dof_idx) = Kg_TBT(dof_idx, dof_idx) + keg_tbt;
            
            %% === PHẦN EULER-BERNOULLI (EBT) ===
            ke_ebt = zeros(6,6); detJ = le / 2; 
            for i = 1:3 
                xi = xi_L(i); wt = wi_L(i);
                Bm = zeros(1, 6); Bm(1, [1, 4]) = [-1/le, 1/le];    
                Bb = zeros(1, 6);
                Bb(1, 2) = - (6 / le^2) * xi; Bb(1, 3) = (1 / le) * (1 - 3*xi);       
                Bb(1, 5) = (6 / le^2) * xi;   Bb(1, 6) = - (1 / le) * (1 + 3*xi);       
                ke_ebt = ke_ebt + (A*(Bm'*Bm) + B*(Bm'*Bb + Bb'*Bm) + D*(Bb'*Bb)) * detJ * wt;
            end
            keg_ebt = zeros(6,6);
            kg_w = (1/(30*le)) * [ 36, 3*le, -36, 3*le; 3*le, 4*le^2, -3*le, -le^2;
                                  -36, -3*le, 36, -3*le; 3*le, -le^2, -3*le, 4*le^2 ];
            keg_ebt([2,3,5,6], [2,3,5,6]) = kg_w; 
            K_EBT(dof_idx, dof_idx)  = K_EBT(dof_idx, dof_idx) + ke_ebt;
            Kg_EBT(dof_idx, dof_idx) = Kg_EBT(dof_idx, dof_idx) + keg_ebt;
        end
        
        % Giải bài toán trị riêng TBT
        [~, D_eig_TBT] = eig(K_TBT(active, active), Kg_TBT(active, active));
        eigenvalues_TBT = diag(D_eig_TBT);
        P_cr_TBT = min(eigenvalues_TBT(eigenvalues_TBT > 1e-5 & isreal(eigenvalues_TBT)));
        P_bar_FEM_TBT = (P_cr_TBT * 12 * L^2) / (E_m * b * h^3);
        
        % Giải bài toán trị riêng EBT FEM
        [~, D_eig_EBT] = eig(K_EBT(active, active), Kg_EBT(active, active));
        eigenvalues_EBT = diag(D_eig_EBT);
        P_cr_EBT_FEM = min(eigenvalues_EBT(eigenvalues_EBT > 1e-5 & isreal(eigenvalues_EBT)));
        P_bar_FEM_EBT = (P_cr_EBT_FEM * 12 * L^2) / (E_m * b * h^3);
        
        % Giải tích EBT (Exact Analytical)
        D_eff = D - (B^2)/A; 
        if BC_type == 1 
            P_cr_EBT_Exact = (4 * pi^2 * D_eff) / (L^2);
        else            
            P_cr_EBT_Exact = (pi^2 * D_eff) / (L^2);
        end
        P_bar_Exact_EBT = (P_cr_EBT_Exact * 12 * L^2) / (E_m * b * h^3);
        
        % Đánh giá sai số
        P_bar_Ref = Ref_Data(idx);
        Error_pct_TBT = abs(P_bar_FEM_TBT - P_bar_Ref) / P_bar_Ref * 100;
        Error_pct_EBT = abs(P_bar_FEM_EBT - P_bar_Exact_EBT) / P_bar_Exact_EBT * 100;
        
        P_FEM_TBT_plot(idx) = P_bar_FEM_TBT;
        P_FEM_EBT_plot(idx) = P_bar_FEM_EBT;
        P_Exact_EBT_plot(idx) = P_bar_Exact_EBT;
        
        % In kết quả theo đúng các cột yêu cầu
        fprintf('| %-8d | %-16.3f | %-12.3f | %-13.2f | %-15.3f | %-12.3f | %-13.2f |\n', ...
                  p, P_bar_Ref, P_bar_FEM_TBT, Error_pct_TBT, P_bar_Exact_EBT, P_bar_FEM_EBT, Error_pct_EBT);
    end
    fprintf('--------------------------------------------------------------------------------------------------------------------\n');
    
    % Vẽ đồ thị kết quả chính
    subplot(1, 2, BC_type);
    plot(p_vals, Ref_Data, '-ko', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k'); hold on;
    plot(p_vals, P_FEM_TBT_plot, '-bs', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    plot(p_vals, P_FEM_EBT_plot, '-r^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    plot(p_vals, P_Exact_EBT_plot, '--m*', 'LineWidth', 1.2, 'MarkerSize', 6);
    grid on;
    xlabel('Chỉ số gradient vật liệu p', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Lực tới hạn không thứ nguyên P_{bar}', 'FontSize', 12);
    if BC_type == 1
        title('Ngàm - Ngàm (C-C)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        title('Khớp - Khớp (P-P)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    legend('HOBT (Trinh 2016)', 'FEM (TBT)', 'FEM (EBT)', 'Giải tích (EBT)', 'Location', 'northeast');
end