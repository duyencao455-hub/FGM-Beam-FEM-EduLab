%% ========================================================================
% BÀI TOÁN 2: KHẢO SÁT HIỆN TƯỢNG KHÓA CẮT THEO TỶ SỐ L/h
% Chứng minh tính hiệu quả của Tích phân giảm chọn lọc (SRI)
% ========================================================================
clc; clear; close all;

% --- 1. KHAI BÁO THÔNG SỐ CƠ BẢN ---
E1 = 14.4e9; E2 = 1.44e9; nu = 0.38; n_val = 5.0;
h = 88e-6; b = 2 * h; 
q0 = 1; 
N_fixed = 20; % Cố định 20 phần tử
BC=2;
% Dải tỷ số mảnh cần khảo sát
L_h_array = [4, 5, 10, 20, 50, 100]; 
w_Euler_ref = 0.5935; % Nghiệm Benchmark của dầm Euler (n=5)

% Khởi tạo mảng lưu kết quả để vẽ đồ thị
w_FI_arr = zeros(length(L_h_array), 1);
w_SRI_arr = zeros(length(L_h_array), 1);

% --- IN TIÊU ĐỀ BẢNG ---
fprintf('\n=======================================================================================\n');
fprintf(' Bảng 1.6: So sánh sai số của mô hình FI và SRI so với nghiệm lý thuyết Euler (BET)\n');
fprintf('=======================================================================================\n');
fprintf(' %-5s | %-10s | %-10s | %-10s | %-15s | %-15s\n', 'L/h', 'FI', 'SRI', 'BET', 'e_FI-BET (%)', 'e_SRI-BET (%)');
fprintf('-------|------------|------------|------------|-----------------|-----------------\n');

% --- 2. VÒNG LẶP TÍNH TOÁN VÀ XUẤT BẢNG ---
for i = 1:length(L_h_array)
    L_h = L_h_array(i);
    L_current = L_h * h;
    
    % a. Tính độ võng bằng thuật toán SRI (Tích phân giảm)
    [~, w_SRI] = solve_Timoshenko_FGM(L_current, h, b, E1, E2, nu, n_val, N_fixed, q0,1,BC);
    w_SRI_arr(i) = w_SRI;
    
    % b. Tính độ võng bằng thuật toán FI (Tích phân đầy đủ)
    [~, w_FI] = solve_Timoshenko_FGM_FI(L_current, h, b, E1, E2, nu, n_val, N_fixed, q0, 'uniform');
    w_FI_arr(i) = w_FI;
    
    % c. Tính toán sai số tương đối (%) so với nghiệm Euler (BET)
    err_FI  = abs(w_FI - w_Euler_ref) / w_Euler_ref * 100;
    err_SRI = abs(w_SRI - w_Euler_ref) / w_Euler_ref * 100;
    
    % d. Định dạng hiển thị: Chuyển dấu chấm thành dấu phẩy giống hình ảnh
    str_FI  = strrep(sprintf('%.4f', w_FI), '.', ',');
    str_SRI = strrep(sprintf('%.4f', w_SRI), '.', ',');
    str_BET = strrep(sprintf('%.4f', w_Euler_ref), '.', ',');
    str_eFI = strrep(sprintf('%.2f', err_FI), '.', ',');
    str_eSRI= strrep(sprintf('%.2f', err_SRI), '.', ',');
    
    % In từng dòng của bảng ra Command Window
    fprintf(' %-5d | %-10s | %-10s | %-10s | %-15s | %-15s\n', L_h, str_FI, str_SRI, str_BET, str_eFI, str_eSRI);
end
fprintf('=======================================================================================\n');
% --- 3. VẼ ĐỒ THỊ TRỰC QUAN HÓA (SHEAR LOCKING) - CHUẨN IN ĐEN TRẮNG ---
figure('Name', 'Shear Locking Phenomenon', 'Color', 'w', 'Position', [100, 100, 800, 550]);
hold on; grid on;

% Làm mờ lưới đồ thị để các nét vẽ đen nổi bật hơn
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12);

% 1. Đường FI (Bị Khóa cắt): Nét liền đen (-), Tam giác đen đặc (^k)
semilogx(L_h_array, w_FI_arr, '-k^', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'k'); 

% 2. Đường SRI (Hội tụ tốt): Nét liền đen (-), Tròn rỗng ruột (ok, nền trắng 'w')
semilogx(L_h_array, w_SRI_arr, '-ko', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'w');

% 3. Đường chuẩn Euler (Benchmark): Nét đứt màu đen (--k)
yline(w_Euler_ref, '--k', 'LineWidth', 1.5); 

% Định dạng trục và Nhãn
xlabel('Tỷ số độ mảnh L/h (Thang Logarit)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Độ võng chuẩn hóa w^* \times 10^2', 'FontSize', 12, 'FontWeight', 'bold');
title('Minh chứng Hiện tượng Khóa biến dạng cắt (Shear Locking)', 'FontSize', 14);

% Tùy chỉnh hiển thị trục X và Y
xticks(L_h_array);
xticklabels(string(L_h_array));
ylim([0, 0.75]); % Nới rộng trục Y một chút để chú thích không đè lên hình

% Chú thích (Legend)
legend({'Timoshenko (FI) - Bị Khóa cắt', 'Timoshenko (SRI) - Hội tụ tốt', 'Euler (BET) - Benchmark'}, ...
    'Location', 'southwest', 'FontSize', 11);
box on;

% Tự động xuất file ảnh độ phân giải cao (300 dpi) để chèn vào Word/LaTeX
print(gcf, 'shear_locking_bw.png', '-dpng', '-r300');

% HÀM 2: THUẬT TOÁN TÍCH PHÂN ĐẦY ĐỦ (FI) - CỐ TÌNH GÂY KHÓA CẮT
% ========================================================================
function [w_max, w_star] = solve_Timoshenko_FGM_FI(L, h, b, E1, E2, nu, n, numElem, q0, load_type)
    numNodes = numElem + 1;
    ks = 5/6; 
    A = 0; B = 0; D = 0; S_val = 0;
    
    % Vẫn dùng Hybrid Integration để tính vật liệu
    if n <= 5
        xi_z = [-0.861136, -0.339981, 0.339981, 0.861136];
        wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];
        for i = 1:4
            z = (h/2) * xi_z(i); 
            J = h/2; 
            Vc = (z/h + 0.5)^n;
            Ez = (E1 - E2) * Vc + E2;
            Gz = Ez / (2 * (1 + nu));
            A = A + wi_z(i) * (Ez) * J;
            B = B + wi_z(i) * (z * Ez) * J;
            D = D + wi_z(i) * (z^2 * Ez) * J;
            S_val = S_val + wi_z(i) * (ks * Gz) * J;
        end
        Cb = [A*b, B*b; B*b, D*b]; 
        Cs = S_val * b;
    else
        E_diff = E1 - E2;
        A_a = b * h * (E_diff / (n + 1) + E2);
        B_a = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        D_a = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
        S_a = ks * A_a / (2 * (1 + nu));
        Cb = [A_a, B_a; B_a, D_a]; 
        Cs = S_a;
    end

    GDof = 3 * numNodes; K = zeros(GDof); F = zeros(GDof, 1);
    le = L / numElem;

    for e = 1:numElem
        idx_e = (3*e-2):(3*e+3); 
        ke = zeros(6);
        
        % DÙNG 2 ĐIỂM GAUSS CHO CẢ UỐN (BEND) VÀ CẮT (SHEAR)
        for k = [-0.57735, 0.57735]
            
            % 1. Nội suy Uốn
            dN_dx = [-0.5, 0.5] * (2/le);
            Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 dN_dx(1) 0 0 dN_dx(2)];
            
            % 2. Nội suy Cắt (Dùng hàm dạng nội suy tại xi = k)
            N1 = 0.5 * (1 - k);
            N2 = 0.5 * (1 + k);
            Bs = [0, dN_dx(1), -N1, 0, dN_dx(2), -N2];
            
            % Cộng dồn cả Uốn và Cắt tại cùng điểm Gauss k
            ke = ke + (Bb' * Cb * Bb + Bs' * Cs * Bs) * (le/2);
        end
        
        fe = zeros(6,1);
        if strcmp(load_type, 'uniform')
            fe(2) = q0 * le / 2; fe(5) = q0 * le / 2; 
        end
        K(idx_e, idx_e) = K(idx_e, idx_e) + ke; 
        F(idx_e) = F(idx_e) + fe;
    end
    
    if strcmp(load_type, 'point')
        center_node = floor(numElem / 2) + 1;
        F(3 * center_node - 1) = q0;
    end

    fixedNodes = [1, 2, 3*numNodes-1]; 
    activeDof = setdiff(1:GDof, fixedNodes);
    U = zeros(GDof, 1); U(activeDof) = K(activeDof, activeDof) \ F(activeDof);
    w_max = max(abs(U(2:3:end)));
    
    I0 = (b * h^3) / 12;
    if strcmp(load_type, 'uniform')
        w_star = (w_max * E2 * I0) / (abs(q0) * L^4) * 100; 
    else
        w_star = (w_max * E2 * I0) / (abs(q0) * L^3) * 100; 
    end
end