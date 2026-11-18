% =========================================================================
% BÀI TOÁN 5: KHẢO SÁT THAY ĐỔI TẢI TRỌNG (PHÂN BỐ ĐỀU q & TẬP TRUNG P)
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
load_val = 1;           % Giá trị tải
BC = 2;                 % Điều kiện biên: 2 = Tựa đơn (S-S)

% Mảng hệ số vật liệu n
n_array = [0; 1; 5; 10; 100]; 

% -------------------------------------------------------------------------
% NGHIỆM THAM CHIẾU (Được trích xuất trực tiếp từ Bảng 1.10 trong ảnh)
% -------------------------------------------------------------------------
w_Ref_BET_q = [0.1302; 0.3047; 0.5935; 0.6526; 1.0534];
w_Ref_TBT_q = [0.1310; 0.3062; 0.5968; 0.6571; 1.0610];
w_Ref_BET_P = [0.2083; 0.4876; 0.9496; 1.0442; 1.6855];
w_Ref_TBT_P = [0.2100; 0.4906; 0.9562; 1.0532; 1.7006];

% Khởi tạo mảng lưu nghiệm PTHH sẽ tính
w_FEM_BET_q = zeros(5, 1); w_FEM_TBT_q = zeros(5, 1);
w_FEM_BET_P = zeros(5, 1); w_FEM_TBT_P = zeros(5, 1);

% =========================================================================
% 2. VÒNG LẶP TÍNH TOÁN (GỌI HÀM PTHH)
% =========================================================================
disp('Đang tính toán Nghiệm Phần tử hữu hạn, vui lòng đợi...');
for i = 1:length(n_array)
    n = n_array(i);
    
    % --- NGHIỆM PTHH: TẢI PHÂN BỐ (load_type = 1) ---
    [~, w_FEM_BET_q(i), ~] = solve_Euler_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 1, BC);
    [~, w_FEM_TBT_q(i), ~] = solve_Timoshenko_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 1, BC);
    
    % --- NGHIỆM PTHH: TẢI TẬP TRUNG (load_type = 2) ---
    [~, w_FEM_BET_P(i), ~] = solve_Euler_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 2, BC);
    [~, w_FEM_TBT_P(i), ~] = solve_Timoshenko_FGM(L, h, b, Ec, Em, nu, n, N_elem, load_val, 2, BC);
end
disp('Hoàn tất! Đang in kết quả và vẽ đồ thị...');

% Hàm định dạng số hiển thị dấu phẩy (chuẩn tiếng Việt)
fmt = @(num) strrep(sprintf('%.4f', num), '.', ',');

% =========================================================================
% 3. IN BẢNG KẾT QUẢ RA COMMAND WINDOW
% =========================================================================
fprintf('\nBảng: Độ võng chuẩn hóa dầm vi mô FGM (L/h = 20, thuật toán Hybrid, N = 100 phần tử)\n');
fprintf('--------------------------------------------------------------------------------\n');
fprintf(' %-3s | %-16s | %-10s | %-13s | %-12s | %-10s\n', 'n', 'Tải trọng', 'Lý thuyết', 'Nghiệm PTHH', 'Tham chiếu', 'Sai số (%)');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:length(n_array)
    % Khối tải phân bố q
    err_BET_q = abs(w_FEM_BET_q(i) - w_Ref_BET_q(i)) / w_Ref_BET_q(i) * 100;
    err_TBT_q = abs(w_FEM_TBT_q(i) - w_Ref_TBT_q(i)) / w_Ref_TBT_q(i) * 100;
    fprintf(' %-3d | Phân bố đều (q)  | %-10s | %-13s | %-12s | %-10s\n', n_array(i), 'BET', fmt(w_FEM_BET_q(i)), fmt(w_Ref_BET_q(i)), strrep(sprintf('%.2f', err_BET_q), '.', ','));
    fprintf('     |                  | %-10s | %-13s | %-12s | %-10s\n', 'TBT', fmt(w_FEM_TBT_q(i)), fmt(w_Ref_TBT_q(i)), strrep(sprintf('%.2f', err_TBT_q), '.', ','));
    
    % Khối tải tập trung P
    err_BET_P = abs(w_FEM_BET_P(i) - w_Ref_BET_P(i)) / w_Ref_BET_P(i) * 100;
    err_TBT_P = abs(w_FEM_TBT_P(i) - w_Ref_TBT_P(i)) / w_Ref_TBT_P(i) * 100;
    fprintf('     |------------------|------------|---------------|--------------|-----------\n');
    fprintf('     | Tập trung (P)    | %-10s | %-13s | %-12s | %-10s\n', 'BET', fmt(w_FEM_BET_P(i)), fmt(w_Ref_BET_P(i)), strrep(sprintf('%.2f', err_BET_P), '.', ','));
    fprintf('     |                  | %-10s | %-13s | %-12s | %-10s\n', 'TBT', fmt(w_FEM_TBT_P(i)), fmt(w_Ref_TBT_P(i)), strrep(sprintf('%.2f', err_TBT_P), '.', ','));
    fprintf('--------------------------------------------------------------------------------\n');
end

% =========================================================================
% 4. TRỰC QUAN HOÁ KẾT QUẢ (VẼ ĐỒ THỊ CHUẨN ĐEN TRẮNG)
% =========================================================================

% --- HÌNH 1.12: BIỂU ĐỒ CỘT (BAR CHART) ĐEN TRẮNG ---
figure('Name', 'Hình 1.12: So sánh độ võng chuẩn hóa', 'Color', 'w', 'Position', [100, 100, 800, 500]);
data_bar = [w_FEM_BET_q, w_FEM_TBT_q, w_FEM_BET_P, w_FEM_TBT_P];
h_bar = bar(data_bar, 'grouped');

% Áp dụng dải màu Grayscale kết hợp với viền đen (EdgeColor)
h_bar(1).FaceColor = [0.2 0.2 0.2]; % Xám đậm gần đen (BET q)
h_bar(2).FaceColor = [0.5 0.5 0.5]; % Xám vừa (TBT q)
h_bar(3).FaceColor = [0.8 0.8 0.8]; % Xám nhạt (BET P)
h_bar(4).FaceColor = [1.0 1.0 1.0]; % Trắng (TBT P)

for i = 1:4
    h_bar(i).EdgeColor = 'k';       % Viền màu đen cho tất cả
    h_bar(i).LineWidth = 1.0;       % Độ dày viền
end

% Hiển thị giá trị text trên đầu cột (Màu đen)
for i = 1:length(h_bar)
    xtips = h_bar(i).XEndPoints; ytips = h_bar(i).YEndPoints;
    labels = string(round(h_bar(i).YData, 2));
    % Cộng thêm một chút vào ytips để chữ không chạm dính vào đỉnh cột
    text(xtips, ytips + 0.05, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'Rotation', 90, 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
end

set(gca, 'XTickLabel', {'n = 0', 'n = 1','n = 5', 'n = 10', 'n = 100'}, 'FontSize', 11);
ylabel('Độ võng chuẩn hóa w^*', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Chỉ số vật liệu n', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 2.2]);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
legend('BET - Phân bố (q)', 'TBT - Phân bố (q)', 'BET - Tập trung (P)', 'TBT - Tập trung (P)', ...
    'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 11);
box on;
print(gcf, 'barchart_bt5_bw.png', '-dpng', '-r300'); % Xuất ảnh

% --- HÌNH 1.13: BIỂU ĐỒ ĐƯỜNG (LINE PLOT) ĐEN TRẮNG ---
figure('Name', 'Hình 1.13: Phần trăm chênh lệch', 'Color', 'w', 'Position', [150, 150, 600, 400]);
diff_q = (w_FEM_TBT_q - w_FEM_BET_q) ./ w_FEM_BET_q * 100;
diff_P = (w_FEM_TBT_P - w_FEM_BET_P) ./ w_FEM_BET_P * 100;

x_idx = 1:5; 
% Đường tải phân bố: Nét liền, Vuông đen đặc
plot(x_idx, diff_q, '-ks', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
% Đường tải tập trung: Nét đứt, Tròn trắng rỗng
plot(x_idx, diff_P, '--ko', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');

set(gca, 'XTick', 1:5, 'XTickLabel', {'0', '1','5', '10', '100'}, 'FontSize', 11);
ylabel('Chênh lệch TBT so với BET (%)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Chỉ số vật liệu n', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 1.2]);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
legend('Tải phân bố (q)', 'Tải tập trung (P)', 'Location', 'northwest', 'FontSize', 11);
box on;
print(gcf, 'lineplot_bt5_bw.png', '-dpng', '-r300'); % Xuất ảnh
