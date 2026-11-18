% ========================================================================
% BÀI TOÁN 3: KHẢO SÁT HỘI TỤ, THAY ĐỔI VẬT LIỆU (n) VÀ HÌNH DÁNG BIẾN DẠNG
% ĐIỀU KIỆN BIÊN: NGÀM 2 ĐẦU (C-C) ĐỂ KHỚP VỚI BENCHMARK CỦA REDDY
% ========================================================================
clc; clear; close all;

% --- 1. THÔNG SỐ ĐẦU VÀO ---
E1 = 200e9;         % Gốm (Pa)
E2 = 20e9;          % Kim loại (Pa)
nu = 0.25;          % Hệ số Poisson
L = 1.0;            % Chiều dài dầm (m)
b = 0.01;           % Chiều rộng (m)
h = 0.01;           % Chiều cao (m)
q0 = 2.5;           % Tải trọng phân bố đều (N/m)

load_type = 1;      % 1: Tải phân bố (q) | 2: Tải tập trung (P)
BC = 1;             % BẮT BUỘC = 1: Ngàm 2 đầu (C-C) để so sánh với Reddy

% Nghiệm giải tích chuẩn của Reddy (C-C)
% Reddy_BET/Reddy_TBT = [n=0, n=1]
Reddy_BET = [3.9062, 9.1104];
Reddy_TBT = [3.9109, 9.1508];


% ========================================================================
% PHẦN A: KHẢO SÁT HỘI TỤ LƯỚI (Với FGM n = 1)
% ========================================================================
n_conv = 1;
IDX_N1 = 2;   % SỬA: đặt tên biến rõ nghĩa thay vì hard-code "Reddy_BET(2)"
              % rải rác trong code — tránh nhầm lẫn nếu sau này đổi n_conv
N_elements = [2, 4, 8, 10, 20, 30, 40, 50, 100];
num_tests = length(N_elements);

w_conv_BET   = zeros(num_tests, 1);
w_conv_TBT   = zeros(num_tests, 1);
err_conv_BET = zeros(num_tests, 1);   % THÊM: lưu sai số để vẽ đồ thị (bản gốc chỉ in ra console)
err_conv_TBT = zeros(num_tests, 1);   % THÊM

fprintf('\n===============================================================================\n');
fprintf(' PHẦN A: KHẢO SÁT HỘI TỤ LƯỚI VỚI DẦM FGM (n = 1), NGÀM C-C\n');
fprintf('===============================================================================\n');
fprintf('  N    |  w*10^5 (Euler) | Sai số Euler | w*10^5 (Timo) | Sai số Timo  \n');
fprintf('-------|-----------------|--------------|---------------|--------------\n');

for i = 1:num_tests
    ne = N_elements(i);

    % Chạy Euler
    [w_max_E, ~, ~] = solve_Euler_FGM(L, h, b, E1, E2, nu, n_conv, ne, q0, load_type, BC);
    w_conv_BET(i)   = w_max_E * 1e5;
    err_conv_BET(i) = abs(w_conv_BET(i) - Reddy_BET(IDX_N1)) / Reddy_BET(IDX_N1) * 100;

    % Chạy Timoshenko
    [w_max_T, ~, ~] = solve_Timoshenko_FGM(L, h, b, E1, E2, nu, n_conv, ne, q0, load_type, BC);
    w_conv_TBT(i)   = w_max_T * 1e5;
    err_conv_TBT(i) = abs(w_conv_TBT(i) - Reddy_TBT(IDX_N1)) / Reddy_TBT(IDX_N1) * 100;

    fprintf(' %3d   |     %8.4f    |    %6.4f%%  |   %8.4f    |   %6.4f%% \n', ...
            ne, w_conv_BET(i), err_conv_BET(i), w_conv_TBT(i), err_conv_TBT(i));
end
fprintf('-------|-----------------|--------------|---------------|--------------\n');
fprintf(' Bench.|     %8.4f    |      -       |   %8.4f    |     -        \n', ...
        Reddy_BET(IDX_N1), Reddy_TBT(IDX_N1));
fprintf('===============================================================================\n\n');

% ------------------------------------------------------------------------
% HÌNH A1: ĐỘ VÕNG w THEO N — bản toàn cảnh (giữ nguyên tinh thần bản gốc)
% Mục đích: cho thấy trực quan mức sụp đổ ở lưới rất thô (N=2, Timo gần 0).
% LƯU Ý: vì trục w tuyến tính 0→10 nuốt mất chi tiết hội tụ ở N lớn (toàn
% bộ N>=8 dồn vào dải 8.5-9.2), hình này KHÔNG dùng để phân tích hội tụ chi
% tiết — chỉ dùng minh họa tổng quan. Phân tích chi tiết dùng Hình A2.
% ------------------------------------------------------------------------
figure('Name', 'Hoi tu luoi - BT3 (do vong, toan canh)', 'Color', 'w', 'Position', [100, 100, 750, 450]);
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12);

plot(N_elements, w_conv_TBT, '-ks', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 7);
plot(N_elements, w_conv_BET, '--ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'MarkerSize', 7);
yline(Reddy_TBT(IDX_N1), 'k:',  'LineWidth', 2);
yline(Reddy_BET(IDX_N1), 'k-.', 'LineWidth', 2);

xlabel('Số phần tử (N)', 'FontSize', 12);
ylabel('Độ võng tại tâm w \times 10^5 (m)', 'FontSize', 12);
xticks(N_elements); xlim([0, 105]); ylim([0, 10]);
legend({'Timoshenko', 'Euler-Bernoulli', 'Benchmark Timo', 'Benchmark Euler'}, ...
       'Location', 'SouthEast', 'FontSize', 11);
box on;
print(gcf, 'hoitu_bt3_dovong_toancanh_bw.png', '-dpng', '-r300');

% ------------------------------------------------------------------------
% HÌNH A2 (KHUYẾN NGHỊ DÙNG CHÍNH TRONG LUẬN VĂN): SAI SỐ (%) THEO N
% Thang log-log cho cả hai trục: xử lý được đồng thời sai số cực lớn ở N=2
% (~100%) và sai số cực nhỏ ở N=100 (~0.01-0.3%) trên cùng một hình, đồng
% thời làm lộ rõ đặc điểm hội tụ KHÔNG đơn điệu của Euler-Bernoulli (giảm
% dần rồi hơi nhích lên do sai số mô hình/phi tuyến hình học, xem thảo
% luận trong Mục 4.4.3), khác hẳn Timoshenko hội tụ đơn điệu liên tục.
% ------------------------------------------------------------------------
figure('Name', 'Hoi tu luoi - BT3 (sai so %)', 'Color', 'w', 'Position', [100, 100, 750, 450]);
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, ...
         'XScale', 'log', 'YScale', 'log');

loglog(N_elements, err_conv_TBT, '-ks',  'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 7);
loglog(N_elements, err_conv_BET, '--ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'MarkerSize', 7);

xlabel('Số phần tử (N)', 'FontSize', 12);
ylabel('Sai số so với nghiệm Reddy (%)', 'FontSize', 12);
xticks(N_elements); xlim([1.5, 120]);
legend({'Timoshenko', 'Euler-Bernoulli'}, 'Location', 'NorthEast', 'FontSize', 11);
box on;
print(gcf, 'hoitu_bt3_saiso_bw.png', '-dpng', '-r300');


% ========================================================================
% PHẦN B: SO SÁNH VẬT LIỆU TẠI LƯỚI ĐÃ HỘI TỤ VÀ VẼ HÌNH DÁNG BIẾN DẠNG
% (Không phát hiện lỗi ở phần này — giữ nguyên logic bản gốc)
% ========================================================================
N_elem = 100; % Lưới đã hội tụ từ Phần A
n_array = [0, 1];
w_plot_array = zeros(N_elem + 1, 1);

fprintf('=================================================================================================\n');
fprintf(' PHẦN B: BẢNG SO SÁNH ĐỘ VÕNG TẠI TÂM (w x 10^5 m) TẠI N = 100\n');
fprintf('=================================================================================================\n');
fprintf(' %-10s | %-28s | %-28s\n', 'Hệ số n', 'Euler (BET)', 'Timoshenko (TBT)');
fprintf(' %-10s | %-9s | %-8s | %-7s | %-9s | %-8s | %-7s\n', '', 'Tính toán', 'Reddy(C-C)', 'Sai số', 'Tính toán', 'Reddy(C-C)', 'Sai số');
fprintf('------------|-----------|----------|---------|-----------|----------|---------\n');

for i = 1:length(n_array)
    n = n_array(i);

    [w_max_BET, ~, ~] = solve_Euler_FGM(L, h, b, E1, E2, nu, n, N_elem, q0, load_type, BC);
    w_FEM_BET = w_max_BET * 1e5;

    [w_max_TBT, ~, U_TBT] = solve_Timoshenko_FGM(L, h, b, E1, E2, nu, n, N_elem, q0, load_type, BC);
    w_FEM_TBT = w_max_TBT * 1e5;

    % Trích xuất mảng độ võng của Timoshenko (n=1) để vẽ đồ thị
    if n == 1
        w_plot_array = U_TBT(2:3:end);
    end

    w_Exact_BET = Reddy_BET(i);
    w_Exact_TBT = Reddy_TBT(i);

    err_BET = abs(w_FEM_BET - w_Exact_BET) / w_Exact_BET * 100;
    err_TBT = abs(w_FEM_TBT - w_Exact_TBT) / w_Exact_TBT * 100;

    if n == 0
        name = '0 (Gốm)';
    else
        name = '1 (FGM)';
    end
    fprintf(' %-10s | %-9.4f | %-8.4f | %5.3f%% | %-9.4f | %-8.4f | %5.3f%%\n', ...
            name, w_FEM_BET, w_Exact_BET, err_BET, w_FEM_TBT, w_Exact_TBT, err_TBT);
end
fprintf('=================================================================================================\n\n');

% --- VẼ ĐỒ THỊ HÌNH DÁNG BIẾN DẠNG (ĐEN TRẮNG) ---
x_coords = linspace(0, 1, N_elem + 1);
w_normalized = -(w_plot_array / max(abs(w_plot_array)));

figure('Name', 'Hinh dang dam bien dang', 'Color', 'w', 'Position', [150, 150, 700, 350]);
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12);

plot(x_coords, zeros(size(x_coords)), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
plot(x_coords, w_normalized, '-k', 'LineWidth', 2.0);

xlabel('Tọa độ dọc trục x/L', 'FontSize', 12);
ylabel('Độ võng w (Phóng đại)', 'FontSize', 12);
xlim([0 1]);
ylim([-1.2 0.2]);
yticks([-1, -0.5, 0]);

legend({'Trước biến dạng', 'Sau biến dạng'}, 'Location', 'SouthEast', 'FontSize', 11);
box on;
print(gcf, 'hinhdang_biendang_bw.png', '-dpng', '-r300');