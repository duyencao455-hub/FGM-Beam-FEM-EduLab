%BÀI TOÁN 1
clc; clear; close all;

% =========================================================================
% 1. KHAI BÁO THÔNG SỐ BÀI TOÁN FGM
% =========================================================================
h = 88e-6;              % Chiều cao tiết diện (m)
L = 20 * h;             % Chiều dài dầm (L/h = 20)
b = 2 * h;              % Bề rộng tiết diện (m)
E1 = 14.4e9;            % Mô-đun pha Gốm (Pa)
E2 = 1.44e9;            % Mô-đun pha Kim loại (Pa)
nu = 0.38;              % Hệ số Poisson
n = 5.0;                % Chỉ số FGM
q0 = 1;                 % Tải trọng phân bố 1 N/m

% [ĐÃ THÊM] - Khai báo Loại tải trọng và Điều kiện biên để khớp với Benchmark
load_type = 1;          % 1: Tải phân bố đều (q)
BC = 2;                 % 2: Tựa đơn hai đầu (S-S)

% =========================================================================
% 2. THIẾT LẬP MẢNG KHẢO SÁT VÀ BENCHMARK
% =========================================================================
% [ĐÃ SỬA] - Thêm 100 vào mảng số phần tử
N_elements = [2, 5, 10, 20, 30, 40, 50, 100]; 
num_tests = length(N_elements);

% Mảng lưu kết quả Euler
w_max_euler = zeros(num_tests, 1);
w_star_euler = zeros(num_tests, 1);
err_euler = zeros(num_tests, 1);
benchmark_euler = 0.5935;

% Mảng lưu kết quả Timoshenko
w_max_timo = zeros(num_tests, 1);
w_star_timo = zeros(num_tests, 1);
err_timo = zeros(num_tests, 1);
benchmark_timo = 0.5968;

% =========================================================================
% 3. CHẠY VÒNG LẶP VÀ TÍNH TOÁN
% =========================================================================
disp('Đang tính toán hội tụ lưới cho cả hai mô hình...');
for i = 1:num_tests
    ne = N_elements(i);
    
    % --- Chạy Euler ---
    [w_max_E, w_star_E, ~] = solve_Euler_FGM(L, h, b, E1, E2, nu, n, ne, q0, load_type, BC);
    w_max_euler(i) = w_max_E;
    w_star_euler(i) = w_star_E;
    err_euler(i) = abs(w_star_E - benchmark_euler) / benchmark_euler * 100;
    
    % --- Chạy Timoshenko ---
    [w_max_T, w_star_T, ~] = solve_Timoshenko_FGM(L, h, b, E1, E2, nu, n, ne, q0, load_type, BC);
    w_max_timo(i) = w_max_T;
    w_star_timo(i) = w_star_T;
    err_timo(i) = abs(w_star_T - benchmark_timo) / benchmark_timo * 100;
end

% =========================================================================
% 4. IN BẢNG KẾT QUẢ
% =========================================================================
fprintf('\n=================================================================================\n');
fprintf('  BẢNG 1.4: ĐÁNH GIÁ SAI SỐ TUYỆT ĐỐI VÀ SỰ HỘI TỤ CỦA PHƯƠNG PHÁP PTHH\n');
fprintf('=================================================================================\n');
fprintf('  N    |  w*10^2 (Euler) | Sai số Euler | w*10^2 (Timo) | Sai số Timo  \n');
fprintf('---------------------------------------------------------------------------------\n');
for i = 1:length(N_elements)
    % Định dạng %3d đã đủ khoảng trống để hiển thị đẹp số 100
    fprintf(' %3d   |     %1.4f      |    %2.4f%%   |    %1.4f     |   %2.4f%% \n', ...
            N_elements(i), w_star_euler(i), err_euler(i), w_star_timo(i), err_timo(i));
end
fprintf('---------------------------------------------------------------------------------\n');
fprintf(' Bench.|     %1.4f      |       -      |    %1.4f     |     -       \n', ...
        benchmark_euler, benchmark_timo);
fprintf('=================================================================================\n\n');

% =========================================================================
% 5. VẼ ĐỒ THỊ HỘI TỤ (TỐI ƯU CHO IN ĐEN TRẮNG VÀ SỬA LỖI FONT)
% =========================================================================
figure('Name', 'Do thi hoi tu', 'NumberTitle', 'off', 'Position', [100, 100, 700, 450]);
hold on; grid on;
% Làm mờ lưới (Grid) để các đường đồ thị màu đen nổi bật hơn
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontSize', 12);

% --- VẼ ĐƯỜNG DỮ LIỆU ---
% 1. Timoshenko: Nét liền (-), Marker vuông (s), tô đen đặc
plot(N_elements, w_star_timo, '-ks', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 7);

% 2. Euler: Nét đứt (--), Marker tròn (o), rỗng ruột (màu trắng)
plot(N_elements, w_star_euler, '--ko', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 7);

% --- VẼ ĐƯỜNG BENCHMARK ---
% 3. Benchmark Timo: Nét chấm chấm (:) dày
yline(benchmark_timo, 'k:', 'LineWidth', 2);

% 4. Benchmark Euler: Nét gạch-chấm (-.) dày
yline(benchmark_euler, 'k-.', 'LineWidth', 2);

% --- ĐỊNH DẠNG TRỤC ---
xlabel('Số phần tử (N)', 'FontSize', 12, 'Interpreter', 'none');

% Đã đổi \bar{w} thành w^* và \cdot thành x để không bị lỗi với Tiếng Việt
ylabel('Độ võng chuẩn hóa w^* x 10^2', 'FontSize', 12, 'Interpreter', 'none');

% --- TÙY CHỈNH TRỤC X, Y ---
xticks(N_elements); 
xlim([0, 105]); 
ylim([min(min(w_star_timo, w_star_euler))*0.95, max(benchmark_timo, benchmark_euler)*1.05]);

% --- CHÚ THÍCH (LEGEND) ---
legend({'Timoshenko', 'Euler-Bernoulli', 'Benchmark Timo', 'Benchmark Euler'}, ...
       'Location', 'SouthEast', 'FontSize', 11);
box on;

% Lưu ảnh chất lượng cao (dpi 300) để in ấn sắc nét hơn
print(gcf, 'hoitu_luoi.png', '-dpng', '-r300');
% =========================================================================
% 6. VẼ HÌNH DẠNG DẦM TRƯỚC VÀ SAU BIẾN DẠNG (PHÓNG ĐẠI) - IN ĐEN TRẮNG
% =========================================================================
figure('Name', 'Hinh dang bien dang', 'NumberTitle', 'off', 'Position', [150, 150, 700, 300]);
hold on; grid on;

% Làm mờ lưới để đường đồ thị nổi bật hơn
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12);

% 1. Tạo trục tọa độ x chuẩn hóa (x/L) từ 0 đến 1
xi = linspace(0, 1, 100);

% 2. Phương trình đường cong uốn (dùng hàm sin tương tự như trong file LaTeX)
% Hàm sin xấp xỉ rất tốt đường đàn hồi của dầm tựa đơn và tạo đường cong mượt
w_deformed = -1 * sin(pi * xi);

% 3. Hình dạng ban đầu (đường thẳng y = 0)
w_initial = zeros(1, 100);

% 4. VẼ ĐỒ THỊ BẰNG NÉT ĐEN TRẮNG
% - Trước biến dạng: Nét đứt (--), màu xám đậm (gray), dày vừa
plot(xi, w_initial, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); 

% - Sau biến dạng: Nét liền (-), màu đen thui (black), thật dày
plot(xi, w_deformed, '-k', 'LineWidth', 2.0); 

% 5. Định dạng trục
xlabel('Tọa độ dọc trục x/L', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Độ võng w', 'FontSize', 12, 'Interpreter', 'none');
xlim([0, 1]);
ylim([-1.2, 0.2]);

% Đặt các mốc chia trên trục x, y cho gọn gàng giống mẫu
xticks(0:0.1:1);
yticks(-1:0.5:0);

% 6. Chú thích
legend({'Trước biến dạng', 'Sau biến dạng'}, 'Location', 'SouthEast', 'FontSize', 11);
box on;

% 7. Xuất ảnh độ phân giải cao 300dpi
print(gcf, 'hinhdang_biendang.png', '-dpng', '-r300');