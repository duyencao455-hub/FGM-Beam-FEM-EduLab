% =========================================================================
% CHƯƠNG TRÌNH KIỂM CHỨNG DAO ĐỘNG TỰ DO DẦM TIMOSHENKO (Bảng 10.2 - Ferreira)
% Đặc điểm: Phần tử 2 nút, Nội suy Lagrange tuyến tính, Tích phân giảm (Reduced)
% =========================================================================
clear; clc;

%% 1. THÔNG SỐ BÀI TOÁN (Theo chuẩn Ferreira: L=1, h=0.001, nu=0.3)
L   = 1.0;          % Chiều dài dầm
h   = 0.001;        % Chiều cao dầm (Dầm rất mảnh h/L = 0.001)
b   = 1.0;          % Chiều rộng dầm
nu  = 0.3;          % Hệ số Poisson
E   = 1e9;          % Giả thiết E (Sẽ triệt tiêu khi chuẩn hóa)
rho = 1000;         % Giả thiết rho (Sẽ triệt tiêu khi chuẩn hóa)
ks  = 5/6;          % Hệ số hiệu chỉnh cắt

% Các đặc trưng mặt cắt
A   = b * h;
I   = (b * h^3) / 12;
G   = E / (2 * (1 + nu));
Cb  = E * I;         % Độ cứng uốn
Cs  = ks * G * A;    % Độ cứng cắt

ne   = 50;           % Số phần tử (Theo đúng Ferreira dùng 50 phần tử)
nn   = ne + 1;
le   = L / ne;
GDof = 2 * nn;       % 2 Bậc tự do mỗi nút: [w, theta]

K = zeros(GDof, GDof);
M = zeros(GDof, GDof);

%% 2. LẬP MA TRẬN ĐỘ CỨNG (K) VÀ KHỐI LƯỢNG (M) TỪNG PHẦN TỬ
for e = 1:ne
    ke = zeros(4,4);
    me = zeros(4,4);
    
    % --- 2.1 TÍCH PHÂN GAUSS 2 ĐIỂM (Cho Uốn và Khối lượng) ---
    gauss_points = [-0.577350269189626, 0.577350269189626];
    gauss_weights = [1, 1];
    
    for i = 1:2
        xi = gauss_points(i);
        wt = gauss_weights(i);
        
        % Hàm dạng Lagrange tuyến tính (Linear Lagrange)
        N1 = 0.5 * (1 - xi);
        N2 = 0.5 * (1 + xi);
        dN_dxi = [-0.5, 0.5];
        dN_dx = dN_dxi * (2/le); % Chuyển đạo hàm sang tọa độ x thực
        
        % Ma trận đạo hàm góc xoay B_bending = [0, dN1/dx, 0, dN2/dx]
        Bb = [0, dN_dx(1), 0, dN_dx(2)];
        ke = ke + (Bb' * Cb * Bb) * (le/2) * wt;
        
        % Ma trận khối lượng (Consistent Mass): Quán tính tịnh tiến + Quán tính quay
        Nw = [N1, 0, N2, 0];  % Cho độ võng w
        Nt = [0, N1, 0, N2];  % Cho góc xoay theta
        me = me + (rho*A*(Nw'*Nw) + rho*I*(Nt'*Nt)) * (le/2) * wt;
    end
    
    % --- 2.2 TÍCH PHÂN GIẢM (REDUCED INTEGRATION - 1 ĐIỂM) Cho Cắt ---
    % Giúp tránh hiện tượng khóa cắt (Shear Locking) cho dầm mảnh
    xi_s = 0;
    wt_s = 2;
    N1_s = 0.5 * (1 - xi_s);
    N2_s = 0.5 * (1 + xi_s);
    dN_dx_s = [-0.5, 0.5] * (2/le);
    
    % Ma trận biến dạng cắt B_shear = [dN/dx, -N] (vì gamma = dw/dx - theta)
    Bs = [dN_dx_s(1), -N1_s, dN_dx_s(2), -N2_s];
    ke = ke + (Bs' * Cs * Bs) * (le/2) * wt_s;
    
    % --- 2.3 LẮP RÁP VÀO MA TRẬN TỔNG ---
    idx = (2*e-1):(2*e+2);
    K(idx, idx) = K(idx, idx) + ke;
    M(idx, idx) = M(idx, idx) + me;
end

%% 3. ĐIỀU KIỆN BIÊN: NGÀM (CANTILEVER)
% Nút 1 ngàm chặt -> Khóa w1 (DOF 1) và theta1 (DOF 2)
fixedDofs = [1, 2];
activeDofs = setdiff(1:GDof, fixedDofs);

K_act = K(activeDofs, activeDofs);
M_act = M(activeDofs, activeDofs);

%% 4. GIẢI BÀI TOÁN TRỊ RIÊNG VÀ CHUẨN HÓA KẾT QUẢ
[V, D] = eig(K_act, M_act);
omega_rad = sqrt(sort(diag(D))); % Tần số góc rad/s

% Công thức chuẩn hóa tần số không thứ nguyên (lambda^2) theo Ferreira
lambda_sq = omega_rad .* (L^2) .* sqrt((rho*A) / (E*I));

%% 5. IN KẾT QUẢ SO SÁNH (Check với Ferreira Table 10.2)
fprintf('\n=================================================================\n');
fprintf(' KIỂM CHỨNG DAO ĐỘNG TỰ DO: DẦM NGÀM TIMOSHENKO (Bảng 10.2)\n');
fprintf(' (L = 1, h = 0.001, nu = 0.3, Lưới = 50 phần tử)\n');
fprintf('=================================================================\n');
fprintf(' Mode |  FEM (Chuẩn hóa)  |  Sách Ferreira  |   Sai số (%%)\n');
fprintf('------|-------------------|-----------------|--------------------\n');

exact_ref = [3.5160, 22.0345]; % Nghiệm chuẩn trong sách

for i = 1:2
    fem_val = lambda_sq(i);
    ref_val = exact_ref(i);
    error = abs(fem_val - ref_val) / ref_val * 100;
    fprintf('  %d   |      %8.4f     |     %8.4f    |      %6.4f %%\n', i, fem_val, ref_val, error);
end
fprintf('-----------------------------------------------------------------\n');