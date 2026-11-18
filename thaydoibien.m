function thaydoibien()
    % =========================================================================
    % PHẦN 1: TẠO GIAO DIỆN NGƯỜI DÙNG (UI)
    % =========================================================================
    fig = uifigure('Name', 'Mô phỏng Khóa Cắt - Dầm FGM Timoshenko', 'Position', [100, 100, 950, 600]);
    
    % Không gian vẽ đồ thị
    ax = uiaxes(fig, 'Position', [50, 50, 600, 500]);
    xlabel(ax, 'Chiều dài dầm x (m)', 'FontWeight', 'bold');
    ylabel(ax, 'Độ võng w (m)', 'FontWeight', 'bold');
    grid(ax, 'on');
    
    % Bảng điều khiển (Panel)
    pnl = uipanel(fig, 'Title', 'Thông số & Điều kiện biên', 'Position', [670, 20, 260, 560]);
    
    % Thanh trượt 1: Hệ số vật liệu n
    uilabel(pnl, 'Position', [10, 500, 200, 22], 'Text', 'Hệ số vật liệu n (0 - 100):');
    sld_n = uislider(pnl, 'Position', [20, 490, 220, 3]);
    sld_n.Limits = [0, 100]; sld_n.Value = 4.5; 
    
    % Thanh trượt 2: Tỷ số L/h (Độ mảnh của dầm)
    uilabel(pnl, 'Position', [10, 440, 200, 22], 'Text', 'Tỷ số L/h (Độ mảnh):');
    sld_Lh = uislider(pnl, 'Position', [20, 430, 220, 3]);
    sld_Lh.Limits = [5, 100]; sld_Lh.Value = 5; 
    
    % Thanh trượt 3: Số lượng phần tử FEM
    uilabel(pnl, 'Position', [10, 380, 200, 22], 'Text', 'Số phần tử FEM (Mesh):');
    sld_elem = uislider(pnl, 'Position', [20, 370, 220, 3]);
    sld_elem.Limits = [2, 100]; sld_elem.Value = 41; 
    
    % Dropdown: Loại tải trọng
    uilabel(pnl, 'Position', [10, 320, 200, 22], 'Text', 'Loại tải trọng:');
    dd_load = uidropdown(pnl, 'Position', [10, 295, 220, 22], ...
        'Items', {'Tải phân bố đều (uniform)', 'Tải tập trung giữa dầm (point)'}, ...
        'ItemsData', {'uniform', 'point'});
        
    % Dropdown: Điều kiện biên TRÁI
    uilabel(pnl, 'Position', [10, 255, 200, 22], 'Text', 'ĐK biên Trái:');
    dd_bc_left = uidropdown(pnl, 'Position', [10, 230, 220, 22], ...
        'Items', {'Ngàm (Clamped)', 'Khớp (Pinned)', 'Gối di động (Roller)', 'Tự do (Free)'}, ...
        'ItemsData', {'clamped', 'pinned', 'roller', 'free'}, 'Value', 'pinned');

    % Dropdown: Điều kiện biên PHẢI
    uilabel(pnl, 'Position', [10, 190, 200, 22], 'Text', 'ĐK biên Phải:');
    dd_bc_right = uidropdown(pnl, 'Position', [10, 165, 220, 22], ...
        'Items', {'Ngàm (Clamped)', 'Khớp (Pinned)', 'Gối di động (Roller)', 'Tự do (Free)'}, ...
        'ItemsData', {'clamped', 'pinned', 'roller', 'free'}, 'Value', 'roller');

    % Bảng hiển thị kết quả
    lbl_res_unlocked = uilabel(pnl, 'Position', [10, 105, 240, 40], ...
        'Text', 'w_star (Giảm): ...', 'FontWeight', 'bold', 'FontColor', 'blue');
    lbl_res_locked = uilabel(pnl, 'Position', [10, 65, 240, 40], ...
        'Text', 'w_star (Đầy đủ): ...', 'FontWeight', 'bold', 'FontColor', 'red');
        
    % Cảnh báo nếu chọn sai điều kiện biên (Hệ biến hình)
    lbl_warning = uilabel(pnl, 'Position', [10, 20, 240, 40], ...
        'Text', '', 'FontWeight', 'bold', 'FontColor', '#D95319');

    % Cài đặt sự kiện
    sld_n.ValueChangedFcn = @(src, event) updateSimulation();
    sld_Lh.ValueChangedFcn = @(src, event) updateSimulation();
    sld_elem.ValueChangedFcn = @(src, event) updateSimulation();
    dd_load.ValueChangedFcn = @(src, event) updateSimulation();
    dd_bc_left.ValueChangedFcn = @(src, event) updateSimulation();
    dd_bc_right.ValueChangedFcn = @(src, event) updateSimulation();
    
    % Chạy mô phỏng lần đầu tiên
    updateSimulation();
    
    % =========================================================================
    % PHẦN 2: HÀM CẬP NHẬT GIAO DIỆN
    % =========================================================================
    function updateSimulation()
        n = sld_n.Value;
        L_over_h = sld_Lh.Value;
        numElem = round(sld_elem.Value); 
        load_type = dd_load.Value;
        bc_left = dd_bc_left.Value;
        bc_right = dd_bc_right.Value;
        
        L = 1.0; b = 0.1; h = L / L_over_h;
        E1 = 380e9; E2 = 70e9; nu = 0.3; q0 = -10000; 
        
        % Truyền thêm 2 biến bc_left, bc_right vào hàm tính toán
        [x_nodes, w_unlocked, w_locked, star_unlocked, star_locked, err_flag] = ...
            solve_Compare_Locking(L, h, b, E1, E2, nu, n, numElem, q0, load_type, bc_left, bc_right);
        
        cla(ax);
        
        if err_flag
            % Báo lỗi biến hình học
            lbl_warning.Text = 'LỖI: Hệ bất biến hình! Hãy chọn lại Gối tựa.';
            lbl_res_unlocked.Text = 'w_star: NaN';
            lbl_res_locked.Text = 'w_star: NaN';
        else
            lbl_warning.Text = ''; % Tắt cảnh báo
            hold(ax, 'on');
            plot(ax, x_nodes, w_unlocked, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'Tích phân giảm (Chính xác)');
            plot(ax, x_nodes, w_locked, 'r--', 'LineWidth', 2, 'DisplayName', 'Tích phân đầy đủ (Bị khóa cắt)');
            hold(ax, 'off');
            
            title(ax, sprintf('Dầm FGM (n=%.1f, L/h=%.1f) - Trái: %s, Phải: %s', n, L_over_h, upper(bc_left), upper(bc_right)));
            legend(ax, 'Location', 'southwest', 'FontSize', 11);
            
            lbl_res_unlocked.Text = sprintf('w_star (Chính xác): %.4f', star_unlocked);
            lbl_res_locked.Text = sprintf('w_star (Bị khóa): %.4f', star_locked);
        end
    end
end

% =========================================================================
% PHẦN 3: CODE CORE FEM 
% =========================================================================
function [x_nodes, w_unlocked, w_locked, star_unlocked, star_locked, err_flag] = solve_Compare_Locking(L, h, b, E1, E2, nu, n, numElem, q0, load_type, bc_left, bc_right)
    numNodes = numElem + 1;
    ks = 5/6; 
    A = 0; B = 0; D = 0; S_val = 0;
    err_flag = false;
    
    % --- Tính ma trận độ cứng vật liệu ---
    if n <= 5
        xi_z = [-0.861136, -0.339981, 0.339981, 0.861136]; wi_z = [0.347855, 0.652145, 0.652145, 0.347855];
        for i = 1:4
            z = (h/2) * xi_z(i); J = h/2; Vc = (z/h + 0.5)^n;
            Ez = (E1 - E2) * Vc + E2; Gz = Ez / (2 * (1 + nu));
            A = A + wi_z(i) * (Ez) * J; B = B + wi_z(i) * (z * Ez) * J;
            D = D + wi_z(i) * (z^2 * Ez) * J; S_val = S_val + wi_z(i) * (ks * Gz) * J;
        end
        Cb = [A*b, B*b; B*b, D*b]; Cs = S_val * b;
    else
        E_diff = E1 - E2;
        A_a = b * h * (E_diff / (n + 1) + E2);
        B_a = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        D_a = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
        S_a = ks * A_a / (2 * (1 + nu));
        Cb = [A_a, B_a; B_a, D_a]; Cs = S_a;
    end
    
    % --- Khởi tạo hệ phương trình ---
    GDof = 3 * numNodes; 
    K_unlocked = zeros(GDof); 
    K_locked = zeros(GDof);   
    F = zeros(GDof, 1);
    le = L / numElem;
    
    % --- Lắp ráp ma trận phần tử ---
    for e = 1:numElem
        idx_e = (3*e-2):(3*e+3); 
        ke_b = zeros(6);       
        ke_s_full = zeros(6);  
        
        for k = [-0.57735, 0.57735]
            dN_dx = [-0.5, 0.5] * (2/le);
            Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 dN_dx(1) 0 0 dN_dx(2)];
            ke_b = ke_b + Bb' * Cb * Bb * (le/2);
            
            N1 = 0.5 * (1 - k); N2 = 0.5 * (1 + k);
            Bs_full = [0, -1/le, -N1, 0, 1/le, -N2];
            ke_s_full = ke_s_full + Bs_full' * Cs * Bs_full * (le/2);
        end
        
        Bs_reduced = [0, -1/le, -0.5, 0, 1/le, -0.5];       
        ke_s_reduced = Bs_reduced' * Cs * Bs_reduced * le; 
        
        fe = zeros(6,1);
        if strcmp(load_type, 'uniform')
            fe(2) = q0 * le / 2; fe(5) = q0 * le / 2; 
        end
        F(idx_e) = F(idx_e) + fe;
        
        K_unlocked(idx_e, idx_e) = K_unlocked(idx_e, idx_e) + ke_b + ke_s_reduced;
        K_locked(idx_e, idx_e)   = K_locked(idx_e, idx_e)   + ke_b + ke_s_full;
    end
    
    if strcmp(load_type, 'point')
        center_node = floor(numElem / 2) + 1; 
        F(3 * center_node - 1) = q0; 
    end
    
    % =====================================================================
    % THIẾT LẬP ĐIỀU KIỆN BIÊN THEO LỰA CHỌN CỦA NGƯỜI DÙNG
    % Bậc tự do mỗi Nút: [u, w, theta]
    % =====================================================================
    fixedNodes = [];
    
    % Xử lý Gối Trái (Nút số 1 -> Bậc tự do: 1, 2, 3)
    switch bc_left
        case 'clamped' % Ngàm (Khóa u, w, theta)
            fixedNodes = [fixedNodes, 1, 2, 3];
        case 'pinned'  % Khớp cố định (Khóa u, w)
            fixedNodes = [fixedNodes, 1, 2];
        case 'roller'  % Gối di động (Chỉ khóa w)
            fixedNodes = [fixedNodes, 2];
        case 'free'    % Tự do (Không khóa gì)
            % Bỏ trống
    end
    
    % Xử lý Gối Phải (Nút cuối cùng -> Bậc tự do: 3N-2, 3N-1, 3N)
    N = numNodes;
    switch bc_right
        case 'clamped' % Ngàm (Khóa u, w, theta)
            fixedNodes = [fixedNodes, 3*N-2, 3*N-1, 3*N];
        case 'pinned'  % Khớp cố định (Khóa u, w)
            fixedNodes = [fixedNodes, 3*N-2, 3*N-1];
        case 'roller'  % Gối di động (Chỉ khóa w)
            fixedNodes = [fixedNodes, 3*N-1];
        case 'free'    % Tự do (Không khóa gì)
            % Bỏ trống
    end
    
    activeDof = setdiff(1:GDof, fixedNodes);
    
    % --- Kiểm tra sự ổn định của hệ (Ngăn lỗi Ma trận suy biến) ---
    % Nếu người dùng chọn Tự Do - Tự Do, K sẽ bằng 0 gây lỗi
    if isempty(activeDof) || rcond(K_unlocked(activeDof, activeDof)) < 1e-12
        err_flag = true;
        x_nodes = linspace(0, L, numNodes); 
        w_unlocked = zeros(numNodes, 1); w_locked = zeros(numNodes, 1);
        star_unlocked = NaN; star_locked = NaN;
        return;
    end
    
    % --- Giải hệ ---
    U_unlocked = zeros(GDof, 1); U_locked = zeros(GDof, 1); 
    
    U_unlocked(activeDof) = K_unlocked(activeDof, activeDof) \ F(activeDof);
    U_locked(activeDof)   = K_locked(activeDof, activeDof) \ F(activeDof);
    
    % --- Trích xuất dữ liệu ---
    x_nodes = linspace(0, L, numNodes); 
    w_unlocked = U_unlocked(2:3:end);               
    w_locked   = U_locked(2:3:end);
    
    I0 = (b * h^3) / 12;
    if strcmp(load_type, 'uniform')
        star_unlocked = (max(abs(w_unlocked)) * E2 * I0) / (abs(q0) * L^4) * 100; 
        star_locked   = (max(abs(w_locked)) * E2 * I0) / (abs(q0) * L^4) * 100; 
    else
        star_unlocked = (max(abs(w_unlocked)) * E2 * I0) / (abs(q0) * L^3) * 100; 
        star_locked   = (max(abs(w_locked)) * E2 * I0) / (abs(q0) * L^3) * 100; 
    end
end