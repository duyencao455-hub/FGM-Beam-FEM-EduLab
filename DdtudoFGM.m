% =========================================================================
% CHUONG TRINH PHAN TICH DAO DONG TU DO DAM FGM TIMOSHENKO (ALL-IN-ONE)
% DOI CHIEU DONG THOI: NGHIEM GIAI TICH NAVIER VA FEM FERREIRA (2020)
% =========================================================================
clear; clc; close all;

% --- 1. THONG SO VAT LIEU VA HINH HOC CHUAN FERREIRA (2020) ---
Ec = 14.4e9;       % Gom: 14.4 GPa
rho_c = 12.2e3;    % Gom: 12.2 x 10^3 kg/m3

Em = 1.44e9;       % Kim loai: 1.44 GPa
rho_m = 1.22e3;    % Kim loai: 1.22 x 10^3 kg/m3

nu = 0.38;         % He so Poisson
h = 17.6e-6;       % Chieu day h = 17.6 micro-meters
b = 2 * h;         % Be rong b = 2h
L = 20 * h;        % Chieu dai nhip L = 20h

% He so chuan hoa tan so theo pha kim loai (Ferreira & Reddy)
A_base = b * h;
I_base = (b * h^3) / 12;
Norm_Factor = L^2 * sqrt((rho_m * A_base) / (Em * I_base));


%% ========================================================================
% PHAN 1: KHAO SAT HOI TU LUOI CHO 3 MODES (N = 2 -> 100)
% ========================================================================
N_elements = [2, 4, 10, 20, 30, 50, 100];

% Nghiem tham chieu Navier tu sach Ferreira (2020) / Reddy (2017)
Ref_M1_n1  = 8.6700;   % Mode 1, n = 1
Ref_M2_n1  = 34.2900;  % Mode 2, n = 1
Ref_M3_n10 = 88.8000;  % Mode 3, n = 10

lam_M1 = zeros(length(N_elements), 1);
lam_M2 = zeros(length(N_elements), 1);
lam_M3 = zeros(length(N_elements), 1);

fprintf('\n===================================================================================================================\n');
fprintf(' PHAN 1: KHAO SAT HOI TU LUOI CHO CA 3 MODES DAO DONG (DOI CHIEU FERREIRA 2020)\n');
fprintf('===================================================================================================================\n');
fprintf('   N    | Mode 1 (n=1) | Sai so (%%) | Mode 2 (n=1) | Sai so (%%) | Mode 3 (n=10) | Sai so (%%)\n');
fprintf('--------|--------------|------------|--------------|------------|---------------|------------\n');

for i = 1:length(N_elements)
    Ne = N_elements(i);
    
    % --- Tinh n = 1 (Lay Mode 1 va Mode 2) ---
    w_n1 = solve_Vibration_Timo_FGM(L, h, b, Ec, Em, rho_c, rho_m, nu, 1, Ne);
    lam_M1(i) = w_n1(1) * Norm_Factor;
    lam_M2(i) = w_n1(2) * Norm_Factor;
    err_1 = abs(lam_M1(i) - Ref_M1_n1) / Ref_M1_n1 * 100;
    err_2 = abs(lam_M2(i) - Ref_M2_n1) / Ref_M2_n1 * 100;
    
    % --- Tinh n = 10 (Lay Mode 3) ---
    w_n10 = solve_Vibration_Timo_FGM(L, h, b, Ec, Em, rho_c, rho_m, nu, 10, Ne);
    lam_M3(i) = w_n10(3) * Norm_Factor; 
    err_3 = abs(lam_M3(i) - Ref_M3_n10) / Ref_M3_n10 * 100;
    
    fprintf('  %-4d  |   %8.4f   |   %6.3f   |   %8.4f   |   %6.3f   |   %8.4f    |   %6.3f\n', ...
            Ne, lam_M1(i), err_1, lam_M2(i), err_2, lam_M3(i), err_3);
end
fprintf('-------------------------------------------------------------------------------------------------------------------\n');
fprintf(' Tham chieu:   %8.4f   |     -      |   %8.4f   |     -      |   %8.4f    |     -    \n', Ref_M1_n1, Ref_M2_n1, Ref_M3_n10);
fprintf('===================================================================================================================\n\n');

% --- VE DO THI HOI TU 3 SUBPLOTS ---
figure('Name', 'Hoi tu 3 Modes', 'Color', 'w', 'Position', [50, 100, 1100, 350]);

subplot(1, 3, 1); hold on; grid on; set(gca, 'GridLineStyle', ':', 'FontSize', 10);
plot(N_elements, lam_M1, '-ks', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 6);
yline(Ref_M1_n1, '--k', 'LineWidth', 1.5);
xlabel('So phan tu (N)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Tan so chuan hoa \lambda', 'FontSize', 10, 'FontWeight', 'bold');
title('Mode 1 (n = 1)', 'FontSize', 11);
xticks([2, 10, 30, 50, 100]); xlim([0, 105]); legend({'Timoshenko', 'Benchmark'}, 'Location', 'NorthEast'); box on;

subplot(1, 3, 2); hold on; grid on; set(gca, 'GridLineStyle', ':', 'FontSize', 10);
plot(N_elements, lam_M2, '-k^', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'MarkerSize', 6);
yline(Ref_M2_n1, '--k', 'LineWidth', 1.5);
xlabel('So phan tu (N)', 'FontSize', 10, 'FontWeight', 'bold');
title('Mode 2 (n = 1)', 'FontSize', 11);
xticks([2, 10, 30, 50, 100]); xlim([0, 105]); legend({'Timoshenko', 'Benchmark'}, 'Location', 'NorthEast'); box on;

subplot(1, 3, 3); hold on; grid on; set(gca, 'GridLineStyle', ':', 'FontSize', 10);
plot(N_elements, lam_M3, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'MarkerSize', 6);
yline(Ref_M3_n10, '--k', 'LineWidth', 1.5);
xlabel('So phan tu (N)', 'FontSize', 10, 'FontWeight', 'bold');
title('Mode 3 (n = 10)', 'FontSize', 11);
xticks([2, 10, 30, 50, 100]); xlim([0, 105]); legend({'Timoshenko', 'Benchmark'}, 'Location', 'NorthEast'); box on;


%% ========================================================================
% PHAN 2: DOI CHIEU TOAN DIEN (BANG 3.9 TRONG LUAN VAN)
% ========================================================================
numElem = 100;
n_list = [0, 1, 10];

% Du lieu chuan tu Bang 15.2 (Ferreira 2020, trang 325):
Data_Ferreira_Navier = [
    9.83,  38.82, 85.63;   % n = 0
    8.67,  34.29, 75.79;   % n = 1
   10.28,  40.47, 88.80    % n = 10
];

Data_Ferreira_FEM40 = [
    9.8353, 38.9412, 86.2024;   % n = 0
    8.6730, 34.3826, 76.0537;   % n = 1
   10.2898, 40.5756, 88.2976    % n = 10
];

fprintf('========================================================================================================================\n');
fprintf(' PHAN 2: DOI CHIEU KET QUA VOI SACH CUA A. J. M. FERREIRA (2020) - BANG 15.2, TRANG 325\n');
fprintf('========================================================================================================================\n');
fprintf('  n  | Mode | Nghiem PTHH (N=100) | Giai tich Navier | Sai so vs Navier | FEM Ferreira (N=40) | Sai so Ferr vs Navier\n');
fprintf('-----|------|---------------------|------------------|------------------|---------------------|----------------------\n');

for i = 1:length(n_list)
    n_val = n_list(i);
    omega_n = solve_Vibration_Timo_FGM(L, h, b, Ec, Em, rho_c, rho_m, nu, n_val, numElem);
    lambda = omega_n * Norm_Factor;
    
    for m = 1:3
        my_fem  = lambda(m);
        nav_ref = Data_Ferreira_Navier(i, m);
        fem_ref = Data_Ferreira_FEM40(i, m);
        
        % Tinh sai so deu so voi nghiem Navier lam chuan tuyet doi:
        err_my_vs_nav   = abs(my_fem - nav_ref) / nav_ref * 100;
        err_ferr_vs_nav = abs(fem_ref - nav_ref) / nav_ref * 100; % <-- DA SUA CHUAN DONG NAY
        
        if m == 2
            fprintf(' %-3d |  %d   |      %8.4f       |     %8.2f     |     %6.3f %%     |       %8.4f      |       %6.3f %%\n', ...
                    n_val, m, my_fem, nav_ref, err_my_vs_nav, fem_ref, err_ferr_vs_nav);
        else
            fprintf('     |  %d   |      %8.4f       |     %8.2f     |     %6.3f %%     |       %8.4f      |       %6.3f %%\n', ...
                    m, my_fem, nav_ref, err_my_vs_nav, fem_ref, err_ferr_vs_nav);
        end
    end
    fprintf('-----|------|---------------------|------------------|------------------|---------------------|----------------------\n');
end


%% ========================================================================
% PHAN 3: HAM GIAI TIMOSHENKO FGM VOI SRI VA TICH PHAN LAI (HYBRID)
% ========================================================================
function omega_n = solve_Vibration_Timo_FGM(L, h, b, E1, E2, rho1, rho2, nu, n, numElem)
    numNodes = numElem + 1;
    ks = 5 * (1 + nu) / (6 + 5 * nu); % He so Cowper

    A = 0; B = 0; D = 0; S_val = 0;
    I0 = 0; I1 = 0; I2 = 0; 
    
    % THUAT TOAN TICH PHAN LAI (HYBRID INTEGRATION)
    if n <= 5
        xi_z = [-0.861136, -0.339981, 0.339981, 0.861136];
        wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];
        for i = 1:4
            z = (h/2) * xi_z(i); 
            J = h/2; wt = wi_z(i) * J * b; 
            
            Vc = (z/h + 0.5)^n;
            Ez = (E1 - E2) * Vc + E2;
            Gz = Ez / (2 * (1 + nu));
            rho_z = (rho1 - rho2) * Vc + rho2; 
            
            A = A + (Ez) * wt;
            B = B + (z * Ez) * wt;
            D = D + (z^2 * Ez) * wt;
            S_val = S_val + (ks * Gz) * wt;
            
            I0 = I0 + rho_z * wt;
            I1 = I1 + rho_z * z * wt;
            I2 = I2 + rho_z * z^2 * wt;
        end
        Cb = [A, B; B, D]; Cs = S_val;
    else
        % Cong thuc giai tich dang dong khi n > 5
        E_diff = E1 - E2;
        A_a = b * h * (E_diff / (n + 1) + E2);
        B_a = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        D_a = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
        S_a = ks * A_a / (2 * (1 + nu));
        Cb = [A_a, B_a; B_a, D_a]; Cs = S_a;
        
        rho_diff = rho1 - rho2;
        I0 = b * h * (rho_diff / (n + 1) + rho2);
        I1 = b * h^2 * rho_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        I2 = b * h^3 * (rho_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + rho2 / 12);
    end

    GDof = 3 * numNodes; 
    K = zeros(GDof); 
    M = zeros(GDof); 
    le = L / numElem;

    for e = 1:numElem
        idx_e = (3*e-2):(3*e+3); 
        ke = zeros(6);
        me = zeros(6); 
        
        % Tich phan 2 diem Gauss cho Ma tran uon & Ma tran khoi luong
        for k = [-0.577350269189626, 0.577350269189626]
            dN_dx = [-0.5, 0.5] * (2/le);
            Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 -dN_dx(1) 0 0 -dN_dx(2)];
            ke = ke + Bb' * Cb * Bb * (le/2);
            
            xi = k;
            N1 = 0.5*(1-xi); N2 = 0.5*(1+xi);
            Nu = [N1, 0, 0, N2, 0, 0];
            Nw = [0, N1, 0, 0, N2, 0];
            Nt = [0, 0, N1, 0, 0, N2];
            me = me + (I0*(Nu'*Nu + Nw'*Nw) - I1*(Nu'*Nt + Nt'*Nu) + I2*(Nt'*Nt)) * (le/2);
        end
        
        % Tich phan giam SRI (1 diem Gauss) cho Ma tran cat
        Bs = [0, -1/le, -0.5, 0, 1/le, -0.5];       
        ke = ke + Bs' * Cs * Bs * le; 
        
        K(idx_e, idx_e) = K(idx_e, idx_e) + ke; 
        M(idx_e, idx_e) = M(idx_e, idx_e) + me; 
    end

    % Dieu kien bien Tuong thich Tua don (Pin-Roller S-S)
    fixedNodes = [1, 2, 3*numNodes-1];
    activeDof = setdiff(1:GDof, fixedNodes);
    
    K_act = K(activeDof, activeDof);
    M_act = M(activeDof, activeDof);

    [~, D_mat] = eig(K_act, M_act);
    omega_rad = sqrt(sort(diag(D_mat))); 
    omega_n = omega_rad(1:3); 
end


%% ========================================================================
% PHAN 4: TRICH XUAT VEC-TO RIENG VA VE DANG DAO DONG (MODE SHAPES)
% ========================================================================
numElem = 100;
numNodes = numElem + 1;
x_coords = linspace(0, 1, numNodes); 

% Chay giai bai toan tri rieng cho truong hop n = 1
[~, Modes_w] = solve_Vibration_With_ModeShapes(L, h, b, Ec, Em, rho_c, rho_m, nu, 1, numElem);

% Ve do thi 3 Mode shapes
figure('Name', 'Dang dao dong Mode Shapes', 'Color', 'w', 'Position', [100, 100, 800, 450]);
hold on; grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11);

plot(x_coords, Modes_w(:,1), '-b', 'LineWidth', 2.2, 'DisplayName', 'Mode 1');
plot(x_coords, Modes_w(:,2), '--r', 'LineWidth', 2.2, 'DisplayName', 'Mode 2');
plot(x_coords, Modes_w(:,3), ':g', 'Color', [0 0.5 0], 'LineWidth', 2.2, 'DisplayName', 'Mode 3');

% Danh dau cac diem nut tinh (w = 0)
plot([0.5, 1/3, 2/3], [0, 0, 0], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
text(0.52, -0.08, 'Nut tinh (w = 0)', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Toa do chuan hoa doc truc dam (x/L)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Bien do dao dong chuan hoa', 'FontSize', 11, 'FontWeight', 'bold');
title('Hinh dang bien dang uon (Mode shapes) bac 1, 2 va 3 cua dam FGM (S-S)', 'FontSize', 12);
ylim([-1.2, 1.2]);
legend('Location', 'SouthOutside', 'Orientation', 'horizontal', 'FontSize', 11);


%% ========================================================================
% HAM PHU TRICH XUAT VEC-TO RIENG CUA TIMOSHENKO
% ========================================================================
function [omega_n, Modes_w] = solve_Vibration_With_ModeShapes(L, h, b, E1, E2, rho1, rho2, nu, n, numElem)
    numNodes = numElem + 1;
    ks = 5 * (1 + nu) / (6 + 5 * nu);
    
    A = 0; B = 0; D = 0; S_val = 0;
    I0 = 0; I1 = 0; I2 = 0; 
    
    if n <= 5
        xi_z = [-0.861136, -0.339981, 0.339981, 0.861136];
        wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];
        for i = 1:4
            z = (h/2) * xi_z(i); 
            J = h/2; wt = wi_z(i) * J * b; 
            
            Vc = (z/h + 0.5)^n;
            Ez = (E1 - E2) * Vc + E2;
            Gz = Ez / (2 * (1 + nu));
            rho_z = (rho1 - rho2) * Vc + rho2; 
            
            A = A + (Ez) * wt;
            B = B + (z * Ez) * wt;
            D = D + (z^2 * Ez) * wt;
            S_val = S_val + (ks * Gz) * wt;
            
            I0 = I0 + rho_z * wt;
            I1 = I1 + rho_z * z * wt;
            I2 = I2 + rho_z * z^2 * wt;
        end
        Cb = [A, B; B, D]; Cs = S_val;
    else
        E_diff = E1 - E2;
        A_a = b * h * (E_diff / (n + 1) + E2);
        B_a = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        D_a = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
        S_a = ks * A_a / (2 * (1 + nu));
        Cb = [A_a, B_a; B_a, D_a]; Cs = S_a;
        
        rho_diff = rho1 - rho2;
        I0 = b * h * (rho_diff / (n + 1) + rho2);
        I1 = b * h^2 * rho_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        I2 = b * h^3 * (rho_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + rho2 / 12);
    end

    GDof = 3 * numNodes; 
    K = zeros(GDof); M = zeros(GDof); 
    le = L / numElem;

    for e = 1:numElem
        idx_e = (3*e-2):(3*e+3); 
        ke = zeros(6); me = zeros(6); 
        for k = [-0.577350269189626, 0.577350269189626]
            dN_dx = [-0.5, 0.5] * (2/le);
            Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 -dN_dx(1) 0 0 -dN_dx(2)];
            ke = ke + Bb' * Cb * Bb * (le/2);
            
            xi = k;
            N1 = 0.5*(1-xi); N2 = 0.5*(1+xi);
            Nu = [N1, 0, 0, N2, 0, 0];
            Nw = [0, N1, 0, 0, N2, 0];
            Nt = [0, 0, N1, 0, 0, N2];
            me = me + (I0*(Nu'*Nu + Nw'*Nw) - I1*(Nu'*Nt + Nt'*Nu) + I2*(Nt'*Nt)) * (le/2);
        end
        Bs = [0, -1/le, -0.5, 0, 1/le, -0.5];       
        ke = ke + Bs' * Cs * Bs * le; 
        K(idx_e, idx_e) = K(idx_e, idx_e) + ke; 
        M(idx_e, idx_e) = M(idx_e, idx_e) + me; 
    end

    fixedNodes = [1, 2, 3*numNodes-1];
    activeDof = setdiff(1:GDof, fixedNodes);
    
    K_act = K(activeDof, activeDof);
    M_act = M(activeDof, activeDof);

    [V_act, D_mat] = eig(K_act, M_act);
    [omega_sq, sort_idx] = sort(diag(D_mat));
    omega_n = sqrt(omega_sq(1:3));
    
    Modes_w = zeros(numNodes, 3);
    for m = 1:3
        V_full = zeros(GDof, 1);
        V_full(activeDof) = V_act(:, sort_idx(m));
        w_nodes = V_full(2:3:end);
        
        [max_val, max_idx] = max(abs(w_nodes));
        if w_nodes(max_idx) < 0
            w_nodes = -w_nodes; 
        end
        Modes_w(:, m) = w_nodes / max_val;
    end
end