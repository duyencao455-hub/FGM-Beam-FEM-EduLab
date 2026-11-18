%% ========================================================================
% HÀM CON 2: GIẢI DẦM TIMOSHENKO (TBT)
% ========================================================================
function [w_max, w_star, U] = solve_Timoshenko_FGM(L, h, b, E1, E2, nu, n, ne, q0, load_type, BC)
    nn = ne + 1;
    ks = 5 * (1 + nu) / (6 + 5 * nu); % He so Cowper theo Reddy (2017); 

    A = 0; B = 0; D = 0; S_val = 0;
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

    GDof = 3 * nn; 
    K = zeros(GDof); 
    F = zeros(GDof, 1);
    le = L / ne;

    for e = 1:ne
        idx_e = (3*e-2):(3*e+3); 
        ke = zeros(6);
        
        for k = [-0.57735, 0.57735]
            dN_dx = [-0.5, 0.5] * (2/le);
            Bb = [dN_dx(1) 0 0 dN_dx(2) 0 0; 0 0 dN_dx(1) 0 0 dN_dx(2)];
            ke = ke + Bb' * Cb * Bb * (le/2);
        end
        Bs = [0 -1/le -0.5 0 1/le -0.5];       
        ke = ke + Bs' * Cs * Bs * le; 
        
        fe = zeros(6,1);
        
        % TỰ ĐỘNG CHỌN TẢI
        if load_type == 1 % Phân bố đều
            fe(2) = q0 * le / 2; 
            fe(5) = q0 * le / 2; 
        end
        
        K(idx_e, idx_e) = K(idx_e, idx_e) + ke; 
        F(idx_e) = F(idx_e) + fe;
    end
    
    % NẾU LÀ TẢI TẬP TRUNG P
    if load_type == 2
        center_node = floor(ne / 2) + 1;
        F(3 * center_node - 1) = q0;
    end
    
    % TỰ ĐỘNG CHỌN ĐIỀU KIỆN BIÊN
    if BC == 1     % Ngàm 2 đầu (C-C)
        fixedNodes = [1, 2, 3, GDof-2, GDof-1, GDof]; 
    elseif BC == 2 % Tựa đơn (S-S)
        fixedNodes = [1, 2, 3*nn-1]; 
    elseif BC == 3 % Ngàm 1 đầu (C-F)
        fixedNodes = [1, 2, 3]; 
    end

    activeDof = setdiff(1:GDof, fixedNodes);
    U = zeros(GDof, 1); 
    U(activeDof) = K(activeDof, activeDof) \ F(activeDof);
    
    w_max = max(abs(U(2:3:end)));
    I0 = (b * h^3) / 12;
    
    % Chuẩn hóa (L^4 cho tải q, L^3 cho tải P)
    if load_type == 1
        w_star = (w_max * E2 * I0) / (abs(q0) * L^4) * 100; 
    else
        w_star = (w_max * E2 * I0) / (abs(q0) * L^3) * 100; 
    end
end