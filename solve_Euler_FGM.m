%% ========================================================================
% HÀM CON 1: GIẢI DẦM EULER-BERNOULLI (BET)
% Đã thêm input: load_type (1: q, 2: P) và BC (1: CC, 2: SS, 3: CF)
% ========================================================================
function [w_max, w_star, U] = solve_Euler_FGM(L, h, b, E1, E2, ~, n, ne, q0, load_type, BC)
    nn = ne + 1;         
    le = L / ne;         
    
    % Dùng lại Hybrid cho n<=5 để khớp tuyệt đối với số liệu bảng
    A = 0; B = 0; D = 0;
    if n <= 5
        xi_z = [-0.861136, -0.339981, 0.339981, 0.861136];
        wi_z = [ 0.347855,  0.652145, 0.652145, 0.347855];
        for i = 1:4
            z = (h/2) * xi_z(i); 
            wt = (h/2) * wi_z(i); 
            Vc = (z/h + 0.5)^n;
            Ez = (E1 - E2) * Vc + E2;
            A = A + Ez * wt * b; 
            B = B + Ez * z * wt * b; 
            D = D + Ez * z^2 * wt * b;
        end
    else
        E_diff = E1 - E2;
        A = b * h * (E_diff / (n + 1) + E2);
        B = b * h^2 * E_diff * (1 / (n + 2) - 1 / (2 * (n + 1)));
        D = b * h^3 * (E_diff * (1 / (n + 3) - 1 / (n + 2) + 1 / (4 * (n + 1))) + E2 / 12);
    end

    GDof = 3 * nn; 
    K = zeros(GDof, GDof); 
    F = zeros(GDof, 1);
    
    xi_L = [-0.7745966692, 0, 0.7745966692]; 
    wi_L = [5/9, 8/9, 5/9];

    for e = 1:ne
        ke = zeros(6,6); 
        fe = zeros(6,1);
        detJ = le / 2; 
        
        for i = 1:3
            xi = xi_L(i); 
            wt = wi_L(i);
            
            Bm = zeros(1, 6); 
            Bm(1, [1, 4]) = [-1/le, 1/le];    
            
            Bb = zeros(1, 6);
            Bb(1, 2) = - (6 / le^2) * xi;             
            Bb(1, 3) =   (1 / le) * (1 - 3*xi);       
            Bb(1, 5) =   (6 / le^2) * xi;             
            Bb(1, 6) = - (1 / le) * (1 + 3*xi);       
            
            ke = ke + (A*(Bm'*Bm) + B*(Bm'*Bb + Bb'*Bm) + D*(Bb'*Bb)) * detJ * wt;
        end
        
        % TỰ ĐỘNG CHỌN TẢI
        if load_type == 1 % Tải phân bố đều
            fe([2,3,5,6]) = [q0*le/2; q0*le^2/12; q0*le/2; -q0*le^2/12];
        end
        
        idx = (3*e-2):(3*e+3);
        K(idx, idx) = K(idx, idx) + ke;
        F(idx) = F(idx) + fe;
    end
    
    % NẾU LÀ TẢI TẬP TRUNG P
    if load_type == 2
        center_node = floor(ne / 2) + 1;
        F(3 * center_node - 1) = q0;
    end
    
    % TỰ ĐỘNG CHỌN ĐIỀU KIỆN BIÊN
    if BC == 1     % Ngàm 2 đầu (C-C)
        fixedDofs = [1, 2, 3, GDof-2, GDof-1, GDof]; 
    elseif BC == 2 % Tựa đơn (S-S)
        fixedDofs = [1, 2, 3*nn-1]; 
    elseif BC == 3 % Ngàm 1 đầu (C-F)
        fixedDofs = [1, 2, 3]; 
    end

    activeDofs = setdiff(1:GDof, fixedDofs);
    U = zeros(GDof, 1); 
    U(activeDofs) = K(activeDofs, activeDofs) \ F(activeDofs);
    
    w_max = max(abs(U(2:3:end)));
    I0_ref = (b * h^3) / 12;
    
    % Chuẩn hóa (L^4 cho tải q, L^3 cho tải P)
    if load_type == 1
        w_star = (w_max * E2 * I0_ref) / (abs(q0) * L^4) * 100;
    else
        w_star = (w_max * E2 * I0_ref) / (abs(q0) * L^3) * 100;
    end
end