% code-btl-dstt.m
%% BÀI TẬP LỚN ĐẠI SỐ TUYẾN TÍNH
clc; clear; close all;

%% BƯỚC 1 (Mạnh): Ma trận đầu vào



                                   Pseudo code
// BƯỚC 1: Nhập và kiểm tra cấp ma trận (n)
LẶP VÔ HẠN:
    THỬ:
        IN "Nhập cấp ma trận vuông n: "
        NHẬP n
        NẾU n là số VÀ n là số đơn VÀ n >= 1 VÀ n là số nguyên THÌ
            THOÁT VÒNG LẶP
        NGƯỢC LẠI:
            IN "Lỗi: Kích thước phải là số nguyên dương! Nhập lại."
    BẮT LỖI:
        IN "Lỗi: Cú pháp không hợp lệ! Nhập lại."
    KẾT THÚC
KẾT THÚC LẶP

// BƯỚC 2: Khởi tạo ma trận A kích thước n x n
Khởi tạo ma trận A gồm n hàng và n cột (tất cả phần tử bằng 0)
IN "Nhập các phần tử từng hàng (cách nhau bằng dấu cách, ví dụ: 1 5 7):"

// BƯỚC 3: Nhập dữ liệu cho từng hàng của ma trận
VỚI i TỪ 1 ĐẾN n LẦN:
    LẶP VÔ HẠN:
        IN "Hàng i: "
        NHẬP chuỗi_hàng
        mảng_số = CHUYỂN_CHUỖI_THÀNH_MẢNG_SỐ(chuỗi_hàng)

        NẾU số_lượng_phần_tử(mảng_số) == n VÀ mảng_số LÀ VECTOR HÀNG THÌ
            GÁN A[i, :] = mảng_số
            THOÁT VÒNG LẶP
        NGƯỢC LẠI:
            IN "Lỗi: Cần nhập đúng n số trên hàng này! Hãy nhập lại."
    KẾT THÚC LẶP
KẾT THÚC FOR

// BƯỚC 4: Hiển thị ma trận A
IN "Ma trận đầu vào A ="
IN ma_trận A

                                    Code

format short g; % Tự động ẩn bớt các số 0 thừa (hiển thị 1 thay vì 1.0000)
% 1. Nhập cấp của ma trận
while true
    try
        n = input('Nhập cấp ma trận vuông n: ');
        if isnumeric(n) && isscalar(n) && n >= 1 && floor(n) == n
            break;
        else
            fprintf('Lỗi: Kích thước phải là số nguyên dương! Nhập lại.\n');
        end
    catch
        fprintf('Lỗi: Cú pháp không hợp lệ! Nhập lại.\n');
    end
end

% 2. Nhập từng hàng
A = zeros(n, n);
fprintf('\nNhập các phần tử từng hàng (cách nhau bằng dấu cách, ví dụ: 1 5 7):\n');

for i = 1:n
    while true
        row_str = input(sprintf('Hàng %d: ', i), 's');
        row_vec = str2num(row_str); 
        
        if numel(row_vec) == n && isrow(row_vec)
            A(i, :) = row_vec;
            break;
        else
            fprintf('  -> Lỗi: Cần nhập đúng %d số trên hàng này! Hãy nhập lại.\n', n);
        end
    end
end

fprintf('\nMa trận đầu vào A =\n');
disp(A);



%% BƯỚC 2 (Mạnh): Đưa về dạng Hessenberg (bằng phương pháp Householder)




                              Pseudocode
// BƯỚC 1: Khởi tạo bản sao ma trận
GÁN H = A  // Giữ nguyên ma trận A gốc

// BƯỚC 2: Thực hiện vòng lặp biến đổi Householder
VỚI k TỪ 1 ĐẾN n - 2 LÀM:
    m = n - k
    x = H[k+1:n, k]  // Trích xuất đoạn cột x từ hàng k+1 đến n tại cột k
    norm_x = chuẩn_Euclid(x)
    NẾU norm_x == 0 THÌ
        BỎ QUA VÀ TIẾP TỤC VÒNG LẶP  // Chuyển sang bước k tiếp theo
    KẾT THÚC NẾU
    NẾU x(1) khác 0 THÌ
        phase = x(1)/(module x(1))
    KẾT THÚC NẾU
    alpha = phase * norm_x
    u = x
    u[1] = u[1] + alpha
    v = u / chuẩn_Euclid(u)  // Chuẩn hóa u thành vectơ đơn vị v

    // Tạo ma trận Householder P_k
    P_prime = ma_trận_đơn_vị(m) - 2 * (v * chuyển_vị(v))
    P_k = ma_trận_đơn_vị(n)
    P_k[k+1:n, k+1:n] = P_prime  // Gán P_prime vào góc dưới bên phải

    // Biến đổi tương tự: H = P_k * H * P_k'
    H = P_k * H * P_k'
KẾT THÚC FOR

// BƯỚC 3: In ma trận kết quả
IN "----------------------"
IN "Ma trận Hessenberg của A là:"
IN "H ="
IN ma_trận H

                                      Code
                                      
H = A; % Tạo bản sao H để giữ nguyên ma trận A gốc
W = eye(n); % Tích lũy các phép biến đổi Householder 
for k = 1:(n - 2)
    m = n - k;
    x = H(k+1:n, k);
    
    norm_x = norm(x);
    if norm_x == 0
        continue;
    end
    
    if x(1) ~= 0
        phase = x(1) / abs(x(1));
    else
        phase = 1; 
    end
    alpha = phase * norm_x;
    
    u = x;
    u(1) = u(1) + alpha;
    v = u / norm(u);
    
    % Tạo ma trận Householder P_k
    P_prime = eye(m) - 2 * (v * v');
    P_k = eye(n);
    P_k(k+1:n, k+1:n) = P_prime;
    
    % Biến đổi tương tự: H = P_k * H * P_k'
    H = P_k * H * P_k';
    % Tích lũy W 
    W = W * P_k; 

end
fprintf('----------------------\n');
fprintf('\nMa trận Hessenberg của A là:\n');
fprintf('H =\n');
disp(H);



%% BƯỚC 3 + 4 (Lạc + Mạnh): Phân rã QR của Hk bằng phép quay Givens + Thuật toán QR lặp



                                    Pseudo code
BẮT ĐẦU
    % 1. Khởi tạo các tham số
    Gán tol = 1e-6
    Gán limit = 1000
    Gán H_k = H
    Gán status = FALSE

    % 2. Vòng lặp thuật toán QR
    VÒNG LẶP i TỪ 1 ĐẾN limit DO:
        R_k = H_k
        Q_k = I_n  (Ma trận đơn vị cấp n)

        % --- Khử Givens để phân tích H_k = Q_k * R_k ---
        VÒNG LẶP j TỪ 1 ĐẾN (n - 1) DO:
            a = R_k[j, j]
            b = R_k[j+1, j]

            NẾU |b| < 1e-12 THÌ:
                Bỏ qua bước này, chuyển sang j tiếp theo
            KẾT THÚC NẾU

            r_val = sqrt(|a|^2 + |b|^2)
            c = a / r_val
            s = b / r_val

            Khởi tạo G_j = I_n
            G_j[j, j]     = conj(c)
            G_j[j, j+1]   = conj(s)
            G_j[j+1, j]   = -s
            G_j[j+1, j+1] = c

            Cập nhật R_k = G_j * R_k
            Cập nhật Q_k = Q_k * G_j'
        KẾT THÚC VÒNG LẶP j

        R_k = Triu(R_k)  (Chỉ giữ lại phần tam giác trên)

        % --- Cập nhật ma trận cho bước lặp tiếp theo ---
        H_k = R_k * Q_k

        % --- Kiểm tra điều kiện hội tụ ---
        Lấy đường chéo phụ bên dưới: subdiag = diag(H_k, -1)
        
        NẾU max(|subdiag|) < tol THÌ:
            status = TRUE
            Thoát vòng lặp i (BREAK)
        KẾT THÚC NẾU
    KẾT THÚC VÒNG LẶP i

    % 3. Xuất kết quả
    NẾU status == TRUE THÌ:
        In "Thuật toán QR hội tụ sau i bước lặp."
    NGƯỢC LẠI:
        In "Cảnh báo: Thuật toán không hội tụ sau limit bước lặp."
    KẾT THÚC NẾU

    In ma trận H_k


                                  Code 

% 4. Thuật toán QR lặp (Viết trực tiếp phép quay Givens bên trong)
Q_tilde = eye(n); % Tích lũy ma trận Givens: Q_tilde = Q_0 * Q_1 * ... * Q_{k-1} 
tol = 1e-6;
limit = 1000;
H_k = H;
status = false;

for i = 1:limit
    % --- Khử Givens để phân tích H_k thành Q_k * R_k ---
    R_k = H_k;
    Q_k = eye(n);
    
    for j = 1:(n - 1)
        a = R_k(j, j);
        b = R_k(j+1, j);
        
        if abs(b) < 1e-12
            continue;
        end
        
        r_val = sqrt(abs(a)^2 + abs(b)^2);
        c = a / r_val;
        s = b / r_val;
        
        G_j = eye(n);
        G_j(j, j)     = conj(c);
        G_j(j, j+1)   = conj(s);
        G_j(j+1, j)   = -s;
        G_j(j+1, j+1) = c;
        
        R_k = G_j * R_k;
        Q_k = Q_k * G_j';
    end
    R_k = triu(R_k);
    
    % Cập nhật H_k = R_k * Q_k cho bước lặp tiếp theo
    H_k = R_k * Q_k;
    % Tích lũy ma trận Givens Q_tilde 
    Q_tilde = Q_tilde * Q_k; 

    
    % Kiểm tra điều kiện dừng
    subdiag = diag(H_k, -1);
    if max(abs(subdiag)) < tol
        status = true;
        break;
    end
end

fprintf('----------------------\n');
if status
    fprintf('Thuật toán QR hội tụ sau %d bước lặp.\n', i);
else
    fprintf('Cảnh báo: Thuật toán không hội tụ sau %d bước lặp.\n', limit);
end

fprintf('\nMa trận H_k cuối cùng:\n');
disp(H_k);




%% BƯỚC 5 (Lạc): Dạng Schur



                                              Pseudo code
GÁN T = H_k 
GÁN Q = W * Q_tilde

 IN_RA "Ma trận dạng Schur T ="
 XUẤT T 

IN_RA "Ma trận trực giao tổng thể Q ="
 XUẤT Q 

GÁN sai_số = chuẩn_Euclid(A - Q * T * chuyển_vị_liên_hợp(Q)) 
IN_RA "Sai số kiểm tra ||A - Q*T*Q*|| =", sai_số 

                                              Code

T = H_k;                 % Ma trận dạng Schur T
Q = W * Q_tilde;         % Ma trận trực giao tổng thể Q = W * Q_tilde

fprintf('----------------------\n');

fprintf('\nMa trận dạng Schur T =\n');
disp(T);

fprintf('\nMa trận trực giao tổng thể Q =\n');
disp(Q);

% Kiểm tra lại phân rã A = Q * T * Q' (Sai số ra 0 thì đúng )
err = norm(A - Q * T * Q');
fprintf('Sai số kiểm tra ||A - Q*T*Q''|| = %.2e\n', err);


%% BƯỚC 6 (Lạc): Trị riêng



                              Pseudo Code

TẠO mảng trị_riêng kích thước n

VÒNG_LẶP i TỪ 1 ĐẾN n:
    GÁN trị_riêng[i] = T[i, i]
KẾT_THÚC_VÒNG_LẶP

IN_RA "Các trị riêng lambda_i = T_ii là:"
XUẤT trị_riêng


                            Code


eigenvalues = diag(T);
fprintf('----------------------\n'); 
fprintf('\nCác trị riêng lambda_i = T_ii:\n'); 
for i = 1:n 
     fprintf(' lambda_%d = %.6f\n', i, eigenvalues(i)); 
end 


%% BƯỚC 7 (Lạc): Vector riêng


                            Pseudo Code

TẠO ma trận Y kích thước (n x n)  // Lưu các vectơ riêng y của T
TẠO ma trận V kích thước (n x n)  // Lưu các vectơ riêng v của A

VÒNG_LẶP i TỪ 1 ĐẾN n:
    GÁN lambda = trị_riêng[i]
    TẠO vectơ y kích thước (n x 1) với các giá trị ban đầu bằng 0
    
    % Chọn phần tử cơ sở y[i] = 1
    GÁN y[i] = 1
    
    % Giải thế ngược từ (i - 1) giảm về 1 cho hệ (T - lambda*I)*y = 0
    VÒNG_LẶP k TỪ (i - 1) GIẢM VỀ 1:
        GÁN tổng = 0
        VÒNG_LẶP j TỪ (k + 1) ĐẾN i:
            GÁN tổng = tổng + T[k, j] * y[j]
        KẾT_THÚC_VÒNG_LẶP
        
        GÁN mẫu_số = T[k, k] - lambda
        NẾU |mẫu_số| < 1e-12 THÌ
            GÁN mẫu_số = 1e-12  // Tránh lỗi chia cho 0 nếu có trị riêng bội
        KẾT_THÚC_NẾU
        
        GÁN y[k] = -tổng / mẫu_số
    KẾT_THÚC_VÒNG_LẶP
    
    % Chuẩn hóa vectơ y
    GÁN y = y / chuẩn_Euclid(y)
    GÁN Y[:, i] = y
    
    % Biến đổi về vectơ riêng của A: v = Q * y
    GÁN v = Q * y
    GÁN v = v / chuẩn_Euclid(v)
    GÁN V[:, i] = v
KẾT_THÚC_VÒNG_LẶP

IN_RA "Ma trận các vectơ riêng V (mỗi cột là một vectơ riêng v_i) ="
XUẤT V

 KIỂM TRA SAI SỐ (||A*v - lambda*v||)

VÒNG_LẶP i TỪ 1 ĐẾN n:
    GÁN sai_số_v = chuẩn_Euclid(A * V[:, i] - trị_riêng[i] * V[:, i])
    IN_RA "Sai số ||A*v_" + i + " - lambda_" + i + "*v_" + i + "|| =", sai_số_v
KẾT_THÚC_VÒNG_LẶP


                                Code

% 1. Giải hệ phương trình T*y = lambda*y <=> (T - lambda*I)*y = 0 để tìm y
% 2. Lấy v = Q * y để thu được vectơ riêng của ma trận A gốc

Y = zeros(n, n); % Lưu các vectơ riêng y của T (mỗi cột là 1 vectơ)
V = zeros(n, n); % Lưu các vectơ riêng v của A (mỗi cột là 1 vectơ)

for i = 1:n
    lambda = eigenvalues(i);
    y = zeros(n, 1);
    
    % Vì T là ma trận tam giác trên, ta chọn phần tử cơ sở y(i) = 1 và y(i+1:n) = 0
    y(i) = 1;
    
    % Thế ngược từ dưới lên để tìm y(1:i-1)
    for k = (i - 1):-1:1
        sum_val = T(k, k+1:i) * y(k+1:i);
        denom = T(k, k) - lambda;
        
        % Tránh lỗi chia cho 0 nếu xuất hiện trị riêng bội
        if abs(denom) < 1e-12
            denom = 1e-12;
        end
        
        y(k) = -sum_val / denom;
    end
    
    % Chuẩn hóa vectơ y
    y = y / norm(y);
    Y(:, i) = y;
    
    % Lấy v = Q * y
    v = Q * y;
    v = v / norm(v); % Chuẩn hóa vectơ riêng v
    V(:, i) = v;
end

fprintf('----------------------\n');
fprintf('\nCác vectơ riêng v_i tương ứng của ma trận A:\n');
for i = 1:n
    fprintf('Vectơ riêng v_%d (ứng với lambda_%d = %.6f):\n', i, i, eigenvalues(i));
    disp(V(:, i));
end

% Kiểm tra lại tính đúng đắn theo phương trình A*v = lambda*v
fprintf('----------------------\n');
fprintf('Kiểm tra sai số ||A*v - lambda*v||:\n');
for i = 1:n
    err_ev = norm(A * V(:, i) - eigenvalues(i) * V(:, i));
    fprintf('  ||A*v_%d - lambda_%d*v_%d|| = %.2e\n', i, i, i, err_ev);
end


