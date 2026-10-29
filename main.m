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
end
fprintf('----------------------\n');
fprintf('\nMa trận Hessenberg của A là:\n');
fprintf('H =\n');
disp(H);



%% BƯỚC 3 (Lạc): Phân rã QR của Hk bằng phép quay Givens



                                Pseudo code

// BƯỚC 1: Khởi tạo ma trận R và Q
GÁN R = H  // Tạo bản sao từ ma trận Hessenberg H
GÁN Q = ma_trận_đơn_vị(n)  // Khởi tạo Q là ma trận đơn vị cấp n


// BƯỚC 2: Thực hiện vòng lặp quay Givens (tương thích số thực & số phức)
VỚI i TỪ 1 ĐẾN n - 1 LÀM:
    a = R[i, i]      // Phần tử trên đường chéo chính
    b = R[i+1, i]    // Phần tử ngay dưới đường chéo chính
    
    NẾU mô_đun(b) == 0 THÌ
        BỎ QUA VÀ TIẾP TỤC VÒNG LẶP  // Đã bằng 0, chuyển sang bước i tiếp theo
    KẾT THÚC NẾU
    
    // Tính độ dài bằng mô-đun tuyệt đối (chuẩn hóa cho cả số phức)
    r = căn_bậc_hai(mô_đun(a)^2 + mô_đun(b)^2)
    c = a / r  // Hệ số c
    s = b / r  // Hệ số s
    
    // Tạo ma trận quay Givens Unitary G_i cấp n x n
    G_i = ma_trận_đơn_vị(n)
    G_i[i, i]     =  liên_hợp_phức(c)
    G_i[i, i+1]   =  liên_hợp_phức(s)
    G_i[i+1, i]   = -s
    G_i[i+1, i+1] =  c
    
    // Cập nhật ma trận R và Q
    R = G_i * R                              // Triệt tiêu phần tử R[i+1, i] về 0
    Q = Q * chuyển_vị_liên_hợp(G_i)         // Nhân dồn ma trận Hermit (G_i^H) vào Q
KẾT THÚC FOR


// BƯỚC 3: Triệt tiêu sai số làm tròn
R = lấy_tam_giác_trên(R)  // Xóa sai số số học bên dưới đường chéo chính


// BƯỚC 4: In kết quả phân rã QR
IN "----------------------"
IN "Ma trận Unitary/Trực giao Q:"
IN ma_trận Q
IN "Ma trận tam giác trên R:"
IN ma_trận R

                                        Code

% Phân rã QR của ma trận Hessenberg H bằng phép quay Givens

R = H;          % Khởi tạo R ban đầu bằng H
Q = eye(n);     % Khởi tạo Q là ma trận đơn vị cấp n

for i = 1:(n - 1)
    % Lấy 2 phần tử trên đường chéo chính và ngay dưới đường chéo chính
    a = R(i, i);
    b = R(i+1, i);
    
    % Nếu phần tử dưới đường chéo chính đã bằng 0 thì không cần khử
    if abs(b) == 0
        continue;
    end
    
    % Tính ma trận quay Givens (cos = c, sin = s)
    r_val = sqrt(abs(a)^2 + abs(b)^2);
    c = a / r_val;
    s = b / r_va;
    
    % Tạo ma trận quay Givens G_i cấp n x n
    G_i = eye(n);
    G_i(i, i)     =conj(c);
    G_i(i, i+1)   = conj(s);
    G_i(i+1, i)   = -s;
    G_i(i+1, i+1) = c;
    
    % Cập nhật R = G_i * R (triệt tiêu phần tử R(i+1, i) về 0)
    R = G_i * R;
    
    % Cập nhật Q = Q * G_i' (vì Q = G1' * G2' * ... * G_{n-1}')
    Q = Q * G_i';
end

% Triệt tiêu các sai số làm tròn nhỏ (gần bằng 0) bên dưới đường chéo chính của R
R = triu(R);

% -------------------------------------------------------------------------
% Hiển thị kết quả
% -------------------------------------------------------------------------
fprintf('----------------------\n');
fprintf('\nMa trận trực giao Q =\n');
disp(Q);

fprintf('Ma trận tam giác trên R =\n');
disp(R);

% Kiểm tra tính đúng đắn: Tích Q * R phải bằng H
fprintf('Kiểm tra tích Q * R (so sánh với H) =\n');
disp(Q * R);



%% BƯỚC 4 (Mạnh): Thuật toán QR lặp


%% BƯỚC 5 (Lạc): Dạng Schur


%% BƯỚC 6 (Lạc): Trị riêng


%% BƯỚC 7 (Lạc): Vectơ riêng
