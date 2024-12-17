function [B, Kcb, L, C, B_tonos, K_tonos, Kb, Zc, size_of_K, size_of_N, filler, ils, sheet, H_new, K, m, n, k] = ldpc_find_parameters(info_bits, BG)

A = [2 4 8 16 32 64 128 256;
    3 6 12 24 48 69 192 384;
    5 10 20 40 80 160 320 0;
    7 14 28 56 112 224 0 0;
    9 18 36 72 144 288 0 0;
    11 22 44 88 176 352 0 0;
    13 26 52 104 208 0 0 0;
    15 30 60 120 240 0 0 0];

Z = [2 3 4 5 6 7 8 9 10	11 12 13 14 15 16 18 20 22 24 26 28 30 32 36 40 44 48 52 56 60 64 72 80 88 96 104 112 120 128 144 160 176 192 208 224 240 256 288 320 352 384];

B = length(info_bits);
% B=info_bits; %info_bits is a number (the range of the segment)

if BG == 1
    Kcb=8448;
    Kb = 22;
else
    Kcb=3840;
    if B>640
        Kb = 10;
    elseif B>560 && B<640
        Kb = 9;
    elseif B>192 && B<560
        Kb = 8;
    else
        Kb = 6;
    end
end

if B < Kcb || B == Kcb
    L=0;
    C=1;
    B_tonos = B;
    K_tonos = (B_tonos/C);
else
    L=24;
    C=ceil(B/(Kcb-L));
    B_tonos = B+(C*L);
    K_tonos = (B_tonos/C);
end


[M,N]=size(Z);

for i=1:M
    for j=1:N
        a = Z(i,j);
        if Kb*a>=K_tonos
            Zc=a;
            break;
        end
    end
end

if BG==1
    excel_file='R1-1711982_BG1.xlsx';
    [BG_row, BG_col] = size(excel_file);
    Zc=a;
    size_of_K=22*Zc;
    size_of_N=66*Zc;
else
    excel_file='R1-1711982_BG2.xlsx';
    [BG_row, BG_col] = size(excel_file);
    Zc=a;
    size_of_K=10*Zc;
    size_of_N=50*Zc;
end

[rA, cA]=size(A);

for i=1:rA
    for j=1:cA
        if A(i,j) == Zc
            a=i-1;
        end
    end
end

ils=a;

switch ils
    case 0
        sheet = 2;
    case 1
        sheet = 3;
    case 2
        sheet = 4;
    case 3
        sheet = 5;
    case 4
        sheet = 6;
    case 5
        sheet = 7;
    case 6
        sheet = 8;
    case 7
        sheet = 9;
end

f=size_of_K-K_tonos;

filler=zeros(f,1);

K=[info_bits; filler];
K=K';
[a,k]=size(K);

P = readmatrix(excel_file, 'Sheet', sheet);

[P_row, P_col]=size(P);

H_new = cell(P_row*Zc,P_col*Zc);

for i = 1:P_row
    for j = 1:P_col
        if P(i,j)==-1
            H_new{i,j} = zeros(Zc,Zc);
        elseif P(i,j)==0
            I=eye(Zc);
            H_new{i,j} = I;
        else
            val = P(i,j);
            I=eye(Zc);
            H_new{i,j} = circshift(I, [0 val]);
        end
    end
end

resultMatrix = zeros(P_row*Zc, P_col*Zc);

for i=1:P_row
    for j=1:P_col
        row_start=(i-1)*Zc+1;
        row_end=i*Zc;
        col_start=(j-1)*Zc+1;
        col_end=j*Zc;
        resultMatrix(row_start:row_end, col_start:col_end)=H_new{i,j};
    end
end

H_new=resultMatrix;

[m,n]=size(H_new);

% H_new = cell2mat(H_new);
% S = sparse(H_new);
end