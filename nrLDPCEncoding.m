function nrLDPCEncoding(initial_message, BG, rate) %initial_message is a column vector

Kf_range = 512;
num_of_columns = ceil(numel(initial_message)/Kf_range);
matrix_of_Kfs = reshape([initial_message; zeros(Kf_range*num_of_columns - numel(initial_message), 1)], Kf_range, num_of_columns);
excess_bits = (Kf_range*num_of_columns)-length(initial_message);
[Kfs_rows, Kfs_cols] = size(matrix_of_Kfs);

snrdB=1:1:10;

FER=zeros(1,numel(snrdB));
BER=zeros(1,numel(snrdB));

for snrdB=1:1:10
    rx=[];
    errored_frames=0;
    errored_bits=0;

    %   num_of_dif_elements=0;

    for i = 1:Kfs_cols
        kf = matrix_of_Kfs(:, i);
        [B, Kcb, L, C, B_tonos, K_tonos, Kb, Zc, size_of_K, size_of_N, filler, ils, sheet, K, H_new, P_row, P_col, m, n, k] = ldpc_find_parameters(kf, BG);
        [Hsys,G] = makeSystematic(H_new,K);
        c=mod(K*G,2);
        codeword_length=size_of_K/rate;

        %Left (and right) side bit puncturing and bit shortening
        transmitted=c(1,[(2*Zc)+1:K_tonos , size_of_K+1:codeword_length+2*Zc]);
        % transmitted=c(1,(2*Zc)+1:codeword_length+2*Zc);

        modulated_segments=((2*transmitted)-1)*2;


        Eb=sumsqr(kf)/numel(kf);
        Ec=rate*Eb;
        a=sqrt(Ec);
        N0=Eb/(10.^(snrdB/10));

        sigma=sqrt(N0/2);
        noise=sigma*randn(size(modulated_segments));
        r=modulated_segments+noise;

        for j=1:length(r)
            p1(j)=1/(1+exp(-2*a*r(j)/(sigma^2)));
        end

        p0=1-p1;

        a=zeros(1,numel(filler));
        b=ones(1,numel(filler));
        q=0.5*ones(1,(2*Zc));
        d=0.5*ones(1,n-(2*Zc)-codeword_length);

        p0=[p0(1:numel(kf)-2*Zc), b, p0(numel(kf)-2*Zc+1:end)];
        p1=[p1(1:numel(kf)-2*Zc), a, p1(numel(kf)-2*Zc+1:end)];

        p0=[q,p0,d];
        p1=[q,p1,d];

        max_iter=10;

        [decoded_word]=decodeLDPC(p0,p1,H_new,max_iter);
        rx_word=decoded_word(1:numel(kf));

        if i~=Kfs_cols
            rx=cat(2,rx,rx_word); %ΑΛΛΑΓΗ ΕΔΩ ΣΤΟ rx ΜΕ ΑΛΛΗ ΜΕΤΑΒΛΗΤΗ! ΕΣΤΩ f ωστε να ερθω εδω κατω μετα δλδ --> a=cat(2,rx,rx_word);
            rx_word=rx_word';
            dif=sum(kf~=rx_word);
            errored_bits = errored_bits+dif;

            if dif~=0
                errored_frames = errored_frames+1;
            end
        else
            rx_w=rx_word(1:end-excess_bits);
            rx=cat(2,rx,rx_w); %ΑΛΛΑΓΗ ΕΔΩ ΣΤΟ rx ΜΕ ΑΛΛΗ ΜΕΤΑΒΛΗΤΗ! ΕΣΤΩ f ωστε να ερθω εδω κατω μετα δλδ --> a=cat(2,rx,rx_word);
            rx_w=rx_w';
            dif=sum(kf(1:end-excess_bits)~=rx_w);
            errored_bits = errored_bits+dif;

            if dif~=0
                errored_frames = errored_frames+1;
            end
        end

    end

    % rx = rx(1:end-excess_bits); %και στην 2η rx να ερθω και να βαλω εστω f δλδ --> rx=a(1:end-excess_bits)
    % rx=rx(:);
    % errored_bits = sum(initial_message~=rx);
    BER(snrdB) = errored_bits/numel(initial_message);
    FER(snrdB)=errored_frames/num_of_columns;
    disp(BER);
    disp(FER);
end

snrdB=1:1:10;

% BER vs SNR(dB)
figure(1);
semilogy(snrdB, BER, '-o');
grid on;
title('BER vs SNR(dB), 1000000 bits, 128 256 512 bits frame, N=10');
xlabel('SNR(dB)');
ylabel('BER');
hold on;

% FER vs SNR(dB)
figure(2);
semilogy(snrdB, FER, '-o');
grid on;
title('FER vs SNR(dB), 1000000 bits, 128 256 512 bits frame, N=10');
xlabel('SNR(dB)');
ylabel('FER');
hold on;

end


