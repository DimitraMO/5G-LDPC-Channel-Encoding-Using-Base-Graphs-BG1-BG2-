function  [Reliability, bitErr, frameErr, EbN0]=nrLDPCEncoding_3(initial_message,BG,rate)

% for BG=1:2
ber=zeros();
fer=zeros();
Reliability=zeros();
EbN0=0:1:10;   
rep=1;

for repetition=1:rep

    Kf_range=512;
    bitsTotal=length(initial_message);
    frames=ceil(bitsTotal/Kf_range);
    matrixOfKfs=zeros(Kf_range, frames);
    extraZeroBits=numel(matrixOfKfs)-bitsTotal;

    for i=1:frames
        s=(i-1)*Kf_range+1;
        e=min(i*Kf_range,bitsTotal);
        matrixOfKfs(1:e-s+1,i)=initial_message(s:e);
    end

    lastColumn=matrixOfKfs(:,end);

    size=length(initial_message);
    maxIter=20;
    err=zeros(1,length(EbN0));
    ferr=zeros(1,length(EbN0));
    BER=zeros(1,length(EbN0));
    nvar=zeros(1,length(EbN0));
    FER=zeros(1,length(EbN0));


    for frame=1:frames
        ERR=zeros(1,length(EbN0));
        FERR=zeros(1,length(EbN0));

        disp(frame);
        bits=matrixOfKfs(:,frame);
        [B, Kcb, L, C, B_tonos, K_tonos, Kb, Zc, size_of_K, size_of_N, filler, ils, sheet, H_new, K, m, n, k] = ldpc_find_parameters(bits, BG);
        [Hsys,G] = makeSystematic(H_new,K);
        c=mod(K*G,2);

        %parity bit puncturing
        c_1=c(1, 1:ceil(k/(rate*Zc))*Zc);
        %systematic bit puncturing
        c_2=c_1(1, (2*Zc)+1:end);

        %BPSK modulation
        tx=2*c_2-1;

        % signal_mean=0;
        signal_mean=sum(tx(1:Kf_range-(2*Zc)))/(Kf_range-(2*Zc));
        signal_variance=sum((tx(1:(Kf_range-(2*Zc)))-signal_mean).^2)/(Kf_range-(2*Zc));

        if signal_variance==0
            signal_variance=signal_variance+1;
        end

        %the new parity check matrix after rate matching
        H=H_new(1:((ceil(k/(rate*Zc))*Zc)-k), 1:ceil(k/(rate*Zc))*Zc);

        for q=1:length(EbN0)
            erroredBits=0;
            snr=EbN0(q);
            disp(snr);
            noise_variance=signal_variance/(2*rate*10^(0.1*snr));
            % noise_variance=0;
            nvar(q)=noise_variance;
            noisyBits=tx+sqrt(noise_variance)*randn(1,length(tx));
            r=noisyBits;

            firstTwoZcBits=zeros(1, 2*Zc);
            % fillerBits=-Inf(length(filler), 1);

            % v=r(1:(Kf_range-2*Zc));
            % l=r((Kf_range-2*Zc)+1:end);
            r_new=[firstTwoZcBits, r];

            for w=1:length(r_new)
                p1(w)=1/(1+exp((-2*r_new(w))/noise_variance));
            end

            p0=1-p1;

            % [decoded_word]=minSumDecoding(H,r_new,Zc,maxIter);
            % out = decodeLDPC(p0,p1,H,maxIter);
            [out, success, k, prob ] = ldpc_decode(p0,p1,H,maxIter);
            % out=out';
            % out=decoded_word';

            erroredBits=sum(bits~=out(1:Kf_range));
            ERR(q)=ERR(q)+erroredBits;
            disp(ERR);
            if erroredBits~=0
                FERR(q)=FERR(q)+1;
            else
                FERR(q)=FERR(q);
            end

        end
        disp(nvar);
        disp(ERR);
        err=err+ERR;
        ferr=ferr+FERR;

    end

    BER=err/length(initial_message);
    ber=ber+BER;
    FER=ferr/frames;
    fer=fer+FER;
end

%Reliability Calculation
for i=1:length(EbN0)
    Reliability(i)=(1-FER(i))*100;
end

bitErr=ber/rep;
frameErr=fer/rep;


figure(1);
semilogy(EbN0,bitErr);
title("BER vs SNR");
xlabel("EbN0");
ylabel("Bit Error Rate (BER)");
grid on;
hold on;

figure(2);
semilogy(EbN0,frameErr);
title("FER vs SNR");
xlabel("EbN0");
ylabel("Frame Error Rate (FER)");
grid on;
hold on;

end

