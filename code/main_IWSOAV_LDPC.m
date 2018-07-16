%%--------------------------------------------------------------------------------------------
% Simulation for BER performance of IW-SOAV in LDPC coded MIMO
%
% Author:
%   Ryo Hayakawa (e-mail: rhayakawa@sys.i.kyoto-u.ac.jp)
% Article:
%   Ryo Hayakawa and Kazunori Hayashi,
%   "Convex Optimization Based Signal Detection for Massive Overloaded MIMO Systems,"
%   IEEE Transactions on Wireless Communications, vol. 16, no. 11, pp. 7080-7091, Nov. 2017. 
%%--------------------------------------------------------------------------------------------

clear;
addpath myfunctions;
addpath data;

%% simulation setup
% number of transmit antennas
n=100;
% number of receive antennas
m=64;
% number of blocks per channel realization
nBlock=1;
% number of channel realization
nSample=1;
% array of SNR per receive antenna
arrSNR=[6,7,7.5,7.75,8,8.25,8.5,8.75,9,9.25,10,11,11.5,12];

% parameter for SOAV optimization problem
alpha=0.01;
% number of LLR update
nIteration=100;
% number of iteration in Douglas-Rachford algorithm
K=30;
% parameter for Douglas-Rachford algorithm
gamma=1;

% parity check matrix
load('parityCheckMatrix_2000_4000.mat');
% load('parityCheckMatrix_4000_8000.mat');
matParityCheck=matParityCheck_sp;
% size of code block
[nBlockSize,nCodeBlockSize]=size(matParityCheck);
% code rate
R=nBlockSize/nCodeBlockSize;
% number of transmission per information block
nTrans=nCodeBlockSize/2/n;

%% comm object
% QPSK Modulator
hMod=comm.QPSKModulator('BitInput',true);
% Interleaver/Deinterleaver
% interleaverIndices=randperm(nCodeBlockSize)';
interleaverIndices=(1:nCodeBlockSize)'; % no interleave
hInt=comm.BlockInterleaver(interleaverIndices);
hDeint=comm.BlockDeinterleaver(interleaverIndices);
% LDPC Encoder/Decoder
hLDPCEnc=comm.LDPCEncoder('ParityCheckMatrix',matParityCheck);
hLDPCDec=comm.LDPCDecoder('ParityCheckMatrix',matParityCheck,'OutputValue','Whole codeword','DecisionMethod','Soft decision','IterationTerminationCondition','Parity check satisfied','FinalParityChecksOutputPort',true);
% BER measurement
hBER=comm.ErrorRate('ResetInputPort',true);

%% run simulation
matNumError=zeros(nIteration,length(arrSNR));
matNumBit=zeros(nIteration,length(arrSNR));
arrNumError_IWSOAV=zeros(1,length(arrSNR));
arrNumBit_IWSOAV=zeros(1,length(arrSNR));
for sampleIndex=1:nSample
  % channel matrix
  [H_comp,H]=makeChannel(m,n);
  HH=H'*H;
  
  % matrix for prox_g
  ProjMat1_SOAV=(eye(2*n)+alpha*gamma*HH)^(-1);
  ProjMat2_SOAV=alpha*gamma*H';
  
  disp(['sampleIndex=' num2str(sampleIndex)]);
  for SNRIndex=1:length(arrSNR)
    rng(sampleIndex);
    disp(['  SNR=' num2str(arrSNR(SNRIndex))]);
    % noise of additive noise
    sigma=sqrt(n*2/(10^(arrSNR(SNRIndex)/10)));
    sigma_r=sigma/sqrt(2);
    
    for blockIndex=1:nBlock
      % information bits
      data=randi([0,1],nBlockSize,1);
      % code bits
      dataEnc=step(hLDPCEnc,data);
      % interleaved code bits
      dataEncInt=step(hInt,dataEnc);
      % transmitted symbols
      dataSymbol=step(hMod,dataEncInt)*sqrt(2);
      matDataSymbol=reshape(dataSymbol,[n,nTrans]);
      % additive noise vector
      V_comp=randn(m,nTrans)*sigma_r+1j*randn(m,nTrans)*sigma_r;
      % received signal vector
      Y_comp=H_comp*matDataSymbol+V_comp;
      Y=[real(Y_comp); imag(Y_comp)];
      
      %% signal detection using W-SOAV
      LLR_toDetector=zeros(2*n*nTrans,1);
      LLR_fromDetector=zeros(2*n*nTrans,1);
      LLR_fromDecoder=zeros(2*n*nTrans,1);
      bParityCheck=0;
      for iterationIndex=1:nIteration
        if bParityCheck==0
          matLLR_fromDetector=zeros(2*n,nTrans);
          for timeIndex=1:nTrans
            tmpLLR_toDetector=[LLR_toDetector(((timeIndex-1)*2*n+2):2:(timeIndex*2*n));...
                               LLR_toDetector(((timeIndex-1)*2*n+1):2:(timeIndex*2*n));];
            % weight update
            w=1-1./(1+exp(tmpLLR_toDetector));
            % weighted SOAV optimization
            y=Y(:,timeIndex);
            s_hat=WSOAV(y,H,ProjMat1_SOAV,ProjMat2_SOAV,gamma,K,w);
            s_hat=min(1,max(-1,s_hat));
            % calculate LLR
            vec_mu=H*s_hat;
            vec_sigma2=H.^(2)*(1-s_hat.^(2))+sigma_r^(2);
            matNumer=2*H.*(y*ones(1,2*n)-(vec_mu*ones(1,2*n)-H.*(ones(2*m,1)*s_hat')));
            matDenom=vec_sigma2*ones(1,2*n)-H.^(2).*(ones(2*m,2*n)-ones(2*m,1)*(s_hat.^(2))');
            matLLR_fromDetector(:,timeIndex)=sum(matNumer./matDenom)';
          end
          vecS_real_est=reshape(matLLR_fromDetector(1:n,:),[n*nTrans,1]);
          vecS_imag_est=reshape(matLLR_fromDetector((n+1):(2*n),:),[n*nTrans,1]);
          LLR_fromDetector_former=LLR_fromDetector;
          LLR_fromDetector=reshape([vecS_imag_est vecS_real_est]',[2*n*nTrans,1]);
          
          LLR_toDecoder=step(hDeint,LLR_fromDetector);
          LLR_fromDecoder_former=LLR_fromDecoder;
          [LLR_fromDecoder,vecParityCheck]=step(hLDPCDec,LLR_toDecoder);
          if sum(vecParityCheck)==0
            bParityCheck=1;
          end
          LLR_toDetector=step(hInt,LLR_fromDecoder);
          
          % calculate BER in each iteration
          LLR_decData=LLR_fromDecoder(1:nBlockSize);
          decData=double(LLR_decData<0);
          BERtmp=step(hBER,data,decData,1);
          matNumError(iterationIndex,SNRIndex)=matNumError(iterationIndex,SNRIndex)+BERtmp(2);
          matNumBit(iterationIndex,SNRIndex)=matNumBit(iterationIndex,SNRIndex)+BERtmp(3);
        else
          matNumError(iterationIndex,SNRIndex)=matNumError(iterationIndex,SNRIndex)+BERtmp(2);
          matNumBit(iterationIndex,SNRIndex)=matNumBit(iterationIndex,SNRIndex)+BERtmp(3);
        end
      end
      
      %% signal detection using IW-SOAV
      matLLR_IWSOAV=zeros(2*n,nTrans);
      for timeIndex=1:nTrans
        LLR=zeros(2*n,1);
        s_init=zeros(2*n,1);
        y=Y(:,timeIndex);
        for iterationIndex=1:5
          % weight update
          w=1-1./(1+exp(LLR));
          % weighted SOAV optimization
          s_hat=WSOAV(y,H,ProjMat1_SOAV,ProjMat2_SOAV,gamma,K,w);
          s_hat=min(1,max(-1,s_hat));
          % LLR calculation
          vec_mu=H*s_hat;
          vec_sigma2=H.^(2)*(1-s_hat.^(2))+sigma_r^(2);
          matNumer=2*H.*(y*ones(1,2*n)-(vec_mu*ones(1,2*n)-H.*(ones(2*m,1)*s_hat')));
          matDenom=vec_sigma2*ones(1,2*n)-H.^(2).*(ones(2*m,2*n)-ones(2*m,1)*(s_hat.^(2))');
          LLR=sum(matNumer./matDenom)';
        end
        matLLR_IWSOAV(:,timeIndex)=LLR;
      end
      vec_real_LLR_IWSOAV=reshape(matLLR_IWSOAV(1:n,:),[n*nTrans,1]);
      vec_imag_LLR_IWSOAV=reshape(matLLR_IWSOAV((n+1):(2*n),:),[n*nTrans,1]);
      LLR_IWSOAV=reshape([vec_imag_LLR_IWSOAV vec_real_LLR_IWSOAV]',[2*n*nTrans,1]);
      
      LLR_IWSOAV_deint=step(hDeint,LLR_IWSOAV);
      [LLR_IWSOAV_dec,vecParityCheck]=step(hLDPCDec,LLR_IWSOAV_deint);
      
      % calculate BER
      LLR_decData_IWSOAV=LLR_IWSOAV_dec(1:nBlockSize);
      decData_IWSOAV=double(LLR_decData_IWSOAV<0);
      BERtmp=step(hBER,data,decData_IWSOAV,1);
      arrNumError_IWSOAV(SNRIndex)=arrNumError_IWSOAV(SNRIndex)+BERtmp(2);
      arrNumBit_IWSOAV(SNRIndex)=arrNumBit_IWSOAV(SNRIndex)+BERtmp(3);
      
    end
  end
end
matBER=matNumError./matNumBit;
arrBER_IWSOAV=arrNumError_IWSOAV./arrNumBit_IWSOAV;

%% Display results
arrMarker=['o';'^';'x';'>';'+';'v';'*';'<';'s';'d';'p';'h';];
arrIteration=[1:5 10:10:100];
close all;
figure;
setLegend={'Independent det./dec.'};
semilogy(arrSNR,arrBER_IWSOAV,'-k','LineWidth',2,'MarkerSize',8);
hold on;
for iterationIndex=arrIteration
  semilogy(arrSNR,matBER(iterationIndex,:),['-' arrMarker(rem(iterationIndex,9)+1)],'LineWidth',2,'MarkerSize',8);
  setLegend=[setLegend; ['Joint det.dec. ($L_{\mathrm{max}}=' num2str(iterationIndex) '$)'];];
  hold on;
end
grid on;
objLegend=legend(setLegend,'Location','northeast');
xlabel('SNR per receive antenna (dB)');
ylabel('BER');
objLegend.Interpreter='latex';
objLegend.FontSize=14;
fig=gca;
fig.FontSize=18;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([arrSNR(1) arrSNR(length(arrSNR)) 1e-5 1]);

