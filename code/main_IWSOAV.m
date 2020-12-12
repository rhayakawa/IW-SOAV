%%--------------------------------------------------------------------------------------------
% Simulation for BER performance of IW-SOAV
%
% Author:
%   Ryo Hayakawa
% Article:
%   Ryo Hayakawa and Kazunori Hayashi,
%   "Convex Optimization Based Signal Detection for Massive Overloaded MIMO Systems,"
%   IEEE Transactions on Wireless Communications, vol. 16, no. 11, pp. 7080-7091, Nov. 2017. 
%%--------------------------------------------------------------------------------------------

clear;
addpath myfunctions;

%% paremeter settings
% number of transmit antennas
n=100;
% number of receive antennas
m=64;
% number of symbol vectors per channel realization and SNR
nSymbolVector=100;
% number of channel realization
nSample=1;
% SNR
arrSNR=0:2.5:30;

%% parameters of IW-SOAV
% number of weight updates
nIteration=5;
% parameter of optimization problem
arrAlpha=[1e-2,1e-2,1e-2,1e-2,1e-2,1e-1,1e-1,1e-1,1e-1,3e-1,1,1,1];
% parameter for Douglas-Rachford algorithm
gamma=1;
% number of iteration in Douglas-Rachford algorithm
K=50;

%% run simulation
matNumError=zeros(nIteration,length(arrSNR));
matNumBit=zeros(nIteration,length(arrSNR));
for i=1:nSample
  % channel matrix
  [H_comp,H]=makeChannel(m,n);  
  HH=H'*H;
  
  disp(['i=' num2str(i)]);
  for SNRIndex=1:length(arrSNR)
    SNR=arrSNR(SNRIndex);
    % parameter of optimization problem
    alpha=arrAlpha(SNRIndex);
    %     disp(['  SNR=' num2str(SNR) '  alpha=' num2str(alpha)]);
    % variance of additive noise
    sigma=sqrt(n*2/(10^(SNR/10)));
    sigma_r=sigma/sqrt(2);
    
    % matrix for prox_g
    ProjMat1_SOAV=(eye(2*n)+alpha*gamma*HH)^(-1);
    ProjMat2_SOAV=alpha*gamma*H';
    
    for symbolVectorIndex=1:nSymbolVector
      % transmitted signal vector
      data=randi([0,1],2*n,1);
      s=-data*2+ones(2*n,1);
      % additive noise vector
      v=randn(2*m,1)*sigma_r;
      % received signal vector
      y=H*s+v;
      
      %% signal detection via IW-SOAV
      LLR=zeros(2*n,1);
      s_hat=LLR;
      for iterationIndex=1:nIteration
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
        % calculate BER in each iteration
        matNumError(iterationIndex,SNRIndex)=matNumError(iterationIndex,SNRIndex)+nnz(s-sign(LLR));
        matNumBit(iterationIndex,SNRIndex)=matNumBit(iterationIndex,SNRIndex)+2*n;
      end
    end
  end
end
matBER=matNumError./matNumBit;

%% Display results
arrMarker=['o';'^';'x';'>';'+';'v';'*';'<';'s';'d';'p';'h';];
close all;
figure;
setLegend={};
for iterationIndex=1:nIteration
  semilogy(arrSNR,matBER(iterationIndex,:),['-' arrMarker(iterationIndex)],'LineWidth',2,'MarkerSize',8);
  setLegend=[setLegend; ['IW-SOAV ($L= ' num2str(iterationIndex) '$)'];];
  hold on;
end
grid on;
objLegend=legend(setLegend,'Location','southwest');
xlabel('SNR per receive antenna (dB)');
ylabel('BER');
objLegend.Interpreter='latex';
objLegend.FontSize=18;
fig=gca;
fig.FontSize=18;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([arrSNR(1) arrSNR(length(arrSNR)) 1e-5 1]);

