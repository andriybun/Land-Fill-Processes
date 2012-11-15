% Approach should be fit for running in parallel mode
clear
clc

pname = './';
addpath(genpath(pname));

%% Initialize DDREAM

    MCMCPar.ndraw = 500;                  % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    MCMCPar.pJumpRate_one = 0.80;           % Probability of selecting a jumprate of 1 --> jump between modes
 
    
    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 2;                     % Number of steps before calculating convergence diagnostics
       
    ExtraPar.pCR = 'Update';                   % Adaptive tuning of crossover values
   
        % --------------------------------------- Added for reduced sample storage -------------------
    ExtraPar.reduced_sample_collection = 'No'; % Thinned sample collection?
    ExtraPar.T = 1000;                         % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % What type of initial sampling
    ExtraPar.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    %                   k(hyd)  qsmax(meth)
    ParRange.minn =    [0.01    0.1]; 
    ParRange.maxn =    [0.1     0.5]; 
    % Define the boundary handling
    ExtraPar.BoundHandling = 'Reflect';
    % Save in memory or not
    ExtraPar.save_in_memory = 'Yes';

    % Constant parameters for bioreactor
    [ExtraPar.ORI ExtraPar.Meas ExtraPar.Comp ExtraPar.S ExtraPar.Rp ExtraPar.I] = bioreactorP();

    %Measurent Data

    Measurement.MeasData = ExtraPar.Meas.data; %the measurements
    Measurement.N = size(Measurement.MeasData(:),1); %number of measurements
    tm = ExtraPar.Meas.tm; ldata= ExtraPar.Meas.ldata; mdata = ExtraPar.Meas.data; 
%     Measurement.Sigma = zeros(2*tml+tm2l,1); %use individual measurement errors for every measurement, assuming the
%     Measurement.Sigma(1:tml) = 0.4252; % Bio
%     Measurement.Sigma(tml+1:2*tml) = 0.3110;  % pH
%     Measurement.Sigma(2*tml+1:end) = 0.0369; % VFA

    %Parameters are:
    MCMCPar.n = 2;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z 
    ExtraPar.FloorParRange = zeros(1,MCMCPar.n);  %these are used in the Run_ADE function to scale the dream parameters between 0 and 10 
    ExtraPar.FactorParRange = ones(1,MCMCPar.n); %back to the actual parameter values: RealPar = FloorParRange + FactorParRange*DreamPar
                
    ModelName = 'bioreactor';

% Define likelihood function
ExtraPar.Option = 3; %(option 6 implemented by S.Korteland, august 2011)

%Initialize random generator
Restart='No';
[Sequences,Reduced_Seq,X,Z,output] = dream_zs_nocrit(MCMCPar,ParRange,Measurement,ModelName,ExtraPar,ExtraPar.Option,Restart);


%% And plot some results
isub = 2;
 
for ii=1:MCMCPar.n
ParData(:,ii) = ExtraPar.FloorParRange(ii)+ExtraPar.FactorParRange(ii).*Z(:,ii);
end

%plot the evolution of archive Z
ylabelset = {'k(hyd)','qsmax(meth)'};
Zind=max(find(Z(:,1)>0)); %the end-value of Z
 figure;
 for ii=1:MCMCPar.n;
     ind=ii;
     subplot(isub,isub,ii);
     plot(ParData(1:Zind,ind),'k.');
      ylabel(ylabelset{ii},'FontSize',10)
      xlabel('Number of samples stored in Z','FontSize',10)
     set(gca,'FontSize',16);
     if ii ~= 11
         axis([0 Zind ParRange.minn(ii) ParRange.maxn(ii)])
     else
     axis([0 Zind -1*ParRange.minn(ii) -1*ParRange.maxn(ii)])
     end
 end
 %also plot the simulation with the minimum objective function/maximum probability
 [minZ minZidx] = max(Z(MCMCPar.m0+1:Zind,MCMCPar.n+1));
 minZidx=minZidx+MCMCPar.m0;  %the index corresponding to the parameter set with the highest probability
 for ii=1:MCMCPar.n;
     ind=ii;
      subplot(isub,isub,ii);hold on;
      plot(minZidx,ParData(minZidx,ind),'ro');
 end

 %Plot the evolution of the individual sequences (the Markov chains)
 figure;
 for ii=1:MCMCPar.n;
     Sind=max(find(Sequences(:,1,1)>0));
     Seq=reshape(Sequences(1:Sind,ii,:),Sind,3);
     subplot(isub,isub,ii);plot(ExtraPar.FloorParRange(ii)+ExtraPar.FactorParRange(ii).*Seq,'.')
 end
 
 %plot histograms of the last 25% of the data (used to calculate mean and
 %standard deviation
 figure('Position',[188 388 1303 600]);
 for ii=1:MCMCPar.n
 subplot(isub,isub,ii);
 hist(ParData(round(0.75.*Zind):Zind,ii),10)
  xlabel(ylabelset{ii},'FontSize',16)
  set(gca,'FontSize',16);
 h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5])
 end
  
 %plot the Gelman-Rubin Convergence criterion, when this criterion <1.2 for
 %all parameters, convergence generally is declared
 figure;hold on;
 GR_ind= max(find(output.R_stat(:,1)>0));
 plot(output.R_stat(1:GR_ind,1),output.R_stat(1:GR_ind,2:end),'.-');
 x=1:(output.R_stat(GR_ind)+200);
 y = ones(size(x)).*1.2;
 plot(x,y,'k--');
 axis([-Inf Inf 0 10])
 hl=findobj(gca,'Type','Line');
%  legend([hl(4) hl(3) hl(2) hl(1)]...
%          ,{'1','2','3','4'})
%  legend([hl(5) hl(4) hl(3) hl(2) hl(1)]...
%          ,{'1','2','3','4','5'})
 xlabel('number of model evaluations')
 ylabel('Gelman-Rubin convergence criterion')
  

 %Simulate and plot optimum solution
 SimPar_opt=Z(minZidx,1:MCMCPar.n);
 SimData = bioreactor(SimPar_opt,ExtraPar);
N = size(Measurement.MeasData,1);
Residual = Measurement.MeasData-SimData;
RMSE = sqrt(sum(Residual.^2)./(N-MCMCPar.n))
% % WSSR = sum((Residual./Measurement.Sigma).^2) %equal to the probabily density function of dream, which is -WSSR (column n+1 in Z)
% 
figure;

subplot(2,2,1)
plot(tm(1:ldata(1)),mdata(1:ldata(1)),'.b-');
hold on;plot(tm(1:ldata(1)),SimData(1:ldata(1)),'or-')
xlabel('Time (min)')
ylabel('Biogas (mol)')
hl=findobj(gca,'Type','Line');
legend([hl(2) hl(1)],{'Measurements','Simulation'},'Location','NorthWest')

subplot(2,2,2)
plot(tm(ldata(1)+1:ldata(1)+ldata(2)),mdata(ldata(1)+1:ldata(1)+ldata(2)),'.b-');
hold on;plot(tm(ldata(1)+1:ldata(1)+ldata(2)),SimData(ldata(1)+1:ldata(1)+ldata(2)),'or-')
xlabel('Time (min)')
ylabel('pH')
hl=findobj(gca,'Type','Line');
legend([hl(2) hl(1)],{'Measurements','Simulation'},'Location','NorthWest')

subplot(2,2,3)
plot(tm(ldata(1)+ldata(2)+1:end),mdata(ldata(1)+ldata(2)+1:end),'.b-');
hold on;plot(tm(ldata(1)+ldata(2)+1:end),SimData(ldata(1)+ldata(2)+1:end),'or-')
xlabel('Time (min)')
ylabel('Acitic acid (mol/L)')
hl=findobj(gca,'Type','Line');
legend([hl(2) hl(1)],{'Measurements','Simulation'},'Location','NorthWest')
