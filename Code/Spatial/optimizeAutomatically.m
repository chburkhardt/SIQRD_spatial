function varargout = optimizeAutomatically()
 
  %Before the simulation is started, the parameter bounds must be adjusted for the individual steps
  %a folder Name must be defined
  %the model must be defined
  %it must be defined which reduction factors should be fitted 
  folder = "SIQRD1stWaveAuto";
  addpath("Optimize")
  parameterArray = read_parameter('../../Results/parameter.txt');
  parameterArray{1}.model = "SIREDmod";
  parameterArray{1}.folderName = [folder, "betaDF"];
  parameterArray{1}.saveVTK = false;
  parameterArray{1}.saveDiagrams = false;
  parameterArray{1}.showDiagrams = false;
  parameterArray{1}.blowUp = false;
  parameterArray{1}.spatial = "germany";
  parameterArray{1}.initial = "RKIfiles";
  parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
  parameterArray{1}.startDate = datenum([2020, 03, 02]);
  %parameterArray{1}.endDateOpt = datenum([2020, 04, 30]);
  parameterArray{1}.beta_cross_county = 1;
  parameterArray{1}.automatically = true;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Step 1: optimize German-wide beta and DF
    
  parameterArray{1}.wRatioGlobStates = 0; % 1 is fully statewise, 0 is only global
  parameterArray{1}.wRatioID = 0.5; % 1 is only infected, 0 only death
  parameterArray{1}.wISlope = 0.05; % weight for slope of infected last date
  parameterArray{1}.reduceToStates= true;
  parameterArray{1}.optFunGermanSum = true;
  parameterArray{1}.schoolClosing	= 1;
  parameterArray{1}.contactRestrictions = 1;
  parameterArray{1}.exitRestrictions	= 1;		
  parameterArray{1}.majorEvents		= 1;
  parameterArray{1}.easingContactRestrictions	= 1;
  parameterArray{1}.endOfVacation	= 1;
  parameterArray{1}.liftingTravelRestictionsEU	= 1;
  parameterArray{1}.lockdownLight =1;
  parameterArray{1}.showStatePlots = true;
    
  parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
  parameterArray{1}.gamma2 = 0.04; % qarantine -> Recovered
  parameterArray{1}.mortality = 0.006; % 0.006
 
	parameterArray{1}.betaSWscaling = 1;		
  
  %delete betaStatewise and beta_cc_statewise from parameters for first fitting steps
  if isfield(parameterArray{1},'betaStatewise')
    parameterArray{1} = rmfield(parameterArray{1},'betaStatewise');
  elseif isfield(parameterArray{1},'beta_cc_statewise')
    parameterArray{1} = rmfield(parameterArray{1},'beta_cc_statewise');
  end
  
  parameterArray{1}.paraNames = {"beta", "contactRestrictions", "darkFigure"};
  parameterArray{1}.x0Cov = [0.157 0.4 10]; 
  parameterArray{1}.lb = [0 0 8];
  parameterArray{1}.ub = [0.5 1 20];
  
  parameterArray{1}.optAlgorithm = "pso";
  parameterArray{1}.PSOiterations = 250;
  
  %"lastWeek" - calculate factorsDF from lastWeek in fitting period (use when data is 
  %alsost complete)
  %"twoWeeksToSimEnd" - calculate factorsDF from week two weeks before fitting period
  %ends (Use when RKI data is most probably not yet complete)
  parameterArray{1}.timeForFactorsDF = "lastWeek"; %"lastWeek"|"twoWeeksToSimEnd"
  
  %"lastWeek" - fit deaths on last day of fitting period 
  %"twoWeeksToSimEnd" - fit deaths two weeks before fitting period ends
  parameterArray{1}.dayForDeathFit = "lastDay"; %"lastDay"|"twoWeeksToSimEnd"
    
  argOut1 = optimizeStateWise(parameterArray);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Step2: Optimize betaStatewise 
  
  for i=1:length(parameterArray{1}.paraNames)
    parameterArray{1} = setfield(parameterArray{1}, parameterArray{1}.paraNames{i}, argOut1(i));
  end
  
  parameterArray{1}.folderName = [folder, "betaSW"];
  parameterArray{1}.wRatioGlobStates = 1; % 1 is fully statewise, 0 is only global
  parameterArray{1}.wRatioID = 0.5; % 1 is only infected, 0 only death
  parameterArray{1}.wISlope = 0.05; % weight for slope of infected last date
  parameterArray{1}.optFunGermanSum = false;
  parameterArray{1}.beta = 1;
  parameterArray{1}.PSOiterations = 400;
  
  parameterArray{1}.paraNames = {"betaCorr1", "betaCorr2", "betaCorr3", "betaCorr4","betaCorr5",...
  "betaCorr6", "betaCorr7", "betaCorr8", "betaCorr9", "betaCorr10",...
  "betaCorr11", "betaCorr12", "betaCorr13", "betaCorr14", "betaCorr15", "betaCorr16"};
  parameterArray{1}.x0Cov = [0.834401 1.479429 0.622385 1.108438 1.139666 0.500897 0.298197 1.147868 1.612068 1.880193 0.772545 1.191203 0.598741 1.202886 0.56924 0.833693]; 
  parameterArray{1}.lb = zeros(1,16);
  parameterArray{1}.ub = parameterArray{1}.lb+0.5;  
  
  argOut2 = optimizeStateWise(parameterArray);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Step3: Rebalance Network
  
  parameterArray{1}.betaStatewise = argOut2(1:16);
  parameterArray{1}.folderName = [folder, "betaCC"];
  
  argOut3 = naiveGradientSearch(parameterArray);
  
endfunction
