## Copyright (C) 2020 andre
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} optimize (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-24

function varargout = optimizeStateWise (parameterArray) %filelname 1 und filename 2
  tOpt = tic;
  addpath('../', '../Pso');
  %% Example  
  
  if nargin == 0
    close all;
    parameterArray = read_parameter('../../Results/parameter.txt');  
    %% the entries can be directly modified in the struct
    %% eg. set a foldername
    parameterArray{1}.model = "SIREDmod";
    %parameterArray{1}.totalRuntime = 40;
    parameterArray{1}.folderName = "Test2ndWave";
    parameterArray{1}.saveVTK = false;
    parameterArray{1}.saveDiagrams = false;
    parameterArray{1}.showDiagrams = false;
    parameterArray{1}.fullConsoleOut = true;
    parameterArray{1}.blowUp = false;
    parameterArray{1}.spatial = "germany";
    parameterArray{1}.initial = "RKIfiles";%"GitData";
    parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
    parameterArray{1}.startDate = datenum([2020, 03, 02]);
    parameterArray{1}.endDateOpt = datenum([2020, 04, 30]);
    %parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
    %parameterArray{1}.startDate = datenum([2020, 08, 01]);
    %parameterArray{1}.endDateOpt = datenum([2020, 11, 25]);
    parameterArray{1}.reduceToStates= true;
    %parameterArray{1}.betaStateWise= false;
    parameterArray{1}.beta_cross_county = 1;
    
    %parameterArray{1}.slopeFitInterval = "end";
    parameterArray{1}.slopeFitInterval = "lastWeek";
    %parameterArray{1}.dayForDeathFit = "twoWeeksToSimEnd";    
    parameterArray{1}.dayForDeathFit = "twoWeeksToSimEnd";
    %parameterArray{1}.timeForFactorsDF = "twoWeeksToSimEnd";
    parameterArray{1}.timeForFactorsDF = "twoWeeksToSimEnd";
    
    parameterArray{1}.wRatioGlobStates = 0; % 1 is fully statewise, 0 is only global
    parameterArray{1}.wRatioID = 1; % 1 is only infected, 0 only death
    parameterArray{1}.wISlope = 0.05; % weight for slope of infected last date
    parameterArray{1}.dayForDeathFit = "lastDay"; 
    parameterArray{1}.optFunGermanSum = true;
    parameterArray{1}.schoolClosing	= 1;
    parameterArray{1}.contactRestrictions = 1;%0.14;%091810;%0.0757;%1;
    parameterArray{1}.exitRestrictions	= 1;		
    parameterArray{1}.majorEvents		= 1;%0.658481;%0.6615;%1;
    parameterArray{1}.easingContactRestrictions	= 1;
    parameterArray{1}.endOfVacation	= 1;
    parameterArray{1}.liftingTravelRestictionsEU	= 1;
    parameterArray{1}.lockdownLight =1;
    %parameterArray{1}.n_runs = 1;
    parameterArray{1}.showStatePlots = true;
    
    parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
    parameterArray{1}.gamma2 = 0.04; % qarantine -> Recovered
    parameterArray{1}.mortality = 0.006; % 0.006
    parameterArray{1}.darkFigure = 18;%11.2;%9.95; %6.5;
    parameterArray{1}.beta = 0.5;%0.138465;%0.1413;%1.842096; %0.196693; %0.17;
    %parameterArray{1}.betaStatewise = [0.148106,0.260469,0.109561,0.195322,0.202218,0.086421,0.052591,0.212167,0.281259,0.329717,0.144268,0.217264,0.106012,0.210734,0.097104,0.155813];
	
	%parameterArray{1}.beta =  1;
	parameterArray{1}.betaSWscaling = 1;		
##	parameterArray{1}.majorEvents	= 0.318781; %0.302252; %0.2092;
##	parameterArray{1}.contactRestrictions = 0.446369;%0.499724; %0.5956;
##	

##    parameterArray{1}.showStatePlots = false;
    ##  ##  % For the countywise model
    ##    parameterArray{1}.reduceToStates= false;
    ##    parameterArray{1}.beta_cross_county =  0.87869;%0.436106125; % 1.08 ohne ss to cs anpassung
  endif
  
  %days können wir später dann in das Parameterarry mit aufnehmen
  %days beschreibt die Zeitspanne in der die Optimierung laufen soll
  %außerdem brauchen wir noch einen Parameter, der bestimmt, ob wir 
  %statewise oder für ganz Deutschland fitten
  
  RKIread = read_case_history_RKIfiles();
  
  %rewrite the cell, since states are arranged differently in RKIread and statewiseSIR
  RKIdata = cell(16,1);
  RKIdata{1} = RKIread{15};
  RKIdata{2} = RKIread{6};
  RKIdata{3} = RKIread{9};
  RKIdata{4} = RKIread{5};
  RKIdata{5} = RKIread{10};
  RKIdata{6} = RKIread{7};
  RKIdata{7} = RKIread{11};
  RKIdata{8} = RKIread{1};
  RKIdata{9} = RKIread{2};
  RKIdata{10} = RKIread{12};
  RKIdata{11} = RKIread{3};
  RKIdata{12} = RKIread{4};
  RKIdata{13} = RKIread{8};
  RKIdata{14} = RKIread{13};
  RKIdata{15} = RKIread{14};
  RKIdata{16} = RKIread{16};
  
  %cumulate infection numbers 
  %create splines
  timesRKI = cell2mat(RKIdata{1}.time);
  parameterArray{1}.dataStatus = datestr(timesRKI(end));
  for i=1:length(RKIdata)
    %RKIinf = RKIdata{i}.infected;
    RKIdata{i}.infected = cell2mat(RKIdata{i}.infected);
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end
  
##    darkFiguresStates = zeros(1, length(RKIdata));
##    for i=1:length(RKIdata)
##      darkFiguresStates(i) = ppval(RKIdata{i}.splineDead, parameterArray{1}.endDateOpt)/...
##      ppval(RKIdata{i}.splineInfected, parameterArray{1}.endDateOpt - 15);
##    end
##    darkFiguresStates /= parameterArray{1}.mortality;
  
##  darkFiguresStates = zeros(1, length(RKIdata));
##  for i=1:length(RKIdata)
##    darkFiguresStates(i) = ppval(RKIdata{i}.splineDead, parameterArray{1}.endDateOpt)/...
##    (ppval(RKIdata{i}.splineInfected, (parameterArray{1}.endDateOpt - 15))*...
##    parameterArray{1}.mortality);
##  end
##  darkFiguresStates /= mean(darkFiguresStates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %calculation of DF factors for states by end of fitting period
##  darkFiguresStates = zeros(1, length(RKIdata));
##  for i=1:length(RKIdata)
##    darkFiguresStates(i) = ppval(RKIdata{i}.splineDead, parameterArray{1}.endDateOpt)/...
##    (ppval(RKIdata{i}.splineInfected, (parameterArray{1}.endDateOpt - 15))*...
##    parameterArray{1}.mortality);
##  end
##  parameterArray{1}.darkFigure = mean(darkFiguresStates);
##  darkFiguresStates /= mean(darkFiguresStates);
##  parameterArray{1}.factorsDF = darkFiguresStates;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %This parameter determines the period of time over which the scaling factors 
  %of the estimated dark figure are to be determined for the individual federal 
  %states. If you fit the first wave, you can use the last week of the fitting 
  %period. So you avoid that deviations on the last day of the fitting period are 
  %too significant. If you are fitting data that is not yet complete 
  %(RKI supplements are very likely), you should determine the factors 2 weeks 
  %before the end of the fitting
  startIndex = parameterArray{1}.startDate - cell2mat(RKIdata{1}.time(1)) + 1;
  endIndex = parameterArray{1}.endDateOpt - cell2mat(RKIdata{1}.time(1)) + 1;
  RT = endIndex - startIndex;
  dfState = zeros(1,16);
  if isfield(parameterArray{1},'timeForFactorsDF')
    if strcmp(parameterArray{1}.timeForFactorsDF, 'lastWeek')
      interval = endIndex-6:endIndex;
    elseif strcmp(parameterArray{1}.timeForFactorsDF, 'twoWeeksToSimEnd')
      interval = endIndex-20:endIndex-14;
    end
  else
    interval = endIndex-6:endIndex;
  end
  for j=1:16
    meanState = 0;
    counter = 0;
    %for i=endIndex-6:endIndex
    for i=interval(1):interval(end)
      meanState += cell2mat(RKIdata{j}.df(i));
      counter += 1;
    end
    RT = counter; 
    meanState /= RT;
    dfState(j) = meanState;
  end
  %parameterArray{1}.darkFigure = mean(dfState); 
  %calculate FactorsDF for sir_eqn_spatial
  factorsDF = dfState/mean(dfState);
  parameterArray{1}.factorsDF = factorsDF;
 
  ##  figure;
  ##  hold on;
  ##  population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
  ##  4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
  ##  1609675, 4077937, 2208321, 2143145];
  ##  for i=1:length(RKIdata)
  ##    semilogy(cell2mat(RKIdata{i}.time), (RKIdata{i}.infected/RKIdata{i}.infected(1))); 
  ##  end
  ##  a=[0.942866 1.256358 0.742975 1.030929 0.987411 0.636118 0.500751 1.119750 1.473729 1.659163 0.889826 0.977812 0.727012 1.210213 0.660300 0.838075];
  ##  figure;
  ##  hist(a);
  
  ##  times = cell2mat(RKIdata{1}.time);
  ##  rkiD = zeros(1,length(RKIdata{1}.time)*16);
  ##  rkiI = zeros(1,length(RKIdata{1}.time)*16);
  ##  %rkiI = ppval(RKIdata{1}.splineInfected, times);
  ##  for i=1:length(RKIdata)
  ##    rkiI(1, (i-1)*length(times)+1:i*length(times)) = ppval(RKIdata{i}.splineInfected, times);
  ##    rkiD(1, (i-1)*length(times)+1:i*length(times)) = ppval(RKIdata{i}.splineDead, times);
  ##  end

  if isfield(parameterArray{1},'paraNames')
    paraNames = parameterArray{1}.paraNames;
    x0Cov = parameterArray{1}.x0Cov;
    lb = parameterArray{1}.lb;
    ub = parameterArray{1}.ub;
  else 
  %setting paramters
##  paraNames = {"betaCorr1", "betaCorr2", "betaCorr3", "betaCorr4","betaCorr5",...
##  "betaCorr6", "betaCorr7", "betaCorr8", "betaCorr9", "betaCorr10",...
##  "betaCorr11", "betaCorr12", "betaCorr13", "betaCorr14", "betaCorr15", "betaCorr16"};
##  x0Cov = [0.834401 1.479429 0.622385 1.108438 1.139666 0.500897 0.298197 1.147868 1.612068 1.880193 0.772545 1.191203 0.598741 1.202886 0.56924 0.833693]; 
##  lb = zeros(16,1);
##  ub = lb+0.5;  

    paraNames = {"beta", "majorEvents", "contactRestrictions", "darkFigure"};
    x0Cov = [0.157 0.61 0.3 10]; 
    lb = [0 0 0 8];
    ub = [0.5 1 1 20];
  end
  %beta: 0.122582 | majorEvents: 0.811800 | contactRestrictions: 0.089347 | darkFigure: 12.323605 
  %0.122583 | majorEvents: 0.811847 | contactRestrictions: 0.089465 | darkFigure: 12.323155 | 

##  parameterArray{1}.beta =  0.122583;
##  parameterArray{1}.majorEvents = 0.811847;
##  parameterArray{1}.contactRestrictions = 0.089465;
##  parameterArray{1}.darkFigure = 12.323155;
##  parameterArray{1}.fitNewInfections = true;
  if isfield(parameterArray{1},'fitNewInfections')
  else 
    parameterArray{1}.fitNewInfections = false;
  end
  if isfield(parameterArray{1},'optAlgorithm')
    optAlgorithm = parameterArray{1}.optAlgorithm;
  else
    optAlgorithm = "justVisualize"; %lsqnonlin|pso
    parameterArray{1}.optAlgorithm = optAlgorithm;
  end
##  if isfield(parameterArray{1}, "paraNamesOpt")
##    paraNames = getfield(parameterArray{1}, "paraNamesOpt");
##    x0Cov = getfield(parameterArray{1}, "x0CovOpt");
##    lb = getfield(parameterArray{1}, "lbOpt");
##    ub = getfield(parameterArray{1}, "ubOpt");  
##    optAlgorithm = getfield(parameterArray{1}, "optAlgorithmOpt");  
##  endif 
  switch optAlgorithm
    case "lsqnonlin"
      parameterArray{1}.fullConsoleOut = false;
      [paramCovBest, resMin] = runLsqNonlin(x0Cov, lb, ub,...
      paraNames, parameterArray, RKIdata);
    case "pso"
      opt.parallel = true;
      opt.visu = false;
      opt.n_particles = 200;
      if isfield(parameterArray{1},'PSOiterations')
        opt.n_iter = parameterArray{1}.PSOiterations;
      else
        opt.n_iter = 150;
      end
      opt.coupling = 5;
      opt.parameter_names = paraNames;
      
      parameterArray{1}.fullConsoleOut = false;
      args = {paraNames, parameterArray, RKIdata};
      [paramCovBest, resMin] = pso(@evalSirLocalOneArgument, args, lb, ub, opt);
    case "justVisualize"
      resMin = inf;
      %parameterArray{1}.saveVTK = true;
      paraNames = {};
      x0Cov = [];
      paramCovBest = x0Cov;
  endswitch
  
  if or(nargin == 0, strcmp(optAlgorithm, "justVisualize"))
    %%%%%%%% Show the Results %%%%%%%%
    optCovid(paramCovBest, paraNames, parameterArray, RKIdata, true);
    ind=find(ismember(paraNames,'startDate'));
    dateInitial = paramCovBest(ind);
    date = paramCovBest(ind);
       
    %print final results console 
    for i=1:length(paraNames)
      if i==ind
        fprintf("%sInitial: %s %s: %s \n",  paraNames{i},...
        datestr(dateInitial), paraNames{i}, datestr(date));
      else
        fprintf("%sInitial: %f %s: %f \n",  paraNames{i},...
        x0Cov(i), paraNames{i}, paramCovBest(i));
      end
    end
    fprintf("\nresnormCovid: %f\nCalculation took %f seconds\n", resMin, toc(tOpt)); 
  end
  
  if nargout == 1 
    varargout = {paramCovBest};
  elseif nargout == 2
    varargout = {paramCovBest, resMin};
  end  
  
  %write final results to .txt file
  folder = ["../../Results/", parameterArray{1}.folderName];
  fid = fopen([folder, "/protokoll_final.txt"], "a");
  fprintf(fid, "\nresnormCovid: %f\nCalculation took %f seconds\n", resMin, toc(tOpt)); 
  fprintf(fid, "\n");
  for i=1:length(paraNames)
##    if i==ind
##      fprintf(fid, "%sInitial: %s %s: %s \n",  paraNames{i},...
##      datestr(dateInitial), paraNames{i}, datestr(date));
##    else
      fprintf(fid, "%sInitial: %f %s: %f \n",  paraNames{i},...
      x0Cov(i), paraNames{i}, paramCovBest(i));
##    end
  end
  fprintf(fid, "\n");
  fclose(fid);   
  
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paramCovBest, resMin] = runLsqNonlin(x0Cov, lb, ub,...
  paraNames, parameterArray, RKIdata)  
  %%%%%%%% Opt package and settings %%%%%%%%
  pkg load optim;
  opts = optimset();%("Algorithm", "levenberg-marquardt");
  ##  opts.TolFun = 1e-10;
  
  % Try to load parallel package but only if its not 4.0.0
  try
    [dummy, info] = pkg('list');
    findParallel = strcmp( cellfun( @(ind) ind.name, info, 'uni', false), {'parallel'} );
    indParallel = find(findParallel==1);
    version = info{indParallel}.version;
    pkg load parallel;
    if version == "4.0.0"
      pkg unload parallel;
      fprintf("unloaded parallel:4.0.0 package, because it's not supported yet.\n");
    end
  catch
    fprintf("No parallel package found or used, code will run serial.\n");
  end_try_catch
  
  % Allocate memory for the parallel runs
  nruns = parameterArray{1}.nRuns;
  nparallel = nproc();
  paramCovid = zeros(nruns, length(x0Cov));
  resNormCovid = inf(nruns, 1);
  runOptimOnce = @(x)lsqnonlinTwoRetvals(x, lb, ub, opts,...
  paraNames, parameterArray, RKIdata, false);  
  
  % Best and initial solution
  paramCovBest = x0Cov;
  resMin = inf(1);
  x0CovInit = x0Cov;
  runs = 0;
  nChunksPerCore = 4;
  while runs < nruns
    x0Array = {};
    for i=1:nChunksPerCore*nparallel
      if runs + i > nruns
        break;
      end    
      % first run with initial paramters, following runs with varied parameters
      x0Array{end+1} = x0Cov;
      x0Cov = paramCovBest .* ((rand(size(paramCovBest)) - 0.5) * 0.5 + 1);
      x0Cov(x0Cov>ub) = ub(x0Cov>ub);
      x0Cov(x0Cov<lb) = lb(x0Cov<lb);
    endfor
    if exist("parcellfun.m", "file")
      [a,b] = parcellfun (nparallel, runOptimOnce , x0Array, "UniformOutput", false);
    else
      [a,b] = arrayfun (runOptimOnce , x0Array, "UniformOutput", false);
    endif
    %% Save parameters to file
    folder = ["../../Results/", parameterArray{1}.folderName];
    fid = fopen([folder, "/protokoll.txt"], "a");
    for i=1:length(a)
      paramCovid(runs+i, :) = cell2mat(a(i));
      resNormCovid(runs+i) = cell2mat(b(i));    
      fprintf(fid, "Res: %f ", resNormCovid(runs+i));
      for j= 1:length(x0Cov)
        fprintf(fid, "%s: %f | ", paraNames{j}, paramCovid(runs+i, j));
      end
      fprintf(fid, "\n");
    end  
    fclose(fid);
    
    runs += length(b);
    indexBest = find((min(resNormCovid) == resNormCovid), 1);
    resMin = resNormCovid(indexBest);
    paramCovBest = paramCovid(indexBest, :);
  endwhile
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paraRun, resRun] = lsqnonlinTwoRetvals(x0Cov, lb, ub, opts,...
  paraNames, parameterArray, RKIdata, showPlots)  
  opt = @(x)optCovid(x, paraNames, parameterArray, RKIdata, showPlots); 
  try
    if iscell(x0Cov)
      x0Cov = cell2mat(x0Cov);
    end     
    [paraRun, resRun] = lsqnonlin(opt, x0Cov, lb, ub, opts);
  catch
    folder = ["../../Results/", parameterArray{1}.folderName];
    fid = fopen([folder, "/protokoll.txt"], "a");
    fprintf(fid, "Error ");
    for j= 1:length(x0Cov)
      fprintf(fid, "%s: %f | ", paraNames{j}, x0Cov(j));
    end
    fprintf(fid, "\n");
    fclose(fid);
    paraRun = inf(size(x0Cov));
    resRun = inf(1);
  end_try_catch  
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
