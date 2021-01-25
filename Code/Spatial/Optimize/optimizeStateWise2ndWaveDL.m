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

function varargout = optimizeStateWise2ndWaveDL (parameterArray) %filelname 1 und filename 2
  tOpt = tic;
  addpath('../', '../Pso');
  %% Example  
  
  if nargin == 0
    close all;
    parameterArray = read_parameter('../../Results/parameter.txt');  
    %% the entries can be directly modified in the struct
    %% eg. set a foldername
    parameterArray{1}.model = "SIRH";
    %parameterArray{1}.totalRuntime = 40;
    parameterArray{1}.folderName = "SIREDmodFitting_Beta_Test3";
    parameterArray{1}.saveVTK = false;
    parameterArray{1}.saveDiagrams = false;
    parameterArray{1}.showDiagrams = false;
    parameterArray{1}.fullConsoleOut = true;
    parameterArray{1}.blowUp = false;
    parameterArray{1}.spatial = "germany";
    parameterArray{1}.initial = "RKIfiles";
    parameterArray{1}.initalDistributionDate = datenum([2020, 10, 09]);
    parameterArray{1}.startDate = datenum([2020, 10, 08]);
    parameterArray{1}.endDateOpt = datenum([2020, 10, 29]);    
    parameterArray{1}.reduceToStates= true;
    parameterArray{1}.betaStateWise= false;
    parameterArray{1}.exitRestrictions= 1;
    parameterArray{1}.schoolClosing= 1;
    parameterArray{1}.beta_cross_county = 1;
    
    parameterArray{1}.wRatioGlobStates = 0; % 1 is fully statewise, 0 is only global
    parameterArray{1}.wRatioID = 0.5; % 1 is only infected, 0 only death
    parameterArray{1}.wISlope = 0.05; % weight for slope of infected last date
    parameterArray{1}.optFunGermanSum = false;
    parameterArray{1}.showStatePlots = false;
    
    parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
    parameterArray{1}.gamma2 = 0.04; % Quarantine -> Recovered
    parameterArray{1}.mortality = 0.006; % 0.006
    parameterArray{1}.darkFigure = 6.5; %9.95; (new)
    parameterArray{1}.beta =  0.17;
      
    
    parameterArray{1}.beta = 2.5;%0.629793;%2;%2.5;%(SIRHtesting), 0.2;%SIREDmod testing
    parameterArray{1}.betaSWscaling = 1;    
    %parameterArray{1}.majorEvents = 0.2102;
    %parameterArray{1}.contactRestrictions = 0.6096;
    
    %parameterArray{1}.reduceToStates = false;
    parameterArray{1}.showStatePlots = true;
    
##    parameterArray{1}.showStatePlots = false;
    ##  ##  % For the countywise model
    ##    parameterArray{1}.reduceToStates= false;
    ##    parameterArray{1}.beta_cross_county =  0.87869;%0.436106125; % 1.08 ohne ss to cs anpassung
  endif
 
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
  for i=1:length(RKIdata)
    %RKIinf = RKIdata{i}.infected;
    RKIdata{i}.infected = cell2mat(RKIdata{i}.infected);
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
    
    %calculate average number of new infections over the last seven days
    RKIdata{i}.newInfections = [diff(RKIdata{i}.infected(1:2),1,2),...
    diff(RKIdata{i}.infected,1,2)];
    %RKIdata{i}.newInfections = movmean(RKIdata{i}.newInfections, [3 3]);
    RKIdata{i}.splineNewInfections = spline(timesRKI, RKIdata{i}.newInfections);
  end
  
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



  %calculate mean of darkFiguresStates over fitting period
  startIndex = parameterArray{1}.startDate - cell2mat(RKIdata{1}.time(1)) + 1;
  endIndex = parameterArray{1}.endDateOpt - cell2mat(RKIdata{1}.time(1)) + 1;
  RT = endIndex - startIndex;
  dfState = zeros(1,16);
  for j=1:16
    meanState = 0;
    for i=startIndex:endIndex
      meanState += cell2mat(RKIdata{j}.df(i));
    end
    meanState /= RT;
    dfState(j) = meanState;
  end
  parameterArray{1}.darkFigure = mean(dfState); 
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
  
  %setting paramters
##  paraNames = {"betaCorr1", "betaCorr2", "betaCorr3", "betaCorr4","betaCorr5",...
##  "betaCorr6", "betaCorr7", "betaCorr8", "betaCorr9", "betaCorr10",...
##  "betaCorr11", "betaCorr12", "betaCorr13", "betaCorr14", "betaCorr15", "betaCorr16"};
##  x0Cov = [0.834401 1.479429 0.622385 1.108438 1.139666 0.500897 0.298197 1.147868 1.612068 1.880193 0.772545 1.191203 0.598741 1.202886 0.56924 0.833693]; 
  
  paraNames = {"beta"};%, "easingContactRestrictions", "endOfVacation", "liftingTravelRestictionsEU"};
  x0Cov = [0.1];% 1 1 1]; 
##  variations = ones(1, length(paraNames)) * 0.2;
##  variations = [0.2 0.2 0.3];
##  lb = x0Cov .* (1 - variations);
##  ub = x0Cov .* (1 + variations);
  lb = [1.5];% 0.8 0.8 0.8];
  ub = [5.5];% 5 5 5];
  
  parameterArray{1}.fitNewInfections = false; %newInfections|Cumulated
  optAlgorithm = "pso"; %lsqnonlin|pso
  if isfield(parameterArray{1}, "paraNamesOpt")
    paraNames = getfield(parameterArray{1}, "paraNamesOpt");
    x0Cov = getfield(parameterArray{1}, "x0CovOpt");
    lb = getfield(parameterArray{1}, "lbOpt");
    ub = getfield(parameterArray{1}, "ubOpt");  
    optAlgorithm = getfield(parameterArray{1}, "optAlgorithmOpt");  
  endif 
  switch optAlgorithm
    case "lsqnonlin"
      parameterArray{1}.fullConsoleOut = false;
      [paramCovBest, resMin] = runLsqNonlin(x0Cov, lb, ub,...
      paraNames, parameterArray, RKIdata);
    case "pso"
      opt.parallel = true;
      opt.visu = true;
      opt.n_particles = 20;
      opt.n_iter = 100;
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
    varargout = paramCovBest;
  elseif nargout == 2
    varargout = {paramCovBest, resMin};
  end  
 
  %write final results to .txt file
  folder = ["../../Results/", parameterArray{1}.folderName];
  fid = fopen([folder, "/protokoll_final.txt"], "a");
  fprintf(fid, "\nresnormCovid: %f\nCalculation took %f seconds\n", resMin, toc(tOpt)); 
  fprintf(fid, "\n");
  for i=1:length(paraNames)
    if i==ind
      fprintf(fid, "%sInitial: %s %s: %s \n",  paraNames{i},...
      datestr(dateInitial), paraNames{i}, datestr(date));
    else
      fprintf(fid, "%sInitial: %f %s: %f \n",  paraNames{i},...
      x0Cov(i), paraNames{i}, paramCovBest(i));
    end
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
