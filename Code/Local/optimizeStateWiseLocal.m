function retval = optimizeStateWiseLocal (parameterArray) %filelname 1 und filename 2
  tOpt = tic;
  addpath('../Spatial', '../Pso', '../Spatial/Optimize');
  %% Example  
  
  if nargin == 0
    close all;
    parameterArray = read_parameter("../../Results/parameter_local.txt");  
    %% the entries can be directly modified in the struct
    %% eg. set a foldername
      %parameter = read_parameter("../../Results/parameter_local.txt"){1};   
    parameterArray{1}.showDiagrams = false;
    parameterArray{1}.model = "SIRH"; %"SIREDmod" | "SIRH" 
    parameterArray{1}.optimizationMode = "Germany"; %"Counties"|"Germany"|"Germany-and-Counties"
    parameterArray{1}.folderName = "SIRH_Local_Fit-20200430_3";
    parameterArray{1}.startDate = datenum([2020, 03, 02]);
    parameterArray{1}.endDateOpt = datenum([2020, 04, 30]);
    parameterArray{1}.betaStatewise	= true;		
    parameterArray{1}.wRatioGlobStates = 0;%.5; % 1 is fully statewise, 0 is only global
    parameterArray{1}.wRatioID = 0.5; % 1 is only infected, 0 only death
    parameterArray{1}.wISlope = 0.1; %0.002 % weight for slope of infected last date
    parameterArray{1}.dayForDeathFit = "lastDay";
    parameterArray{1}.timeForFactorsDF = "lastWeek";

    
    %%%%%%%% parameters from lsqnonlin25runs (14.05.2020) %%%%%%%% 
    %%%%%%%% median from good residuals at similat start dates %%%%%%%%
    %parameter.gamma1 = 0.095196; 
    %parameter.gamma2 = 0.041280;
    %parameter.mortality = 0.0068125;
    %parameter.darkFigure = 7.6285;
    %parameter.startDate = 737839.81073;
    
    %%%%%%%% parameters from literature %%%%%%%%
    parameterArray{1}.mortality = 0.006;
    parameterArray{1}.gamma1 = 0.067; 
    parameterArray{1}.gamma2 = 0.04;
    parameterArray{1}.beta = 0.43;%1.4;
    parameterArray{1}.majorEvents = 0.64;
  endif
  
  %days können wir später dann in das Parameterarry mit aufnehmen
  %days beschreibt die Zeitspanne in der die Optimierung laufen soll
  %außerdem brauchen wir noch einen Parameter, der bestimmt, ob wir 
  %statewise oder für ganz Deutschland fitten
  
  mkdir(["../../Results/", parameterArray{1}.folderName]);
  
  RKIread = read_case_history_RKIfiles();
  
  %rewrite the cell, since states are arranged differently in RKIread and statewiseSIR
  RKIdata = cell(17,1);
  RKIdata{1} = RKIread{15};
  RKIdata{1}.population = 2896712;
  RKIdata{2} = RKIread{6};
  RKIdata{2}.population = 1841179;
  RKIdata{3} = RKIread{9};
  RKIdata{3}.population = 7982448;
  RKIdata{4} = RKIread{5};
  RKIdata{4}.population = 682986; 
  RKIdata{5} = RKIread{10};
  RKIdata{5}.population = 17932651;
  RKIdata{6} = RKIread{7};
  RKIdata{6}.population = 6265809;
  RKIdata{7} = RKIread{11};
  RKIdata{7}.population = 4084844;
  RKIdata{8} = RKIread{1};
  RKIdata{8}.population = 11069533;
  RKIdata{9} = RKIread{2};
  RKIdata{9}.population = 13076721; 
  RKIdata{10} = RKIread{12};
  RKIdata{10}.population = 990509;
  RKIdata{11} = RKIread{3};
  RKIdata{11}.population = 3644826;
  RKIdata{12} = RKIread{4};
  RKIdata{12}.population = 2511917;
  RKIdata{13} = RKIread{8};
  RKIdata{13}.population = 1609675;
  RKIdata{14} = RKIread{13};
  RKIdata{14}.population = 4077937;
  RKIdata{15} = RKIread{14};
  RKIdata{15}.population = 2208321; 
  RKIdata{16} = RKIread{16};
  RKIdata{16}.population = 2143145;
  
  timesRKI = cell2mat(RKIread{1}.time);
  
  for i=1:length(RKIdata)-1
    %RKIinf = RKIdata{i}.infected;
    RKIdata{i}.infected = cell2mat(RKIdata{i}.infected);
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end
  
  RKIinfGermany = zeros(size(timesRKI));
  RKIdeadGermany = zeros(size(timesRKI));
  for i=1:length(RKIread)
    RKIinfGermany += cell2mat(RKIread{i}.infected);
    RKIdeadGermany += cell2mat(RKIread{i}.dead);
  end
  RKIdata{17}.infected = RKIinfGermany;
  RKIdata{17}.dead = RKIdeadGermany;
  RKIdata{17}.time = mat2cell(timesRKI, 1);
  RKIdata{17}.population = 83019213;
  RKIdata{17}.name = 'Germany';
  RKIdata{17}.splineInfected = spline(timesRKI, RKIdata{17}.infected);
  RKIdata{17}.splineDead = spline(timesRKI, RKIdata{17}.dead);
  
  %cumulate infection numbers 
  %create splines
  
  ##    darkFiguresStates = zeros(1, length(RKIdata));
  ##    for i=1:length(RKIdata)
  ##      darkFiguresStates(i) = ppval(RKIdata{i}.splineDead, parameterArray{1}.endDateOpt)/...
  ##      ppval(RKIdata{i}.splineInfected, parameterArray{1}.endDateOpt - 15);
  ##    end
  ##    darkFiguresStates /= parameterArray{1}.mortality;
  
  %calculation of DF factors for states
##  darkFiguresStates = zeros(1, length(RKIdata));
##  for i=1:length(RKIdata)
##    darkFiguresStates(i) = RKIdata{i}.dead(parameterArray{1}.endDateOpt-parameterArray{1}.startDate+1)/...
##    (RKIdata{i}.infected(parameterArray{1}.endDateOpt-parameterArray{1}.startDate-12)*...
##    parameterArray{1}.mortality);
##  end
##  darkFiguresStates /= mean(darkFiguresStates);
##  
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
  
  %setting paramters
  ##  paraNames = {"betaCorr1", "betaCorr2", "betaCorr3", "betaCorr4","betaCorr5",...
  ##  "betaCorr6", "betaCorr7", "betaCorr8", "betaCorr9", "betaCorr10",...
  ##  "betaCorr11", "betaCorr12", "betaCorr13", "betaCorr14", "betaCorr15", "betaCorr16"};
  ##  x0Cov = [0.834401 1.479429 0.622385 1.108438 1.139666 0.500897 0.298197 1.147868 1.612068 1.880193 0.772545 1.191203 0.598741 1.202886 0.56924 0.833693]; 
  
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
  parameterArray{1}.darkFigure = mean(dfState); 
  %calculate FactorsDF for sir_eqn_spatial
  factorsDF = dfState/mean(dfState);
  parameterArray{1}.factorsDF = factorsDF;
  
  paraNames = {"beta", "majorEvents", "contactRestrictions","darkFigure"};
  lb = [1 0 0.49 7]; %0.3
  ub = [6 1 1 16];  
  x0Cov = (lb+ub)/2;  
  
  %%%%%%%% optimizationModes %%%%%%%%
  switch parameterArray{1}.optimizationMode
    case "Germany"
      optMode = 17:17;
    case "Counties"
      optMode = 1:16;
    case "Germany-and-Counties"
      optMode = 1:17;
  endswitch
  
  parameterArray{1}.endoptMode = optMode(end);
  
  %%%%%%%% fit data for each state %%%%%%%%
  for k=optMode(1):optMode(end)  
    parameterArray{1}.AGS_state = k;  
    N = RKIdata{k}.population;
    if parameterArray{1}.AGS_state == 17
      parameterArray{1}.betaStatewise = false; 
    endif
    
    
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
        paraNames, parameterArray, RKIdata{k});
      case "pso"
        opt.parallel = true;
        opt.visu = true;
        opt.n_particles = 200;
        opt.n_iter = 100;
        opt.coupling = 5;
        opt.parameter_names = paraNames;
        
        parameterArray{1}.fullConsoleOut = false;
        args = {paraNames, parameterArray, RKIdata{k}};
        [paramCovBest, resMin] = pso(@evalSirLocalOneArgumentLocal, args, lb, ub, opt);
      case "justVisualize"
        resMin = inf;
        %parameterArray{1}.saveVTK = true;
        paraNames = {};
        x0Cov = [];
        paramCovBest = x0Cov;
    endswitch
    
    if or(nargin == 0, strcmp(optAlgorithm, "justVisualize"))
      %%%%%%%% Show the Results %%%%%%%%
      optCovidLocal(paramCovBest, paraNames, parameterArray, RKIdata{k}, true);
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
  endfor
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paramCovBest, resMin] = runLsqNonlin(x0Cov, lb, ub,...
  paraNames, parameterArray, RKIdata)  
  %%%%%%%% Opt package and settings %%%%%%%%
  pkg load optim;
  opts = optimset();%("Algorithm", "levenberg-marquardt");
  ##  opts.TolFun = 1e-10;
  paramCovBest = x0Cov;
  x0CovInit = x0Cov;
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
##          runOptimOnce = @(x)lsqnonlinTwoRetvals(x, lb, ub, opts,...
##        paraNames, parameter, RKIdata{k}, false, false);  
  
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
    fprintf(fid, "%scounty: %s", RKIdata.name);
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
  opt = @(x)optCovidLocal(x, paraNames, parameterArray, RKIdata, showPlots); 
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