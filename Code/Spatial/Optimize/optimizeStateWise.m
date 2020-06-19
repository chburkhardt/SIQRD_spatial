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
	if nargin == 0
  	close all;
		parameterArray = read_parameter('../../Results/parameter.txt');  
		%% the entries can be directly modified in the struct
		%% eg. set a foldername
		parameterArray{1}.model = "SIREDmod";
		%parameterArray{1}.totalRuntime = 40;
		parameterArray{1}.folderName = "Test_Plots";
		parameterArray{1}.saveVTK = false;
		parameterArray{1}.showDiagrams = false;
		parameterArray{1}.fullConsoleOut = true;
		parameterArray{1}.blowUp = false;
		parameterArray{1}.spatial	= "germany";
		parameterArray{1}.initial = "GitData";
		parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
		parameterArray{1}.startDate = datenum([2020, 03, 02]);
		parameterArray{1}.endDateOpt = datenum([2020, 04, 25]);
		parameterArray{1}.reduceToStates= true;
		parameterArray{1}.betaStateWise= false;
		parameterArray{1}.exitRestrictions= 1;
		parameterArray{1}.schoolClosing= 1;
		parameterArray{1}.beta_cross_county	= 1;
		
		parameterArray{1}.wRatioGlobStates = 0;%.5; % 1 is fully statewise, 0 is only global
		parameterArray{1}.wRatioID = 1; % 1 is only infected, 0 only death
    parameterArray{1}.wISlope = 0.002; % weight for slope of infected last date
    parameterArray{1}.optFunGermanSum = true;
		parameterArray{1}.showStatePlots = false;
		
		parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
		parameterArray{1}.gamma2 = 0.04; % Quarantine -> Recovered
		parameterArray{1}.mortality = 0.006; % 0.006
		parameterArray{1}.darkFigure = 6;
		parameterArray{1}.beta =  0.8;
		parameterArray{1}.betaSWscaling =  2;
		
		parameterArray{1}.majorEvents	= 0.37;
		parameterArray{1}.contactRestrictions = 0.188669;  
		
	endif
	  
  RKIread = read_case_history();
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
    RKIdata{i}.infected = cumsum(cell2mat(RKIdata{i}.infected));
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end
	  
  %setting paramters
  paraNames = {"betaCorr1", "betaCorr2", "betaCorr3", "betaCorr4","betaCorr5",...
	"betaCorr6", "betaCorr7", "betaCorr8", "betaCorr9", "betaCorr10",...
	"betaCorr11", "betaCorr12", "betaCorr13", "betaCorr14", "betaCorr15", "betaCorr16"};
  x0Cov = [0.834401 1.479429 0.622385 1.108438 1.139666 0.500897 0.298197 1.147868 1.612068 1.880193 0.772545 1.191203 0.598741 1.202886 0.56924 0.833693];	
	
  paraNames = {"beta", "majorEvents", "contactRestrictions"};
  x0Cov = [0.157 0.61 0.5]; 
  variations = ones(1, length(paraNames)) * 0.2;
  variations = [0.2 0.2 0.7];
  lb = x0Cov .* (1 - variations);
  ub = x0Cov .* (1 + variations);
	
  optAlgorithm = "justVisualize"; %lsqnonlin|pso
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
      opt.n_particles = 300;
      opt.n_iter = 300;
      opt.coupling = 5;
			opt.parameter_names = paraNames;
			
			parameterArray{1}.fullConsoleOut = false;
      evalSirLocalOneArgument= @(x) norm(optCovid(x, paraNames,...
      parameterArray, RKIdata, false));
      [paramCovBest, resMin] = pso(evalSirLocalOneArgument, lb, ub, opt);
    case "justVisualize"
      resMin = inf;
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
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paramCovBest, resMin] = runLsqNonlin(x0Cov, lb, ub,...
	paraNames, parameterArray, RKIdata)  
	%%%%%%%% Opt package and settings %%%%%%%%
	pkg load optim;
	opts = optimset();%("Algorithm", "levenberg-marquardt");
	
	% Try to load parallel package but only if its not 4.0.0
	try
		[dummy, info] = pkg('list');
		findParallel = strcmp( cellfun( @(ind) ind.name, info, 'uni', false), {'parallel'} );
		indParallel = find(findParallel==1);
		version = info{indParallel}.version;
		pkg load parallel;
		if version == "4.0.0"
			pkg unload parallel;
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
function retval = optCovid(xCov, paraNames, parameterArray, RKIdata, showPlots)
	%% this is how a parameter struct can be loaded from a parameterfile
	% alle Pfade sind so gesetzt das man im "Spatial" ordner starten muss,
	% daher auch hier  
	if iscell(xCov)
		xCov = xCov{1};
	endif
	%vary paramters
	for i=1:length(paraNames)
		parameterArray{1} = setfield(parameterArray{1}, paraNames{i}, xCov(i));
	end
	
	% Modell soll bis zum xx.yy.2020 laufen
	parameterArray{1}.totalRuntime =...
	parameterArray{1}.endDateOpt - parameterArray{1}.startDate;
	startDateSim = parameterArray{1}.startDate; 
	
	%% run code and store the workspace to avoid io
	workspaceCalculation = sir_spatial(parameterArray);
	
	%% load results statewise with ne new option to avoid io
	stateWiseSIR = extractStatewiseResults (workspaceCalculation, "workspace",...
	["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
	
	##  for i=1:length(paraNames)
	##    fprintf("%s= %f\n", paraNames{i}, xCov(i));
	##  end
	
	timespan = [max([min(cell2mat(RKIdata{1}.time)), min(stateWiseSIR{1}.time + startDateSim)]),...
	min([max(cell2mat(RKIdata{1}.time)), max(stateWiseSIR{1}.time + startDateSim)])];
	times = linspace(timespan(1), timespan(2), 100); 
	
  darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1});	
  for i=1:length(stateWiseSIR)
    %%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
    % discovered = \int_0^t alpha * E * dT 
    alpha = parameterArray{1}.gamma1 / (darkFigures(i) - 1);
    splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
    % Integriere E*alpha to get discovered (realy infected)
    get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
    [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
    max(stateWiseSIR{i}.time)], stateWiseSIR{i}.SIR_vs_time(1,2));
    % get solutionDiscovered at the same timesteps as the other solutions
    solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);
    
    stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
    solutionDiscovered);    
    stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,5));
    stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,2));
    stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,3));
    stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,4));
  end
	
	%write vector that contains data for each state 
	simdataI = zeros(1,length(times)*length(stateWiseSIR));
	simdataIslope = zeros(1,length(stateWiseSIR));
	simdataD = zeros(1,length(times(end))*length(stateWiseSIR));
	
	wIstatewise = zeros(1,length(times)*length(stateWiseSIR));
	wDstatewise = zeros(1,length(times(end))*length(stateWiseSIR));
	
	rkiI = zeros(1,length(times)*length(stateWiseSIR));
	rkiIslope = zeros(1,length(stateWiseSIR));
	rkiD = zeros(1,length(times(end))*length(stateWiseSIR));  
	
	wwI = parameterArray{1}.wRatioID;
	wwD = 1 - parameterArray{1}.wRatioID;
  wwISlope = parameterArray{1}.wISlope;
	for i=1:length(stateWiseSIR)
		simdataIstate = ppval(stateWiseSIR{i}.splineInfected, times);
		simdataI(1, (i-1)*length(times)+1:i*length(times)) = simdataIstate;
    simdataIslope(1, i) = diff(simdataIstate)(end);
		
		simdataDstate = ppval(stateWiseSIR{i}.splineDead, times(end));
		simdataD(1, (i-1)*length(times(end))+1:i*length(times(end))) = simdataDstate;
		
		rkiIstate = ppval(RKIdata{i}.splineInfected, times);
		rkiI(1, (i-1)*length(times)+1:i*length(times)) = rkiIstate;
    rkiIslope(1, i) = diff(rkiIstate)(end);
		wIstatewise(1, (i-1)*length(times)+1:i*length(times)) = wwI./max(rkiIstate);
		
		rkiDstate = ppval(RKIdata{i}.splineDead, times(end));
		rkiD(1, (i-1)*length(times(end))+1:i*length(times(end))) = rkiDstate;
		wDstatewise(1, (i-1)*length(times(end))+1:i*length(times(end))) = wwD./max(rkiDstate); 
	end
	
	wIglob = wwI/max(rkiI);
	wDglob = wwD/max(rkiD) * length(rkiI) / length(rkiD);
  wIslope = wwISlope/max(rkiIslope) * length(rkiI) / length(rkiIslope);
	
	wDstatewise *= wwD * length(rkiI) / length(rkiD);
	
	ratio = parameterArray{1}.wRatioGlobStates;
	wI = wIglob * (1 - ratio) + wIstatewise * ratio;
	wD = wDglob * (1 - ratio) + wDstatewise * ratio;
	
	%opt function
  if parameterArray{1}.optFunGermanSum 
    % sum germany
    retval = [sum(reshape((simdataI - rkiI) .* wI, length(simdataI)/length(stateWiseSIR),...
    length(stateWiseSIR))', 1),...
    sum(reshape((simdataD - rkiD) .* wD, length(simdataD)/length(stateWiseSIR),...
    length(stateWiseSIR)), 2),...
    sum(reshape((simdataIslope -rkiIslope) .* wIslope,...
    length(simdataIslope)/length(stateWiseSIR), length(stateWiseSIR)), 2)];
  else
    retval = [(simdataI - rkiI) .* wI, (simdataD - rkiD) .* wD,...
    (simdataIslope -rkiIslope) .* wIslope];
	end
	
	folder = ["../../Results/", parameterArray{1}.folderName];
	fid = fopen([folder, "/protokoll_lsqnonlin.txt"], "a");
	fprintf(fid, "Res: %f ", norm(retval));
	for j= 1:length(xCov)
		fprintf(fid, "%s: %f | ", paraNames{j}, xCov(j));
	end
	fprintf(fid, "\n");
	fclose(fid);  
	
	if showPlots
		%sum up all states
		rkiIsum =  sum(reshape(rkiI, length(times), length(stateWiseSIR)), 2);
		rkiDsum =  sum(reshape(rkiD, length(times(end)), length(stateWiseSIR)), 2);
		figure(17);
		clf;
		subplot(1,2,1), hold on;
		subplot(1,2,2), hold on;
		subplot(1,2,1), plot(times, rkiIsum, "-x");
		subplot(1,2,2), plot(times(end), rkiDsum, "-x");
		namesI = {"RKI infected (cumulated)"};
		namesD = {"RKI dead"};
		
    simdataIsum =  sum(reshape(simdataI, length(times), length(stateWiseSIR)), 2);
    subplot(1,2,1), plot(times, simdataIsum, "-o");
    namesI{end + 1} = "Cummulated Infected";
    simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2);
    subplot(1,2,2), plot(times(end), simdataDsum, "-o");
    namesD{end + 1} = "Dead";
       
		subplot(1,2,1), xlim([min(times), max(times)]);
		subplot(1,2,1), ylim([0, max(rkiIsum)*1.5]);
		subplot(1,2,1), legend(namesI);
		subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
		subplot(1,2,2), xlim([min(times), max(times)]);
		subplot(1,2,2), ylim([0, max(rkiDsum)*1.5]);
		subplot(1,2,2), legend(namesD);
		subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
		title("Germany");
		fname = ["../../Results/", parameterArray{1}.folderName];
		filename = "plot_Germany";
		drawnow;
		pause(0.3);
		saveas(gcf, fullfile(fname, filename), 'jpeg');
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		simdataI2 = zeros(1,length(times)*length(stateWiseSIR));
		simdataR = zeros(1,length(times)*length(stateWiseSIR));
		simdataE = zeros(1,length(times)*length(stateWiseSIR));
		for i=1:length(stateWiseSIR)
			simdataI2(1, (i-1)*length(times)+1:i*length(times)) = ppval(stateWiseSIR{i}.splineI, times);
			simdataR(1, (i-1)*length(times)+1:i*length(times)) =...
			ppval(stateWiseSIR{i}.splineR, times);
			simdataE(1, (i-1)*length(times)+1:i*length(times)) =...
			ppval(stateWiseSIR{i}.splineE, times);
		end
		simdataIplot =  reshape(simdataI, length(times), length(stateWiseSIR));
		simdataI2plot = reshape(simdataI2, length(times), length(stateWiseSIR));
		simdataRplot = reshape(simdataR, length(times), length(stateWiseSIR));
		simdataEplot = reshape(simdataE, length(times), length(stateWiseSIR));
		simdataDplot =  reshape(simdataD, length(times(end)), length(stateWiseSIR));
		
		%statewise plots
		for i=1:length(RKIdata)
			##    population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
			##    4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
			##    1609675, 4077937, 2208321, 2143145];
			if !parameterArray{1}.showStatePlots
				break;
			end
			rkiIplot =  reshape(rkiI, length(times), length(stateWiseSIR));
			rkiDplot =  reshape(rkiD, length(times(end)), length(stateWiseSIR));
			
			caseFatilityRateRKI = rkiDplot(end, i) / rkiIplot(end, i);
			
			figNum = figure;
			subplot(1,2,1), hold on;
			subplot(1,2,2), hold on;
			subplot(1,2,1), plot(times, rkiIplot(:,i), "-x");
			subplot(1,2,2), plot(times(end), rkiDplot(:,i), "-x");
			namesI = {"RKI infected (cumulated)"};
			namesD = {"RKI dead"};

      subplot(1,2,1), plot(times, simdataIplot(:,i), "-o");
      namesI{end + 1} = "Cummulated Infected";
      subplot(1,2,1), plot(times, simdataI2plot(:,i));
      namesI{end + 1} = "Infected (with symptoms)";
      subplot(1,2,1), plot(times, simdataRplot(:,i));
      namesI{end + 1} = "Recovered";
      subplot(1,2,1), plot(times, simdataEplot(:,i));
      namesI{end + 1} = "Exposed + no symptoms";     
      subplot(1,2,2), plot(times(end), simdataDplot(:,i), "-o");
      namesD{end + 1} = "Dead";

			subplot(1,2,1), xlim([min(times), max(times)]);
			subplot(1,2,1), ylim([0, max(rkiIplot(:,i))*1.5]);
			subplot(1,2,1), legend(namesI);
			subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
			subplot(1,2,2), xlim([min(times), max(times)]);
			subplot(1,2,2), ylim([0, max(max(rkiDplot(:,i))*1.5 ,1)]);
			subplot(1,2,2), legend(namesD);
			subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
			
			caseFatilityRateSim = simdataDplot(end, i)/simdataIplot(end, i);
			fprintf("%f, ", rkiIplot(end, i)/simdataIplot(end, i));
			
			name = [RKIdata{i}.name, " df_error: ",...
			num2str(caseFatilityRateRKI / caseFatilityRateSim),...
			" Ratio rkiI/simI: ", num2str(rkiIplot(end, i)/simdataIplot(end, i))];
			set(figNum, 'Name', name);
			fname = ["../../Results/", parameterArray{1}.folderName];
			filename = ["plot_", RKIdata{i}.name];
			drawnow;
			pause(0.3);
			saveas(gcf, fullfile(fname, filename), 'jpeg');
		endfor
	endif
endfunction
