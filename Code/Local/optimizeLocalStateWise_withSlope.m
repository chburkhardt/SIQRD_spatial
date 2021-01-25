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
## Created: 2020-05-04

function retval = optimizeLocalStateWise (input1, input2)
  addpath("../Spatial", "../Spatial/Optimize", "../Pso");
  tic; 
  %%%%%%%% Parameter %%%%%%%%
  parameter = read_parameter("../../Results/parameter_local.txt"){1};   
  parameter.showDiagrams = false;
  parameter.model = "SIR"; %"SIREDmod" | "SIRH" 
  parameter.optimizationMode = "Germany"; %"Counties"|"Germany"|"Germany-and-Counties"
  parameter.folderName = "2020-06-04_PSO_slopeFit";
  parameter.startDate = datenum([2020, 03, 02]);
  parameter.betaStatewise	= true;		
  parameter.wRatioGlobStates = 0;%.5; % 1 is fully statewise, 0 is only global
  parameter.wRatioID = 1; % 1 is only infected, 0 only death
  parameter.wISlope = 0.05; %0.002 % weight for slope of infected last date
  
  %%%%%%%% parameters from lsqnonlin25runs (14.05.2020) %%%%%%%% 
  %%%%%%%% median from good residuals at similat start dates %%%%%%%%
  %parameter.gamma1 = 0.095196; 
  %parameter.gamma2 = 0.041280;
  %parameter.mortality = 0.0068125;
  %parameter.darkFigure = 7.6285;
  %parameter.startDate = 737839.81073;
  
  %%%%%%%% parameters from literature %%%%%%%%
  parameter.mortality = 0.006;
  parameter.gamma1 = 0.067; 
  parameter.gamma2 = 0.04;
  mkdir(["../../Results/", parameter.folderName]);
  
  %%%%%%%% RKI data %%%%%%%%
  RKIread = read_case_history_RKIfiles();
  
  %%%%%%%% rearrange RKI data so that they are sorted by AGS %%%%%%%%
  %%%%%%%% source for population: Wikipedia, 2018 %%%%%%%%
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

  %%%%%%%% germany %%%%%%%%
  timesRKI = cell2mat(RKIread{1}.time);
  RKIinfGermany = zeros(size(timesRKI));
  RKIdeadGermany = zeros(size(timesRKI));
  
  % cumulate RKI data from states and build splines
  for i=1:length(RKIread)
    RKIinfGermany += cell2mat(RKIread{i}.infected);
    RKIdeadGermany += cell2mat(RKIread{i}.dead);
  end
  RKIdata{17}.infected = RKIinfGermany;
  RKIdata{17}.dead = RKIdeadGermany;
  RKIdata{17}.time = mat2cell(timesRKI, 1);
  RKIdata{17}.population = 8e7;
  RKIdata{17}.name = 'Germany';
  RKIdata{17}.infected = cumsum(RKIdata{17}.infected);
  RKIdata{17}.splineInfected = spline(timesRKI, RKIdata{17}.infected);
  RKIdata{17}.splineDead = spline(timesRKI, RKIdata{17}.dead);

  %%%%%%%% build splines for each state %%%%%%%%
  for i=1:length(RKIdata)-1
    %RKIdata{i}.infected = cumsum(cell2mat(RKIdata{i}.infected));
    RKIdata{i}.infected = cell2mat(RKIdata{i}.infected);
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end   
  
  maxRKIIglob = 0; 
  maxRKII = 0;
  maxRKID = 0;
  maxRKIDglob = 0;
  maxRKIIslope = 0;
  maxRKIIslopeglob = 0;

  for i=1:length(RKIdata)-1
    maxRKII = max(RKIdata{i}.infected);
    maxRKID = max(RKIdata{i}.dead);
    maxRKIIslope = diff(RKIdata{i}.infected)(end);
    if maxRKII > maxRKIIglob
      maxRKIIglob = maxRKII;
    end
    if maxRKID > maxRKIDglob
      maxRKIDglob = maxRKID;
    end
    if maxRKIIslope > maxRKIIslopeglob
      maxRKIIslopeglob = maxRKIIslope;
    end
  end
  
  wwI = parameter.wRatioID;
	wwD = 1 - parameter.wRatioID;
  wwISlope = parameter.wISlope;
  
  wIglob = wwI/maxRKIIglob;
	wDglob = wwD/maxRKIDglob * 100;
  
  parameter.wIglob = wIglob;
  parameter.wDglob = wDglob;
  parameter.wISlopeGlob = maxRKIIslopeglob;
 
  %%%%%%%% optimizationModes %%%%%%%%
  switch parameter.optimizationMode
    case "Germany"
      optMode = 17:17;
    case "Counties"
      optMode = 1:16;
    case "Germany-and-Counties"
      optMode = 1:17;
  endswitch
  
  parameter.endoptMode = optMode(end);
  
  %%%%%%%% fit data for each state %%%%%%%%
  for k=optMode(1):optMode(end)  
    parameter.AGS_state = k;  
    N = RKIdata{k}.population;
    if parameter.AGS_state == 17
      parameter.betaStatewise = false; 
    endif
    
##    switch parameter.model
##      case "SIRH"
##        I0 = RKIdata{k}.infected(1);      
##        R0 = 0; 
##        S0 = N - I0 - R0;   
##        X0sir = S0 / N; 
##        
##        paraNames = {"beta", "startDate", "majorEvents", "schoolClosing", "contactRestrictions", "exitRestrictions"};
##        x0Cov = [0.47, datenum([2020, 02, 19]), 0.47, 0.66, 0.38, 1.52];
##        lb = [0.4, datenum([2020, 02, 1]), 0.1, 0.1, 0.1, 0.1];
##        ub = [0.65, datenum([2020, 03, 1]), 1, 2, 2, 4];   
##      case "SIREDmod"
##				I0 = ppval(RKIdata{k}.splineInfected, parameter.startDate);
##				D0 = 0;
##				E0 = ppval(RKIdata{k}.splineInfected, parameter.startDate)*...
##				(parameter.beta*(parameter.darkFigure - 1)/parameter.gamma1 - parameter.darkFigure);
##				R0 = 0; 
##				S0 = N - I0 - R0 - E0 - D0;
##				X0sir = [S0, I0, R0, E0, D0] / N;     
        
        ##      paraNames = {"beta", "gamma1", "gamma2", "mortality", "darkFigure", "startDate", "majorEvents", "schoolClosing", "contactRestrictions", "exitRestrictions"};
        ##      x0Cov = [0.441, 0.0236, 0.0028, 0.008, 8, datenum([2020, 02, 24]), 0.58, 0.64, 0.77, 0.3];
        ##      lb = [0.3, 0.01, 0.002, 0.003, 3, datenum([2020, 02, 15]), 0.1, 0.5, 0.1, 0.001];
        ##      ub = [0.5, 0.05, 0.01, 0.015, 30, datenum([2020, 03, 2]), 1, 2, 0.89, 0.45];      
        
##        paraNames = {"beta", "gamma1", "gamma2", "mortality", "darkFigure", "startDate",...
##        "majorEvents", "schoolClosing", "contactRestrictions"};
##        ##      x0Cov = [0.411, 0.069, 0.0023, 0.01, 37.1, datenum([2020, 02, 13]),...
##        ##      0.5, 0.48, 0.35]; % funktioniert!!! alt andi
##        x0Cov = [0.411061, 0.069801, 0.02302, 0.009986, 5.7, datenum([2020, 02, 13]),...
##        0.497353, 0.488051, 0.345821]; % funktioniert!!! dominik
##        lb = [0.36, 0.001, 0.0001, 0.005, 5, datenum([2020, 02, 15]),...
##        0.3, 0.25, 0.3];
##        ub = [0.50, 0.12, 0.12, 0.01, 10, datenum([2020, 02, 23]),...
##        1, 1, 1];  

        %gamma1, gamma2, mortality, darkFigure, startDate festgesetzt
##        paraNames = {"beta", "gamma1", "gamma2", "mortality", "darkFigure", "majorEvents", "schoolClosing", "contactRestrictions"};
##        x0Cov = [0.411, 0.069801, 0.02302, 0.009986, 5.7, 0.497353, 0.488051, 0.345821]; 
##        lb = [0.3, 0.001, 0.0001, 0.005, 5, 0.3, 0.25, 0.3];
##        ub = [0.5,  0.12, 0.12, 0.01, 10, 1, 1, 1];  
##        x0Cov = (lb+ub)/2; 

%%%%% Parameters from local PSO %%%%%%%%%
   
##        paraNames = {"beta", "darkFigure", "majorEvents", "schoolClosing", "contactRestrictions"};
##        x0Cov = [0.411, 5.7, 0.497353, 0.488051, 0.345821]; 
##        lb = [0.15, 5, 0.3, 0.25, 0.3]; %0.3
##        ub = [0.5, 10, 1, 1, 1];  
##        x0Cov = (lb+ub)/2;  

%%%%% Parameters for Counties - changed %%%%%%%%%
##        paraNames = {"beta", "darkFigure", "majorEvents", "schoolClosing", "contactRestrictions"};
##        x0Cov = [0.411, 5.7, 0.497353, 0.488051, 0.345821]; 
##        lb = [0.15, 4, 0.3, 0.2, 0.25]; %0.3
##        ub = [0.5, 15, 1, 1, 1];  
##        x0Cov = (lb+ub)/2;  

        paraNames = {"beta", "darkFigure", "majorEvents", "schoolClosing", "contactRestrictions"};
        x0Cov = [0.411, 5.7, 0.497353, 0.488051, 0.345821]; 
        lb = [0.15, 4, 0, 0, 0]; %0.3
        ub = [0.5, 15, 1, 1, 1];  
        x0Cov = (lb+ub)/2;  
       
    %endswitch
    
    optAlgorithm = "lsqnonlin"; %lsqnonlin|pso
    switch optAlgorithm
      case "lsqnonlin"
        %%%%%%%% Optfunction %%%%%%%%
        pkg load optim;
        opts = optimset();%("Algorithm", "levenberg-marquardt");
        ##  opts.TolFun = 1e-10;
        paramCovBest = x0Cov;
        x0CovInit = x0Cov;
        
        % use parallel package if available
        [dummy, info] = pkg('list');
        findParallel = strcmp( cellfun( @(ind) ind.name, info, 'uni', false), {'parallel'} );
        indParallel = find(findParallel==1);
        version = info{indParallel}.version;
        
        try
          pkg load parallel;
          if version == "4.0.0"
            pkg unload parallel;
          end    
        catch
          fprintf("No parallel package installed, code will run serial.\n");
        end_try_catch
        
        % change start values for every run randomly
        nruns = parameter.nRuns;
        nparallel = nproc();
        paramCovid = zeros(nruns, length(x0Cov));
        resNormCovid = inf(nruns, 1);
        runOptimOnce = @(x)lsqnonlinTwoRetvals(x, lb, ub, opts,...
        paraNames, parameter, RKIdata{k}, false, false);  
        
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
          % Save parameters to file
          folder = ["../../Results/", parameter.folderName];
          fid = fopen([folder, "/protokoll.txt"], "a");
          fprintf(fid, "%scounty: %s", RKIdata{k}.name);
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
          paramCovBest = paramCovid(indexBest, :);
        endwhile
        
        %%%%%%%% Results %%%%%%%%
        evalSirLocal(paramCovBest, paraNames, parameter, RKIdata{k}, true, true);
        ind=find(ismember(paraNames,'startDate'));
        dateInitial = paramCovBest(ind);
        date = paramCovBest(ind);
        for i=1:length(paraNames)
          if i==ind
            fprintf("%sInitial: %s %s: %s \n",  paraNames{i},...
            datestr(dateInitial), paraNames{i}, datestr(date));
          else
            fprintf("%sInitial: %f %s: %f \n",  paraNames{i},...
            x0CovInit(i), paraNames{i}, paramCovBest(i));
          end
        end
        fprintf("\n");  
        fprintf("resnormCovid: %f\n ", min(resNormCovid));
        fprintf("Calculation took %f seconds\n", toc); 
      case "pso"
      opt.parallel = false;
      opt.visu = true;
      opt.n_particles = 250;
      opt.n_iter = 200;
      opt.coupling = 5;
      opt.parameter_names = paraNames;
        
        evalSirLocalOneArgument= @(x) norm(evalSirLocal(x, paraNames,...
        parameter, RKIdata{k}, false, false));
        [paramCovBest, resMin] = pso(evalSirLocalOneArgument, lb, ub, opt);
        
        % Save parameters to file
        folder = ["../../Results/", parameter.folderName];
        fid = fopen([folder, "/protokoll.txt"], "a");
        fprintf(fid, "%s:\n", RKIdata{k}.name);   
        fprintf(fid, "Res: %f ", resMin);
        for j= 1:length(paramCovBest)
          fprintf(fid, "%s: %f | ", paraNames{j}, paramCovBest(j));
        end
        fprintf(fid, "\n");  
        fclose(fid);
        
        %%%%%%%% Results %%%%%%%%
        evalSirLocal(paramCovBest, paraNames, parameter, RKIdata{k}, true, true);
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
        fprintf("\n");  
        fprintf("resnormCovid: %f\n ", resMin);
        fprintf("Calculation took %f seconds\n", toc);  
    endswitch
		
  end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paraRun, resRun] = lsqnonlinTwoRetvals(x0Cov, lb, ub, opts,...
  paraNames, parameter, rkiData, showPlots, writeTexTable) 
  % saves fitted parameters for every single run 
  sirLocal = @(x)evalSirLocal(x, paraNames, parameter, rkiData, showPlots, writeTexTable); 
  try
    if iscell(x0Cov)
      x0Cov = cell2mat(x0Cov);
    end     
    [paraRun, resRun] = lsqnonlin(sirLocal, x0Cov, lb, ub, opts);
  catch
    folder = ["../../Results/", parameter.folderName];
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
function retval = evalSirLocal(Xpara, paraNames, parameter, rkiData, showPlots, writeTexTable)  
  % set parameters 
  if iscell(Xpara)
    Xpara = Xpara{1};
  endif
  for i=1:length(paraNames)
    parameter = setfield(parameter, paraNames{i}, Xpara(i));
  end
  
  % Modell soll bis zum xx.yy.2020 laufen
  parameter.totalRuntime = parameter.endDateOpt - parameter.startDate;
  startDateSim = parameter.startDate;  

  %set initial values für I and E
 shiftDays = 4; 
 sumIsimStartP1 = rkiData.infected(1);
 sumIsimStartP2 = rkiData.infected(1 + shiftDays);
 
 base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
  
  N = rkiData.population;
  I0 = ppval(rkiData.splineInfected, parameter.startDate);
  D0 = 0;
  %E0 = I0*(base_Q - 1)*((parameter.darkFigure * 2 - 1) / parameter.gamma1);	
  %E0 =  (parameter.beta*(parameter.darkFigure - 1)/parameter.gamma1 - parameter.darkFigure) * I0;
  E0 = I0*log(base_Q)*((parameter.darkFigure - 1) / parameter.gamma1)*exp(1);
  R0 = 0; 
  S0 = N - I0 - R0 - E0 - D0;
  X0sir = [S0, I0, R0, E0, D0] / N;
  
  % calculate model
  tic1 = tic;
  xSim = SIR_run(parameter, X0sir); 
  
  if !exist("parcellfun.m", "file")
    fprintf("Calculation starting at %s with beta = %d took %fs\n",...
    datestr(startDateSim), Xpara(1), toc(tic1));
  end
  
  % Spline für die Simulationsergebnisse im Zeitformat der RKIdaten
  % xSim.x (Zeit) laeuft von 0 bis totalRuntime -> +startDateSim
  switch parameter.model
    case "SIRH"
      % weights
      wI = 1;
      wD = 50;
      solutionDiscovered = xSim.y(2,:); % discovered
      solutionDead = xSim.y(3,:); % dead
    case "SIREDmod"
      % weights
      wI = 1;
      wD = 1;
      % discovered = \int_0^t alpha * E * dT       
      splineE = spline(xSim.x, xSim.y(4,:));
      alpha = parameter.gamma1 / (parameter.darkFigure - 1);
      get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
      [t2, discovered] = ode23(get_dDiscovering_dT, [min(xSim.x), max(xSim.x)], xSim.y(2,1));
      solutionDiscovered = spline(t2, discovered, xSim.x);
      solutionDead = xSim.y(5,:); % dead
  endswitch  
  splineSimDiscovered = spline(xSim.x + startDateSim, solutionDiscovered);
  splineSimDead = spline(xSim.x + startDateSim, solutionDead);
  
  % relevante Zeitspanne für den Fit ist die Ueberschneidung von Simulation
  % und RKIdaten
  timespan = [max([min(cell2mat(rkiData.time)), min(xSim.x + startDateSim)]),...
  min([max(cell2mat(rkiData.time)), max(xSim.x + startDateSim)])];
  times = linspace(timespan(1), timespan(2), 100);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %write vector that contains data for each state 
	simdataI = zeros(1,length(times)*length(rkiData));
	simdataIslope = zeros(1,length(rkiData));
	simdataD = zeros(1,length(times(end))*length(rkiData));
	
	wIstatewise = zeros(1,length(times)*length(rkiData));
##	wDstatewise = zeros(1,length(times(end))*length(rkiData));
	
	rkiI = zeros(1,length(times));
	rkiIslope = zeros(1,length(rkiData));
	rkiD = zeros(1,length(times(end)));  
	
	wwI = parameter.wRatioID;
	wwD = 1 - parameter.wRatioID;
  wwISlope = parameter.wISlope;
  
  simdataI = ppval(splineSimDiscovered, times)*rkiData.population;
  simdataIslope = diff(simdataI)(end);
  simdataD = ppval(splineSimDead, times(end))*rkiData.population;
  
  rkiI = ppval(rkiData.splineInfected, times);
  rkiIslope = diff(rkiI)(end);
  rkiD = ppval(rkiData.splineDead, times(end));
   
  wIstatewise(1, 1:length(times)) = wwI/max(rkiI);
  wDstatewise = wwD/max(rkiD); 

  wIslope = wwISlope/parameter.wISlopeGlob * length(rkiI) / length(rkiIslope);
	
	wDstatewise *= wwD * length(rkiI) / length(rkiD);
	
	ratio = parameter.wRatioGlobStates;
	wI = parameter.wIglob * (1 - ratio) + wIstatewise * ratio;
	wD = parameter.wDglob * (1 - ratio) + wDstatewise * ratio;

  retval = [(simdataI - rkiI) .* wI, (simdataD - rkiD) .* wD,...
  (simdataIslope - rkiIslope) .* wIslope];
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  % Differenzenvektor berechnen [Infected, Dead]
  % all for Infected, last point for dead
##  simdataI = ppval(splineSimDiscovered, times)*rkiData.population;
##  simdataD = ppval(splineSimDead, times(end))*rkiData.population;
##  rkiI = ppval(rkiData.splineInfected, times);
##  rkiD = ppval(rkiData.splineDead, times(end));
##  
##  ##  wI *= 1/max(rkiI);
##  ##  wD *= 1/max(rkiD);  
##  ##  retval = [(simdataI - rkiI) * wI, (simdataD - rkiD) * wD];
##  
##  wI *= 1/max(rkiI);
##  wD *= 0.1/max(rkiD)*length(rkiI)/length(rkiD); % "normalize" different length of I and D 
##  retval = [(simdataI - rkiI) * wI, (simdataD - rkiD) * wD];
  
  % evtl Plot anzeigen
  if showPlots
    close all;
    name = rkiData.name;
    figure('Name',name);
    subplot(1,2,1), hold on;
    subplot(1,2,2), hold on;
    subplot(1,2,1), plot(cell2mat(rkiData.time), rkiData.infected, "-x");%, times, retval(1:end/2));
    subplot(1,2,2), plot(cell2mat(rkiData.time), rkiData.dead, "-x");%, times, retval(end/2 + 1:end));
    namesI = {"RKI infected (cumulated)"};
    namesD = {"RKI dead"};
    switch parameter.model
      case "SIRH"
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*8e7, "-o");
        namesI{end + 1} = "discovered";
        subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(3,:)*8e7, "-o");
        namesD{end + 1} = "Dead";
      case "SIREDmod"
        subplot(1,2,1), plot(times, simdataI, "-o");
        namesI{end + 1} = "Cummulated Infected";             
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*rkiData.population);
        namesI{end + 1} = "Infected (with symptoms)";
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(3,:)*rkiData.population);
        namesI{end + 1} = "Recovered";
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(4,:)*rkiData.population);
        namesI{end + 1} = "Exposed + no symptoms";
        subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(5,:)*rkiData.population, "-o");
        namesD{end + 1} = "Dead";
    endswitch
    subplot(1,2,1), xlim([min(times), max(times)]);
    subplot(1,2,1), ylim(ppval(rkiData.splineInfected,[min(times), max(times)]).*[1, 1.5]);
    subplot(1,2,1), legend(namesI);
    subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
    subplot(1,2,2), xlim([min(times), max(times)]);
    subplot(1,2,2), ylim(ppval(rkiData.splineDead,[min(times), max(times)]).*[1, 1.5]);
    subplot(1,2,2), legend(namesD);
    subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
    
    %title(rkiData.name);
    fname = ["../../Results/", parameter.folderName];
    filename = ["plot_", rkiData.name];
    saveas(gcf, fullfile(fname, filename), 'png');
  end  
  
  if writeTexTable 
    
    stateNames= {"Schleswig-Holstein", "Freie und Hansestadt Hamburg",...
    "Niedersachsen", "Freie Hansestadt Bremen", "Nordrhein-Westfalen",...
    "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Freistaat Bayern",...
    "Saarland", "Berlin", "Brandenburg", "Mecklenburg-Vorpommern"...
    "Freistaat Sachsen", "Sachsen-Anhalt", "Freistaat Thueringen", "Germany"};
    
    k = parameter.AGS_state;
    switch parameter.model
      case {"SIRED","SIREDmod"}
        titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
        "Cumsum Infected", "RKI Infected", "RKI Dead"};
        ##    case "SIRH"
        ##      titles = {"Susceptible", "Exposed", "Removed", "Infectious",...
        ##      "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
    endswitch
    
    nVals = size(xSim.y, 1) + 3;
    persistent out = zeros(length(times), 1+(nVals)*length(stateNames));
    out(:,1) = times;
    
    splineS = spline(xSim.x + startDateSim, xSim.y(1,:));
    splineI2 = spline(xSim.x + startDateSim, xSim.y(2,:));
    splineR = spline(xSim.x + startDateSim, xSim.y(3,:));
    splineE = spline(xSim.x + startDateSim, xSim.y(4,:));
    simdataS = ppval(splineS, times)*rkiData.population;
    simdataI2 = ppval(splineI2, times)*rkiData.population;
    simdataR = ppval(splineR, times)*rkiData.population;
    simdataE = ppval(splineE, times)*rkiData.population;
    
    out(:,(k-1)*nVals+2) = simdataS;
    out(:,(k-1)*nVals+3) = simdataI2;
    out(:,(k-1)*nVals+4) = simdataR;
    out(:,(k-1)*nVals+5) = simdataE;
    out(:,(k-1)*nVals+6) = 0;
    out(size(out,1),(k-1)*nVals+6) = simdataD;
    out(:,(k-1)*nVals+7) = simdataI;
    out(:,(k-1)*nVals+8) = rkiI;
    out(:,(k-1)*nVals+9) = 0;
    out(size(out,1),(k-1)*nVals+9) = rkiD;
    
    if k==parameter.endoptMode
      folder = ["../../Results/", parameter.folderName];
      fid = fopen([folder, "/textable"], 'w'); 
      fprintf(fid, "%s ", "time");
      for i=1:k
        for j=1:length(titles)
          fprintf(fid, "%s ", [strrep(titles{j}," ","_"), "_",...
          strrep(stateNames{i}," ","_")]);
        endfor
      endfor
      fprintf(fid, "\n");
      for i=1:size(out, 1)
        fprintf(fid, "%f ", out(i,:));
        fprintf(fid, "\n");
      endfor  
      fclose(fid);
    end
  end
  
  end