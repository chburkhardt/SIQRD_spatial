function retval = optCovidLocal(xCov, paraNames, parameterArray, RKIdata, showPlots)
  %% this is how a parameter struct can be loaded from a parameterfile
  % alle Pfade sind so gesetzt das man im "Spatial" ordner starten muss,
  % daher auch hier

  addpath('../Spatial', '../Pso', '../Spatial/Optimize');
  
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
##  
##  shiftDays = 4; 
##  sumIsimStartP1 = RKIdata.infected(1);
##  sumIsimStartP2 = RKIdata.infected(1 + shiftDays);
  
##  base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
## switch parameterArray{1}.model
##   case "SIREDmod"
##    N = RKIdata.population;
##    I0 = ppval(RKIdata.splineInfected, parameterArray{1}.startDate);
##    D0 = 0;
##    %E0 = I0*(base_Q - 1)*((parameterArray{1}.darkFigure * 2 - 1) / parameterArray{1}.gamma1);	
##    %E0 =  (parameterArray{1}.beta*(parameterArray{1}.darkFigure - 1)/parameterArray{1}.gamma1 - parameterArray{1}.darkFigure) * I0;
##    E0 = I0*log(base_Q)*((parameterArray{1}.darkFigure - 1) / parameterArray{1}.gamma1)*exp(1);
##    R0 = 0; 
##    S0 = N - I0 - R0 - E0 - D0;
##    X0sir = [S0, I0, R0, E0, D0] / N;
##  case "SIR"
##    N = RKIdata.population;
##    I0 = ppval(RKIdata.splineInfected, parameterArray{1}.startDate);
##    D0 = 0;
##    %E0 = I0*(base_Q - 1)*((parameterArray{1}.darkFigure * 2 - 1) / parameterArray{1}.gamma1);	
##    %E0 =  (parameterArray{1}.beta*(parameterArray{1}.darkFigure - 1)/parameterArray{1}.gamma1 - parameterArray{1}.darkFigure) * I0;
##    E0 = I0*log(base_Q)*((parameterArray{1}.darkFigure - 1) / parameterArray{1}.gamma1)*exp(1);
##    R0 = 0; 
##    S0 = N - I0 - R0 - E0 - D0;
##    X0sir = [S0, I0, R0, 0, 0] / N;
##case "SIRH"
##  N = RKIdata.population;
##  I0 = ppval(RKIdata.splineInfected, parameterArray{1}.startDate);
##  D0 = 0;
##  %E0 = I0*(base_Q - 1)*((parameterArray{1}.darkFigure * 2 - 1) / parameterArray{1}.gamma1);	
##  %E0 =  (parameterArray{1}.beta*(parameterArray{1}.darkFigure - 1)/parameterArray{1}.gamma1 - parameterArray{1}.darkFigure) * I0;
##  E0 = I0*log(base_Q)*((parameterArray{1}.darkFigure - 1) / parameterArray{1}.gamma1)*exp(1);
##  R0 = 0; 
##  S0 = N - I0 - R0 - E0 - D0;
##  X0sir = [S0, I0, R0, 0, 0] / N;
##endswitch
 
  dataStart = datenum(2020,03,02);
  index = parameterArray{1}.startDate - dataStart+1;
  darkFigures = sir_eqn_spatial("getDarkFigureStatewise", parameterArray{1});
  k = parameterArray{1}.AGS_state;
  parameterArray{1}.population = RKIdata.population;
  if k ==17
    darkFigure = parameterArray{1}.darkFigure;
  else 
    darkFigure = darkFigures(k);
  end
  shiftDays = 4; 
  sumIsimStartP1 = RKIdata.infected(1);
  sumIsimStartP2 = RKIdata.infected(1 + shiftDays);
  
  base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
  lags = sir_eqn_spatial("lags");  
  if parameterArray{1}.startDate<737866
    i0 = RKIdata.infected(index)/RKIdata.population;
    %handle states that don't have infections yet (otherwise there would be a division by zero)
    if i0 == 0
      i0 = 1/RKIdata.population;
    end
    r0 = 0;
    d0 = RKIdata.dead(index)/RKIdata.population;
    e0 = i0*log(base_Q)*((parameterArray{1}.darkFigure - 1) / parameterArray{1}.gamma1)*exp(1);;
    e0l = RKIdata.infected(index+3)/RKIdata.population;
    s0 = 1-i0-r0-d0-e0;
  else 
    %r0 = (datavec.recovered(index)*darkFigures(k)-datavec.dead(index))/n_city;
    r0 = (RKIdata.infected(index-14)*darkFigure - RKIdata.dead(index-14))...
    /RKIdata.population;
    d0 = RKIdata.dead(index)/RKIdata.population;
    i0 = (RKIdata.infected(index)-RKIdata.infected(index-14))/RKIdata.population;
    e0l = (RKIdata.infected(index+3)-RKIdata.infected(index))/RKIdata.population;
    e0 = i0*darkFigures(k);  
    s0 = 1-i0-d0-e0-r0;
  end
  switch parameterArray{1}.model
     case {"SIRED", "SIREDmod", "SIREDYO"}
        X0sir = [s0, i0, r0, e0, d0];
        %darkFigure = darkFigures(k);
    case {"SIREDLiterature"}
        X0sir = [s0, i0, r0, e0l, d0];
        %darkFigure = darkFigures(k);
      case {"SIR"}
        X0sir = [1-i0-r0, i0, r0]';%, 0, 0];%
      case {"SIRH"}
        %cities{j}.SIR = [1-i0*parameter.darkFigure, i0, r0, 0, 0];
        %[1-i0-d0-r0, i0, r0, 0, 0];
        %cities{j}.baseQ = base_Q;
        if (parameterArray{1}.startDate>737908)
          %sum up population of all cities in county otherwise history would
          %be negative in case of small cities in big counties            
          infhistory = RKIdata.infected(index-28:index)
          i0h = infhistory/RKIdata.population;
          history = 1-i0h*darkFigure;%*factor;%*factor;%+i0h*parameter.darkFigure; % +r0h;%+d0h
          splineH = pchip(1:29, history, 29 - linspace(lags(1), lags(end), 40));
          %cities{j}.history = splineH;
          %cities{j}.SIR = [1-i0*parameter.darkFigure, i0, r0, 0, 0];
          X0sir = [1-i0*darkFigure, i0, r0, 0, 0];
          deathStart = d0;
          %darkFigure = darkFigures(k);
        else
          X0sir = [1-i0, i0, r0, 0, 0];
          %cities{j}.deathStart = d0;
        end
  endswitch
  
  %X0sir = spatialDistinction{parameterArray{1}.AGS_state}.SIR;
 
  %% run code and store the workspace to avoid io
  if (parameterArray{1}.startDate>737908)&&(strcmp("RKIfiles",parameterArray{1}.initial))&&(strcmp("SIRH",parameterArray{1}.model))
    xSim = SIR_run(parameterArray{1}, X0sir, splineH, deathStart, darkFigure, RKIdata.population);
  else% (parameterArray{1}.startDate>737908)&&(strcmp("RKIfiles",parameterArray{1}.initial))&&(strcmp("SIREDmod",parameterArray{1}.model))
    xSim = SIR_run(parameterArray{1}, X0sir);
  end
  for i = 1:length(paraNames)
    fprintf("%s: %d ", paraNames{i}, xCov(i));
  end
  fprintf("\n");
  
  %% load results statewise with ne new option to avoid io
##  stateWiseSIR = extractStatewiseResults (workspaceCalculation, "workspace",...
##  ["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
##  
  ##  for i=1:length(paraNames)
  ##    fprintf("%s= %f\n", paraNames{i}, xCov(i));
  ##  end

  
  switch parameterArray{1}.model
    case "SIRH"
##      % weights
##      wI = 1;
##      wD = 50;
      solutionDiscovered = xSim.y(2,:); % discovered
      solutionDead = xSim.y(3,:); % dead
    case "SIREDmod"
##      % weights
##      wI = 1;
##      wD = 1;
      % discovered = \int_0^t alpha * E * dT       
      splineE = spline(xSim.x, xSim.y(4,:));
      alpha = parameterArray{1}.gamma1 / (darkFigure - 1);
      get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
      %[t2, discovered] = ode23(get_dDiscovering_dT, [min(xSim.x), max(xSim.x)], xSim.y(2,1));
      [t2, discovered] = ode15s(get_dDiscovering_dT, [min(xSim.x),max(xSim.x)],...
      xSim.y(3,1)/darkFigure+xSim.y(2,1));
      solutionDiscovered = spline(t2, discovered, xSim.x);
      solutionDead = xSim.y(5,:); % dead
      
##      alpha = parameterArray{1}.gamma1 / (darkFigure - 1);
##      splineE = spline(x.x, x.y(4,:));
##      % Integriere E*alpha to get discovered (realy infected)
##      get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
####      [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
####      max(stateWiseSIR{i}.time)],...
####      stateWiseSIR{i}.SIR_vs_time(1,3)/darkFigures(i)+stateWiseSIR{i}.SIR_vs_time(1,2));      [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
##      [t2, discovered] = ode15s(get_dDiscovering_dT, [min(x.x),...
##      max(x.x)],...
##      x.y(3,1)/darkFigure+x.y(2,1));%-...
##      %stateWiseSIR{i}.SIR_vs_time(1,5)+stateWiseSIR{i}.SIR_vs_time(1,4)/darkFigures(i));
##      % get solutionDiscovered at the same timesteps as the other solutions
##      solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time); 
##      x.y(6,:) = solutionDiscovered;
  case "SIREDLiterature"
  
        darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1});  
        splineE = spline(xSim.x, xSim.y(4,:));
        alpha = 1/3;
        get_dDiscovering_dTL = @(t_eq, x0) [alpha*ppval(splineE, t_eq)];
        [t2, discovered] = ode15s(get_dDiscovering_dTL, [min(xSim.x),max(xSim.x)],xSim.y(2,1));
      %xSim.y(3,1)/darkFigure+xSim.y(2,1));
        %ode23(get_dDiscovering_dTL, [min(stateWiseSIR{i}.time),...
        %max(stateWiseSIR{i}.time)],stateWiseSIR{i}.SIR_vs_time(1,2));
        solutionDiscovered = spline(t2, discovered, xSim.x);        
        %stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
        %stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
        %solutionDiscovered);   
##        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,2));
##        stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,4));
##        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,5));
##        stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,3));
        solutionDead = xSim.y(5,:); 

   case "SIR"
     %solutionDiscovered = xSim.y(2,:);
     solutionDiscovered = 1-xSim.y(1,:);
     solutionDead = zeros(1,length(xSim.x));
  endswitch  
  splineSimDiscovered = spline(xSim.x + startDateSim, solutionDiscovered);
  splineSimDead = spline(xSim.x + startDateSim, solutionDead);
  
##  
##  switch parameterArray{1}.model
##    case "SIREDmod"
##      darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1});  
##      for i=1:length(stateWiseSIR)
##        %%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
##        % discovered = \int_0^t alpha * E * dT 
##        alpha = parameterArray{1}.gamma1 / (darkFigures(i) - 1);
##        splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
##        % Integriere E*alpha to get discovered (realy infected)
##        get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
##        [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
##        max(stateWiseSIR{i}.time)], stateWiseSIR{i}.SIR_vs_time(1,2));
##        % get solutionDiscovered at the same timesteps as the other solutions
##        solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);
##        
##        stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
##        solutionDiscovered);    
##        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,5));
##        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,2));
##        stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,3));
##        stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,4));
##      end
##    case "SIRH"
##      for i=1:length(stateWiseSIR)
##        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim,...
##        stateWiseSIR{i}.SIR_vs_time(:,8));
##        stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
##        stateWiseSIR{i}.SIR_vs_time(:,9));  
##        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim,...
##        stateWiseSIR{i}.SIR_vs_time(:,3));
##        stateWiseSIR{i}.splineR = spline([0 1], [0 0]);
##        stateWiseSIR{i}.splineE = spline([0 1], [0 0]);    
##      endfor
##  endswitch
  
  timespan = [max([min(cell2mat(RKIdata.time)), min(xSim.x + startDateSim)]),...
  min([max(cell2mat(RKIdata.time)), max(xSim.x + startDateSim)])];
  times = linspace(timespan(1), timespan(2), 100);
  
  
  %write vector that contains data for each state 
  simdataI = zeros(1,length(times)*length(RKIdata));
  simdataIslope = zeros(1,length(RKIdata));
  simdataD = zeros(1,length(times(end))*length(RKIdata));
  
  wIstatewise = zeros(1,length(times)*length(RKIdata));
  wDstatewise = zeros(1,length(times(end))*length(RKIdata));
  
  rkiI = zeros(1,length(times));
  rkiIslope = zeros(1,length(RKIdata));
  rkiD = zeros(1,length(times(end)));  
  
  wwI = parameterArray{1}.wRatioID;
  wwD = 1 - parameterArray{1}.wRatioID;
  wwISlope = parameterArray{1}.wISlope;
    simdataIstate = ppval(splineSimDiscovered, times)*RKIdata.population;%...
##    +(ppval(RKIdata{i}.splineInfected, times(1))-...
##    ppval(stateWiseSIR{i}.splineInfected, times(1)));
    simdataI(1, :) = simdataIstate;
    simdataIslope(1, 1) = ppval(splineSimDiscovered, times(end))-ppval(splineSimDiscovered, times(end-2));
    
    simdataDstate = ppval(splineSimDead, times(end))*RKIdata.population;%+...
##    (ppval(RKIdata{i}.splineDead, times(1))-...
##    ppval(stateWiseSIR{i}.splineDead, times(1)));
    simdataD(1, :) = simdataDstate;
    
    rkiIstate = ppval(RKIdata.splineInfected, times);
    rkiI(1, :) = rkiIstate;
    rkiIslope(1, 1) = ppval(RKIdata.splineInfected, times(end))-ppval(RKIdata.splineInfected, times(end-2));;
    wIstatewise = wwI./max(rkiIstate);
    
    rkiDstate = ppval(RKIdata.splineDead, times(end));
    rkiD(1, :) = rkiDstate;
    
    wIstatewise(1, 1:length(times)) = wwI/max(rkiI);
    wDstatewise = wwD./max(rkiDstate); 

  
  wIglob = wwI/max(rkiI);
  wDglob = wwD/max(rkiD) * length(rkiI) / length(rkiD);
  wIslope = wwISlope/max(rkiIstate);
  %wIslope = wwISlope/max(rkiIslope) * length(rkiI) / length(rkiIslope);
  
  wDstatewise *= wwD * length(rkiI) / length(rkiD);
  
  %ratio = parameterArray{1}.wRatioGlobStates;
  wI = wIglob;
  wD = wDglob;
  
  %opt function
##  if parameterArray{1}.optFunGermanSum 
##    % sum germany
##    retval = [sum(reshape((simdataI - rkiI) .* wI, length(simdataI)/length(xSim.y),...
##    length(xSim.y))', 1),...
##    sum(reshape((simdataD - rkiD) .* wD, length(simdataD)/length(xSim.y),...
##    length(xSim.y)), 2),...
##    sum(reshape((simdataIslope -rkiIslope) .* wIslope,...
##    length(xSim.y)/length(stateWiseSIR), length(xSim.y)), 2)];
##  else
  switch parameterArray{1}.model
    case "SIR"
    retval = [(simdataI - rkiI) .* wI,...
    (simdataIslope -rkiIslope) .* wIslope];
  case {"SIRED", "SIREDmod", "SIREDLiterature", "SIRH"}
     retval = [(simdataI - rkiI) .* wI, (simdataD - rkiD) .* wD,...
    (simdataIslope -rkiIslope) .* wIslope];
  endswitch
##  end
  
  folder = ["../../Results/", parameterArray{1}.folderName];
  fid = fopen([folder, "/protokoll_lsqnonlin.txt"], "a");
  fprintf(fid, "Res: %f ", norm(retval));
  for j= 1:length(xCov)
    fprintf(fid, "%s: %f | ", paraNames{j}, xCov(j));
  end
  fprintf(fid, "\n");
  fclose(fid);  
  
  % evtl Plot anzeigen
if showPlots
    close all;
    name = RKIdata.name;
    figure('Name',name);
    subplot(1,2,1), hold on;
    subplot(1,2,2), hold on;
    subplot(1,2,1), plot(cell2mat(RKIdata.time), RKIdata.infected, "-x");%, times, retval(1:end/2));
    subplot(1,2,2), plot(cell2mat(RKIdata.time), RKIdata.dead, "-x");%, times, retval(end/2 + 1:end));
    namesI = {"RKI infected (cumulated)"};
    namesD = {"RKI dead"};
    switch parameterArray{1}.model
      case "SIRH"
        %subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*8e7, "-o");
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*RKIdata.population, "-o");
        namesI{end + 1} = "discovered";
        %subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(3,:)*8e7, "-o");
        subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(3,:)*RKIdata.population, "-o");
        namesD{end + 1} = "Dead";
      case {"SIREDmod", "SIREDLiterature"}
        subplot(1,2,1), plot(times, simdataI, "-o");
        namesI{end + 1} = "Cummulated Infected";             
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*RKIdata.population);
        namesI{end + 1} = "Infected (with symptoms)";
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(3,:)*RKIdata.population);
        namesI{end + 1} = "Recovered";
        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(4,:)*RKIdata.population);
        namesI{end + 1} = "Exposed + no symptoms";
        subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(5,:)*RKIdata.population, "-o");
        namesD{end + 1} = "Dead";
      case "SIR"
        %subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*8e7, "-o");
        subplot(1,2,1), plot(xSim.x + startDateSim, (1-xSim.y(1,:))*RKIdata.population, "-o");
        namesI{end + 1} = "discovered";
        %subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(3,:)*8e7, "-o");
        %subplot(1,2,2), plot(xSim.x + startDateSim, xSim.y(3,:)*RKIdata.population, "-o");
        %namesD{end + 1} = "Dead";
    endswitch
    subplot(1,2,1), xlim([min(times), max(times)]);
    subplot(1,2,1), ylim(ppval(RKIdata.splineInfected,[min(times), max(times)]).*[1, 1.5]);
    subplot(1,2,1), legend(namesI);
    subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
    subplot(1,2,2), xlim([min(times), max(times)]);
    subplot(1,2,2), ylim(ppval(RKIdata.splineDead,[min(times), max(times)]).*[1, 1.5]);
    subplot(1,2,2), legend(namesD);
    subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
    
    %title(RKIdata.name);
    fname = ["../../Results/", parameterArray{1}.folderName];
    filename = ["plot_", RKIdata.name];
    saveas(gcf, fullfile(fname, filename), 'png');
  end  
  
##  if writeTexTable 
##    
##    stateNames= {"Schleswig-Holstein", "Freie und Hansestadt Hamburg",...
##    "Niedersachsen", "Freie Hansestadt Bremen", "Nordrhein-Westfalen",...
##    "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Freistaat Bayern",...
##    "Saarland", "Berlin", "Brandenburg", "Mecklenburg-Vorpommern"...
##    "Freistaat Sachsen", "Sachsen-Anhalt", "Freistaat Thueringen", "Germany"};
##    
##    k = parameterArray{1}.AGS_state;
##    switch parameterArray{1}.model
##      case {"SIRED","SIREDmod"}
##        titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
##        "Cumsum Infected", "RKI Infected", "RKI Dead"};
##        ##    case "SIRH"
##        ##      titles = {"Susceptible", "Exposed", "Removed", "Infectious",...
##        ##      "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
##    endswitch
##    
##    nVals = size(xSim.y, 1) + 3;
##    persistent out = zeros(length(times), 1+(nVals)*length(stateNames));
##    out(:,1) = times;
##    
##    splineS = spline(xSim.x + startDateSim, xSim.y(1,:));
##    splineI2 = spline(xSim.x + startDateSim, xSim.y(2,:));
##    splineR = spline(xSim.x + startDateSim, xSim.y(3,:));
##    splineE = spline(xSim.x + startDateSim, xSim.y(4,:));
##    simdataS = ppval(splineS, times)*RKIdata.population;
##    simdataI2 = ppval(splineI2, times)*RKIdata.population;
##    simdataR = ppval(splineR, times)*RKIdata.population;
##    simdataE = ppval(splineE, times)*RKIdata.population;
##    
##    out(:,(k-1)*nVals+2) = simdataS;
##    out(:,(k-1)*nVals+3) = simdataI2;
##    out(:,(k-1)*nVals+4) = simdataR;
##    out(:,(k-1)*nVals+5) = simdataE;
##    out(:,(k-1)*nVals+6) = 0;
##    out(size(out,1),(k-1)*nVals+6) = simdataD;
##    out(:,(k-1)*nVals+7) = simdataI;
##    out(:,(k-1)*nVals+8) = rkiI;
##    out(:,(k-1)*nVals+9) = 0;
##    out(size(out,1),(k-1)*nVals+9) = rkiD;
##    
##    if k==parameterArray{1}.endoptMode
##      folder = ["../../Results/", parameterArray{1}.folderName];
##      fid = fopen([folder, "/textable"], 'w'); 
##      fprintf(fid, "%s ", "time");
##      for i=1:k
##        for j=1:length(titles)
##          fprintf(fid, "%s ", [strrep(titles{j}," ","_"), "_",...
##          strrep(stateNames{i}," ","_")]);
##        endfor
##      endfor
##      fprintf(fid, "\n");
##      for i=1:size(out, 1)
##        fprintf(fid, "%f ", out(i,:));
##        fprintf(fid, "\n");
##      endfor  
##      fclose(fid);
##    end
##  end
##  
endfunction