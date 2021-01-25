% optCovid for statewise optimization
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
  
  if parameterArray{1}.fitNewInfections
    times = timespan(1):timespan(2);
  else
    times = linspace(timespan(1), timespan(2), 100); 
  end
  
  time2weeksToSimEnd = timespan(2)-14;
  indTwoWeeksToSimEnd = find(time2weeksToSimEnd==round(times))(1);
  time1WeekToSimEnd = timespan(2)-6;
  ind1WeekToSimEnd = find(time1WeekToSimEnd==round(times))(1);
  
  %Idee: addiere einfach die Anzahl der kumulierten Infizierten vom RKI am ersten
  %Tag der Optimierung drauf (siehe Bild plot_Germany)
  switch parameterArray{1}.model
    case "SIREDmod"
      darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1});  
      for i=1:length(stateWiseSIR)
        %%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
        % discovered = \int_0^t alpha * E * dT 
        alpha = parameterArray{1}.gamma1 / (darkFigures(i) - 1);
        splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
        % Integriere E*alpha to get discovered (realy infected)
        get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
        [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
        max(stateWiseSIR{i}.time)],...
        stateWiseSIR{i}.SIR_vs_time(1,3)/darkFigures(i)+stateWiseSIR{i}.SIR_vs_time(1,2));%-...
        %stateWiseSIR{i}.SIR_vs_time(1,5)+stateWiseSIR{i}.SIR_vs_time(1,4)/darkFigures(i));
        % get solutionDiscovered at the same timesteps as the other solutions
        solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);        
        
        stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
        solutionDiscovered);    
        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,5));
        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,2));
        stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,3));
        stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,4));
      end
    case "SIRH"
      for i=1:length(stateWiseSIR)
        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim,...
        stateWiseSIR{i}.SIR_vs_time(:,8));
        stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
        stateWiseSIR{i}.SIR_vs_time(:,9));  
        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim,...
        stateWiseSIR{i}.SIR_vs_time(:,3));
        stateWiseSIR{i}.splineR = spline([0 1], [0 0]);
        stateWiseSIR{i}.splineE = spline([0 1], [0 0]); 
        stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim,...
        stateWiseSIR{i}.SIR_vs_time(:,2));   
      endfor
 case "SIREDLiterature"
   for i=1:length(stateWiseSIR)
        darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1});  
        splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
        alpha = 1/3;
        get_dDiscovering_dTL = @(t_eq, x0) [alpha*ppval(splineE, t_eq)];
        if (parameterArray{1}.startDate<737866)
          [t2, discovered] = ode23(get_dDiscovering_dTL, [min(stateWiseSIR{i}.time),...
          max(stateWiseSIR{i}.time)],stateWiseSIR{i}.SIR_vs_time(1,2));
        else
          [t2, discovered] = ode23(get_dDiscovering_dTL, [min(stateWiseSIR{i}.time),...
          max(stateWiseSIR{i}.time)],(stateWiseSIR{i}.SIR_vs_time(1,3)/darkFigures(i)...
          +stateWiseSIR{i}.SIR_vs_time(1,2)));
        end
        solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);        
        stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
        solutionDiscovered);   
        stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,2));
        stateWiseSIR{i}.splineE = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,4));
        stateWiseSIR{i}.splineDead = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,5));
        stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim, stateWiseSIR{i}.SIR_vs_time(:,3));
        
   endfor
  case "SIR"
      for i=1:length(stateWiseSIR)
         darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameterArray{1}); 
       if parameterArray{1}.startDate<737866
         stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
         stateWiseSIR{i}.SIR_vs_time(:,4)); 
       else
         stateWiseSIR{i}.splineInfected = spline(stateWiseSIR{i}.time + startDateSim,...
         stateWiseSIR{i}.SIR_vs_time(:,4)./darkFigures(i)); 
       end
         stateWiseSIR{i}.splineI = spline(stateWiseSIR{i}.time + startDateSim,...
         stateWiseSIR{i}.SIR_vs_time(:,2));
         stateWiseSIR{i}.splineR = spline(stateWiseSIR{i}.time + startDateSim,...
         stateWiseSIR{i}.SIR_vs_time(:,3));
      end
  endswitch
  
  if isfield(parameterArray{1},'slopeFitInterval')
    if strcmp(parameterArray{1}.slopeFitInterval, 'end')
      timeS = times(end-2);
    elseif strcmp(parameterArray{1}.slopeFitInterval, 'lastWeek')
      timeS = times(ind1WeekToSimEnd);
    end
  else 
    timeS = times(end-2);
  end
  
  switch parameterArray{1}.model
    case {"SIRH", "SIREDmod", "SIREDLiterature"}
      if isfield(parameterArray{1},'dayForDeathFit')
        if strcmp(parameterArray{1}.dayForDeathFit, 'lastDay')
          timeD = times(end);
        elseif strcmp(parameterArray{1}.dayForDeathFit, 'twoWeeksToSimEnd')
          timeD = times(indTwoWeeksToSimEnd);
        end
      else 
        timeD = times(end);
      end
      if parameterArray{1}.fitNewInfections;%and(parameter.initial = "RKIfiles", parameterArray{1}.startDate>737908)
        %pkg load data-smoothing
          simdataI = zeros(1,(length(times)-7)*length(stateWiseSIR));
          simdataD = zeros(1,length(times(end))*length(stateWiseSIR));
          simdataDComplete = zeros(1,length(times)*length(stateWiseSIR));
          
          wIstatewise = zeros(1,length(times)-7*length(stateWiseSIR));
          wDstatewise = zeros(1,length(times(end))*length(stateWiseSIR));
          
          rkiI = zeros(1,(length(times)-7)*length(stateWiseSIR));
          rkiD = zeros(1,length(times(end))*length(stateWiseSIR)); 
          rkiDComplete = zeros(1,length(times)*length(stateWiseSIR));
          
          wwI = parameterArray{1}.wRatioID;
          wwD = 1 - parameterArray{1}.wRatioID;

        for i=1:length(stateWiseSIR)
          diffSim = diff(ppval(stateWiseSIR{i}.splineInfected, times));
          simdataIstate = movmean(diffSim, [6 0], "Endpoints", "discard");
          diffRKI = diff(ppval(RKIdata{i}.splineInfected, times));
          rkiIstate = movmean(diffRKI, [6 0], "Endpoints", "discard");
    
          simdataI(1, (i-1)*(length(times)-7)+1:i*(length(times)-7)) = simdataIstate;  
          rkiI(1, (i-1)*(length(times)-7)+1:i*(length(times)-7)) = rkiIstate;
          
          wIstatewise(1, (i-1)*(length(times)-7)+1:i*(length(times)-7)) = wwI./max(rkiIstate);
      
          simdataDstate = ppval(stateWiseSIR{i}.splineDead, timeD)-ppval(stateWiseSIR{i}.splineDead, times(1));
          simdataD(1, (i-1)*length(times(end))+1:i*length(times(end))) = simdataDstate;
          simdataDstateComplete =  ppval(stateWiseSIR{i}.splineDead, times);
          simdataDComplete(1, (i-1)*length(times)+1:i*length(times)) = simdataDstateComplete;
          
          rkiDstate = ppval(RKIdata{i}.splineDead, timeD)-ppval(RKIdata{i}.splineDead, times(1));
          rkiD(1, (i-1)*length(times(end))+1:i*length(times(end))) = rkiDstate;
          rkiDstateComplete =  ppval(RKIdata{i}.splineDead, times);
          rkiDComplete(1, (i-1)*length(times)+1:i*length(times)) = rkiDstateComplete;
          
          wDstatewise(1, (i-1)*length(times(end))+1:i*length(times(end))) = wwD./max(rkiDstate);  
        end
        wIglob = wwI/max(rkiI);
        wDglob = wwD/max(rkiD) * length(rkiI) / length(rkiD);
        
        wDstatewise *= wwD * length(rkiI) / length(rkiD);
        
        ratio = parameterArray{1}.wRatioGlobStates;
        wI = wIglob * (1 - ratio) + wIstatewise * ratio;
        wD = wDglob * (1 - ratio) + wDstatewise * ratio;     
      else
        %write vector that contains data for each state 
        simdataI = zeros(1,length(times)*length(stateWiseSIR));
        simdataIslope = zeros(1,length(stateWiseSIR));
        simdataD = zeros(1,length(times(end))*length(stateWiseSIR));
        simdataDComplete = zeros(1,length(times)*length(stateWiseSIR));
        
        wIstatewise = zeros(1,length(times)*length(stateWiseSIR));
        wDstatewise = zeros(1,length(times(end))*length(stateWiseSIR));
        
        rkiI = zeros(1,length(times)*length(stateWiseSIR));
        rkiIslope = zeros(1,length(stateWiseSIR));
        rkiD = zeros(1,length(times(end))*length(stateWiseSIR)); 
        rkiDComplete = zeros(1,length(times)*length(stateWiseSIR));
        
        wwI = parameterArray{1}.wRatioID;
        wwD = 1 - parameterArray{1}.wRatioID;
        wwISlope = parameterArray{1}.wISlope;
        for i=1:length(stateWiseSIR)
          simdataIstate = ppval(stateWiseSIR{i}.splineInfected, times);
          simdataI(1, (i-1)*length(times)+1:i*length(times)) = simdataIstate;
          simdataIslope(1, i) = ppval(stateWiseSIR{i}.splineInfected, times(end))-ppval(stateWiseSIR{i}.splineInfected, timeS);

          simdataDstate = ppval(stateWiseSIR{i}.splineDead, timeD);
          simdataD(1, (i-1)*length(times(end))+1:i*length(times(end))) = simdataDstate;
          simdataDstateComplete =  ppval(stateWiseSIR{i}.splineDead, times);
          simdataDComplete(1, (i-1)*length(times)+1:i*length(times)) = simdataDstateComplete;
          rkiIstate = ppval(RKIdata{i}.splineInfected, times);
          rkiI(1, (i-1)*length(times)+1:i*length(times)) = rkiIstate;
          rkiIslope(1, i) = ppval(RKIdata{i}.splineInfected, times(end))-ppval(RKIdata{i}.splineInfected, timeS);

          wIstatewise(1, (i-1)*length(times)+1:i*length(times)) = wwI./max(rkiIstate);
          
          rkiDstate = ppval(RKIdata{i}.splineDead, timeD);
          rkiD(1, (i-1)*length(times(end))+1:i*length(times(end))) = rkiDstate;
          rkiDstateComplete =  ppval(RKIdata{i}.splineDead, times);
          rkiDComplete(1, (i-1)*length(times)+1:i*length(times)) = rkiDstateComplete;
          wDstatewise(1, (i-1)*length(times(end))+1:i*length(times(end))) = wwD./max(rkiDstate); 
        end
        
        wIglob = wwI/max(rkiI);
        wDglob = wwD/max(rkiD) * length(rkiI) / length(rkiD);
        wIslope = wwISlope/max(rkiIslope) * length(rkiI) / length(rkiIslope);
        
        wDstatewise *= wwD * length(rkiI) / length(rkiD);
        
        ratio = parameterArray{1}.wRatioGlobStates;
        wI = wIglob * (1 - ratio) + wIstatewise * ratio;
        wD = wDglob * (1 - ratio) + wDstatewise * ratio;
        
      end
      
      %opt function
      if parameterArray{1}.fitNewInfections
        if parameterArray{1}.optFunGermanSum 
          % sum germany
          retval = [sum(reshape((simdataI - rkiI) .* wI, length(simdataI)/length(stateWiseSIR),...
          length(stateWiseSIR))', 1),...
          sum(reshape((simdataD - rkiD) .* wD, length(simdataD)/length(stateWiseSIR),...
          length(stateWiseSIR)), 2)];
        else
          retval = [(simdataI - rkiI) .* wI, (simdataD - rkiD) .* wD];
        end
      else
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
      end
      
    case "SIR"
      for i=1:length(stateWiseSIR)
        simdataIstate = ppval(stateWiseSIR{i}.splineInfected, times);
        simdataI(1, (i-1)*length(times)+1:i*length(times)) = simdataIstate;
        simdataIslope(1, i) = diff(simdataIstate)(end);
        
        rkiIstate = ppval(RKIdata{i}.splineInfected, times);
        rkiI(1, (i-1)*length(times)+1:i*length(times)) = rkiIstate;
        rkiIslope(1, i) = diff(rkiIstate)(end);
        wIstatewise(1, (i-1)*length(times)+1:i*length(times)) = wwI./max(rkiIstate);
      end
      wIglob = wwI/max(rkiI);
      wIslope = wwISlope/max(rkiIslope) * length(rkiI) / length(rkiIslope);
      ratio = parameterArray{1}.wRatioGlobStates;
      wI = wIglob * (1 - ratio) + wIstatewise * ratio;
      
      if parameterArray{1}.optFunGermanSum 
        
        retval = [sum(reshape((simdataI - rkiI) .* wI, length(simdataI)/length(stateWiseSIR),...
        length(stateWiseSIR))', 1),...
        sum(reshape((simdataIslope -rkiIslope) .* wIslope,...
        length(simdataIslope)/length(stateWiseSIR), length(stateWiseSIR)), 2)];
      else
        retval = [(simdataI - rkiI) .* wI, (simdataIslope -rkiIslope) .* wIslope];
      end    
  endswitch
  
  if strcmp(parameterArray{1}.optAlgorithm,"lsqnonlin")
    folder = ["../../Results/", parameterArray{1}.folderName];
    fid = fopen([folder, "/protokoll_lsqnonlin.txt"], "a");
    fprintf(fid, "Res: %f ", norm(retval));
    for j= 1:length(xCov)
      fprintf(fid, "%s: %f | ", paraNames{j}, xCov(j));
    end
    fprintf(fid, "\n");
    fclose(fid);  
  end
  
  % evtl Plot anzeigen
  if showPlots
    %sum up all states
    if parameterArray{1}.fitNewInfections
      reshapeSize = length(times)-7;
      plotSize = times(8:end);
    else
      reshapeSize = length(times);
      plotSize = times;
    end
    rkiIsum =  sum(reshape(rkiI, reshapeSize, length(stateWiseSIR)), 2);
    rkiDsumComplete = sum(reshape(rkiDComplete, length(times), length(stateWiseSIR)), 2);
    if parameterArray{1}.fitNewInfections
      rkiDsum =  sum(reshape(rkiD, length(times(end)), length(stateWiseSIR)), 2)+...
      rkiDsumComplete(1);
    else
      rkiDsum =  sum(reshape(rkiD, length(times(end)), length(stateWiseSIR)), 2);
    end
    ##    rkiI =  sum(reshape(rkiI, length(times), length(stateWiseSIR)), 2);
    ##    rkiI =  sum(reshape(rkiI, length(times), length(stateWiseSIR)), 2);
    figure(17);
    clf;
    subplot(1,2,1), hold on;
    subplot(1,2,2), hold on;
    
    subplot(1,2,1), plot(plotSize, rkiIsum, "-x");
    
    subplot(1,2,2), plot(timeD, rkiDsum, "-x");%, times, retval(end/2 + 1:end));
    subplot(1,2,2), plot(times, rkiDsumComplete, ":");
    if parameterArray{1}.fitNewInfections;
      namesI = {"7-day average of new infections"};
    else
      namesI = {"RKI infected (cumulated)"};
    end
    namesD = {"RKI dead"};
    
    switch parameterArray{1}.model
      case {"SIRH", "SIREDmod"}
        for i=1:length(stateWiseSIR)
          simdataEstate = ppval(stateWiseSIR{i}.splineE, times);
          simdataE(1, (i-1)*length(times)+1:i*length(times)) = simdataEstate;
        end
    endswitch
    
    
    switch parameterArray{1}.model
      case "SIRH"
        if parameterArray{1}.fitNewInfections; %objectiveFunction = "newInfections";%and(parameter.initial = "RKIfiles", parameterArray{1}.startDate>737908)
          simdataIsum =  sum(reshape(simdataI, reshapeSize, length(stateWiseSIR)), 2);
          subplot(1,2,1), plot(plotSize, simdataIsum, "-o");
          namesI{end + 1} = "New Infections";
        else
          simdataIsum =  sum(reshape(simdataI, reshapeSize, length(stateWiseSIR)), 2);
          subplot(1,2,1), plot(times, simdataIsum, "-o");
          namesI{end + 1} = "Cummulated Infected";
        end
        simdataDsumComplete = sum(reshape(simdataDComplete, length(times), length(stateWiseSIR)), 2);
        subplot(1,2,2), plot(times, simdataDsumComplete, "--");
        namesD{end + 1} = "Dead Complete";
        if parameterArray{1}.fitNewInfections
          simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2)+...
          simdataDsumComplete(1);
        else
          simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2);
        end
        subplot(1,2,2), plot(timeD, simdataDsum, "-o");
        namesD{end + 1} = "Dead";
        simdataEsum =  sum(reshape(simdataE, length(times), length(stateWiseSIR)), 2);
        subplot(1,2,1), plot(times, simdataEsum, "-o");
        namesI{end + 1} = "Exposed";
        subplot(1,2,2), xlim([min(times), max(times)]);
        
        if parameterArray{1}.startDate < 737900
          subplot(1,2,2), ylim([0, max(rkiDsum)*1.5]);
        else
          subplot(1,2,2), ylim([min(rkiDsumComplete), max(rkiDsum)*1.2]);
        end
        subplot(1,2,2), legend(namesD);
        subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
    case {"SIREDmod", "SIREDLiterature"}
      if parameterArray{1}.fitNewInfections;
        simdataIsum =  sum(reshape(simdataI, reshapeSize, length(stateWiseSIR)), 2);
        subplot(1,2,1), plot(plotSize, simdataIsum, "-o");
        namesI{end + 1} = "New Infections";
      else
        simdataIsum =  sum(reshape(simdataI, length(times), length(stateWiseSIR)), 2);
        subplot(1,2,1), plot(times, simdataIsum, "-o");
        namesI{end + 1} = "Cummulated Infected";
      end
        ##        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(2,:)*8e7);
        ##        namesI{end + 1} = "Infected (with symptoms)";
        ##        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(3,:)*8e7);
        ##        namesI{end + 1} = "Recovered";
        ##        subplot(1,2,1), plot(xSim.x + startDateSim, xSim.y(4,:)*8e7);
        ##        namesI{end + 1} = "Exposed + no symptoms";
        simdataDsumComplete = sum(reshape(simdataDComplete, length(times), length(stateWiseSIR)), 2);
        subplot(1,2,2), plot(times, simdataDsumComplete, "--");
        namesD{end + 1} = "Dead Complete";
        if parameterArray{1}.fitNewInfections
          simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2)+...
          simdataDsumComplete(1);
        else
          simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2);
        end
        subplot(1,2,2), plot(timeD, simdataDsum, "-o");
        namesD{end + 1} = "Dead";
        subplot(1,2,2), xlim([min(times), max(times)]);
        subplot(1,2,2), ylim([0, max(rkiDsum)*1.5]);
        subplot(1,2,2), legend(namesD);
        subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
   case "SIR"
     if parameterArray{1}.fitNewInfections;
        simdataIsum =  sum(reshape(simdataI, reshapeSize, length(stateWiseSIR)), 2);
        subplot(1,2,1), plot(plotSize, simdataIsum, "-o");
        namesI{end + 1} = "New Infections";
     else
        simdataIsum =  sum(reshape(simdataI, length(times), length(stateWiseSIR)), 2);
        subplot(1,2,1), plot(times, simdataIsum, "-o");
        namesI{end + 1} = "Cummulated Infected";
     end
    endswitch
    subplot(1,2,1), xlim([min(times), max(times)]);
    if parameterArray{1}.startDate < 737900
        subplot(1,2,1), ylim([0, max(rkiIsum)*1.5]);
    else
        subplot(1,2,1), ylim([min(rkiIsum), max(rkiIsum)*1.1]);
    end
    subplot(1,2,1), legend(namesI);
    subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
    title("Germany");
    fname = ["../../Results/", parameterArray{1}.folderName];
    filename = "plot_Germany";
    drawnow;
    pause(0.3);
    saveas(gcf, fullfile(fname, filename), 'jpeg');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch parameterArray{1}.model
      case {"SIRH", "SIREDmod", "SIREDLiterature"}
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
        simdataIplot =  reshape(simdataI, reshapeSize, length(stateWiseSIR));
        simdataI2plot = reshape(simdataI2, length(times), length(stateWiseSIR));
        simdataRplot = reshape(simdataR, length(times), length(stateWiseSIR));
        simdataEplot = reshape(simdataE, length(times), length(stateWiseSIR));
        %simdataDplot =  reshape(simdataD, length(times(end)), length(stateWiseSIR));
        simdataDplotComplete = reshape(simdataDComplete, length(times), length(stateWiseSIR));
        if parameterArray{1}.fitNewInfections
          simdataDplot =  reshape(simdataD, length(times(end)), length(stateWiseSIR))+...
          simdataDplotComplete(1,:);
        else
          simdataDplot =  reshape(simdataD, length(times(end)), length(stateWiseSIR));
        end
      case "SIR"
        simdataR = zeros(1,length(times)*length(stateWiseSIR));
        for i=1:length(stateWiseSIR)
          simdataR(1, (i-1)*length(times)+1:i*length(times)) =...
          ppval(stateWiseSIR{i}.splineR, times);
        end
        simdataIplot =  reshape(simdataI, reshapeSize, length(stateWiseSIR));
        simdataRplot = reshape(simdataR, length(times), length(stateWiseSIR));
    endswitch
    ##    end
    %statewise plots
    for i=1:length(RKIdata)
      ##    population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
      ##    4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
      ##    1609675, 4077937, 2208321, 2143145];
      
      if !parameterArray{1}.showStatePlots
        break;
      end
      ##      if parameter.objectiveFunction = "newInfections";
      ##        rkiIplot = reshape(rkiI, length(times)-1, length(stateWiseSIR));
      ##        rkiDplot =  reshape(rkiD, length(times(end)), length(stateWiseSIR));
      ##      else
      rkiIplot =  reshape(rkiI, reshapeSize, length(stateWiseSIR));
      rkiDplotComplete = reshape(rkiDComplete, length(times), length(stateWiseSIR));
      if parameterArray{1}.fitNewInfections
          rkiDplot =  reshape(rkiD, length(times(end)), length(stateWiseSIR))+...
          rkiDplotComplete(1,:);
      else
        rkiDplot =  reshape(rkiD, length(times(end)), length(stateWiseSIR));
      end
      
      ##      end
      
      caseFatilityRateRKI = rkiDplot(end, i) / rkiIplot(end, i);
      
      figNum = figure;
      subplot(1,2,1), hold on;
      subplot(1,2,2), hold on;
      subplot(1,2,1), plot(plotSize, rkiIplot(:,i), "-x");%, times, retval(1:end/2));
      subplot(1,2,2), plot(timeD, rkiDplot(:,i), "-x");%, times, retval(end/2 + 1:end));
      subplot(1,2,2), plot(times, rkiDplotComplete(:,i), ":");
      if parameterArray{1}.fitNewInfections;
        namesI = {"7-day average of new infections"};
      else
        namesI = {"RKI infected (cumulated)"};
      end
      namesD = {"RKI dead"};
      ##      fname = ["../../Results/", parameterArray{1}.folderName];
      ##      filename = ["plot_Germany"];
      ##      saveas(gcf, fullfile(fname, filename), 'jpeg');;
      
      switch parameterArray{1}.model
        case "SIRH"
          if parameterArray{1}.fitNewInfections;%and(parameter.initial = "RKIfiles", parameterArray{1}.startDate>737908)
            ##        simdataIsum =  sum(reshape(simdataI, length(times)-1, length(stateWiseSIR)), 2);
            ##        subplot(1,2,1), plot(times(2:end), simdataIsum, "-o");
            ##        namesI{end + 1} = "Cummulated Infected";
            ##        simdataDsum =  sum(reshape(simdataD, length(times(end)), length(stateWiseSIR)), 2);
            ##        subplot(1,2,2), plot(times(end), simdataDsum, "-o");
            ##        namesD{end + 1} = "Dead";
            ##        simdataEsum =  sum(reshape(simdataE, length(times), length(stateWiseSIR)), 2);
            ##        subplot(1,2,1), plot(times, simdataEsum, "-o");
            ##        namesI{end + 1} = "Exposed";
            subplot(1,2,1), plot(plotSize, simdataIplot(:,i), "-o");
            namesI{end + 1} = "New Infections";
          else
            subplot(1,2,1), plot(times, simdataIplot(:,i), "-o");
            namesI{end + 1} = "discovered";
          end
          subplot(1,2,2), plot(timeD, simdataDplot(:,i), "-o");
          namesD{end + 1} = "Dead";   
          subplot(1,2,2), plot(times, simdataDplotComplete(:,i), "--");
          namesD{end + 1} = "Dead Complete";
          subplot(1,2,1), plot(times, simdataEplot(:,i));
          namesI{end + 1} = "Exposed";  
          subplot(1,2,2), hold on;
          subplot(1,2,2), plot(timeD, rkiDplot(:,i), "-x");%, times, retval(end/2 + 1:end));
          subplot(1,2,2), xlim([min(times), max(times)]);
          if parameterArray{1}.startDate < 737900
            subplot(1,2,2), ylim([0, max(max(rkiDplot(:,i))*1.5 ,1)]);
          else
            subplot(1,2,2), ylim([min(rkiDplotComplete(:,i)), max(rkiDplot(:,i))*1.2]);
          end
          subplot(1,2,2), legend(namesD);
          subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
##          namesI = {"RKI infected (cumulated)"};
##          namesD = {"RKI dead"};
##                        
          caseFatilityRateSim = simdataDplot(end, i)/simdataIplot(end, i);
          %fprintf("%f\n", caseFatilityRateRKI / caseFatilityRateSim);
          fprintf("%f, ", rkiIplot(end, i)/simdataIplot(end, i));
          
          name = [RKIdata{i}.name, " df_error: ",...
          num2str(caseFatilityRateRKI / caseFatilityRateSim),...
          " Ratio rkiI/simI: ", num2str(rkiIplot(end, i)/simdataIplot(end, i))];
      
      case {"SIREDmod", "SIREDLiterature"}    
        if parameterArray{1}.fitNewInfections;
          subplot(1,2,1), plot(plotSize, simdataIplot(:,i), "-o");
          namesI{end + 1} = "New Infections";
        else
          subplot(1,2,1), plot(times, simdataIplot(:,i), "-o");
          namesI{end + 1} = "Cummulated Infected";
        end
          subplot(1,2,1), plot(times, simdataI2plot(:,i));
          namesI{end + 1} = "Infected (with symptoms)";
          subplot(1,2,1), plot(times, simdataRplot(:,i));
          namesI{end + 1} = "Recovered";
          subplot(1,2,1), plot(times, simdataEplot(:,i));
          namesI{end + 1} = "Exposed + no symptoms";     
          subplot(1,2,2), plot(timeD, simdataDplot(:,i), "-o");
          namesD{end + 1} = "Dead";
          subplot(1,2,2), plot(times, simdataDplotComplete(:,i), "--");
          namesD{end + 1} = "Dead Complete";
          subplot(1,2,2), hold on;
          subplot(1,2,2), plot(timeD, rkiDplot(:,i), "-x");%, times, retval(end/2 + 1:end));
##          namesI = {"RKI infected (cumulated)"};
##          namesD = {"RKI dead"};
          subplot(1,2,2), xlim([min(times), max(times)]);
          if parameterArray{1}.startDate < 737900
            subplot(1,2,2), ylim([0, max(max(rkiDplot(:,i))*1.5 ,1)]);
          else
            subplot(1,2,2), ylim([min(rkiDplotComplete(:,i)), max(rkiDplot(:,i))*1.2]);
          end
          %subplot(1,2,2), ylim([0, max(max(rkiDplot(:,i))*1.5 ,1)]);
          subplot(1,2,2), legend(namesD);
          subplot(1,2,2), datetick ("x", "dd.mm.yyyy", "keeplimits");
                
          caseFatilityRateSim = simdataDplot(end, i)/simdataIplot(end, i);
          %fprintf("%f\n", caseFatilityRateRKI / caseFatilityRateSim);
          fprintf("%f, ", rkiIplot(end, i)/simdataIplot(end, i));
          
          name = [RKIdata{i}.name, " df_error: ",...
          num2str(caseFatilityRateRKI / caseFatilityRateSim),...
          " Ratio rkiI/simI: ", num2str(rkiIplot(end, i)/simdataIplot(end, i))];       
         case"SIR"
          subplot(1,2,1), plot(times, simdataIplot(:,i), "-o");
          namesI{end + 1} = "Cummulated Infected";
          subplot(1,2,1), plot(times, simdataRplot(:,i));
          namesI{end + 1} = "Recovered";
          name = [RKIdata{i}.name];
      endswitch
      subplot(1,2,1), xlim([min(times), max(times)]);
      if parameterArray{1}.startDate < 737900
        subplot(1,2,1), ylim([0, max(rkiIplot(:,i))*1.5]);
      else
        subplot(1,2,1), ylim([min(rkiIplot(:,i)), max(rkiIplot(:,i))*1.1]);
      end
      %subplot(1,2,1), ylim([0, max(rkiIplot(:,i))*1.5]);
      subplot(1,2,1), legend(namesI);
      subplot(1,2,1), datetick ("x", "dd.mm.yyyy", "keeplimits");
      
      %write ratios in protokoll
      folder = ["../../Results/", parameterArray{1}.folderName];
      fid = fopen([folder, "/protokoll_final.txt"], "a");
      switch parameterArray{1}.model
      case {"SIRH", "SIREDmod"}
        fprintf(fid, "\county: %s\ndf_error: %f\nRatio rkiI/simI: %f ", RKIdata{i}.name, (caseFatilityRateRKI / caseFatilityRateSim), (rkiIplot(end, i)/simdataIplot(end, i))); 
        fprintf(fid, "\n");
      endswitch
      fclose(fid);  
      
      set(figNum, 'Name', name);
      fname = ["../../Results/", parameterArray{1}.folderName];
      filename = ["plot_", RKIdata{i}.name];
      drawnow;
      pause(0.3);
      saveas(gcf, fullfile(fname, filename), 'jpeg');
    endfor
  endif
endfunction
