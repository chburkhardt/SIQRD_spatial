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
## @deftypefn {} {@var{retval} =} standAlonePostProcessing (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-13

function retval = standAlonePostProcessing (folder, filename,...
  workspace)
  tstart = tic;
  if nargin == 0
    error("No filename given for standAlonePostProcessing");
  end
  
  if nargin < 3
    load([folder, "/", filename]);
    workspace.x = x;
    t = x.x';
    x = x.y'; 
    workspace.parameter = parameter;
    if exist("initialHistory", "var")
      workspace.initialHistory = initialHistory;
    end
    workspace.cities = cities;
    workspace.distances = distances;
  else
    t = workspace.x.x';
    x = workspace.x.y';
    parameter = workspace.parameter;
    if isfield(workspace , "initialHistory")
      initialHistory = workspace.initialHistory;
    end    
    cities = workspace.cities;
    filename = "notValid";
  end
  
  
  if parameter.vtkStyleMap    
    load("../../Daten/shapeData.mat");
    agsVecShape = zeros(length(shapeData) ,1);
    for i=1:length(agsVecShape)
      agsVecShape(i) = shapeData{i}.ags;
      for j=1:length(shapeData{i}.coordinates)
        [coordsX, coordsY] = latLon2XY(shapeData{i}.coordinates{j}(:,2),...
        shapeData{i}.coordinates{j}(:,1));
        shapeData{i}.coordinates{j} = [coordsX, coordsY];
      end
    end 
    
    if isfield(parameter, "keepAGS")
      load("../../Daten/GemeindeshapeData.mat");
      agsVecShapeGemeinden = zeros(length(shapeDataGemeinden) ,1);
      for i=1:length(agsVecShapeGemeinden)
        agsVecShapeGemeinden(i) = shapeDataGemeinden{i}.ags;
        for j=1:length(shapeDataGemeinden{i}.coordinates)
          [coordsX, coordsY] = latLon2XY(shapeDataGemeinden{i}.coordinates{j}(:,2),...
          shapeDataGemeinden{i}.coordinates{j}(:,1));
          shapeDataGemeinden{i}.coordinates{j} = [coordsX, coordsY];
        end
      end 
    end
    
    for i=1:length(cities)    
      if isfield(cities{i}, "verbandschluessel")
        index = find(agsVecShapeGemeinden ==...
        cities{i}.ags*1e7+cities{i}.verbandschluessel*1e3+cities{i}.gemeindekennzahl);
        if index
          shapeDataGemeinde = shapeDataGemeinden{index};
          cities{i}.shapeData = shapeDataGemeinde;
        end   
      else
        index = find(agsVecShape == cities{i}.ags);
        if index
          cities{i}.shapeData = shapeData{index};
        end   
      end     
    end
  end
  
  
  
  if strcmp(parameter.model, "SIRH")
    x = expandSIRH(x, t, cities, parameter, initialHistory);
    csvwrite(strcat(folder,'/sirh-values.csv'), x);
  endif
  opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep, "AbsTol", 1e-12);
  if and(strcmp(parameter.model, "SIREDmod"), parameter.discoveredInVTK)    
    % calculate discovered
    darkFigures = sir_eqn_spatial ("getDarkFigure", parameter, cities);
    alpha = parameter.gamma1 ./ (darkFigures - 1);
    E = x(:,3:4:end);
    splineE = spline(t, E);
    if (parameter.startDate<737866)
      % SQID frueher( SIED )
      Q0 = x(1, 2:4:end);
      %%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
      % discovered = \int_0^t alpha * E * dT 
      % Integriere E*alpha to get discovered (realy infected)
      get_dDiscovering_dT = @(t_eq, x0) [alpha .* ppval(splineE, t_eq)];
      [t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)], Q0, opt);
    else
      Q0 = x(1, 2:4:end);
      E0 = x(1, 3:4:end);
      S0 = x(1, 1:4:end);
      D0 = x(1, 4:4:end);
      R0 = 1-(Q0+E0+S0+D0)+Q0;
      get_dDiscovering_dT = @(t_eq, x0) [alpha .* ppval(splineE, t_eq)];
      [t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)], R0, opt);
    end
    % get solutionDiscovered at the same timesteps as the other solutions
    solutionDiscoveredCW = spline(t2, discovered, t);
  else
    solutionDiscoveredCW = zeros(length(cities), length(t));
  endif
  
  % write out parameters
  names = fieldnames(parameter);
  fid = fopen([folder, "/parameters_out.txt"], 'w'); 
  for i=1:length(names)
    if ischar(parameter.(names(i){1}))
      fprintf(fid, "set %30s = %s\n", names(i){1}, parameter.(names(i){1}));
    elseif iscell(parameter.(names(i){1}))
      %nothing to do
    else  
      fprintf(fid, "set %30s = %s\n", names(i){1}, mat2str(parameter.(names(i){1})));
    end
  end
  fclose(fid);
  
  % count citizens
  N = 0;
  N_old = 0;
  for i=1:length(cities)
    N +=cities{i}.population * (1 - cities{i}.fracOld);
    N_old +=cities{i}.population * cities{i}.fracOld;
  end  
  
  switch parameter.model
    case "SIR"
      k=2;
    case {"SIRED","SIREDLiterature"}
      k=4;
    case "SIREDmod"
      k=4;
    case "SIREDYO"
      k=8;
    case "SIRH"
      %reshaped SIR contains R as well
      k=9;
  endswitch
  
  % build a matrix that contains the sums of all variables
  population = zeros(length(cities), 1);  
  for i=1:length(cities)
    population(i) = cities{i}.population;
  end
  sumMatrix = zeros(length(t), k);
  for i = 1:k
    sumMatrix(:, i) = x *   reshape([zeros(length(cities), i-1),...
    population, zeros(length(cities), k-i)]', size(x,2),1);
  endfor
  
  if parameter.showDiagrams
    ##    colorvector = ["-g";"-r";"-c";"-m";"-g";"-r";"-c";"-m"];
    ##    figure; 
    ##    hold on;
    ##    for i=1:length(cities)
    ##      for j = 1:k
    ##        plot(t, x(:,k*i-(k-j)), colorvector(j,:));
    ##      end
    ##    end
    ##    hold off;
    
    % statewise
    
    % total sum
    legendvector = {"S","I","R","E","D","S_o","I_o","R_o","E_o","D_o"};
    colorvector2 = ["-g";"-r";"-b";"-k";"-y";"-m";"-c";"-r";"-g";"-b"];
    switch parameter.model
      case "SIR"
        datavector = [sumMatrix, N-sum(sumMatrix(:,1:2), 2)];
      case "SIREDLiterature"
        datavector = [sumMatrix(:,1:2), N-sum(sumMatrix(:,1:2), 2), sumMatrix(:,3:4)];
      case {"SIRED","SIREDmod"}
        
        %%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
        % discovered = \int_0^t alpha * E * dT 
        alpha = parameter.gamma1 / (parameter.darkFigure - 1);
        splineE = spline(t, sumMatrix(:,3));
        % Integriere E*alpha to get discovered (realy infected)
        get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
        [t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)],...
        sumMatrix(1,2));
        % get solutionDiscovered at the same timesteps as the other solutions
        solutionDiscovered = spline(t2, discovered, t);
        
        datavector = [sumMatrix(:,1:2), N-sum(sumMatrix(:,1:4), 2),...
        sumMatrix(:,3:4), solutionDiscovered];
        legendvector = {"S","I","R","E","D","Discovered"};
      case "SIREDYO"
        % N = N_young, Nold
        datavector = [sumMatrix(:,1:2), N-sum(sumMatrix(:,1:4), 2),...
        sumMatrix(:,3:6), N_old - sum(sumMatrix(:,5:8), 2), sumMatrix(:, 7:8)];
      case "SIRH"
        datavector = sumMatrix;
        legendvector = {"Susceptible", "Exposed", "Removed", "Infectious",...
        "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
    endswitch
    
    
    figure();
    hold on;
    for i =1:size(datavector,2)
      plot(t(datavector(:,i) > 0),...
      datavector((datavector(:,i) > 0),i),colorvector2(i,:))
    endfor
    xlim([min(t), max(t)]);
    ##    ylim(max(max(datavector)) *[1e-6, 1]);
    ylim([0, 150000]);
    xlabel("Time in days","fontweight","bold");
    ylabel("Number","fontweight","bold");
    h = legend(legendvector{[,1:size(datavector, 2)]});
    legend(h,"show");
    hold off;    
  endif
  
  
  if 0
    if parameter.fullConsoleOut
      fprintf("Doing extra processing which has something to do with 2020-03-28.\n");
    end
    a=1;
    ##    solutionDiscoveredCW, 401 x timestep -> spline fuer Zeit, auswerten am passenden tag
    dayOfInterest = datenum(2020,03,28);
    t_daywise = linspace(min(t),max(t),(max(t)-min(t))+1);
    solutionDiscoveredCWtmp = spline(t, sumMatrix(:,9), t_daywise); 
    solutionDiscovereddayofinterest = solutionDiscoveredCWtmp(:,(dayOfInterest - parameter.startDate)+1);
    solutionDiscovereddayofinterestPopultaion =solutionDiscovereddayofinterest .*population;
    InfectedRKIdataAGSwise = read_data_rki(dayOfInterest); 
    infectedRKI = zeros(length(cities), 1);   
    agsRKI = zeros(length(InfectedRKIdataAGSwise),1);
    for k=1:length(InfectedRKIdataAGSwise)
      infectedRKI(k) = InfectedRKIdataAGSwise{k}.Infected;
      agsRKI(k) = InfectedRKIdataAGSwise{k}.ags;
    end
    agsCities = zeros(length(cities), 1);
    for k=1:length(cities)
      agsCities(k) = cities{k}.ags;
    end
    [SortedRKI, IndicesRKI] = sort (agsRKI);
    [SortedCities, IndicesCities] = sort (agsCities);
    
    
    %loop over all states
    for i=1:max(floor(SortedCities/1000))
      fprintf("State Number %i\n", i);
      IndicesCurrentStateRKI = and(SortedRKI>=(i*1000),SortedRKI<((i+1)*1000));
      IndicesCurrentStateCities = and(SortedCities>=(i*1000),SortedCities<((i+1)*1000));
      [r1, p1] = corrcoef(solutionDiscovereddayofinterestPopultaion(IndicesCurrentStateCities),...
      infectedRKI(IndicesCurrentStateRKI));
      j = length(r1);
      correlationCoeff = r1(1,j)
      pValue = p1(1,j)
    endfor
    
    %subplot(1,2,1), plot(solutionDiscoveredCWtmp(IndicesCities).*population(IndicesCities) ,infectedRKI(IndicesRKI), '*');
    %xlabel("model");
    %ylabel("RKI");
    %[r1, p1] = corrcoef(solutionDiscoveredCWtmp(IndicesCities).*population(IndicesCities) ,infectedRKI(IndicesRKI));
    %correlationCoeff = r1(1,2)
    %pValue = p1(1,2)
    
    ##    tmp = (solutionDiscoveredCWtmp(IndicesCities).*population - infectedRKI(IndicesRKI)) ./ infectedRKI(IndicesRKI);
    ##    [Sorted, Indices] = sort (tmp);
    ##    [r2, p2] = corrcoef(solutionDiscoveredCWtmp(...
    ##    IndicesCities(Indices(10:end-10))).*population(IndicesCities(Indices(10:end-10))) ,infectedRKI(IndicesRKI(Indices(10:end-10))));
    ##    correlationCoeff2 = r2(1,2);
    ##    p2 = p2(1,2);
    ##    subplot(1,2,2), plot(solutionDiscoveredCWtmp(...
    ##    IndicesCities(Indices(10:end-10))).*population(IndicesCities(Indices(10:end-10))) ,infectedRKI(IndicesRKI(Indices(10:end-10))), '*');
  end
  
  if parameter.saveDiagrams    
    addpath("Optimize");
    extractStatewiseResults(workspace, "workspace", "saveFigure",...
    [folder, "/", "dataTable.dat"]);
  end
  
  if parameter.saveVTK
    if parameter.OutputComparisonRKI
      error("Copmarable RKI data is always included now")
      ComparisonDates = [datenum(2020,03,28),datenum(2020,04,11),datenum(2020,04,25)];
      % assuming the minimum and maximum entries are integers
      n_days = parameter.totalRuntime+1;
      tdaywise = linspace(min(t), max(t), n_days);
      xUniform = spline(t, x, tdaywise )';
      discoveredUniform = spline(t, solutionDiscoveredCW, tdaywise)';
      
      if parameter.fullConsoleOut
        fprintf("\rPostprocessing 0/%i ", length(ComparisonDates));
      end
      for i=1:length(ComparisonDates)
        position = ComparisonDates(i)-parameter.startDate + 1;
        cities = setSIR(cities, xUniform(position,:), parameter);
        cities = setDiscovered(cities, discoveredUniform(position,:));
        InfectedRKIdataAGSwise = read_data_rki(ComparisonDates(i));
        %%
        %% convert infectedRKIdataAGSwise to vector
        InfectedRKI = zeros(length(InfectedRKIdataAGSwise),1);
        agsRKI = zeros(length(InfectedRKIdataAGSwise),1);
        for k=1:length(InfectedRKIdataAGSwise)
          InfectedRKI(k) = InfectedRKIdataAGSwise{k}.Infected;
          agsRKI(k) = InfectedRKIdataAGSwise{k}.ags;
        end
        cities = setInfected(cities, InfectedRKI, agsRKI);
        generateVTK(cities, i, folder, parameter,...
        workspace.distances, tdaywise(position));
        if parameter.fullConsoleOut
          fprintf("\rPostprocessing %i/%i ", i, length(ComparisonDates));
        end
      endfor
      
    else
      %old version
##      rkiDataHistory = read_history_data(parameter.startDate,1,max(t),parameter);
##      rkitime = rkiDataHistory.t;
##      InfectedRKI = rkiDataHistory.Infected;
##      agsRKI = rkiDataHistory.ags;
      AGSPopulation = zeros(2,length(cities));
      for i = 1:length(cities)
        AGSPopulation(1,i) = cities{i}.population;
        AGSPopulation(2,i) = cities{i}.ags;
      end
      
      if parameter.reduceToStates	== false
        filename_infected = '../../Daten/cases-rki-by-ags_current.csv';
        delimiterIn = ',';
        infected = importdata(filename_infected,delimiterIn);
        indexStart = find(infected(1,:) == parameter.startDate);
        indexEnd = find(infected(1,:) == parameter.startDate+parameter.totalRuntime);
      elseif parameter.reduceToStates == true
         dataRKI = read_case_history_RKIfiles();
         infected = zeros(17,length(dataRKI{1}.infected)+1);
         infected(1,2:end) = cell2mat(dataRKI{1}.time(:));
         infected(2:end,1) = AGSPopulation(2,:)';
         for i = 2:17
            infected(i,2:end) = cell2mat(dataRKI{i-1}.infected(:));
         end
      end
      
      agsRKI = zeros(size(infected,1)-1,1);
      InfectedRKI = zeros(size(infected,2)-1,size(infected,1)-1);
      RKI7DayIncidence = zeros(size(infected,2)-1,size(infected,1)-1);
     
      for i = 2:size(infected,1)
        agsRKI(i-1) = infected(i,1);
        RKIInf(i-1,:) = infected(i, 2:end);
        for j = 1:size(RKIInf,2)
          if j<=7
            IncidenceRKI(i-1,j) = (RKIInf(i-1,j)-RKIInf(i-1,1))*100000/AGSPopulation(1,i-1);
          else
            IncidenceRKI(i-1,j) = (RKIInf(i-1,j)-RKIInf(i-1,j-7))*100000/AGSPopulation(1,i-1);
          end
        end
      end
      rkitime = zeros(size(infected,2)-1,1);
      for i = 2:size(infected,2)
        rkitime(i-1) = infected(1,i)-parameter.startDate;
      end
      InfectedRKI = RKIInf';
      RKI7DayIncidence = IncidenceRKI';
      n_pics = parameter.resolutionVTK;  
      % make solution uniform in time
      switch parameter.model
        case "SIREDmod"
          Sim7DayIncidence = zeros(size(solutionDiscoveredCW,1),size(solutionDiscoveredCW,2));
          for i = 1:size(solutionDiscoveredCW,1)
            SimInf(i,:) = solutionDiscoveredCW(i, 1:end);
            for j = 1:size(solutionDiscoveredCW,2)
              if j<=7
                IncidenceSim(i,j) = (SimInf(i,j)-SimInf(i,1))*100000;
              else
                IncidenceSim(i,j) = (SimInf(i,j)-SimInf(i,j-7))*100000;
              end
            end
          end
          Sim7DayIncidence = IncidenceSim';
      case "SIRH"
        Sim7DayIncidence = zeros(size(x,2)/9,length(t));
          for i = 1:size(x,2)/9
            SimInf(i,:) = x(:, i*9);
            for j = 1:length(t)
              if j<=7
                IncidenceSim(i,j) = (SimInf(i,j)-SimInf(i,1))*100000;
              else
                IncidenceSim(i,j) = (SimInf(i,j)-SimInf(i,j-7))*100000;
              end
            end
          end
          Sim7DayIncidence = IncidenceSim';
       endswitch
      
      tUniform = linspace(min(t), max(t), n_pics);
      xUniform = spline(t, x, tUniform )';
      Sim7DayIncidenceUniform = spline(t, Sim7DayIncidence, tUniform)';
      discoveredUniform = spline(t, solutionDiscoveredCW, tUniform)';
     
      if max(t) > max(rkitime)
        for r = max(rkitime):1:max(t)
          rkitime(end+1) = r+1;
          InfectedRKI = [InfectedRKI; InfectedRKI(end,:)];
        endfor
      endif
      %tRKIUniform = linspace(min(rkitime), max(rkitime), n_pics);
      %RKIUniform = zeros(n_pics, size(InfectedRKI, 2));
      RKI7DayIncidenceUniform = spline(rkitime, RKI7DayIncidence, tUniform)';
      RKIUniform = spline(rkitime, InfectedRKI, tUniform )';
      
      ##    for i = 1:size(x, 2)
      ##      xUniform(:, i) = spline(t, x(:,i), tUniform);
      ##    end

      if parameter.fullConsoleOut
        fprintf("\rPostprocessing 0/%i ", length(tUniform));
      end
      for i = 1:length(tUniform)
        cities = setSIR(cities, xUniform(i,:), parameter);
        cities = setInfected(cities, RKIUniform(i,:), agsRKI);
        cities = setDiscovered(cities, discoveredUniform(i,:));
        cities = set7DayIncidenceRKI(cities, RKI7DayIncidenceUniform(i,:), agsRKI);
        cities = set7DayIncidenceSim(cities, Sim7DayIncidenceUniform(i,:));
        %generateGif(cities, false, i);

        generateVTK(cities, i, folder, parameter,...
        workspace.distances, tUniform(i));
        if parameter.fullConsoleOut
          fprintf("\rPostprocessing %i/%i ", i, length(tUniform));
        end
      end
      
    endif
  end
  if parameter.fullConsoleOut
    fprintf("\rPostprocessing took %is\n", toc(tstart));
  end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set sir in all cities for the following postprocessing
##function cities = setSIR(cities, X)
##  % X=[S_1, I_1, S_2, I_2 ,..., S_n, I_n] 
##  for i=1:length(cities)
##    cities{i}.SIR = [X(2*i-1), X(2*i), 1 - X(2*i-1)- X(2*i)];
##  end  
##endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set SIR/SEIDR in all cities for the following postprocessing
function cities = setSIR(cities, X, parameter)
  % X=[S_1, I_1, E_1, D_1, S_2, E_2, I_2, D_2, ..., S_n, E_n, I_n, D_n] 
  % number of unknowns in SIR model (except of R) so 2 for SIR and 4 for SEIDR  
  switch parameter.model
    case "SIR"
      k=2;
    case {"SIRED","SIREDmod","SIREDLiterature"}
      k=4;
    case "SIREDYO"
      k=8;
    case "SIRH"
      k=9;
  endswitch
  solMat = reshape(X, k, length(cities))';
  for i=1:length(cities)
    switch parameter.model
      case "SIR"
        cities{i}.SIR = [solMat(i, 1:2), 1 - sum(solMat(i, 1:2))];
      case {"SIRED","SIREDmod","SIREDLiterature"}
        cities{i}.SIR = [solMat(i, 1:2), 1 - sum(solMat(i, 1:4)),...
        solMat(i,3:4)];
      case "SIREDYO"
        % N = N_young, Nold
        cities{i}.SIR = [solMat(i, 1:2), (1 - cities{i}.fracOld) - sum(solMat(i, 1:4)),...
        solMat(i,3:6), cities{i}.fracOld - sum(solMat(i, 5:8)), solMat(i, 7:8)];
      case "SIRH"
        cities{i}.SIR = [solMat(i, :)];
    endswitch    
  endfor
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set Infected values of RKI in all cities
function cities = setInfected(cities, X, agsVec)
  k = 1;
  for i=1:length(cities)
    index = find(agsVec == cities{i}.ags);
    if index
      cities{i}.RKIInfected = X(index);  
    else 
      cities{i}.RKIInfected = 0;
      if length(cities) == length(agsVec)
        fprintf("ags not found while setting infected from rki: %s\n",...
        cities{i}.name);
      end
      ##      error("AGS not found.");
    end
  endfor
endfunction
function cities = set7DayIncidenceRKI(cities, X, agsVec)
  k = 1;
  for i=1:length(cities)
    index = find(agsVec == cities{i}.ags);
    if index
      cities{i}.SevenDayIncidenceRKI = X(index);  
    else 
      cities{i}.SevenDayIncidenceRKI = 0;
      if length(cities) == length(agsVec)
        fprintf("ags not found while setting infected from rki: %s\n",...
        cities{i}.name);
      end
      ##      error("AGS not found.");
    end
  endfor
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set Infected values of RKI in all cities
function cities = setDiscovered(cities, discovered)
  for i=1:length(cities)    
    cities{i}.discovered = discovered(i); 
  endfor
endfunction
function cities = set7DayIncidenceSim(cities, discovered)
  for i=1:length(cities)    
    cities{i}.SevenDayIncidenceSim = discovered(i); 
  endfor
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = latLon2XY(lat, lon)
  circumferenceAtLat = cos(lat.*0.01745329251).*111.305;
  x = lon.*circumferenceAtLat;
  y = lat.*110.919;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%