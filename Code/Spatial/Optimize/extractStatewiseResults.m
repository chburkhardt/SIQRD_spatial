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
## @deftypefn {} {@var{retval} =} extractStatewiseResults (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-24

function stateWiseSIR = extractStatewiseResults (filenameOrWorkspace, varargin)
  addpath("..");
  if nargin == 0;
    filenameOrWorkspace = "../../../Results/result.mat";
  end
  stateNames= {"Schleswig-Holstein", "Freie und Hansestadt Hamburg",...
  "Niedersachsen", "Freie Hansestadt Bremen", "Nordrhein-Westfalen",...
  "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Freistaat Bayern",...
  "Saarland", "Berlin", "Brandenburg", "Mecklenburg-Vorpommern"...
  "Freistaat Sachsen", "Sachsen-Anhalt", "Freistaat Thueringen"};
  
  if and(nargin > 1, strcmp(varargin{1}, "workspace"))
    t = filenameOrWorkspace.x.x';
    x = filenameOrWorkspace.x.y';
    if isfield(filenameOrWorkspace , "initialHistory")
      initialHistory = filenameOrWorkspace.initialHistory;
    end    
    cities = filenameOrWorkspace.cities;
    parameter = filenameOrWorkspace.parameter;
  else    
    load(filenameOrWorkspace); 
    t = x.x';
    x = x.y'; 
  end 
  if strcmp(parameter.model, "SIRH")
    x = expandSIRH(x, t, cities, parameter, initialHistory);
  endif
  
  AGSvector = zeros(length(cities), 1);
  populationvector = zeros(length(cities), 1);
  for i=1:length(cities)
    AGSvector(i) = cities{i}.ags;
    populationvector(i) = cities{i}.population;
  end
  states = floor(AGSvector / 1000);
  
  switch parameter.model 
    case "SIR"
      k=4;
    case {"SIRED","SIREDmod", "SIREDLiterature"}
      k=5;
    case "SIREDYO"
      k=12;      
    case "SIRH"
      k=9;
  endswitch
  
  stateWiseSIR = cell(16, 1);
  for i=1:length(stateWiseSIR)
    stateWiseSIR{i}.SIR_vs_time = zeros(length(t), k);
    stateWiseSIR{i}.time = t;
    stateWiseSIR{i}.parameter = parameter;
    stateWiseSIR{i}.name = stateNames{i};
  end  
  
  for timestep = 1:length(t)
    SIRmatrix = getSIRmatrix(cities, x(timestep, :), parameter);
    for i=1:length(stateWiseSIR);
      stateWiseSIR{i}.SIR_vs_time(timestep, :) =...
      sum(SIRmatrix(states==i,:) .* populationvector(states==i), 1);
    end
  end 
	
	% integrate discovered for sired and siredmod
	if or(strcmp(parameter.model, "SIRED"), strcmp(parameter.model, "SIREDmod"))
		darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameter);
		opt = odeset();%'NormControl', 'on', 'MaxStep', parameter.maxTimeStep, "AbsTol", 1e-12);
		for i=1:length(stateWiseSIR);		
			%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
			% discovered = \int_0^t alpha * E * dT 
			alpha = parameter.gamma1 / (darkFigures(i) - 1);
			splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
			% Integriere E*alpha to get discovered (realy infected)
			get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
      if (parameter.startDate<737866)
        [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
        max(stateWiseSIR{i}.time)], stateWiseSIR{i}.SIR_vs_time(1,2), opt);
      else 
        [t2, discovered] = ode23(get_dDiscovering_dT, [min(stateWiseSIR{i}.time),...
        max(stateWiseSIR{i}.time)], stateWiseSIR{i}.SIR_vs_time(1,3)/darkFigures(i)...
        +stateWiseSIR{i}.SIR_vs_time(1,2), opt);
      end
			% get solutionDiscovered at the same timesteps as the other solutions
			solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);
			stateWiseSIR{i}.SIR_vs_time = [stateWiseSIR{i}.SIR_vs_time, solutionDiscovered];
		end
  elseif strcmp(parameter.model, "SIREDLiterature")
    darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameter);
    opt = odeset();
    splineE = spline(stateWiseSIR{i}.time, stateWiseSIR{i}.SIR_vs_time(:,4));
    alpha = 1/3;
    get_dDiscovering_dTL = @(t_eq, x0) [alpha*ppval(splineE, t_eq)];
    if (parameter.startDate<737866)
      [t2, discovered] = ode23(get_dDiscovering_dTL, [min(stateWiseSIR{i}.time),...
       max(stateWiseSIR{i}.time)],stateWiseSIR{i}.SIR_vs_time(1,2));
    else
      [t2, discovered] = ode23(get_dDiscovering_dTL, [min(stateWiseSIR{i}.time),...
      max(stateWiseSIR{i}.time)], (stateWiseSIR{i}.SIR_vs_time(1,3)/darkFigures(i)...
          +stateWiseSIR{i}.SIR_vs_time(1,2)));
    end
	  solutionDiscovered = spline(t2, discovered, stateWiseSIR{i}.time);
		stateWiseSIR{i}.SIR_vs_time = [stateWiseSIR{i}.SIR_vs_time, solutionDiscovered];
	end
  
  if nargin > 1
    for i=1:length(varargin)
      switch varargin{i}
        case "plot"
          plotStates(stateWiseSIR, parameter);
        case "saveFigure"
          filename = varargin{i+1};
          writeTexTable(stateWiseSIR, parameter, filename);
      endswitch        
    endfor
  endif
  
  
endfunction

function plotStates(stateWiseSIR, parameter)
  switch parameter.model
    case "SIR"
      titles = {"Susceptible", "Infected", "Removed", "Cumsum Infected"};
    case "SIRED"
      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
      "Discovered"};
    case "SIREDmod"
      titles = {"Susceptible", "Quarantined", "Recovered", "Infected", "Dead",...
      "Discovered"};
    case "SIREDLiterature"
       titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
      "Discovered"};
    case "SIREDYO"
      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
      "Cumsum Infected", "Susceptible risk", "Infected risk", "Recovered risk",...
      "Exposed risk", "Dead risk", "Cumsum Infected risk"};
		case "SIRH"
			%SIRH = [S, X, X, effectiveInfective, symptoms, hospital, icu, dead, discovered];
      titles = {"Susceptible", "Exposed", "Removed", "Infectious",...
      "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
  endswitch
  
  stateNames = cell(length(stateWiseSIR), 1);
  ncols = ceil(length(titles)^0.5);
  nrows = ceil(length(titles)/ncols);
  figure;
  for i=1:length(titles)
    for j=1:length(stateWiseSIR)
      subplot(nrows, ncols, i), hold on;
      plot(stateWiseSIR{j}.time, stateWiseSIR{j}.SIR_vs_time(:,i));
      stateNames{j} = stateWiseSIR{j}.name;
    end
    hold off;
    title(titles{i});
  end
  ##  legend(stateNames, "location", "northeastoutside");
end
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
function solMat = getSIRmatrix(cities, X, parameter)
  switch parameter.model
    case "SIR"
      k=2;
    case {"SIRED","SIREDmod", "SIREDLiterature"}
      k=4;
    case "SIREDYO"
      k=8;
      fracOld = zeros(length(cities), 1);
      for i=1:length(cities)
        fracOld(i) = cities{i}.fracOld;
      end
    case "SIRH"
      k=9;
  endswitch
  
  solMat = reshape(X, k, length(cities))';   
  
  ##    case "SIR"
  ##      titles = {"Susceptible", "Infected", "Removed", "Cumsum Infected"};
  ##    case "SIRED"
  ##      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
  ##      "Cumsum Infected"};
  ##    case "SIREDYO"
  ##      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
  ##      "Cumsum Infected", "Susceptible risk", "Infected risk", "Recovered risk",...
  ##      "Exposed risk", "Dead risk"};
  ##    case "SIRH"
  ##      titles = {"Susceptible", "Exposed", "Removed", "Infectious",...
  ##      "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
  
  switch parameter.model
    case "SIR"
##      solMat = [solMat(:, 1:2), 1 - sum(solMat(:, 1:2), 2), 1 - solMat(:,1)];
      solMat = [solMat(:, 1:2), 1 - solMat(:, 1) - solMat(:,2), 1 - solMat(:,1)];
    case {"SIRED","SIREDmod"}
      solMat = [solMat(:, 1:2), 1 - sum(solMat(:, 1:4), 2),...
      solMat(:,3:4)];
    case {"SIREDLiterature"}
      solMat = [solMat(:, 1:2), 1 - sum(solMat(:, 1:4), 2),...
      solMat(:,3:4)];
      %solMat = [solMat(:, 1:2), 1 - solMat(:, 1) - sum(solMat(:, 2:4), 2),...
    case "SIREDYO"
      % N = N_young, Nold
      solMat = [solMat(:, 1:2), (1 - fracOld) - sum(solMat(i, 1:4), 2),...
      solMat(i,3:4), (1 - fracOld) -solMat(:,1),...
      solMat(i,5:6), cities{i}.fracOld - sum(solMat(i, 5:8), 2),...
      solMat(i, 7:8), fracOld - solMat(:,5)];
    case "SIRH"
      % nothing needs to be done
  endswitch    
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%