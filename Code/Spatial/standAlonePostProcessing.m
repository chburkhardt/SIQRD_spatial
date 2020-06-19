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
		for i=1:length(cities)			
			index =	find(agsVecShape == cities{i}.ags);
			if index
				cities{i}.shapeData = shapeData{index};
			end
		end
	end
	
  
	
	if strcmp(parameter.model, "SIRH")
		x = expandSIRH(x, t, cities, parameter, initialHistory);
	endif
	opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep, "AbsTol", 1e-12);
	if and(strcmp(parameter.model, "SIREDmod"), parameter.discoveredInVTK)		
		% calculate discovered
		darkFigures = sir_eqn_spatial ("getDarkFigure", parameter, cities);
		alpha = parameter.gamma1 ./ (darkFigures - 1);
		% SQID frueher( SIED )
		E = x(:,3:4:end);
		Q0 = x(1, 2:4:end);
		%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
		% discovered = \int_0^t alpha * E * dT 
		splineE = spline(t, E);
		% Integriere E*alpha to get discovered (realy infected)
		get_dDiscovering_dT = @(t_eq, x0) [alpha .* ppval(splineE, t_eq)];
		[t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)], Q0, opt);
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
	
  k=4;
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
		% total sum
		legendvector = {"S","I","R","E","D","S_o","I_o","R_o","E_o","D_o"};
		colorvector2 = ["-g";"-r";"-b";"-k";"-y";"-m";"-c";"-r";"-g";"-b"];
		%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
		% discovered = \int_0^t alpha * E * dT 
		alpha = parameter.gamma1 / (parameter.darkFigure - 1);
		splineE = spline(t, sumMatrix(:,3));
		% integrate E*alpha to get discovered (realy infected)
		get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
		[t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)],...
		sumMatrix(1,2));
		% get solutionDiscovered at the same timesteps as the other solutions
		solutionDiscovered = spline(t2, discovered, t);
				
		datavector = [sumMatrix(:,1:2), N-sum(sumMatrix(:,1:4), 2),...
		sumMatrix(:,3:4), solutionDiscovered];
		legendvector = {"S","I","R","E","D","Discovered"};	
		
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
	
	
	if parameter.saveVTK   
		if parameter.OutputComparisonRKI
			ComparisonDates = [datenum(2020,03,28),datenum(2020,04,11),datenum(2020,04,25)];
			% assuming the minimum and maximum entries are integers
			n_days = parameter.totalRuntime+1;
			tdaywise = linspace(min(t), max(t), n_days);
			xUniform = spline(t, x, tdaywise )';
			discoveredUniform = spline(t, solutionDiscoveredCW, tdaywise)';
			
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
					fprintf("Postprocessing %i/%i\n", i, length(ComparisonDates));
				end
			endfor
			
		else
			rkiDataHistory = read_history_data(parameter.startDate,1,max(t));    
			rkitime = rkiDataHistory.t;
			InfectedRKI = rkiDataHistory.Infected;
			agsRKI = rkiDataHistory.ags;
			
			n_pics = parameter.resolutionVTK;  
			% make solution uniform in time
			tUniform = linspace(min(t), max(t), n_pics);
			xUniform = spline(t, x, tUniform )';
			discoveredUniform = spline(t, solutionDiscoveredCW, tUniform)';
			
			if max(t) > max(rkitime)
				for r = max(rkitime):1:max(t)
					rkitime(end+1) = r+1;
					InfectedRKI = [InfectedRKI; InfectedRKI(end,:)];
				endfor
			endif
			%tRKIUniform = linspace(min(rkitime), max(rkitime), n_pics);
			RKIUniform = zeros(n_pics, size(InfectedRKI, 2));
			RKIUniform = spline(rkitime, InfectedRKI, tUniform )';
			
			##    for i = 1:size(x, 2)
			##      xUniform(:, i) = spline(t, x(:,i), tUniform);
			##    end
			for i = 1:length(tUniform)
				cities = setSIR(cities, xUniform(i,:), parameter);
				cities = setInfected(cities, RKIUniform(i,:), agsRKI);
				cities = setDiscovered(cities, discoveredUniform(i,:));
				generateVTK(cities, i, folder, parameter,...
				workspace.distances, tUniform(i));
				if parameter.fullConsoleOut
					fprintf("Postprocessing %i/%i\n", i, length(tUniform));
				end
			end
			
		endif
	end
	if parameter.fullConsoleOut
		fprintf("Postprocessing took %is\n", toc(tstart));
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
	k = 4;
	solMat = reshape(X, k, length(cities))';
	for i=1:length(cities)
			cities{i}.SIR = [solMat(i, 1:2), 1 - sum(solMat(i, 1:4)),...
				solMat(i,3:4)];
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
			##			error("AGS not found.");
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = latLon2XY(lat, lon)
	circumferenceAtLat = cos(lat.*0.01745329251).*111.305;
	x = lon.*circumferenceAtLat;
	y = lat.*110.919;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%