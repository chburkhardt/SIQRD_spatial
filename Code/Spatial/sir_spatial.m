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
## @deftypefn {} {@var{retval} =} sir_spatial (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-03

function out = sir_spatial (parameterArray)
  if nargin == 0
    close all;
    parameterArray = read_parameter("../../Results/parameter.txt"); 
  end  
  for run = 1:length(parameterArray)
    parameter = parameterArray{run};
		cities = setup_system(parameter); 
    
    % reduction to counties or not
    if parameter.reduce
      [spatialDistinction, reduction] = reduceCities(cities, parameter); 
    else
      spatialDistinction = cities;
    end
    
    % t in days
    t = [0, parameter.totalRuntime];  
    
    % calculate distances between cities and X0
    [distances, population, biggestCityInCounty] =...
    calcDistancesPopulation(spatialDistinction, 100, parameter);  
    X_0 = getSIR(spatialDistinction, parameter);  
    
    %% solve ode and save solution
    t1 = tic;
    opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
		get_dXdT = @(t_eq, x0) sir_eqn_spatial("rates", x0, t_eq, distances,...
      			population, parameter, biggestCityInCounty,spatialDistinction);
		x = ode45(get_dXdT, t, X_0, opt);

  if parameter.fullConsoleOut
    fprintf("ode45 took %fs\n", toc(t1));
  end
  
  % choose if blow up reduction from counties to cities if reduced
  if and(parameter.reduce, parameter.blowUp)
    x.y = blowUpCounties(reduction, cities, x.y, parameter);
    if exist("initialHistory", "var")
      initialHistory = blowUpCounties(reduction, cities, initialHistory, parameter);
    end
  else
    cities = spatialDistinction;
  endif
  
  parent = "../../Results/";
  if strcmp(parameter.folderName, "auto")
    folder = strftime("%y_%m_%d_%H_%M_%S", localtime(time()));  
  else
    folder = parameter.folderName;
  end
  filename = "result.mat";  
  mkdir(parent, folder);
	if exist("initialHistory", "var")
		save([parent, folder, "/", filename],...
		"x", "cities", "initialHistory", "parameter",...
		"distances");
	else
		save([parent, folder, "/", filename], "x", "cities", "parameter",...
		"distances");
	end
  
  % postprocess data
  workspace.x = x;
  workspace.parameter = parameter;
  workspace.distances = distances;
  workspace.biggestCityInCounty = biggestCityInCounty;
  workspace.population = population;
  if exist("initialHistory", "var")
    workspace.initialHistory = initialHistory;
  end
	if exist("reduction", "var")
		workspace.reduction = reduction;
	end	
  workspace.cities = cities;
  standAlonePostProcessing([parent, folder], filename, workspace); 
end  
if nargout == 1
  out = workspace;
  if length(parameterArray) > 1
    warning("For more than one run, the outputfiles should be read from hdd");
  end    
end 
endfunction

% get vector of SIREDmod of all cities
function X = getSIR(cities, parameter)
% number of unknowns in SIREDmod model (except of R) so 4
k=4;
% build matrix with [nCities; k]
Xmatrix = zeros(length(cities), k);
relevantIndices = [1,2,4,5,6,7,9,10](1:k);
for i=1:length(cities)
  Xmatrix(i,:) = cities{i}.SIR(relevantIndices);
end
##The elements of the matrix are accessed in column-major order (like
##Fortran arrays are stored).
X = reshape(Xmatrix', numel(Xmatrix), 1);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the distance matrix with (i,j) the distance from city_i to city_j
function [distances, population, biggestCityInCounty] =...
calcDistancesPopulation(cities, cutoffRadius, parameter)
t1 = tic;
% population vector
population = zeros(length(cities), 1);
biggestCityInCounty = zeros(length(cities), 1);
x = zeros(length(cities), 1);
y = zeros(length(cities), 1);

% distances matrix as sparse matrix
m=[];
n=[];
val=[];

if nargin == 0
  cutoffRadius = 10;
end

for i = 1:length(cities)
  population(i) = cities{i}.population;
  biggestCityInCounty(i) = cities{i}.biggestCityInCounty;
  x(i) = cities{i}.x;
  y(i) = cities{i}.y;
end
largestCity = 3e6;

for i = 1:length(cities)
  distance = ((x(i) - x).^2 + (y(i) - y).^2).^0.5;
  
  [nVector] = or(distance < cutoffRadius,...
  distance ./ max(population(i), population)...
  < (cutoffRadius / (largestCity / 40)));
  
  m = [m; i * ones(sum(nVector), 1)];
  n = [n; find(nVector)];
  val = [val; distance(nVector)];
end

distances = sparse(m,n,val);

if parameter.fullConsoleOut
  fprintf("%i/%i connections, calculation took %is\n", length(m),...
  length(cities)^2, toc(t1));
end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [counties, reduction] = reduceCities(cities, parameter)
%% citiesMatrix(i, j)
% cities(i)
% j= 1: population; 2: x; 3:y; 4: ags, 5: area, 6: fracOld, 7_till_:end : SIR
citiesMatrix = zeros(length(cities), 6 + length(cities{1}.SIR));

for i=1:length(cities)
  citiesMatrix(i,:)=[cities{i}.population, cities{i}.x, cities{i}.y,...
  cities{i}.ags, cities{i}.area, cities{i}.fracOld, cities{i}.SIR];
end  
if parameter.reduceToStates
  ags = unique(floor(citiesMatrix(:,4)/1000));
else
  ags = unique(citiesMatrix(:,4));
end

%% cell array that hold all cities of county i in reduction{i}
reduction = cell(1, length(ags));
counties = cell(1, length(ags));
for i=1:length(reduction)
  if parameter.reduceToStates
    reduction{i} = find(floor(citiesMatrix(:,4)/1000) == ags(i));
  else
    reduction{i} = find(citiesMatrix(:,4) == ags(i));
  end
  indexFirstBiggestCity =  find(citiesMatrix (reduction{i}, 1) ==...
  max(citiesMatrix (reduction{i},1)))(1);
  county.name = cities{reduction{i}(indexFirstBiggestCity)}.name;
  county.population = sum(citiesMatrix(reduction{i}, 1));  
  if length(reduction{i}) > 1
    county.biggestCityInCounty = max(citiesMatrix(reduction{i}, 1));
  else
    county.biggestCityInCounty = cities{reduction{i}}.biggestCityInCounty;
  end  
  county.x = sum(citiesMatrix(reduction{i}, 2).*citiesMatrix(reduction{i}, 1))...
  /county.population;
  county.y = sum(citiesMatrix(reduction{i}, 3).*citiesMatrix(reduction{i}, 1))...
  /county.population;
  county.area = sum(citiesMatrix(reduction{i}, 5));
  SIR = zeros(1,length(cities{1}.SIR));
  for j=1:length(SIR)
    SIR(j) = sum(citiesMatrix(reduction{i}, j+6).*citiesMatrix(reduction{i}, 1))...
    /county.population;
  end
  county.SIR = SIR;
  county.fracOld = sum(citiesMatrix(reduction{i}, 6).*citiesMatrix(reduction{i}, 1))...
  /county.population;
  if parameter.reduceToStates
		county.ags = ags(i) * 1000;
	else
		county.ags = ags(i);
	end
  counties{i} = county;
end  
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = blowUpCounties(reduction, cities, x_red, parameter)
nValuesPerCity = 4;
x = zeros(nValuesPerCity * length(cities), size(x_red, 2));

% add entries countywise from reduced to full solution matrix
for i=1:length(reduction)
  indicesFull = reduction{i};
  tmp = [];
  for j=1:nValuesPerCity
    tmp = [tmp, (nValuesPerCity * (indicesFull - 1) + j)];
  end
  
  indicesFull = reshape(tmp', nValuesPerCity * length(indicesFull), 1);    
  ##    indicesFull = reshape([(2 * indicesFull - 1), (2 * indicesFull)]',...
  ##    2 * length(indicesFull), 1);
  
  indicesReduced = [];
  for j=1:nValuesPerCity
    indicesReduced = [indicesReduced; (nValuesPerCity * (i - 1) + j)];      
  end
  %indicesReduced = [(2 * i - 1); (2 * i)];    
  %indicesReduced = [(2 * (i - 1) + 1); (2 * (i - 1) + 2)];    
  
  ##    for j = 1:length(indicesFull)
  ##      x(indicesFull(j), :) = x_red(indicesReduced( -(mod(j,2)-2)),:);
  ##    end   
  
  for j = 1:length(indicesFull)
    % indicesReduced is length nValuesPerCity
    % indicesFull is numberOfCitiesInCounty x nValuesPerCity
    entryInSIR = mod(j - 1, nValuesPerCity) + 1;
    x(indicesFull(j), :) = x_red(indicesReduced(entryInSIR), :);
  end       
  
end
endfunction















