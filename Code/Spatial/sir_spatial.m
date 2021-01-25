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
    ##    if !exist("citiesIn", "var")
    ##      persistent citiesIn = setup_system(parameter);  
    ##    end
    ##    cities = citiesIn; 
    cities = setup_system(parameter);
    
    parameter.dataStatus = cities{1}.dataStatus;
    
    % reduction to counties or not
    if parameter.reduce
      [spatialDistinction, reduction] = reduceCities(cities, parameter); 
    else
      spatialDistinction = cities;
    end
    
    % t in days
    t = [0, parameter.totalRuntime];  
    
    % calculate distances between cities and X0
    % X=[S_1, I_1, S_2, I_2 ,..., S_n, I_n] 
    [distances, population, biggestCityInCounty] =...
    calcDistancesPopulation(spatialDistinction, 100, parameter);  
    X_0 = getSIR(spatialDistinction, parameter);
    
    %% solve ode and save solution
    t1 = tic;
    opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
    switch parameter.model
      case {"SIR", "SIRED","SIREDmod", "SIREDYO","SIREDLiterature"}
        if parameter.analytically
          x = analyticalModel (spatialDistinction, parameter, distances, X_0, t);
        else
          get_dXdT = @(t_eq, x0) sir_eqn_spatial("rates", x0, t_eq, distances,...
          population, parameter, biggestCityInCounty, spatialDistinction);
          x = ode45(get_dXdT, t, X_0, opt);
        endif
      case "SIRH"
        get_dXdT = @(t_eq, x0, hist) [sir_eqn_spatial("rates", x0, t_eq, distances,...
        population, parameter, biggestCityInCounty, spatialDistinction, hist)];        
        lags = sir_eqn_spatial("lags");
        [initialHistory, X_0] = getHistory(X_0, parameter, 40, spatialDistinction,...
        distances, population, biggestCityInCounty);
        if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate>737881)
          initialHistory = zeros(length(spatialDistinction),40);
          X_0 = zeros(length(spatialDistinction),1);
          
          darkFigures = sir_eqn_spatial ("getDarkFigureStatewise", parameter);
          courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);
          
          for i=1:length(spatialDistinction)
            initialHistory(i,:) = 1-(spatialDistinction{i}.history...
            /spatialDistinction{i}.population);%*darkFigures(i)); %ones(1,40)
          end
          X_0 = initialHistory(:,end);
        end
        x = ode23d_fixed(get_dXdT, t, X_0', lags, initialHistory, opt);
    otherwise
      error("Unknown model");
  endswitch
  if parameter.fullConsoleOut
    fprintf("\rode45 (parameter optimization) took %.2fs\n", toc(t1));
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

% get vector of SIR/ SEID / ... of all cities
function X = getSIR(cities, parameter)
% number of unknowns in SIR model (except of R) so 2 for SIR, 4 for SEIDR
% and 8 for SsEeIiDdRr
switch parameter.model
  case "SIRH"
    k=1;
  case "SIR"
    k=2;
  case {"SIRED","SIREDmod","SIREDLiterature"}
    k=4;
  case "SIREDYO"
    k=8;
endswitch
% build matrix with [nCities; k]
Xmatrix = zeros(length(cities), k);
relevantIndices = [1,2,4,5,6,7,9,10](1:k);
for i=1:length(cities)
  Xmatrix(i,:) = cities{i}.SIR(relevantIndices);
end
##The elements of the matrix are accessed in column-major order (like
##Fortran arrays are stored).
X = reshape(Xmatrix', numel(Xmatrix), 1);
% X=[S_1, I_1, S_2, I_2, ..., S_n, I_n] 
% X=[S_1, I_1, E_1, D_1, S_2, I_2, E_2, D_2, ..., S_n, I_n, E_n, D_n] 
% X=[Sy_1, Iy_1, Ey_1, Dy_1, So_1, Io_1, Eo_1, Do_1, 
%    Sy_2, Ey_2, Io_2, Do_2,  So_2, Io_2, Eo_2, Do_2, 
%     ..., Sy_n, Ey_n, Io_n, Do_n,  So_n, Io_n, Eo_n, Do_n] 
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
function showConnections(cities, distances)
[m, n, val] = find(distances);
population = zeros(length(cities), 1);
for i=1:length(cities)
  population(i) = cities{i}.population;
end

traffic = zeros(length(cities), 1);
for i=1:length(cities)
  traffic(i) = sum(distances(distances(:,i)>0, i).^(-2))*population(i);
end

traffic /= max(traffic);

for i=1:length(cities)
  cities{i}.SIR = [0, traffic(i), 0];
end
generateGif(cities, false, 0);

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [counties, reduction] = reduceCities(cities, parameter)
%% citiesMatrix(i, j)
% cities(i)
% j= 1: population; 2: x; 3:y; 4: ags, 5: area, 6: fracOld, 7_till_:end : SIR
citiesMatrix = zeros(length(cities), 6 + length(cities{1}.SIR));
devAgeAverageVec= zeros(length(cities),1);
if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737881)
  histMatrix = zeros(length(cities), 40);
  deathVec = zeros(length(cities),1);
end
agsVectorCities = zeros(length(cities), 1);
population = zeros(length(cities), 1);        
for i=1:length(cities)
  agsVectorCities(i) = cities{i}.ags;
  population(i) = cities{i}.population;
end
for i=1:length(cities)
  citiesMatrix(i,:)=[cities{i}.population, cities{i}.x, cities{i}.y,...
  cities{i}.ags*1e7+cities{i}.verbandschluessel*1e3+cities{i}.gemeindekennzahl,...
  cities{i}.area, cities{i}.fracOld, cities{i}.SIR];
  devAgeAverageVec(i) = cities{i}.devAgeAverage*cities{i}.population;
##  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737908)
##    n_city = 0;
##    indags = find(agsVectorCities == cities{i}.ags);
##    indags = [indags];
##    for k = indags(1):indags(end)
##      n_city += cities{k}.population;
##    end
##    histMatrix(i,:) = cities{i}.history*n_city;
##    deathVec(i) = cities{i}.deathStart*n_city;
##    for l=1:length(cities{i}.SIR)
##      citiesMatrix(i,6+l) = citiesMatrix(i,6+l)*n_city;
##    end
##  end
  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737881)
    n_city = cities{i}.population;
##    indags = find(agsVectorCities == cities{i}.ags);
##    indags = [indags];
##    for k = indags(1):indags(end)
##      n_city += cities{k}.population;
##    end
    %if strcmp(parameter.model,"SIRH")
      histMatrix(i,:) = cities{i}.history*n_city;
      deathVec(i) = cities{i}.deathStart*n_city;
    %end
    for l=1:length(cities{i}.SIR)
      citiesMatrix(i,6+l) = citiesMatrix(i,6+l)*n_city;
    end
  end
end  
if parameter.reduceToStates
  % ags vector contains values from 1 to 16/17  
  ags = unique(floor(citiesMatrix(:,4)/1e10));
else
  if isfield(parameter, "keepAGS")    
    toKeep = [];
    ags = unique(floor(citiesMatrix(:,4)/1e7));
    for i=1:length(parameter.keepAGS)
      toKeep = [toKeep;...
      citiesMatrix(floor(citiesMatrix(:,4)/1e7)==parameter.keepAGS(i), 4)];
      
      ags = ags(ags!=parameter.keepAGS(i));
    end
  else  
    ags = unique(floor(citiesMatrix(:,4)/1e7));
  end
end

%% cell array that hold all cities of county i in reduction{i}
reduction = cell(1, length(ags));
counties = cell(1, length(ags));
for i=1:length(reduction)
  if parameter.reduceToStates
    reduction{i} = find(floor(citiesMatrix(:,4)/1e10) == ags(i));
  else
    reduction{i} = find(floor(citiesMatrix(:,4)/1e7) == ags(i));
  end
  indexFirstBiggestCity =  find(citiesMatrix (reduction{i}, 1) ==...
  max(citiesMatrix (reduction{i},1)))(1);
  county.name = cities{reduction{i}(indexFirstBiggestCity)}.name;
  county.population = sum(citiesMatrix(reduction{i}, 1));
  county.devAgeAverage = sum(devAgeAverageVec(reduction{i}))/county.population;
  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737881)
    county.history = sum(histMatrix(reduction{i}, :),1);  
    county.deathStart = sum(deathVec(reduction{i}));
    %county.devAgeAverage = sum(devAgeAverageVec(reduction{i}))/county.population;
  end
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
  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737881)
    for j=1:length(SIR)
      %SIR(j) = sum(citiesMatrix(reduction{i}, j+6).*citiesMatrix(reduction{i}, 1))...
      %/county.population;
      SIR(j) = sum(citiesMatrix(reduction{i}, j+6))/county.population;
    end
  else
    for j=1:length(SIR)
      SIR(j) = sum(citiesMatrix(reduction{i}, j+6).*citiesMatrix(reduction{i}, 1))...
      /county.population;
      %SIR(j) = sum(citiesMatrix(reduction{i}, j+6))/county.population;
    end
  end
  county.SIR = SIR;
  county.fracOld = sum(citiesMatrix(reduction{i}, 6).*citiesMatrix(reduction{i}, 1))...
  /county.population;
  if parameter.reduceToStates
    % ags vector is [1; 17]   
    county.ags = ags(i) * 1000;
  else
    % ags vector is [1; 17]*1e3
    county.ags = ags(i);
  end 
  counties{i} = county;
end  
% add cities that are not reduced
if exist("toKeep", "var")
  for i=1:length(toKeep)
    counties{end+1} =  cities{find(citiesMatrix (:,4) == toKeep(i))};
    reduction{end+1} = find(citiesMatrix (:,4) == toKeep(i));
  end
end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = blowUpCounties(reduction, cities, x_red, parameter)
switch parameter.model
  case "SIRH"  
    nValuesPerCity = 1;
  case "SIR"     
    nValuesPerCity = 2;
  case {"SIRED","SIREDmod","SIREDLiterature"}
    nValuesPerCity = 4;
  case "SIREDYO" 
    nValuesPerCity = 8;
endswitch
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [history, X_0] = getHistory(X_0, parameter, size, cities,...
distances, population, biggestCityInCounty)
lags = sir_eqn_spatial("lags");
courses = sir_eqn_spatial("totalCourseOfDisease", parameter);
nu = 0.345;
darkFigure = sir_eqn_spatial ("getDarkFigure", parameter, cities);
betas = sir_eqn_spatial("betas", parameter, 0, cities);

lags2 = (lags(1:end-1)+lags(2:end))/2;

discoveredFun = @(tau) exp(-nu*tau).*interp1(lags2, courses.discovering, tau, "linear");
factor = 1./(integral(discoveredFun, 0.5,27.5) ./ darkFigure * (1-exp(-nu*28)));

% x_0 is equal to the number of discovered infected
history = 1-((1-X_0).*factor).*exp(-nu*linspace(max(lags), min(lags), size));
X_0 = history(:,end);
endfunction














