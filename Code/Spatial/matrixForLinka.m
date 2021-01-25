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
## @deftypefn {} {@var{retval} =} matrixForLinka (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-08-21

function retval = matrixForLinka (input1, input2)
	
	
	parameter = read_parameter("../../Results/parameter.txt"){1}; 
	cities = setup_system(parameter);
	[spatialDistinction, reduction] = reduceCities(cities, parameter);
	
	
	[distances, population, biggestCityInCounty] =...
	calcDistancesPopulation(spatialDistinction, 10000, parameter);  
	
	interactionMatrix = constantPartMatrix(spatialDistinction, distances,...
	population, biggestCityInCounty);
	
	agsVec = zeros(length(cities), 1);
	for i=1:length(cities)		
		agsVec(i) = cities{i}.ags*1e7+cities{i}.verbandschluessel*1e3+cities{i}.gemeindekennzahl;
	endfor
	
	names = {};
	for i=1:length(spatialDistinction)
		names{i} = spatialDistinction{i}.name;
	endfor
	
	fid = fopen("../../Results/names.txt", "w");
	fprintf(fid, "%s\n", names{:});
	fclose(fid);
	
	fid = fopen("../../Results/population.txt", "w");
	fprintf(fid, "%i\n", population);
	fclose(fid);
	
	csvwrite("../../Results/interactionMatrix.txt", interactionMatrix);
	
	fid = fopen("../../Results/reduction.txt", "w");
	for i=1:length(reduction)
		if length(reduction{i}) > 1	
			fprintf(fid, "%i,", agsVec(reduction{i}(1:end-1)));
		endif
		fprintf(fid, "%i\n", agsVec(reduction{i}(end)));
	endfor
	fclose(fid);
	
	
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = constantPartMatrix(spatialDistinction, distances,...
	population, biggestCityInCounty)
	
	retval = eye(size(distances));
	
	exp_alpha = [0.46, 0.35];
	exp_gamma = [0.64, 0.37];
	factor_r = [82, 1];
	rcutoff = 300;
	maxBiggestCityInCounty = 3e6;
	
	
	for i = 1:size(distances, 1) - 1 
		
		tmp_sparse = [spdiags(distances, i)(i+1:end);...
		spdiags(distances, i - length(distances))(1:i)];
		indices = tmp_sparse > 0;   
		
		betaCrossfull = zeros(length(indices),1);
		index_low = and((tmp_sparse < rcutoff), (tmp_sparse > 0));
		index_up = (tmp_sparse >= rcutoff);
		
		betaCrossfull (index_low) = (population(index_low).^exp_alpha(1))...
		.* (shift(population, -i)(index_low).^exp_gamma(1))...
		./ exp(tmp_sparse(index_low) ./factor_r(1)) /...
		maxBiggestCityInCounty^(exp_gamma(1) + exp_alpha(1));
		
		betaCrossfull (index_up) = (population(index_up).^exp_alpha(2))...
		.* (shift(population, -i)(index_up).^exp_gamma(2))...
		./ exp(tmp_sparse(index_up) ./factor_r(2)) /...
		maxBiggestCityInCounty^(exp_gamma(2) + exp_alpha(2));
		
		% dSmdT abhängig von E in stadt m+i
		tmp = zeros(size(betaCrossfull));
		tmp(indices) = betaCrossfull;
		##							tmp ./= tmp;
		##							tmp *= i;
		##							tmp+=(1:16)';
		for m = 1:length(tmp)
			retval(m, mod(m+i-1, length(spatialDistinction)) + 1) = tmp(m);
		end
	endfor
	
	
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [counties, reduction] = reduceCities(cities, parameter)
	%% citiesMatrix(i, j)
	% cities(i)
	% j= 1: population; 2: x; 3:y; 4: ags, 5: area, 6: fracOld, 7_till_:end : SIR
	citiesMatrix = zeros(length(cities), 6 + length(cities{1}.SIR));
	
	for i=1:length(cities)
		citiesMatrix(i,:)=[cities{i}.population, cities{i}.x, cities{i}.y,...
		cities{i}.ags*1e7+cities{i}.verbandschluessel*1e3+cities{i}.gemeindekennzahl,...
		cities{i}.area, cities{i}.fracOld, cities{i}.SIR];
	end  
  % ags vector contains values from 1 to 16/17  
  agsStates = unique(floor(citiesMatrix(:,4)/1e10));
	
  agsCounties = unique(floor(citiesMatrix(:,4)/1e7));
	agsBayern = agsCounties(floor(agsCounties/1000)==9);
	ags = [agsStates(agsStates!=9); agsBayern];
	
	
	%% cell array that hold all cities of county i in reduction{i}
	reduction = cell(1, length(ags));
	counties = cell(1, length(ags));
	for i=1:length(reduction)
		if ags(i)<100
			reduction{i} = find(floor(citiesMatrix(:,4)/1e10) == ags(i));
		else
			reduction{i} = find(floor(citiesMatrix(:,4)/1e7) == ags(i));
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
