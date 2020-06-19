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
## @deftypefn {} {@var{retval} =} generateVTK (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-17

function generateVTK (cities, index, folder, parameter, distances, time)  
  filename = [folder, "/vtkTest_", sprintf("%i", index), ".vtk"];
  
  coords = zeros(length(cities), 2);
  radius = zeros(length(cities), 1);
  SIR = zeros(length(cities), length(cities{1}.SIR));
  population = zeros(length(cities), 1);
  populationDensity = zeros(length(cities), 1);
  infectedCumRKI = zeros(length(cities),1);
  discovered = zeros(length(cities),1);
  ags = zeros(length(cities), 1);	
	darkFiguresStatewise = sir_eqn_spatial ("getDarkFigureStatewise", parameter);
	darkFiguresCitywise = zeros(length(cities), 1);
	
	% Points per circle
  resolution = 25;
	trisPerCity = 1;	
	nElements = 0;
	nPoints = 0;	
  for i=1:length(cities)
    coords(i,:) = [cities{i}.x, cities{i}.y];
    radius(i) = (cities{i}.area)^0.5;
    population(i) = cities{i}.population;
    populationDensity(i) = cities{i}.population / cities{i}.area;
    ags(i) = cities{i}.ags;
    infectedCumRKI(i) = cities{i}.RKIInfected;
    discovered(i) = cities{i}.discovered;
    SIR(i, :) = cities{i}.SIR;
		darkFiguresCitywise(i) = darkFiguresStatewise(floor(cities{i}.ags/1000));
		
		if parameter.vtkStyleMap
			if isfield(cities{i}, "shapeData")
				for j=1:length(cities{i}.shapeData.coordinates)	
					nElements += 1;
					nPoints += size(cities{i}.shapeData.coordinates{j}, 1) - 1;
				end
			else				
				nElements += 1;
				nPoints += resolution;
			end
		end 
  end
  
	if !parameter.vtkStyleMap
		nElements = trisPerCity * length(cities);
		nPoints = resolution * length(cities);
 	end 
	
  x=zeros(1, nPoints);
  y=zeros(1, nPoints);
  z=zeros(1, nPoints);
  polys=cell(nElements, 1);
  vals=zeros(nPoints, size(SIR, 2)+6);	

	
	iSPoints = 1;
	iSEl = 1;
  for i=1:length(cities)
		cityValues = [darkFiguresCitywise(i), ags(i), populationDensity(i),...
		population(i), infectedCumRKI(i), SIR(i,:)];
		cityValues = [cityValues, discovered(i)];

		if and(parameter.vtkStyleMap, isfield(cities{i}, "shapeData"))
			[xtmp, ytmp, ztmp, polystmp, valstmp] = getPolygons(cities{i}, cityValues);
		else	
			[xtmp, ytmp, ztmp, polystmp, valstmp] = getCircle(coords(i,1), coords(i,2),...
			radius(i) / 2, cityValues, resolution);
			% zValues to avoid bigger circles hiding smaler ones
			ztmp += (max(radius) - radius(i))/max(radius);
		end
		
		% elements
		for j=1:length(polystmp)
			polys{iSEl + (j - 1)} = polystmp{j} + (iSPoints - 1);
		end	
		
    x(iSPoints:iSPoints + length(xtmp) - 1) = xtmp;
    y(iSPoints:iSPoints + length(xtmp) - 1) = ytmp;
    z(iSPoints:iSPoints + length(xtmp) - 1) = ztmp;
    vals(iSPoints:iSPoints + length(xtmp) - 1, :) = valstmp;
		
		iSEl += length(polystmp);
		iSPoints += length(xtmp);
  end 
	
  lines = [];
  
##  if strfind(parameter.spatial, "europeanCapitals")
##    % load countries cell array
##    % cell2mat(countries{i}.boarders{j}) for all i, j are the boarders
##    load("../../Daten/boarders60m.mat");
##    pointsBoarders = [];
##    linesBoarders = [];
##    for i=1:length(countries)
##      for j=1:length(countries{i}.boarders)
##        nPointsOld = size(pointsBoarders, 1);
##        [xLinePoints, yLinePoints] = latLon2XY(cell2mat(countries{i}.boarders{j}));
##        pointsBoarders = [pointsBoarders; [xLinePoints, yLinePoints]];
##        % new indices = [nPointsOld + 1 ; size(pointsBoarders, 1)]
##        nAdditionalPoints = size(cell2mat(countries{i}.boarders{j}), 1);
##        tmp = [1:nAdditionalPoints-1; 2:nAdditionalPoints]' + nPointsOld;
##        linesBoarders = [linesBoarders; tmp];
##      end
##    end
##    
##    linesBoarders += length(x);
##    vals = [vals; zeros(size(pointsBoarders, 1), size(vals, 2))];
##    pointsBoarders = [pointsBoarders, zeros(size(pointsBoarders, 1), 1)];
##    x = [x, pointsBoarders(:,1)'];
##    y = [y, pointsBoarders(:,2)'];
##    z = [z, pointsBoarders(:,3)'];
##    lines = linesBoarders;
##  end
  
	names = {"Darkfigure", "AGS", "Populationdensity","Population","InfectedRKI",...
			"Susceptible", "Quarantined", "Recovered", "Infected", "Dead", "Discovered"};

	
	if parameter.saveConnections
		population = zeros(length(cities), 1);
		biggestCity = zeros(length(cities), 1);
		coords = zeros(2, length(cities));
		for i=1:length(cities)
			population(i) = cities{i}.population;
			biggestCity(i) = cities{i}.biggestCityInCounty;
			coords(:,i) = [cities{i}.x, cities{i}.y];
		end
		maxBiggestCity = 3e6;
		interactionMatrix = sir_eqn_spatial("InteractionMatrix", nan, time, distances,...
		population, parameter, maxBiggestCity, cities);
		
		
		[m,n,val] = find(interactionMatrix);
		
		pointsConnections = zeros(2, 2 * length(m));
		linesConnections = zeros(2, length(m));
		for i=1:length(m)
			pointsIndices = 2*i-1:2*i;
			pointsConnections(:, pointsIndices) = coords(:, [m(i),n(i)]);
			linesConnections(:,i) = pointsIndices;
		end
		linesConnections += length(x);
		lines = [lines; linesConnections'];
		names = {"Traffic", names{:}};
		vals = [zeros(size(vals, 1), 1), vals];
		vals = [vals; [reshape([val, val]', 1, size(val,1) * 2)',...
		zeros(size(pointsConnections, 2), size(vals, 2)-1)]];
		pointsConnections = [pointsConnections; zeros(1, size(pointsConnections, 2))];
		x = [x, pointsConnections(1,:)];
		y = [y, pointsConnections(2,:)];
		z = [z, pointsConnections(3,:)];
		
	end
	
	if (parameter.saveConnections)
		vtkwrite(filename, 'polydata', 'triangle', x ,y, z, polys,...
		"lines", lines, "names", names, "values", vals);
	else
		vtkwrite(filename, 'polydata', 'polygon', x ,y, z, polys,...
		"names", names, "values", vals);
	end
endfunction
################################################################################
function [x, y, z, tri, vals] = getCircle(centerX, centerY, radius, values, resolution)  
	% first one is the center
	x = centerX + [cos(linspace(0, 2 * pi, resolution))] * radius;
	y = centerY + [sin(linspace(0, 2 * pi, resolution))] * radius;
	z = zeros(size(x));  
	tri{1} = 1:length(x);
	vals = ones(length(x), 1) * values;
endfunction
################################################################################
function [x, y, z, polys, vals] = getPolygons(city, values)	
	iSPoints = 1;
	vals = [];
	for i=1:length(city.shapeData.coordinates)
		xtmp = city.shapeData.coordinates{i}(1:end-1,1);
		ytmp = city.shapeData.coordinates{i}(1:end-1,2);
		
		polytmp = [1:length(xtmp)];
		
		polys{i} = polytmp + (iSPoints - 1);
		
    x(iSPoints:iSPoints + length(xtmp) - 1) = xtmp;
    y(iSPoints:iSPoints + length(xtmp) - 1) = ytmp;
    vals = [vals; ones(length(xtmp), 1) * values];
		
		iSPoints += length(xtmp);
	end
	z=zeros(size(x));
endfunction
################################################################################
function [x, y] = latLon2XY(coords)
	lat = coords(:,2);
	lon = coords(:,1);
	circumferenceAtLat = cos(lat*0.01745329251)*111.305;
	x = lon.*circumferenceAtLat;
	y = lat*110.919;
endfunction
################################################################################