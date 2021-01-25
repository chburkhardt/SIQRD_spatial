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
## @deftypefn {} {@var{retval} =} checkInteractionMatrix (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-05-21

function retval = checkInteractionMatrix (input1, input2) 
  parameterArray = read_parameter('../../Results/parameter.txt');
  parameterArray{1}.model = "SIREDmod";
  parameterArray{1}.folderName = "Test_Plots";
  parameterArray{1}.saveVTK = false;
  parameterArray{1}.saveDiagrams = false;
  parameterArray{1}.showDiagrams = false;
  parameterArray{1}.fullConsoleOut = false;
  parameterArray{1}.blowUp = false;
  parameterArray{1}.spatial = "germany";
  parameterArray{1}.initial = "GitData";
  parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
  parameterArray{1}.startDate = datenum([2020, 03, 02]);
  parameterArray{1}.endDateOpt = datenum([2020, 04, 25]);
  parameterArray{1}.betaStateWise= false;
  parameterArray{1}.exitRestrictions= 1;  
  parameterArray{1}.gamma1 = 0.12;
  parameterArray{1}.gamma2 = 0.065;
  parameterArray{1}.mortality = 0.0083;
  parameterArray{1}.darkFigure = 8.5; 
  parameterArray{1}.totalRuntime = 54;  
  parameterArray{1}.beta_cross_county = 1;
  
  parameterArray{1}.beta = 0.163;
  parameterArray{1}.majorEvents   = 0.67;
  parameterArray{1}.schoolClosing = 0.415;
  parameterArray{1}.contactRestrictions = 0.415;
  
  % Change here if states or counties
  parameterArray{1}.reduceToStates= true;   
  workspaceStatewise = sir_spatial(parameterArray); 
  
  parameterArray{1}.reduceToStates= false;
  parameterArray{1}.beta_cross_county = 1.08;
  
  workspaceCountywise = sir_spatial(parameterArray);
  t = min(workspaceStatewise.x.x):max(workspaceStatewise.x.x);
  
  
  xSW = spline(workspaceStatewise.x.x, workspaceStatewise.x.y, t);
  distancesSW = workspaceStatewise.distances;
  spatialDistinctionSW = workspaceStatewise.spatialDistinction;
  populationSW = workspaceStatewise.population;
  biggestCityInCountySW = workspaceStatewise.biggestCityInCounty;
  
  
  xCW = spline(workspaceCountywise.x.x, workspaceCountywise.x.y, t);
  distancesCW = workspaceCountywise.distances;
  spatialDistinctionCW = workspaceCountywise.spatialDistinction;
  populationCW = workspaceCountywise.population;
  biggestCityInCountyCW = workspaceCountywise.biggestCityInCounty;
  
  transformationMatrix = zeros(length(spatialDistinctionSW), length(spatialDistinctionCW));
  for i=1:size(transformationMatrix, 1)
    for j=1:size(transformationMatrix, 2)
      if floor(spatialDistinctionCW{j}.ags/1000) == spatialDistinctionSW{i}.ags/1000
        transformationMatrix(i, j) = spatialDistinctionCW{j}.population /...
        spatialDistinctionSW{i}.population;
      end
    end
  end
  
  correctionFactorsBeta = zeros(length(spatialDistinctionSW), length(t));
  for i=1:size(correctionFactorsBeta, 2)
    interactionMatrixSW = sir_eqn_spatial("InteractionMatrix", xSW(:,i), t(i), distancesSW,...
    populationSW, parameterArray{1}, biggestCityInCountySW, spatialDistinctionSW);
    interactionMatrixCW = sir_eqn_spatial("InteractionMatrix", xCW(:,i), t(i), distancesCW,...
    populationCW, parameterArray{1}, biggestCityInCountyCW, spatialDistinctionCW);
    
    ESW = xSW(3:4:end, i);
    ECW = xCW(3:4:end, i);
    
    % Mit Diagonale
    dISW_dt = interactionMatrixSW * ESW;
    dICW_dt = interactionMatrixCW * ECW;    
    dICS_s = transformationMatrix * dICW_dt;
    correctionFactorsBeta(:, i) = dISW_dt./dICS_s;
    
    % Ohne Diagonale
    dISW_dt2 = (interactionMatrixSW - diag(diag(interactionMatrixSW))) * ESW;
    dICW_dt2 = (interactionMatrixCW - diag(diag(interactionMatrixCW))) * ECW;
    dICS_s2 = transformationMatrix * dICW_dt2;
    correctionFactorsBetaCC(:, i) = dISW_dt2./dICS_s2;
  end
  
  slopeBeginning = (mean(correctionFactorsBeta(:,1:5)') - 1).^3 + 1;
  
endfunction



















