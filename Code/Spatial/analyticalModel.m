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
## @deftypefn {} {@var{retval} =} analyticalModel (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-05-25

function retval = analyticalModel (spatialDistinction, parameter, distances, X0, timespan)
  % addpath pso
  
  population = zeros(length(spatialDistinction), 1);
  biggestCity = zeros(length(spatialDistinction), 1);
  coords = zeros(2, length(spatialDistinction));
  for i=1:length(spatialDistinction)
    population(i) = spatialDistinction{i}.population;
    biggestCity(i) = spatialDistinction{i}.biggestCityInCounty;
    coords(:,i) = [spatialDistinction{i}.x, spatialDistinction{i}.y];
  end
  maxBiggestCity = 3e6;
  interactionMatrix = sir_eqn_spatial("InteractionMatrix", nan, 0, distances,...
  population, parameter, maxBiggestCity, spatialDistinction);
  
  x0State = X0;
  I0 = x0State(3:4:end); %Entries S, Q, I, D
  [V,D] = eig(interactionMatrix);
  %Compute alpha
  alpha = parameter.gamma1 / (parameter.darkFigure - 1);
  t = linspace(min(timespan), max(timespan), 100);
  Q_RKI = zeros(length(I0),length(t));
  for i=1:length(t)
    Q_RKI(:,i) = alpha * V * inv(D) * diag(exp(t(i)*diag(D))) * inv(V) * I0;
  end
  plot(Q_RKI');
  
  %% opt
##  opt.parallel = true;
##      opt.visu = true;
##      opt.n_particles = 300;
##      opt.n_iter = 100;
##      opt.coupling = 5;
##      opt.parameter_names = paraNames;
##      
##      parameterArray{1}.fullConsoleOut = false;
##      evalSirLocalOneArgument= @(x) somethingThatReturnsNormOfError(x);
##      [paramCovBest, resMin] = pso(evalSirLocalOneArgument, lb, ub, opt);
  
##  retval.x = times;
##  retval.y = solutionmatrix;
endfunction


function x = somethingThatReturnsNormOfError()

endfunction
