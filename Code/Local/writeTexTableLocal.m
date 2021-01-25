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
## @deftypefn {} {@var{retval} =} writeTexTable (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-26

function out = writeTexTableLocal(stateWiseSIR, parameter, filename)  
  switch parameter.model
    case "SIR"
      titles = {"Susceptible", "Infected", "Removed", "Cumsum Infected"};
    case {"SIRED","SIREDmod", "SIREDLiterature"}
      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
      "Cumsum Infected"};
    case "SIREDYO"
      titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
      "Cumsum Infected", "Susceptible risk", "Infected risk", "Recovered risk",...
      "Exposed risk", "Dead risk", "Cumsum Infected risk"};
    case "SIRH"
      %SIRH = [S, X, X, effectiveInfective, symptoms, hospital, icu, dead, discovered];
      titles = {"Susceptible", "Exposed", "Removed", "Infectious",...
      "has symptoms", "in hospital", "needs intensive care", "Death", "Discovered"};
  endswitch
  
  nVals = size(stateWiseSIR.y, 1);
  out = zeros(length(stateWiseSIR.x), nVals + 1);
  out(:,1) = stateWiseSIR.x;
  out(:, (1 + 1):(1 + nVals)) = stateWiseSIR.y'*parameter.population;

  
  fid = fopen(filename, 'w'); 
  fprintf(fid, "%s ", "time");
  for j=1:length(titles)
      fprintf(fid, "%s ", [strrep(titles{j}," ", "_")]);
  endfor
  fprintf(fid, "\n");
  for i=1:size(out, 1)
    fprintf(fid, "%f ", out(i,:));
    fprintf(fid, "\n");
  endfor  
  fclose(fid);
endfunction


