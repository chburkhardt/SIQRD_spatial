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

function out = writeTexTable(stateWiseSIR, parameter, filename)  
  switch parameter.model
    case "SIR"
      titles = {"Susceptible", "Infected", "Removed", "Cumsum Infected"};
    case {"SIRED","SIREDmod"}
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
  
  sumGermany.SIR_vs_time = zeros(size(stateWiseSIR{1}.SIR_vs_time));
  sumGermany.name = "Sum Germany";
  for i=1:length(stateWiseSIR)
    sumGermany.SIR_vs_time += stateWiseSIR{i}.SIR_vs_time;
  end
  stateWiseSIR{end+1} = sumGermany;
  
  nVals = size(stateWiseSIR{1}.SIR_vs_time, 2);
  out = zeros(length(stateWiseSIR{1}.time), length(stateWiseSIR) * nVals + 1);
  out(:,1) = stateWiseSIR{1}.time;
  for i=1:length(stateWiseSIR)
    out(:, (1 + nVals*(i-1)+1):(1 + nVals*i)) = stateWiseSIR{i}.SIR_vs_time;
  endfor
  
  fid = fopen(filename, 'w'); 
  fprintf(fid, "%s ", "time");
  for i=1:length(stateWiseSIR)
    for j=1:length(titles)
      fprintf(fid, "%s ", [strrep(titles{j}," ","_"), "_",...
      strrep(stateWiseSIR{i}.name," ","_")]);
    endfor
  endfor
  fprintf(fid, "\n");
  for i=1:size(out, 1)
    fprintf(fid, "%f ", out(i,:));
    fprintf(fid, "\n");
  endfor  
  fclose(fid);
endfunction


