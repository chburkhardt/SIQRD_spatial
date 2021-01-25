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
## @deftypefn {} {@var{retval} =} checkData (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-06-17

function retval = checkData (input1, input2)
  close all;
  
  oldData = read_case_history();
  newDataMeldung = read_case_history ("Meldung");
  newDataErkrankung = read_case_history ("Erkrankungsbeginn");
  
  figure;
  for i=1:length(oldData)
    subplot(4,4,i), plot(cell2mat(oldData{i}.time),...
    cumsum(cell2mat(oldData{i}.infected)));
    hold on;
    subplot(4,4,i), plot(cell2mat(newDataMeldung{i}.time),...
    cumsum(cell2mat(newDataMeldung{i}.infected)));
    subplot(4,4,i), plot(cell2mat(newDataErkrankung{i}.time),...
    cumsum(cell2mat(newDataErkrankung{i}.infected)));
    title(oldData{i}.name);
  end
  
endfunction
