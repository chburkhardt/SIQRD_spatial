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
## @deftypefn {} {@var{retval} =} sir_history (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-22

function xDot = sir_history (t, x, history, parameter)
  courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);
  
  % Parameter values
  beta=parameter.beta;
  
  factorsBeta = sir_eqn_spatial ("getFactorsBeta", parameter, t);
  population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
  4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
  1609675, 4077937, 2208321, 2143145];
  beta *=  sum(factorsBeta(1:end-1).*population')/sum(population);
  
  % Define variables
  S = x;
  S_hist = history;  
  
  oldDeltaS = -diff(S_hist, 1, 2);  
  dSdT=-beta.*S.*(oldDeltaS*courseOfDisease.infective'); 
  xDot = dSdT;
endfunction
