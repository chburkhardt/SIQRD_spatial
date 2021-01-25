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
## @deftypefn {} {@var{retval} =} sir_eqn (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-03


%% x = [S, I, R]
%% t is not used
function xdot = sir_eqn(mode, x, t, parameter)
  % Parameter values
##  if t<=19
##    beta=0.3;
##    gamma=0.05;
##  else
##    beta = 0.15;
##    gamma=0.05;
##  end
##  
  factorsBeta = sir_eqn_spatial ("getFactorsBeta", parameter, t);
  betaHelp = factorsBeta(1:length(factorsBeta)-1);
  betaHelp(17) = 1;
  population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
  4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
  1609675, 4077937, 2208321, 2143145];
  
  beta = parameter.beta;
  gamma = parameter.gamma1;
  
  if parameter.betaStatewise
    beta *= betaHelp(parameter.AGS_state);
  else
    beta *=  sum(factorsBeta(1:end-1).*population')/sum(population);
  end
  
  % Define variables
  s = x(1);
  y = x(2);
  r = x(3);
 if strcmp(mode, "rates")
    
    % Define ODEs
    ds=-beta*s*y;
    dy=beta*s*y-gamma*y;
    dr=gamma*y;
    
    % Return gradients
    xdot = [ds,dy,dr];
  elseif strcmp(mode, "jacobi")
    % Return jacobi
    xdot = [-beta*i, -beta*s, 0; beta*i, beta*s-gamma, 0 ; 0, gamma, 0 ]
   
  elseif 
      error("Unknown mode in sir_eqn");
  endif
endfunction
