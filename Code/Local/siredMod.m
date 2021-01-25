## Copyright (C) 2020 iwtm96_user
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
## @deftypefn {} {@var{retval} =} sir_extendet_eqn (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-04-07

function xdot = siredMod (t, x, parameter)
  % Parameter values
  darkFigure = parameter.darkFigure;
  mortality = parameter.mortality;  
  beta = parameter.beta;
  gamma1 = parameter.gamma1;
  gamma2 = parameter.gamma2;  
  alpha = gamma1 / (darkFigure - 1);
  delta = gamma2 * (darkFigure * mortality)/(1 - darkFigure * mortality);
  
  factorsBeta = sir_eqn_spatial ("getFactorsBeta", parameter, t);
  betaHelp = factorsBeta(1:length(factorsBeta)-1);
  betaHelp(17) = 1;
  population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
  4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
  1609675, 4077937, 2208321, 2143145];
  if parameter.betaStatewise
    beta *= betaHelp(parameter.AGS_state);
  else
    beta *=  sum(factorsBeta(1:end-1).*population')/sum(population);
  end
 
  % Define variables
  s = x(1);
  i = x(2);
  r = x(3);
  e = x(4);
  d = x(5);
  
  % Define ODEs
  ds=-beta*s*(e);
  de=+beta*s*(e) - alpha*e - gamma1*e;
  di=+alpha*e-gamma2*i-delta*i;
  dr=+gamma1*e + gamma2*i;
  dd=+delta*i;
  
  % Return gradients
  xdot = [ds,di,dr,de,dd];
endfunction
