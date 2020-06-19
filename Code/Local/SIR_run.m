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
## @deftypefn {} {@var{retval} =} SIR_test (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-01

function retval = SIR_run (parameter, X0)
  addpath("../Spatial")
  if nargin == 0
    close all;
    parameter = read_parameter("../../Results/parameter_local.txt"){1}; 
    
    N = 8e7;
    I0 = 274; %Q Gruppe
    E0 = 31119; %I Gruppe
    base_Q(1) = (938/274)^(1/2);
    base_Q(2) = (1976/274)^(1/4);
    base_Q(3) = (3717/274)^(1/6);
    E0 = I0 * log(base_Q(2)) * ((parameter.darkFigure - 1) / parameter.gamma1) * exp(1);
    D0 = 0;
    R0 = 0;
    S0 = N - I0 - E0 - D0 -R0;  
    X0=[S0, I0, R0, E0, D0] / N;
  end  
  % t in days
  t = [0, parameter.totalRuntime]; 
  get_dXdT = @ (t, x, hist) [siredMod(t, x, parameter)];   
  vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
  x = ode23(get_dXdT, t, X0, vopt);
  
  if parameter.showDiagrams
    figure();
    colorvector = ["-r";"-g";"-b";"-c";"-m";"-k"];
    legendvector = {"S","Q","R","I","D",};%"CHECK"};
    for i = 1:size(x.y,1)
      plot(x.x(:),x.y(i,:)*N,colorvector(i,:));
      hold on;
    endfor
    
    hold off;
    xlabel("Time","fontweight","bold");
    ylabel("Number","fontweight","bold");
    h = legend(legendvector{[,1:size(x.y,1)]});
    legend(h,"show");
  end
  if nargout > 0
    retval = x;
  end
endfunction


