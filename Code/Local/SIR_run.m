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

function retval = SIR_run (parameter, X0, history, deathStart, darkFigure, population)
  addpath("../Spatial")
  if nargin == 0
    close all;
    parameter = read_parameter("../../Results/parameter_local.txt"){1}; 
    N = 8e7;
    switch parameter.model
      case "SIR"
        I0 = 1000;      
        R0 = 0; 
        S0 = N - I0 - R0;   
        X0=[S0, I0, R0] / N; 
      case "SIRH"
        I0 = 1000;      
        R0 = 0; 
        S0 = N - I0 - R0;   
        X0 = S0 / N; 
      case {"SIRED","SIREDmod"}
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
    endswitch
  end  
  
  % t in days
  t = [0, parameter.totalRuntime];
  
  %setup system
 
  % N S I R
  switch parameter.model
    case "SIR"
      get_dXdT = @(t_eq, x0) sir_eqn("rates", x0, t_eq, parameter);     
      vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
      x = ode23(get_dXdT,t, X0,vopt);
      x.y = [x.y; 1-x.y(1,:)];
    case "SIRED"
      vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
      x = ode23("sir_extendet_eqn", t, X0, vopt);
    case "SIREDmod"
      get_dXdT = @ (t, x, hist) [siredMod(t, x, parameter)];   
      vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
      x = ode23(get_dXdT, t, X0, vopt);
      darkFigure = parameter.darkFigure;  
      splineE = spline(x.x, x.y(4,:));
      alpha = parameter.gamma1 / (darkFigure - 1);
      get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
      %[t2, discovered] = ode23(get_dDiscovering_dT, [min(xSim.x), max(xSim.x)], xSim.y(2,1));
      [t2, discovered] = ode15s(get_dDiscovering_dT, [min(x.x),max(x.x)],...
      x.y(3,1)/darkFigure+x.y(2,1));
      solutionDiscovered = spline(t2, discovered, x.x);
      x.y = [x.y; solutionDiscovered];      
  case "SIREDLiterature"
      get_dXdT = @ (t, x, hist) [siredLiterature(t, x, parameter)];   
      vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
      x = ode23(get_dXdT, t, X0, vopt);
      darkFigure = parameter.darkFigure;  
      splineE = spline(x.x, x.y(4,:));
      alpha = 1/3;
      
      get_dDiscovering_dTL = @(t_eq, x0) [alpha*ppval(splineE, t_eq)];
      [t2, discovered] = ode15s(get_dDiscovering_dTL, [[min(x.x), max(x.x)]], x.y(2,1));

      solutionDiscovered = spline(t2, discovered, x.x);        
      x.y = [x.y; solutionDiscovered];
      
    case "SIRH"    
##      sizeHistory = 29;
      X0 = X0(1);
##      if nargin <= 2
##        history = [linspace(1, 1, floor(sizeHistory*0.60)),...
##        linspace(1, X0, floor(sizeHistory*0.40))];
##      end
      [history, X0] = getHistory(X0, parameter, 29);
      vopt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
      get_dXdT = @ (t, x, hist) [sir_history(t, x, hist, parameter)];    
      lags = sir_eqn_spatial("lags");  
      tic;
      x = ode23d_fixed(get_dXdT, t, X0, lags, history, vopt);
      
      % postprocess for discovered people
      S = x.y;
      S_history = history;
      courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);
      PPS = spline(...
      [linspace(-max(lags), min(x.x), size(S_history, 2)), x.x(2:end)],...
      [S_history, S(2:end)]);
  
    opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
    if nargin <= 2
      get_dDiscovering_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags))*courseOfDisease.discovering']/parameter.darkFigure;
    else
      get_dDiscovering_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags),1,2)*courseOfDisease.discovering']/parameter.darkFigure;
    end
    [t2, discovered] = ode15s(get_dDiscovering_dT, [min(x.x), max(x.x)], ...
    (1-S_history(end))./parameter.darkFigure, vopt);
    discovered = spline(t2, discovered, x.x);%*population;
      
      get_dDying_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags))*courseOfDisease.dying'];
      if nargin <= 2
        [t3, dead] = ode15s(get_dDying_dT, [min(x.x), max(x.x)], 0, vopt);
        dead = spline(t3, dead, x.x);%*population;
      else 
        [t3, dead] = ode15s(get_dDying_dT, [min(x.x), max(x.x)], deathStart, vopt);
        dead = spline(t3, dead, x.x);%*population;
      end
        x.y = [x.y; discovered; dead];
        x.y = [x.y; discovered; dead];
      
      fprintf("ode23d took %is\n", toc);
  endswitch
  
  if isfield(parameter, 'writeTexTableLocal')
    if parameter.writeTexTableLocal
    writeTexTableLocal(x, parameter, ["../../Results/", parameter.folderName, "/", "dataTable.dat"]);
    
    names = fieldnames(parameter);
    fid = fopen(["../../Results/", parameter.folderName, "/parameters_out.txt"], 'w'); 
    for i=1:length(names)
    if ischar(parameter.(names(i){1}))
      fprintf(fid, "set %30s = %s\n", names(i){1}, parameter.(names(i){1}));
    elseif iscell(parameter.(names(i){1}))
      %nothing to do
    else  
      fprintf(fid, "set %30s = %s\n", names(i){1}, mat2str(parameter.(names(i){1})));
    end
  end
  fclose(fid);
  end
  end
  
  if parameter.showDiagrams
    figure();
    colorvector = ["-r";"-g";"-b";"-c";"-m";"-k"];
    legendvector = {"S","Q","R","I","D",};%"CHECK"};
    for i = 1:size(x.y,1)
      plot(x.x(:),x.y(i,:)*N,colorvector(i,:));
      hold on;
    endfor
    %plot(x.x(:),x.y(1,:) + x.y(2,:) + x.y(3,:) + x.y(4,:) + x.y(5,:),colorvector(6,:));
    
    %semilogy(t,x(:,1),"-r",t,x(:,2),"-g",t,x(:,3),"-b",...
    %[min(t), max(t)], [1,1]*n_beds/ratio_needs_bed/dark_ratio/N, "c",...
    %t, x(:,2) + x(:,3), "m",...
    %data(:,1), data(:,2)/N, "k");
    
    hold off;
    %  xlim([min(t), max(t)]);
    xlabel("Time","fontweight","bold");
    ylabel("Number","fontweight","bold");
    %h = legend("S","I","R", "Beds", "I+R", "ger");
    h = legend(legendvector{[,1:size(x.y,1)]});
    legend(h,"show");
  end
  if nargout > 0
    retval = x;
  end
endfunction

function [history, X_0] = getHistory(X_0, parameter, size)
lags = sir_eqn_spatial("lags");
courses = sir_eqn_spatial("totalCourseOfDisease", parameter);
nu = 0.345;
darkFigure = parameter.darkFigure;
beta = parameter.beta;

lags2 = (lags(1:end-1)+lags(2:end))/2;

discoveredFun = @(tau) exp(-nu*tau).*interp1(lags2, courses.discovering, tau, "linear");
factor = 1./(integral(discoveredFun, 0.5,27.5) ./ darkFigure * (1-exp(-nu*28)));

% x_0 is equal to the number of discovered infected
history = 1-((1-X_0).*factor).*exp(-nu*linspace(max(lags), min(lags), size));
X_0 = history(:,end);
endfunction
