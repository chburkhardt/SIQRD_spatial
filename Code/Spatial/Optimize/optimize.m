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
## @deftypefn {} {@var{retval} =} optimize (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-24

function retval = optimize () %filelname 1 und filename 2
  addpath('../');
  %% Example  
  
  %%% initial guess for beta and beta_cross_county
  x0Cov = [0.29586 1e-2];
  x0CovNames = {'beta', 'beta_cross_county'};
  
  parameterArray = read_parameter('../../Results/parameter.txt');
  
  %% the entries can be directly modified in the struct
  %% eg. set a foldername
  parameterArray{1}.model = "SIRED";
  parameterArray{1}.totalRuntime = 40;
  parameterArray{1}.folderName = "Optimize_germany2";
  parameterArray{1}.saveVTK = false;
  parameterArray{1}.saveDiagrams = false;
  parameterArray{1}.showDiagrams = false;
  parameterArray{1}.fullConsoleOut = false;
  
  
  
  %days können wir später dann in das Parameterarry mit aufnehmen
  %days beschreibt die Zeitspanne in der die Optimierung laufen soll
  days = 14;
  
  RKIread = read_case_history("../../Daten/Infectionnumbers.txt", "../../Daten/Deathnumbers.txt");
  
  %rewrite the cell, since states are arranged differently in RKIread and statewiseSIR
  RKIdata = cell(16,1);
  RKIdata{1} = RKIread{15};
  RKIdata{2} = RKIread{6};
  RKIdata{3} = RKIread{9};
  RKIdata{4} = RKIread{4};
  RKIdata{5} = RKIread{10};
  RKIdata{6} = RKIread{7};
  RKIdata{7} = RKIread{11};
  RKIdata{8} = RKIread{1};
  RKIdata{9} = RKIread{2};
  RKIdata{10} = RKIread{12};
  RKIdata{11} = RKIread{3};
  RKIdata{12} = RKIread{4};
  RKIdata{13} = RKIread{8};
  RKIdata{14} = RKIread{13};
  RKIdata{15} = RKIread{14};
  RKIdata{16} = RKIread{16};
  
  %todo: hier brauchen wir dann einen Parameter, der im Parametefile 
  %übergeben wird und die Tage angibt, die wir fitten wollen
  rkihelp = zeros(length(RKIdata)*days,1);
  for i=1:length(RKIdata)
    for j=1:days
      rkihelp((i-1)*days+j) = cell2mat(RKIdata{i,1}.infected(j));
    end
  end
  
  %rki daten kumulieren
  rkidata = zeros(size(rkihelp),1);
  for i=1:size(rkihelp,1)
      if mod(i-1, days) == 0
        sum = rkihelp(i,1);
      elseif mod((i-1),days) == 0 | i == 1
        sum = rkihelp(i,1);
      else
        sum = rkidata(i-1,1) + rkihelp(i,1);
      end
      rkidata(i,1) = sum;
  end

  
  
  runoptCovid = @(xCov) optCovid(xCov, rkidata, parameterArray, days); 
  pkg load optim;
  [paramCovid,resNormCovid] = lsqnonlin(runoptCovid,x0Cov,[0,0],[1,1]);
  resNormCovid
  
  %write fitted parameters in .txt file
  folder = ["../../Results/", parameterArray{1}.folderName];
  fid = fopen([folder, "/fitted_parametrs.txt"], 'w');
  for i=1:length(x0CovNames)
    fprintf(fid, '%s\t %f\t %10f\n', x0CovNames{i}, paramCovid(i,1), resNormCovid);
  end
  fclose(fid);
  paramCovid
end 
function opt = optCovid(xCov, rkidata, parameterArray, days)
  
  %% this is how a parameter struct can be loaded from a parameterfile
  % alle Pfade sind so gesetzt das man im "Spatial" order starten muss,
  % daher auch hier

  parameterArray{1}.beta = xCov(1);
%  parameterArray{1}.gamma = xCov(2);  
  parameterArray{1}.beta_cross_county	= xCov(2);
  disp(xCov);
  
  %% run code and store the workspace to avoid io
  workspaceCalculation = sir_spatial(parameterArray);
  
  %% load results statewise with ne new option to avoid io
  stateWiseSIR = extractStatewiseResults (workspaceCalculation, "workspace",...
  ["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
  
  t = stateWiseSIR{1}.time;
  
##  tUniform = linspace(min(t), max(t), n_pics);
  
  %% stateWiseSIR contains SIR, names, times
  %% see example
##  stateNames = cell(length(stateWiseSIR), 1);
##  figure;
##  hold on;
##  for i=1:length(stateWiseSIR)
##    plot(t, stateWiseSIR{i}.SIR_vs_time(:,2));%, RKIdata{i}.infected);
##    stateNames{i} = stateWiseSIR{i}.name;
##  end
##
##  legend(stateNames, "location", "northeast");

  tUniform = linspace(min(t), max(t), parameterArray{1}.totalRuntime + 1);
  xUniform = zeros(length(tUniform), 16);
  for i = 1:size(xUniform,2)
    xUniform(:, i) = spline(t, stateWiseSIR{i}.SIR_vs_time(:,6), tUniform);
  end
  
  %%%%%%%%
  % write an array with data from 2020-03-02 until 2020-03-13
  % todo: siehe rkidata
  
  simdata = zeros(length(stateWiseSIR)*days,1);

  for i=1:size(xUniform,2);
    for j=28:size(xUniform,1);
      simdata((i-1)*days+(j-27)) = xUniform(j,i);
      %simdata((i-1)*days+(j-(size(xUniform)-days))) = xUniform(j,i);
    end
  end
  
%%optFunktion
  opt = rkidata - simdata;

%  legend(stateNames, "location", "northeastoutside");
    
end