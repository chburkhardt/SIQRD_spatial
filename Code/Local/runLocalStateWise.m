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
## Created: 2020-05-04

function retval = runLocalStateWise (input1, input2)
  addpath("../Spatial", "../Spatial/Optimize", "../Pso");
  tic; 
  %%%%%%%% Parameter %%%%%%%%
  parameter = read_parameter("../../Results/parameter_local.txt"){1};   
  parameter.showDiagrams = false;
  parameter.model = "SIREDmod"; %"SIREDmod" 
  parameter.folderName = "runStateWiseE0";
  parameter.startDate = datenum([2020, 03, 02]);
  betaStateWise	= true;
  
  %%%%%%%% parameters from lsqnonlin25runs (14.05.2020) %%%%%%%% 
  %%%%%%%% median from good residuals at similat start dates %%%%%%%%
  %parameter.gamma1 = 0.095196; 
  %parameter.gamma2 = 0.041280;
  %parameter.mortality = 0.0068125;
  %parameter.darkFigure = 7.6285;
  %parameter.startDate = 737839.81073;
  
  %%%%%%%% parameters from literature %%%%%%%%
  parameter.mortality = 0.006;
  parameter.gamma1 = 0.067; 
  parameter.gamma2 = 0.04;
  mkdir(["../../Results/", parameter.folderName]);
  
  %%%%%%%% read state parameters %%%%%%%
  stateParameter = csvread("../../Daten/PSO_Parameters.csv");
  betaStateWise = stateParameter(2:17,3);
  darkFigureStateWise = stateParameter(2:17,4);
	majorEventsStateWise = stateParameter(2:17,5);
	schoolClosingStateWise = stateParameter(2:17,6);
	contactRestrictionStateWise = stateParameter(2:17,7);
  
  betaStateWise(17) = 0.315259;
  darkFigureStateWise(17) = 9.192634;
  majorEventsStateWise(17) = 0.641772;
  schoolClosingStateWise(17) = 0.493906;
  contactRestrictionStateWise(17) = 0.424934;
  
  paraNames = {"beta", "darkFigure", "majorEvents", "schoolClosing", "contactRestrictions"};
  
  %%%%%%%% RKI data %%%%%%%%
  RKIread = read_case_history("../../Daten/Infectionnumbers.txt",...
  "../../Daten/Deathnumbers.txt");
  
  %%%%%%%% rearrange RKI data so that they are sorted by AGS %%%%%%%%
  %%%%%%%% source for population: Wikipedia, 2018 %%%%%%%%
  RKIdata = cell(17,1);
  RKIdata{1} = RKIread{15};
  RKIdata{1}.population = 2896712;
  RKIdata{2} = RKIread{6};
  RKIdata{2}.population = 1841179;
  RKIdata{3} = RKIread{9};
  RKIdata{3}.population = 7982448;
  RKIdata{4} = RKIread{5};
  RKIdata{4}.population = 682986; 
  RKIdata{5} = RKIread{10};
  RKIdata{5}.population = 17932651;
  RKIdata{6} = RKIread{7};
  RKIdata{6}.population = 6265809;
  RKIdata{7} = RKIread{11};
  RKIdata{7}.population = 4084844;
  RKIdata{8} = RKIread{1};
  RKIdata{8}.population = 11069533;
  RKIdata{9} = RKIread{2};
  RKIdata{9}.population = 13076721; 
  RKIdata{10} = RKIread{12};
  RKIdata{10}.population = 990509;
  RKIdata{11} = RKIread{3};
  RKIdata{11}.population = 3644826;
  RKIdata{12} = RKIread{4};
  RKIdata{12}.population = 2511917;
  RKIdata{13} = RKIread{8};
  RKIdata{13}.population = 1609675;
  RKIdata{14} = RKIread{13};
  RKIdata{14}.population = 4077937;
  RKIdata{15} = RKIread{14};
  RKIdata{15}.population = 2208321; 
  RKIdata{16} = RKIread{16};
  RKIdata{16}.population = 2143145;
  
  %%%%%%%% germany %%%%%%%%
  timesRKI = cell2mat(RKIread{1}.time);
  RKIinfGermany = zeros(size(timesRKI));
  RKIdeadGermany = zeros(size(timesRKI));
  % cumulate RKI data from states and build splines
  for i=1:length(RKIread)
    RKIinfGermany += cell2mat(RKIread{i}.infected);
    RKIdeadGermany += cell2mat(RKIread{i}.dead);
  end
  RKIdata{17}.infected = RKIinfGermany;
  RKIdata{17}.dead = RKIdeadGermany;
  RKIdata{17}.time = mat2cell(timesRKI, 1);
  RKIdata{17}.population = 8e7;
  RKIdata{17}.name = 'Germany';
  RKIdata{17}.infected = cumsum(RKIdata{17}.infected);
  RKIdata{17}.splineInfected = spline(timesRKI, RKIdata{17}.infected);
  RKIdata{17}.splineDead = spline(timesRKI, RKIdata{17}.dead);
  
  %%%%%%%% build splines for each state %%%%%%%%
  for i=1:length(RKIdata)-1
    RKIdata{i}.infected = cumsum(cell2mat(RKIdata{i}.infected));
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end 
  
  %%%%%%%% fit data for each state %%%%%%%%
  for k=1:17  
    parameter.AGS_state = k;  
    N = RKIdata{k}.population;
    if parameter.AGS_state == 17
      parameter.betaStateWise = false; 
    endif
    
    parameter = setfield(parameter, "beta", betaStateWise(k));
    parameter = setfield(parameter, "darkFigure", darkFigureStateWise(k));
    parameter = setfield(parameter, "majorEvents", majorEventsStateWise(k));
    parameter = setfield(parameter, "schoolClosing", schoolClosingStateWise(k));
    parameter = setfield(parameter, "contactRestrictions", contactRestrictionStateWise(k));
    
    shiftDays = 4; 
    sumIsimStartP1 = RKIdata{k}.infected(1);
    sumIsimStartP2 = RKIdata{k}.infected(1 + shiftDays);
    
    base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
    
    N = RKIdata{k}.population;
    I0 = ppval(RKIdata{k}.splineInfected, parameter.startDate);
    D0 = 0;
    E0 = I0*log(base_Q)*((parameter.darkFigure - 1) / parameter.gamma1)*exp(1);
    R0 = 0; 
    S0 = N - I0 - R0 - E0 - D0;
    X0sir = [S0, I0, R0, E0, D0] / N;
    
    % calculate model
    tic1 = tic;
    xSim = SIR_run(parameter, X0sir); 
    
    startDateSim = parameter.startDate;  
    
    % discovered = \int_0^t alpha * E * dT       
    splineE = spline(xSim.x, xSim.y(4,:));
    alpha = parameter.gamma1 / (parameter.darkFigure - 1);
    get_dDiscovering_dT = @(t_eq, x0) [alpha * ppval(splineE, t_eq)];
    [t2, discovered] = ode23(get_dDiscovering_dT, [min(xSim.x), max(xSim.x)], xSim.y(2,1));
    solutionDiscovered = spline(t2, discovered, xSim.x);
    solutionDead = xSim.y(5,:); % dead
    
    splineSimDiscovered = spline(xSim.x + startDateSim, solutionDiscovered);
    splineSimDead = spline(xSim.x + startDateSim, solutionDead);
    
    timespan = [max([min(cell2mat(RKIdata{k}.time)), min(xSim.x + startDateSim)]),...
    min([max(cell2mat(RKIdata{k}.time)), max(xSim.x + startDateSim)])];
    times = linspace(timespan(1), timespan(2), 100);
    
    
    % Differenzenvektor berechnen [Infected, Dead]
    % all for Infected, last point for dead
    simdataI = ppval(splineSimDiscovered, times)*RKIdata{k}.population;
    simdataD = ppval(splineSimDead, times(end))*RKIdata{k}.population;
    rkiI = ppval(RKIdata{k}.splineInfected, times);
    rkiD = ppval(RKIdata{k}.splineDead, times(end));
    
    
    stateNames= {"Schleswig-Holstein", "Freie und Hansestadt Hamburg",...
    "Niedersachsen", "Freie Hansestadt Bremen", "Nordrhein-Westfalen",...
    "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Freistaat Bayern",...
    "Saarland", "Berlin", "Brandenburg", "Mecklenburg-Vorpommern"...
    "Freistaat Sachsen", "Sachsen-Anhalt", "Freistaat Thueringen", "Germany"};
    
    titles = {"Susceptible", "Infected", "Recovered", "Exposed", "Dead",...
        "Cumsum Infected", "RKI Infected", "RKI Dead"};

    
    nVals = size(xSim.y, 1) + 3;
    persistent out = zeros(length(times), 1+(nVals)*length(stateNames));
    out(:,1) = times;
    
    splineS = spline(xSim.x + startDateSim, xSim.y(1,:));
    splineI2 = spline(xSim.x + startDateSim, xSim.y(2,:));
    splineR = spline(xSim.x + startDateSim, xSim.y(3,:));
    splineE = spline(xSim.x + startDateSim, xSim.y(4,:));
    simdataS = ppval(splineS, times)*RKIdata{k}.population;
    simdataI2 = ppval(splineI2, times)*RKIdata{k}.population;
    simdataR = ppval(splineR, times)*RKIdata{k}.population;
    simdataE = ppval(splineE, times)*RKIdata{k}.population;
    
    out(:,(k-1)*nVals+2) = simdataS;
    out(:,(k-1)*nVals+3) = simdataI2;
    out(:,(k-1)*nVals+4) = simdataR;
    out(:,(k-1)*nVals+5) = simdataE;
    out(:,(k-1)*nVals+6) = 0;
    out(size(out,1),(k-1)*nVals+6) = simdataD;
    out(:,(k-1)*nVals+7) = simdataI;
    out(:,(k-1)*nVals+8) = rkiI;
    out(:,(k-1)*nVals+9) = 0;
    out(size(out,1),(k-1)*nVals+9) = rkiD;
    
    if k==length(stateNames)
      folder = ["../../Results/", parameter.folderName];
      fid = fopen([folder, "/textable"], 'w'); 
      fprintf(fid, "%s ", "time");
      for i=1:k
        for j=1:length(titles)
          fprintf(fid, "%s ", [strrep(titles{j}," ","_"), "_",...
          strrep(stateNames{i}," ","_")]);
        endfor
      endfor
      fprintf(fid, "\n");
      for i=1:size(out, 1)
        fprintf(fid, "%f ", out(i,:));
        fprintf(fid, "\n");
      endfor  
      fclose(fid);
    end
    
    persistent saveE0 = zeros(length(stateNames), 1);
    saveE0(k, 1) = E0; 
    if k==length(stateNames)
      folder = ["../../Results/", parameter.folderName];
      fid = fopen([folder, "/textable_E0"], 'w'); 
      for i=1:size(saveE0, 1)
        fprintf(fid, "%f ", saveE0(i,1));
        fprintf(fid, "\n");
      endfor  
      fclose(fid);
    end
    
  end
  
  endfunction