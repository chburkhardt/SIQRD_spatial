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
## @deftypefn {} {@var{retval} =} naiveGradientSearch (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-05-29

function retval = naiveGradientSearch (input1, input2)
	close all;
	% Parameter and rkiData
	parameterArray = getParameterArray();	
	RKIdata = getRKIdata();
	
	% set opt part in parameters
  parameterArray{1}.paraNamesOpt = {"beta_cross_county"};
  parameterArray{1}.x0CovOpt = [1]; 
  variations = ones(1, length(parameterArray{1}.paraNamesOpt)) * 0.2;
  parameterArray{1}.lbOpt = parameterArray{1}.x0CovOpt .* (1 - variations);
  parameterArray{1}.ubOpt = parameterArray{1}.x0CovOpt .* (1 + variations);	
  parameterArray{1}.optAlgorithmOpt = "lsqnonlin"; %lsqnonlin|pso|justVisualize
	
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% calculate reduced to states result	
	parameterArray{1}.reduceToStates= true;
	if !exist("stateWiseSIRSW", "var")
		workspaceCalculationSW = sir_spatial(parameterArray);
		persistent stateWiseSIRSW = extractStatewiseResults (workspaceCalculationSW, "workspace",...
		["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
		for i=1:length(stateWiseSIRSW)
			stateWiseSIRSW{i}.splineDiscovered = spline(stateWiseSIRSW{i}.time +...
			parameterArray{1}.startDate, stateWiseSIRSW{i}.SIR_vs_time(:,6)); 
		end
	end
	timespan = [max([min(cell2mat(RKIdata{1}.time)),...
	min(stateWiseSIRSW{1}.time + parameterArray{1}.startDate)]),...
	min([max(cell2mat(RKIdata{1}.time)),...
	max(stateWiseSIRSW{1}.time + parameterArray{1}.startDate)])];
	time = timespan(2);  
	simDataSW = zeros(length(stateWiseSIRSW) ,1);
	rkiData = zeros(length(stateWiseSIRSW) ,1);
	for i=1:length(stateWiseSIRSW)
		simDataSW(i) = ppval(stateWiseSIRSW{i}.splineDiscovered, time);		
		rkiData(i) = ppval(RKIdata{i}.splineInfected, time);		
	end
	
	%% run code and store the workspace to avoid io
	parameterArray{1}.reduceToStates= false;
  parameterArray{1}.optAlgorithmOpt = "justVisualize"; %lsqnonlin|pso|justVisualize
	parameterArray{1}.showStatePlots = false;
	parameterArray{1}.beta_cross_county	= 0.525; 	
	
	counter = 0;
	optVals = ones(1, length(rkiData)).*parameterArray{1}.beta_cross_county';
	residuum = [];
	while counter < 100
		t1 = tic;
		workspaceCalculation = sir_spatial(parameterArray);
		
		%% load results statewise with ne new option to avoid io
		stateWiseSIR = extractStatewiseResults (workspaceCalculation, "workspace",...
		["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
		for i=1:length(stateWiseSIR)
			stateWiseSIR{i}.splineDiscovered = spline(stateWiseSIR{i}.time +...
			parameterArray{1}.startDate, stateWiseSIR{i}.SIR_vs_time(:,6)); 
		end
		simData = zeros(length(stateWiseSIR) ,1);
		for i=1:length(stateWiseSIR)
			simData(i) = ppval(stateWiseSIR{i}.splineDiscovered, time);		
		end

		ratio = (simDataSW - simData)./simDataSW;
		maxScalePerStep = 0.25;
		exponent = 1.5;		
		scaleFactors = ratio;
		scaleFactors(abs(scaleFactors) > maxScalePerStep^(1/exponent)) =...
		sign(scaleFactors(abs(scaleFactors) > maxScalePerStep^(1/exponent)))*...
		maxScalePerStep^(1/exponent);
		
		residuum(end+1,:) = ratio;
		
		parameterArray{1}.beta_cross_county .*=...
		(1 + sign(scaleFactors).*abs(scaleFactors).^exponent); 
		optVals(end+1, :) = parameterArray{1}.beta_cross_county;
		counter += 1;
		if counter > 1
			figure(20);
			clf;
			plot(optVals);
			pause(0.2);
		end
		fprintf("Iteration: %i took %fs with res: %f Parameters: ", counter, toc(t1), norm(residuum(end)));
		fprintf("%f ", parameterArray{1}.beta_cross_county);
		fprintf("\n");
		if norm(residuum(end)) < 1e-4
			fprintf("Terminate due to small enough residuum");
			break;
		end	
	end
	
	optimizeStateWise (parameterArray);
	
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RKIdata = getRKIdata()
	RKIread = read_case_history("../../Daten/Infectionnumbers.txt", "../../Daten/Deathnumbers.txt");
  
  %rewrite the cell, since states are arranged differently in RKIread and statewiseSIR
  RKIdata = cell(16,1);
  RKIdata{1} = RKIread{15};
  RKIdata{2} = RKIread{6};
  RKIdata{3} = RKIread{9};
  RKIdata{4} = RKIread{5};
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
  
  %cumulate infection numbers 
  %create splines
  timesRKI = cell2mat(RKIdata{1}.time);
  for i=1:length(RKIdata)
    RKIdata{i}.infected = cumsum(cell2mat(RKIdata{i}.infected));
    RKIdata{i}.splineInfected = spline(timesRKI, RKIdata{i}.infected);
    RKIdata{i}.dead = cell2mat(RKIdata{i}.dead);
    RKIdata{i}.splineDead = spline(timesRKI, RKIdata{i}.dead);
  end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameterArray = getParameterArray()
	
	parameterArray = read_parameter('../../Results/parameter.txt');  
	%% the entries can be directly modified in the struct
	%% eg. set a foldername
  parameterArray{1}.model = "SIREDmod";
  parameterArray{1}.folderName = "Test_Plots";
  parameterArray{1}.saveVTK = false;
  parameterArray{1}.saveDiagrams = false;
  parameterArray{1}.showDiagrams = false;
  parameterArray{1}.fullConsoleOut = true;
  parameterArray{1}.blowUp = false;
  parameterArray{1}.spatial	= "germany";
  parameterArray{1}.initial = "GitData";
  parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
  parameterArray{1}.startDate = datenum([2020, 03, 02]);
  parameterArray{1}.endDateOpt = datenum([2020, 04, 25]);
  parameterArray{1}.reduceToStates= true;
  parameterArray{1}.betaStateWise= false;
  parameterArray{1}.exitRestrictions= 1;
  parameterArray{1}.schoolClosing= 1;
  parameterArray{1}.beta_cross_county	= 1;
  
  parameterArray{1}.wRatioGlobStates = 0; % 1 is fully statewise, 0 is only global
  parameterArray{1}.wRatioID = 1; % 1 is only infected, 0 only death
  parameterArray{1}.wISlope = 0.002; % weight for slope of infected last date
  parameterArray{1}.optFunGermanSum = true;
  parameterArray{1}.showStatePlots = true;
  
  parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
  parameterArray{1}.gamma2 = 0.04; % Quarantine -> Recovered
  parameterArray{1}.mortality = 0.006; % 0.006
  parameterArray{1}.darkFigure = 8.6403;
  parameterArray{1}.beta =  0.175564;
  
  parameterArray{1}.majorEvents		= 0.378322;
  parameterArray{1}.contactRestrictions = 0.188669;
  
	parameterArray{1}.fullConsoleOut = false;
endfunction
