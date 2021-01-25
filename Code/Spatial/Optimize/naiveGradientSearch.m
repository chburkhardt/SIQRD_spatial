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
  if nargin == 0
	  parameterArray = getParameterArray();	
  else
    parameterArray = input1;
  end
	RKIdata = getRKIdata();
	
	% set opt part in parameters
  parameterArray{1}.paraNamesOpt = {"beta_cross_county"};
  parameterArray{1}.x0CovOpt = [1]; 
  variations = ones(1, length(parameterArray{1}.paraNamesOpt)) * 0.4;
  parameterArray{1}.lbOpt = parameterArray{1}.x0CovOpt .* (1 - variations);
  parameterArray{1}.ubOpt = parameterArray{1}.x0CovOpt .* (1 + variations);	
  parameterArray{1}.optAlgorithmOpt = "lsqnonlin"; %lsqnonlin|pso|justVisualize
	
  mkdir(["../../Results/", parameterArray{1}.folderName]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% calculate reduced to states result	
	parameterArray{1}.reduceToStates= true;
	if !exist("stateWiseSIRSW", "var")
		workspaceCalculationSW = sir_spatial(parameterArray);
		persistent stateWiseSIRSW = extractStatewiseResults (workspaceCalculationSW, "workspace",...
		["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
		for i=1:length(stateWiseSIRSW)
			switch parameterArray{1}.model
				case "SIREDmod"					
					stateWiseSIRSW{i}.splineDiscovered = spline(stateWiseSIRSW{i}.time +...
					parameterArray{1}.startDate, stateWiseSIRSW{i}.SIR_vs_time(:,6)); 
				case "SIRH"
					stateWiseSIRSW{i}.splineDiscovered = spline(stateWiseSIRSW{i}.time +...
					parameterArray{1}.startDate, stateWiseSIRSW{i}.SIR_vs_time(:,9));  
			endswitch
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
	##		parameterArray{1}.beta_cross_county = rand(16,1)+0.5; %TODO remove this
  parameterArray{1}.optAlgorithmOpt = "lsqnonlin"; %lsqnonlin|pso|justVisualize
	parameterArray{1}.showStatePlots = false;
	##	plotAllStates = true;
	##	##	% For the countywise model
	##	parameterArray{1}.reduceToStates= false;
	parameterArray{1}.beta_cross_county	= 1; % 1.08 ohne ss to cs anpassung
	parameterArray{1}.beta_cc_statewise = ones(16,1);
	##  parameterArray{1}.beta_cross_county	= 0.9; % 0.8 ist noch etwas zu hoch, probiere 0.73
	%optimizeStateWise (parameterArray);
	
	
	counter = 0;
	optVals = ones(1, length(rkiData)).*parameterArray{1}.beta_cross_county';
	residuum = [];
  folder = ["../../Results/", parameterArray{1}.folderName];
  
	while counter < 100
    fid = fopen([folder, "/protokoll_naiveGradientSearch.txt"], "a");
		t1 = tic;
		workspaceCalculation = sir_spatial(parameterArray);
		
		%% load results statewise with ne new option to avoid io
		stateWiseSIR = extractStatewiseResults (workspaceCalculation, "workspace",...
		["../../Results/", parameterArray{1}.folderName, "/result.mat"]);
		for i=1:length(stateWiseSIR)
			switch parameterArray{1}.model
				case "SIREDmod"					
					stateWiseSIR{i}.splineDiscovered = spline(stateWiseSIR{i}.time +...
					parameterArray{1}.startDate, stateWiseSIR{i}.SIR_vs_time(:,6)); 
				case "SIRH"
					stateWiseSIR{i}.splineDiscovered = spline(stateWiseSIR{i}.time +...
					parameterArray{1}.startDate, stateWiseSIR{i}.SIR_vs_time(:,9)); 
			endswitch
		end
		simData = zeros(length(stateWiseSIR) ,1);
		for i=1:length(stateWiseSIR)
			simData(i) = ppval(stateWiseSIR{i}.splineDiscovered, time);		
		end
		##		ratioSW = simDataSW ./ rkiData;
		##		ratioCW = simData ./ rkiData;
		
		ratio = (simDataSW - simData)./simDataSW;
		maxScalePerStep = 0.25;
		exponent = 1.5;		
		scaleFactors = ratio;
		scaleFactors(abs(scaleFactors) > maxScalePerStep^(1/exponent)) =...
		sign(scaleFactors(abs(scaleFactors) > maxScalePerStep^(1/exponent)))*...
		maxScalePerStep^(1/exponent);
		
		residuum(end+1,:) = ratio;
		
		parameterArray{1}.beta_cc_statewise .*=...
		(1 + sign(scaleFactors).*abs(scaleFactors).^exponent); 
		optVals(end+1, :) = parameterArray{1}.beta_cc_statewise;
		counter += 1;
		if counter > 1
			figure(20);
			clf;
			plot(optVals);
			pause(0.2);
		end
    %save .txt file
		fprintf(fid, "Iteration: %i took %fs with res: %f Parameters: ", counter, toc(t1), norm(residuum(end)));
		fprintf(fid, "%f ", parameterArray{1}.beta_cc_statewise);
		fprintf(fid, "\n");
		if norm(residuum(end)) < 1e-4
			fprintf(fid, "Terminate due to small enough residuum");
			break;
		end	
    fprintf(fid, "\n");
    fclose(fid);   
    %console output
    fprintf("Iteration: %i took %fs with res: %f Parameters: ", counter, toc(t1), norm(residuum(end)));
		fprintf("%f ", parameterArray{1}.beta_cc_statewise);
		fprintf("\n");
		if norm(residuum(end)) < 1e-4
			fprintf("Terminate due to small enough residuum");
			break;
		end	
	end
	
	%optimizeStateWise (parameterArray);
	
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RKIdata = getRKIdata()
	RKIread = read_case_history_RKIfiles();
  
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
    %RKIinf = RKIdata{i}.infected;
    %RKIdata{i}.infected = cumsum(cell2mat(RKIdata{i}.infected));
    RKIdata{i}.infected = cell2mat(RKIdata{i}.infected);
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
	%% the entries can be directly modified in the struct
	%% eg. set a foldername
	parameterArray{1}.model = "SIREDmod";
	%parameterArray{1}.totalRuntime = 40;
	parameterArray{1}.folderName = "SIREDmod_1stWave_Paper";
	parameterArray{1}.saveVTK = false;
	parameterArray{1}.saveDiagrams = false;
	parameterArray{1}.showDiagrams = false;
	parameterArray{1}.fullConsoleOut = true;
	parameterArray{1}.blowUp = false;
	parameterArray{1}.spatial	= "germany";
	parameterArray{1}.initial = "RKIfiles";
	parameterArray{1}.initalDistributionDate = datenum([2020, 03, 16]);
	parameterArray{1}.startDate = datenum([2020, 03, 02]);
	parameterArray{1}.endDateOpt = datenum([2020, 04, 30]);
	parameterArray{1}.reduceToStates= true;
	%parameterArray{1}.betaStateWise= false;
	parameterArray{1}.exitRestrictions= 1;
	parameterArray{1}.schoolClosing= 1;
	parameterArray{1}.beta_cross_county	= 1;
	
##	parameterArray{1}.wRatioGlobStates = 0.5; % 1 is fully statewise, 0 is only global
##	parameterArray{1}.wRatioID = 1; % 1 is only infected, 0 only death
##	parameterArray{1}.wISlope = 0.00; % weight for slope of infected last date
	parameterArray{1}.optFunGermanSum = false;
	parameterArray{1}.showStatePlots = false;	
	
	parameterArray{1}.gamma1 = 0.067; % Infected (not detected) -> Recovered
	parameterArray{1}.gamma2 = 0.04; % Quarantine -> Recovered
	parameterArray{1}.mortality = 0.006; % 0.006
	parameterArray{1}.darkFigure = 12.324015;%11.2;%9.95; %6.5;
	%parameterArray{1}.beta =  1.842096; %0.196693; %0.17;
	
	parameterArray{1}.beta =  1;
	parameterArray{1}.betaSWscaling = 1;		
  parameterArray{1}.schoolClosing	= 1;
  parameterArray{1}.contactRestrictions = 0.089893;
  parameterArray{1}.exitRestrictions	= 1;		
  parameterArray{1}.majorEvents		= 0.812462;
  parameterArray{1}.easingContactRestrictions	= 1;
  parameterArray{1}.endOfVacation	= 1;
  parameterArray{1}.liftingTravelRestictionsEU	= 1;
  parameterArray{1}.lockdownLight =1;

##	parameterArray{1}.majorEvents	= 0.318781; %0.302252; %0.2092;
##	parameterArray{1}.contactRestrictions = 0.446369;%0.499724; %0.5956;
  %parameterArray{1}.betaStatewise = [0.074255 0.002135 0.009843 0.163156 0.063369 0.093791 0.023009 0.034313 0.033831 0.147122 0.004086 0.093879 0.087273 0.189055 0.043525 0.030305];
	parameterArray{1}.betaStatewise = [0.093468 0.223027 0.078771 0.107133 0.124950 0.059428 0.046771 0.136943 0.211577 0.242791 0.123193 0.174029 0.055153 0.149495 0.069747 0.139088];
	parameterArray{1}.fullConsoleOut = false;
endfunction
