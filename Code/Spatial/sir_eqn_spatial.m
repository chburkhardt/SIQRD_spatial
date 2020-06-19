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
## @deftypefn {} {@var{retval} =} sir_eqn_spatial (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-08

function out = sir_eqn_spatial (mode, x, t, distances, population, parameter, biggestCityInCounty, cities, history)
  
	calcInteractionMatrix = false;
	switch mode
		case "lags"
			out = linspace(28, 0, 29);  
			return;  
		case "coursesOfDisease"
			%second parameter needs to be parameter
			%third parameter needs to be cities
			parameter = x;
			out = getCoursesOfDisease(parameter);
			return;  
		case "totalCourseOfDisease"    
			%second parameter needs to be parameter
			%third parameter needs to be cities
			parameter = x;
			out = getTotalCourseOfDisease(parameter);
			return;  
		case "showCourses"
			showCoursesOfDisease();
			return;  
		case "texCourses"
			tmp = showCoursesOfDisease();
			fid = fopen("../../Results/courses.dat", "w");
			fprintf(fid, "%s", tmp);
			fclose(fid);
			return;  
		case "getBeta"
			% Second parameter needs to be parameter struct, third is time, fourth is cities
			parameter = x;
			cities = distances;
			out = getBeta(parameter, cities, t);   
			return;  
		case "getDarkFigure"
			% Second parameter needs to be parameter struct, third is cities
			parameter = x;
			cities = t;
			out = getDarkFigures(cities, parameter);
			return;  
		case "getDarkFigureStatewise"
			% Second parameter needs to be parameter struct
			parameter = x;
			out = getDarkFigures(nan, parameter, "statewise");
			return;
		case "getFactorsBeta"
			% Second parameter needs to be parameter struct
			out = getFactorsBeta(x, t);   
			return;  
		case "showFactorsBeta"
			% Second parameter needs to be parameter struct
			showFactorsBeta(x);  
			return;  
		case "InteractionMatrix"
			calcInteractionMatrix = true;
			x = ones(length(cities) * 4, 1);
		endswitch 
		
		if parameter.fullConsoleOut
			fprintf("Day %f\n", t);
		endif
		
		% Parameter values
		darkFigureCityWise = getDarkFigures(cities, parameter);
		mortality = parameter.mortality;  
		beta = getBeta(parameter,cities,t);		
		gamma1 = parameter.gamma1;
		gamma2 = parameter.gamma2; 
		alpha = gamma1 ./ (darkFigureCityWise - 1);
		delta = gamma2 .* (darkFigureCityWise * mortality)./...
		(1 - darkFigureCityWise * mortality);
		n_traffic = parameter.n_traffic;
		
		beta_cross_county = getBetaCC(parameter, cities, t);
		
		model = parameter.model;
		
		if nargin < 7
			biggestCityInCounty = population;
		end
		
		if calcInteractionMatrix
			interactionMatrix = zeros(length(cities));
		end
		
		% Define variables
		S = x(1:4:end);
		I = x(2:4:end);
		E = x(3:4:end);
		D = x(4:4:end);
		
		% Define ODEs without interaction
		dSdT= -beta.*S.*E;
		dEdT= +beta.*S.*E -alpha .* E - gamma1 *E ;
		dIdT= +alpha.*E - gamma2.*I -delta.*I;
		dDdT= +delta.*I;				
		
		% add interaction term
		maxBiggestCityInCounty = 3e6;
		for i = 1:size(distances, 1) - 1    
			tmp_sparse = [spdiags(distances, i)(i+1:end);...
			spdiags(distances, i - length(distances))(1:i)];
			indices = tmp_sparse > 0;    
			
			if beta_cross_county > 0
				% differentiate if beta is a single scalar (for constant beta) or is a vector (spatially varying beta)
				if !parameter.MobNetworkModel
					if length(beta)==1
						error("This point should not be reached");
						betaCross =  beta_cross_county * beta   ./...
						(tmp_sparse(indices)).^n_traffic ...
						.* biggestCityInCounty(indices).^2 / maxBiggestCityInCounty^2;
					else
						betaCross =  beta_cross_county * beta(indices)   ./...
						(tmp_sparse(indices)).^n_traffic ...
						.* biggestCityInCounty(indices).*shift(biggestCityInCounty, -i)(indices) /...
						maxBiggestCityInCounty^2;
					endif
				else
					exp_alpha = [0.46, 0.35];
					exp_gamma = [0.64, 0.37];
					factor_r = [82, 1];
					rcutoff = 300;
					betaCross = zeros(sum(indices),1);
					betaCrossfull = zeros(length(indices),1);
					index_low = and((tmp_sparse < rcutoff), (tmp_sparse > 0));
					index_up = (tmp_sparse >= rcutoff);
					betaCrossfull (index_low) = beta_cross_county(index_low) .*...
					(population(index_low).^exp_alpha(1))...
					.* (shift(population, -i)(index_low).^exp_gamma(1))...
					./ exp(tmp_sparse(index_low) ./factor_r(1)) /...
					maxBiggestCityInCounty^(exp_gamma(1) + exp_alpha(1));
					betaCrossfull (index_up) = beta_cross_county(index_up) .*...
					(population(index_up).^exp_alpha(2))...
					.* (shift(population, -i)(index_up).^exp_gamma(2))...
					./ exp(tmp_sparse(index_up) ./factor_r(2)) /...
					maxBiggestCityInCounty^(exp_gamma(2) + exp_alpha(2));
					betaGeoMean = (beta .* shift(beta, -i)).^0.5;
					betaCross = betaCrossfull(indices).*betaGeoMean(indices);
				end
        dSdT(indices) -= betaCross ...
        .* S(indices) .* shift(E, -i)(indices);
        dEdT(indices) += betaCross ...
        .* S(indices) .* shift(E, -i)(indices);		
        
        if calcInteractionMatrix
          tmp = zeros(size(dSdT));
          tmp(indices) = betaCross;
          for m = 1:length(tmp)
            interactionMatrix(m, mod(m+i-1, length(cities)) + 1) = tmp(m);
          end
        end
			endif
		end
		% Return gradients
    out = reshape([dSdT, dIdT, dEdT, dDdT]', length(x), 1);
		
		if calcInteractionMatrix
			out = interactionMatrix;
		end 
		
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function betaCityWise = getBeta(parameter,cities,t)
		if isfield(parameter, "betaCorr1")
			factorsBetaCorrection = [parameter.betaCorr1, parameter.betaCorr2,...
			parameter.betaCorr3, parameter.betaCorr4, parameter.betaCorr5,...
			parameter.betaCorr6, parameter.betaCorr7, parameter.betaCorr8,...
			parameter.betaCorr9, parameter.betaCorr10, parameter.betaCorr11,...
			parameter.betaCorr12, parameter.betaCorr13, parameter.betaCorr14,...
			parameter.betaCorr15, parameter.betaCorr16, 1];
		elseif isfield(parameter, "betaStatewise")
			factorsBetaCorrection = parameter.betaStatewise ./ parameter.beta;
			if isfield(parameter, "betaSWscaling")
				factorsBetaCorrection *= parameter.betaSWscaling;
			end	
		end
		factorsBetaRestrictions = getFactorsBeta (parameter, t);
		betaCityWise = ones(length(cities), 1) * parameter.beta;
		
		for i=1:length(cities)
			betaCityWise(i) *= factorsBetaRestrictions(floor(cities{i}.ags/1000));
			if exist("factorsBetaCorrection", "var")
				betaCityWise(i) *= factorsBetaCorrection(floor(cities{i}.ags/1000));
			end
		endfor  
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function datesSchoolClosing = getDatesSchoolClosing ()
		datesSchoolClosing = inf(17,2);
		datesSchoolClosing(1,1) = datenum(2020,03,16);
		datesSchoolClosing(2,1) = datenum(2020,03,16);
		datesSchoolClosing(3,1) = datenum(2020,03,16);
		datesSchoolClosing(4,1) = datenum(2020,03,16);
		datesSchoolClosing(5,1) = datenum(2020,03,16);
		datesSchoolClosing(6,1) = datenum(2020,03,16);
		datesSchoolClosing(7,1) = datenum(2020,03,16);
		datesSchoolClosing(8,1) = datenum(2020,03,17);
		datesSchoolClosing(9,1) = datenum(2020,03,16);
		datesSchoolClosing(10,1) = datenum(2020,03,16);
		datesSchoolClosing(11,1) = datenum(2020,03,17);
		datesSchoolClosing(12,1) = datenum(2020,03,18);
		datesSchoolClosing(13,1) = datenum(2020,03,16);
		datesSchoolClosing(14,1) = datenum(2020,03,16);
		datesSchoolClosing(15,1) = datenum(2020,03,16);
		datesSchoolClosing(16,1) = datenum(2020,03,17);
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function datesContactRestrictions = getDatesContactRestrictions ()
		datesContactRestrictions = inf(17,2);
		datesContactRestrictions(1,1) = datenum(2020,03,22);
		datesContactRestrictions(2,1) = datenum(2020,03,22);
		datesContactRestrictions(3,1) = datenum(2020,03,22);
		datesContactRestrictions(4,1) = datenum(2020,03,22);
		datesContactRestrictions(5,1) = datenum(2020,03,22);
		datesContactRestrictions(6,1) = datenum(2020,03,22);
		datesContactRestrictions(7,1) = datenum(2020,03,22);
		datesContactRestrictions(8,1) = datenum(2020,03,22);
		datesContactRestrictions(9,1) = datenum(2020,03,22);
		datesContactRestrictions(10,1) = datenum(2020,03,22);
		datesContactRestrictions(11,1) = datenum(2020,03,22);
		datesContactRestrictions(12,1) = datenum(2020,03,22);
		datesContactRestrictions(13,1) = datenum(2020,03,22);
		datesContactRestrictions(14,1) = datenum(2020,03,22);
		datesContactRestrictions(15,1) = datenum(2020,03,22);
		datesContactRestrictions(16,1) = datenum(2020,03,22);
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	##  startDatesExitRestrictions(1) = 0-calculationStartDate;
	function datesMajorEventsRestrictions = getDatesMajorEventsRestrictions ()
		% values equal to inf mean that there have been no exit restrictions
		datesMajorEventsRestrictions = [datenum(2020,03,8) * ones(17,1) ,inf(17,1)];
		##  datesMajorEventsRestrictions(10,1) = inf;
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	##  startDatesExitRestrictions(1) = 0-calculationStartDate;
	function datesExitRestrictions = getDatesExitRestrictions ()
		% values equal to inf mean that there have been no exit restrictions
		datesExitRestrictions = inf(17,2);
		datesExitRestrictions(9,1) = datenum(2020,03,22);
		datesExitRestrictions(10,1) = datenum(2020,03,20);
		datesExitRestrictions(10,2) = datenum(2020,04,28);
		datesExitRestrictions(11,1) = datenum(2020,03,22);
		datesExitRestrictions(11,2) = datenum(2020,04,22);
		datesExitRestrictions(12,1) = datenum(2020,03,22);
		datesExitRestrictions(14,1) = datenum(2020,03,22);
		datesExitRestrictions(14,2) = datenum(2020,04,20);
		datesExitRestrictions(15,1) = datenum(2020,03,22);
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function factorsBeta = getFactorsBeta (parameter, t)
		if !exist("datesSchoolClosing", "var")
			persistent datesSchoolClosing = getDatesSchoolClosing (parameter.startDate);
			persistent datesContactRestrictions = getDatesContactRestrictions (parameter.startDate);
			persistent datesExitRestrictions = getDatesExitRestrictions (parameter.startDate);
			persistent datesMajorEventsRestrictions = getDatesMajorEventsRestrictions (parameter.startDate);
		end
		
		t+=parameter.startDate;
		
		factorsBeta = ones(size(datesSchoolClosing, 1), 1);
		factorsBeta(and(datesSchoolClosing(:,1) < t,...
		t < datesSchoolClosing(:,2))) *= parameter.schoolClosing;
		
		factorsBeta(and(datesContactRestrictions(:,1) < t,...
		t < datesContactRestrictions(:,2))) *= parameter.contactRestrictions;
		
		factorsBeta(and(datesExitRestrictions(:,1) < t,...
		t < datesExitRestrictions(:,2))) *= parameter.exitRestrictions;
		
		factorsBeta(and(datesMajorEventsRestrictions(:,1) < t,...
		t < datesMajorEventsRestrictions(:,2))) *= parameter.majorEvents;    
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function showFactorsBeta (parameter)
		datesSchoolClosing = getDatesSchoolClosing (parameter);
		datesContactRestrictions = getDatesContactRestrictions (parameter);
		datesExitRestrictions = getDatesExitRestrictions (parameter);
		datesMajorEventsRestrictions = getDatesMajorEventsRestrictions (parameter);
		
		t= -10+min([datesExitRestrictions(!isinf(datesExitRestrictions));...
		datesContactRestrictions(!isinf(datesContactRestrictions));...
		datesSchoolClosing(!isinf(datesSchoolClosing));...
		datesMajorEventsRestrictions(!isinf(datesMajorEventsRestrictions))]):...  
		max([datesExitRestrictions(!isinf(datesExitRestrictions));...
		datesContactRestrictions(!isinf(datesContactRestrictions));...
		datesSchoolClosing(!isinf(datesSchoolClosing));...
		datesMajorEventsRestrictions(!isinf(datesMajorEventsRestrictions))])+10;
		
		population = [2896712, 1841179, 7982448, 682986, 17932651, 6265809,...
		4084844, 11069533, 13076721, 990509, 3644826, 2511917,...
		1609675, 4077937, 2208321, 2143145];
		close all;
		figure;
		hold on;
		factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesSchoolClosing(:,1) < t,...
		t < datesSchoolClosing(:,2))) *= parameter.schoolClosing;
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
		
		factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesContactRestrictions(:,1) < t,...
		t < datesContactRestrictions(:,2))) *= parameter.contactRestrictions;
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
		
		factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesExitRestrictions(:,1) < t,...
		t < datesExitRestrictions(:,2))) *= parameter.exitRestrictions;
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
		
		factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesMajorEventsRestrictions(:,1) < t,...
		t < datesMajorEventsRestrictions(:,2))) *= parameter.majorEvents; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
		
		legend("Schools closing", "Contact restrictions", "Exit restrictions",...
		"Major events restrictions", "location", "northeastoutside");
		datetick ("x", "dd.mm.yyyy", "keeplimits");
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function darkFigureCityWise = getDarkFigures(cities, parameter, statewise)			
		factorsDfCorrection = [ 1.06263   0.92143   1.13461   1.39972   1.01612   1.18285 ...
    0.68441   1.15415   1.13741   1.23239   0.57678 1.11515   0.64394   0.91610   0.71415   1.10816 1];
		factorsDfCorrection /= mean(factorsDfCorrection);
	
		if nargin == 3
			if strcmp(statewise, "statewise")
				darkFigureCityWise = parameter.darkFigure * factorsDfCorrection;
				return;
			end
		end
    
		darkFigureCityWise = ones(length(cities), 1) * parameter.darkFigure;	
		for i=1:length(cities)
			darkFigureCityWise(i) *= factorsDfCorrection(floor(cities{i}.ags/1000));
		endfor  	
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function beta_cross_county = getBetaCC(parameter, cities, t)
		beta_cross_county = ones(length(cities), 1) * parameter.beta_cross_county;
		
		if isfield(parameter, "beta_cc_statewise")
			beta_cc_statewise	= parameter.beta_cc_statewise;
			if length(beta_cc_statewise) == 16
				beta_cc_statewise(end+1) = 1; % other countries
			endif			
			for i=1:length(cities)
				beta_cross_county(i) = beta_cc_statewise(floor(cities{i}.ags/1000));
			end
		end
		
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	
