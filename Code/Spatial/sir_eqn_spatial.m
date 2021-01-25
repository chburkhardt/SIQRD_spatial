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

function out = sir_eqn_spatial (mode, x, t, distances, population, parameter,...
	biggestCityInCounty, cities, history)
  
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
			%second parameter needs to be parameter
			parameter = x;
			showCoursesOfDisease(parameter);
			return;  
		case "texCourses"
			%second parameter needs to be parameter
			parameter = x;
			tmp = showCoursesOfDisease(parameter);
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
		case "getBetaCC"
			% Second parameter needs to be parameter struct, third is time, fourth is cities
			parameter = x;
			cities = distances;
			out = getBetaCC(parameter, cities, t);   
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
    case "getDelta"
      parameter = x;
			cities = t;
			out = getDelta(cities, parameter);
		case "getFactorsBeta"
			% Second parameter needs to be parameter struct
			out = getFactorsBeta(x, t);   
			return;  
		case "showFactorsBeta"
			% Second parameter needs to be parameter struct
			showFactorsBeta(x);  
			return;  
		case "betas"	
			%Second parameter needs to be parameter, third distances	
			parameter = x;	
			cities = distances;
			out = getBeta(parameter,cities,t);	
			return;	
		case "InteractionMatrix"
			calcInteractionMatrix = true;
			switch parameter.model
				case "SIR"
					x = ones(length(cities) * 2, 1);	
				case "SIRH"
					x = ones(length(cities), 1);
					lags = sir_eqn_spatial("lags");
					history = ones(length(cities), length(lags));
				case {"SIRED","SIREDmod","SIREDLiterature"}
					x = ones(length(cities) * 4, 1);
				case "SIREDYO"
					x = ones(length(cities) * 8, 1);
			endswitch
     case "mortalityFactor"
			% Second parameter needs to be parameter struct
      parameter = x; 
      cities = t;
      t = distances;
			out = getMortalityCourse(parameter, cities, t);  
			return;  
		endswitch 
    
		
		% update status message
		if parameter.fullConsoleOut
			fprintf("\rDay %.1f ", t);
		endif
		% x=[S_1, I_1, S_2, I_2 ,..., S_n, I_n]
		
		% Parameter values
		% to do: extend the Beta model to be temporally dependend
		
		darkFigureCityWise = getDarkFigures(cities, parameter);
    mortality = parameter.mortality;
    if isfield(parameter,"considereAgeOfNewlyInfected")
      if parameter.considereAgeOfNewlyInfected
        factorMortality = getMortalityCourse(parameter, cities, t);
        mortality *= factorMortality;
      end  
    end
		beta = getBeta(parameter,cities,t);		
		gamma1 = parameter.gamma1;
		gamma2 = parameter.gamma2; 
		%	gamma2 = 0.045 - 0.00027*darkFigureCityWise; %!!!!!!!!!!!!
		alpha = gamma1 ./ (darkFigureCityWise - 1);
    
    if isfield(parameter,"considereAgeOfNewlyInfected")
      if parameter.considereAgeOfNewlyInfected
        delta = gamma2 .* (darkFigureCityWise .* mortality)./...
		    (1 - darkFigureCityWise .* mortality);
      else
        delta = gamma2 .* (darkFigureCityWise * mortality)./...
		    (1 - darkFigureCityWise * mortality);
      end 
    else
      delta = gamma2 .* (darkFigureCityWise * mortality)./...
		  (1 - darkFigureCityWise * mortality);
    end
    
    %SIREDLiterature Model
    if strcmp(parameter.model, "SIREDLiterature")
      %delta = 1/17*0.015;%*(1/17);%CFR * mean duration to die after entering I compartment      
      if isfield(parameter,"considereAgeOfNewlyInfected")
        delta = parameter.delta;
      else
        parameter.delta = 0.005;
      end
      delta = getDelta(cities, parameter);
      alpha = 1/3; %mean latency period = 3 days
    end
    
		c_s = parameter.c_s;
		c_i = parameter.c_i;
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
		% for SIR model
		switch model
			case "SIR"
				S = x(1:2:end);
				I = x(2:2:end);
			case "SIRH"
				S = x;
				S_hist = history;
			case {"SIRED","SIREDmod","SIREDLiterature"}
				S = x(1:4:end);
				I = x(2:4:end);
				E = x(3:4:end);
				D = x(4:4:end);
			case "SIREDYO"
				S = x(1:8:end);
				I = x(2:8:end);
				E = x(3:8:end);
				D = x(4:8:end);
				Sold = x(5:8:end);
				Iold = x(6:8:end);
				Eold = x(7:8:end);
				Dold = x(8:8:end);
		endswitch
		
		% Define ODEs without diffusion
		% for SIR model
		switch model
			case "SIR"
				dSdT=-beta.*S.*I;
				dIdT=beta.*S.*I-gamma1*I;
			case "SIRH"      
				courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);        
				oldDeltaS = -diff(S_hist, 1, 2);
				
				##        dSdT=-beta*S.*I;   
				dSdT=-beta.*S.*(oldDeltaS*courseOfDisease.infective');  
				if calcInteractionMatrix
					interactionMatrix += diag(beta);
				end    
			case "SIRED"
				dSdT= -beta.*S.*(I+E);
				dEdT= +beta.*S.*(I+E) -alpha.*E;
				dIdT= +alpha.*E - gamma1*I -delta.*I;
				dDdT= +delta.*I;
      case "SIREDLiterature"
				dSdT= -beta.*S.*I;
				dEdT= +beta.*S.*I - alpha.*E;
				dIdT= +alpha.*E - gamma1*I -delta.*I;
				dDdT= +delta.*I;
			case "SIREDmod"
				dSdT= -beta.*S.*E;
				dEdT= +beta.*S.*E -alpha .* E - gamma1 *E ;
				dIdT= +alpha.*E - gamma2.*I -delta.*I;
				dDdT= +delta.*I;				
				
				if calcInteractionMatrix
					interactionMatrix += diag(beta - alpha - gamma1);
				end
			case "SIREDYO"
				dSdT= -beta.*S.*(E+I) - beta.*S.*(Eold+Iold);
				dEdT= +beta.*S.*(E+I) + beta.*S.*(Eold+Iold) -alpha .* E;
				dIdT= +alpha.*E - gamma1*I -delta.*I;
				dDdT= +delta.*I;
				dSolddT= -beta.*Sold.*(I+E) -beta.*Sold.*(Iold+Eold);
				dEolddT= +beta.*Sold.*(I+E) +beta.*Sold.*(Iold+Eold) -alpha .* Eold;
				dIolddT= +alpha.*Eold - gamma1*Iold -delta.*Iold;
				dDolddT= +delta.*Iold;
		endswitch
		
		% add diffusion term
		maxBiggestCityInCounty = 3e6;
		for i = 1:size(distances, 1) - 1    
			tmp_sparse = [spdiags(distances, i)(i+1:end);...
			spdiags(distances, i - length(distances))(1:i)];
			indices = tmp_sparse > 0;    
			
			if or(c_s > 0, c_i > 0)
				error("Diffusion not implemented anymore, set c_i = c_s = 0");
				dSdT(indices) -= (c_s ./ tmp_sparse(indices)).^n_traffic...
				.* (S(indices) - shift(S, -i)(indices)) .* shift(population, -i)(indices);
				dIdT(indices) -= (c_i ./ tmp_sparse(indices)).^n_traffic...
				.*(I(indices) - shift(I, -i)(indices)) .* shift(population, -i)(indices);
			endif
			
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
						##            .* biggestCityInCounty(indices).^2 / maxBiggestCityInCounty^2;
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
				switch model
					case "SIR"          
						dSdT(indices) -= betaCross ...
						.* S(indices) .* shift(I, -i)(indices);
						dIdT(indices) += betaCross ...
						.* S(indices) .* shift(I, -i)(indices);
					case "SIRH"         
						courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);        
						oldDeltaS = -diff(S_hist, 1, 2);            
						effectiveInfective = oldDeltaS*courseOfDisease.infective';
						
						dSdT(indices) -= betaCross ...
						.* S(indices) .* shift(effectiveInfective, -i)(indices);
					case "SIRED"
						% define a factor to scale cross county infections for infected down
						% because exposed people carry around the disease (without knowing)
						factorCrossInfected = 0.1;
						dSdT(indices) -= factorCrossInfected * betaCross ...
						.* S(indices) .* shift(I, -i)(indices);
						dSdT(indices) -= betaCross ...
						.* S(indices) .* shift(E, -i)(indices);
						dEdT(indices) += factorCrossInfected * betaCross ...
						.* S(indices) .* shift(I, -i)(indices);
						dEdT(indices) += betaCross ...
						.* S(indices) .* shift(E, -i)(indices);
          case "SIREDLiterature"
						dSdT(indices) -= betaCross ...
						.* S(indices) .* shift(I, -i)(indices);
						dEdT(indices) += betaCross ...
						.* S(indices) .* shift(I, -i)(indices);	
					case "SIREDmod"
						dSdT(indices) -= betaCross ...
						.* S(indices) .* shift(E, -i)(indices);
						dEdT(indices) += betaCross ...
						.* S(indices) .* shift(E, -i)(indices);		
						
					case "SIREDYO"
						% define a factor to scale cross county infections for infected down
						% because exposed people carry around the disease (without knowing)
						factorCrossInfected = 0.1;
						dSdT(indices) -= factorCrossInfected * betaCross ...
						.* S(indices) .* (shift(I, -i)(indices) + shift(Iold, -i)(indices));
						dSdT(indices) -= betaCross ...
						.* S(indices) .* (shift(E, -i)(indices) + shift(Eold, -i)(indices));
						dEdT(indices) += factorCrossInfected * betaCross ...
						.* S(indices) .* (shift(I, -i)(indices) + shift(Iold, -i)(indices));
						dEdT(indices) += betaCross ...
						.* S(indices) .* (shift(E, -i)(indices) + shift(Eold, -i)(indices));
						dSolddT(indices) -= factorCrossInfected * betaCross ...
						.* Sold(indices) .* (shift(I, -i)(indices) + shift(Iold, -i)(indices));
						dSolddT(indices) -= betaCross ...
						.* Sold(indices) .* (shift(E, -i)(indices) + shift(Eold, -i)(indices));
						dEolddT(indices) += factorCrossInfected * betaCross ...
						.* Sold(indices) .* (shift(I, -i)(indices) + shift(Iold, -i)(indices));
						dEolddT(indices) += betaCross ...
						.* Sold(indices) .* (shift(E, -i)(indices) + shift(Eold, -i)(indices));
				endswitch
				
				if calcInteractionMatrix
					% dSmdT abhängig von E in stadt m+i
					tmp = zeros(size(dSdT));
					tmp(indices) = betaCross;
					##							tmp ./= tmp;
					##							tmp *= i;
					##							tmp+=(1:16)';
					for m = 1:length(tmp)
						interactionMatrix(m, mod(m+i-1, length(cities)) + 1) = tmp(m);
					end
				end
			endif
		end
		
		% Return gradients
		switch model
			case "SIR"
				out = reshape([dSdT, dIdT]', length(x), 1);
			case "SIRH"
				out = dSdT;
			case {"SIRED","SIREDmod","SIREDLiterature"}
				out = reshape([dSdT, dIdT, dEdT, dDdT]', length(x), 1);
			case "SIREDYO"
				out = reshape([dSdT, dIdT, dEdT, dDdT, dSolddT, dIolddT, dEolddT, dDolddT]', length(x), 1);
		endswitch
		
		if calcInteractionMatrix
			out = interactionMatrix;
		end 
		
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function out = getCoursesOfDisease(parameter)
		% all vectors must have length length(lags)-1
		% case 1
		out{1}.infective = [ linspace(0,0,19), linspace(0,1,4), linspace(1,1,3), linspace(1,0,1), linspace(0,0,1)]; 
		out{1}.infective /=  sum(filter([1,1]/2,1, out{1}.infective));
		out{1}.symptoms = [ linspace(0,0,15), linspace(1,1,9), linspace(0,0,4)]; 
		out{1}.hospital = linspace(0,0,28);
		out{1}.icu = linspace(0,0,28);
		out{1}.recovering = diff([linspace(0,0,11), linspace(0,1,4), linspace(1,1,14)]);
		out{1}.dying = diff(linspace(0,0,29));
		out{1}.discovering = diff([linspace(0,0,19), linspace(0,1,8), linspace(1,1,2)]);
		out{1}.share = 0.95;
		% case 2 Hospital
		out{2}.infective = [ linspace(0,0,20), linspace(0,1,3), linspace(1,1,3), linspace(1,0,1), linspace(0,0,1)]; 
		out{2}.infective /=  sum(filter([1,1]/2,1, out{2}.infective));
		out{2}.symptoms = [ linspace(0,0,6), linspace(1,1,18), linspace(0,0,4)]; 
		out{2}.hospital = [ linspace(0,0,6), linspace(1,1,14), linspace(0,0,8)];
		out{2}.icu = linspace(0,0,28);
		out{2}.recovering = diff([linspace(0,0,2), linspace(0,1,4), linspace(1,1,23)]);
		out{2}.dying = diff(linspace(0,0,29));
		out{2}.discovering = diff([linspace(0,0,19), linspace(0,1,8), linspace(1,1,2)]);
		out{2}.share = 0.04;
		% case 3 ICO survive
		out{3}.infective = [ linspace(0,0,20), linspace(0,1,3), linspace(1,1,3), linspace(1,0,1), linspace(0,0,1)]; 
		out{3}.infective /=  sum(filter([1,1]/2,1, out{3}.infective));
		out{3}.symptoms = [ linspace(0,0,9), linspace(1,1,15), linspace(0,0,4)]; 
		out{3}.hospital = [ linspace(0,0,6), linspace(1,1,14), linspace(0,0,8)];
		out{3}.icu = [ linspace(0,0,9), linspace(1,1,10), linspace(0,0,9)];
		out{3}.recovering = diff([linspace(0,0,2), linspace(0,1,4), linspace(1,1,23)]);
		out{3}.dying = diff(linspace(0,0,29));
		out{3}.discovering = diff([linspace(0,0,19), linspace(0,1,8), linspace(1,1,2)]);
		out{3}.share = 0.01 - parameter.mortality;
		% case 4 ICU die
		out{4}.infective = [ linspace(0,0,20), linspace(0,1,3), linspace(1,1,3), linspace(1,0,1), linspace(0,0,1)]; 
		out{4}.infective /=  sum(filter([1,1]/2,1, out{4}.infective));
		out{4}.symptoms = [ linspace(0,0,9), linspace(1,1,15), linspace(0,0,4)]; 
		out{4}.hospital = [ linspace(0,0,9), linspace(1,1,11), linspace(0,0,8)];
		out{4}.icu = [ linspace(0,0,9), linspace(1,1,10), linspace(0,0,9)];
		out{4}.recovering = diff(linspace(0,0,29));
		out{4}.dying = diff([linspace(0,0,5), linspace(0,1,4), linspace(1,1,20)]);
		out{4}.discovering = diff([linspace(0,0,19), linspace(0,1,8), linspace(1,1,2)]);
		out{4}.share = parameter.mortality;
		if parameter.mortality > 0.01
			error("Max mortality 0.01 exeeded");	
		end	
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function out = getTotalCourseOfDisease(parameter)	
    %check if totalCourseOfDisease should be recalculated
		if ~isfield(parameter,"courseNew")
      if exist("tmp", "var")
        out = tmp;
        return;
      end	
    end
		coursesOfDisease = getCoursesOfDisease(parameter);
		persistent tmp
    %delete persistent tmp if age dependency of mortality should be considered
    %Only called if considereAgeOfNewlyInfected true
    if isfield(parameter,"courseNew")
      if parameter.courseNew
          clear("tmp");
       end
    end
		tmp.infective = coursesOfDisease{1}.share * coursesOfDisease{1}.infective;
		tmp.symptoms = coursesOfDisease{1}.share * coursesOfDisease{1}.symptoms;
		tmp.hospital = coursesOfDisease{1}.share * coursesOfDisease{1}.hospital;
		tmp.icu = coursesOfDisease{1}.share * coursesOfDisease{1}.icu;
		tmp.recovering = coursesOfDisease{1}.share * coursesOfDisease{1}.recovering;
		tmp.dying = coursesOfDisease{1}.share * coursesOfDisease{1}.dying;
		tmp.discovering = coursesOfDisease{1}.share * coursesOfDisease{1}.discovering;
		for i=2:length(coursesOfDisease)
			tmp.infective += coursesOfDisease{i}.share * coursesOfDisease{i}.infective;
			tmp.symptoms += coursesOfDisease{i}.share * coursesOfDisease{i}.symptoms;
			tmp.hospital += coursesOfDisease{i}.share * coursesOfDisease{i}.hospital;
			tmp.icu += coursesOfDisease{i}.share * coursesOfDisease{i}.icu;
			tmp.recovering += coursesOfDisease{i}.share * coursesOfDisease{i}.recovering;
			tmp.dying += coursesOfDisease{i}.share * coursesOfDisease{i}.dying;
			tmp.discovering += coursesOfDisease{i}.share * coursesOfDisease{i}.discovering;
		end 
    	out = tmp;
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function out = showCoursesOfDisease(parameter)
		lags = sir_eqn_spatial("lags");
		x = sir_eqn_spatial("coursesOfDisease", parameter);
		texOut = "lag ";
		texVals = lags(1:end-1)';%length(x) * (length(fieldnames(x{1})) - 1));
		colors = ["-xk";"-xr";"-xg";"-xb";"-xy";"-xm";"-xc"];
		nrows = ceil(length(x)^0.5);
		ncols = ceil(length(x)/nrows);
		figure;
		for i=1:length(x)    
			% write out parameters
			names = fieldnames(x{i});    
			for j=1:length(names)-1
				texOut = sprintf("%scase_%i_%s ", texOut, i, names{j});
				subplot(2,2,i), hold on;
				if strfind(names{j}, "ing")
					plot(lags(1:end-1), sum(x{i}.(names(j){1}))-...
					cumsum(x{i}.(names(j){1})), colors(j,:));
					texVals (:, end+1) = (sum(x{i}.(names(j){1}))-...
					cumsum(x{i}.(names(j){1})));
				else
					plot(lags(1:end-1), x{i}.(names(j){1}), colors(j,:));
					texVals (:, end+1) = x{i}.(names(j){1});
				end
			end    
			hold off;
			legend(names{1:end-1}, "location", "northeastoutside");
			title(["course of disease ", num2str(i), " concerns ", num2str(x{i}.share), "%"]);
		end  
		if nargout == 1
			for i=1:size(texVals, 1)
				nextLine = sprintf("%f ", texVals(i, :));
				texOut = sprintf("%s\n%s", texOut, nextLine);
			end
			out = texOut;
		end  
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function betaCityWise = getBeta(parameter,cities,t)
		
		##		%% BIP
		##		bipPerPerson = [33.712,66.879,38.423,49.215,39.678,46.923,35.457,...
		##		47.290,48.323,36.684,41.967,29.541,28.940,31.453,28.880,29.883, 40];	
		##		factorsBetaBIP = exp((bipPerPerson / max(bipPerPerson)-1)/5);
		##		factorsBetaBIP = factorsBetaBIP + mean(1 - factorsBetaBIP);
		##		
		##		%% Population density
		##		factorsBetaPopulationDensity = zeros(length(cities), 1);
		##		for i=1:length(cities)
		##			factorsBetaPopulationDensity(i) = cities{i}.population/cities{i}.area;
		##		end
		##		factorsBetaPopulationDensity = exp((factorsBetaPopulationDensity /...
		##		max(factorsBetaPopulationDensity )-1)/5);
		##		factorsBetaPopulationDensity = factorsBetaPopulationDensity +...
		##		mean(1 - factorsBetaPopulationDensity);
		##		
		##		%% Bevölkerungswachstum
		##		populationGrowth = [6,3.9,4.1,-4.2,0.3,3.1,4.4,5.7,8,-7.7,...
		##		-2.1,-3.7,-15.4,-13.4,-20,-15.6, 0];
		##		factorsPopulationGrowth = exp((populationGrowth / max(populationGrowth)-1)/8);
		##		factorsPopulationGrowth = factorsPopulationGrowth + mean(1 - factorsPopulationGrowth);
		
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
			##    betaCityWise(i) *= factorsBetaBIP(floor(cities{i}.ags/1000));
			##    betaCityWise(i) *= factorsBetaPopulationDensity(cities{i}.ags);
			##    betaCityWise(i) *= factorsPopulationGrowth(floor(cities{i}.ags/1000));
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
  %dates where meetings with up to 10 persons are allowed again
  function datesEasingContactRestrictions = getDatesEasingContactRestrictions ()
		datesEasingContactRestrictions = inf(17,2);
		datesEasingContactRestrictions(1,1) = datenum(2020,06,08);
		datesEasingContactRestrictions(2,1) = datenum(2020,07,01);
		datesEasingContactRestrictions(3,1) = datenum(2020,06,22);
		datesEasingContactRestrictions(4,1) = datenum(2020,06,26);
		datesEasingContactRestrictions(5,1) = datenum(2020,05,30);
		datesEasingContactRestrictions(6,1) = datenum(2020,06,11);
		datesEasingContactRestrictions(7,1) = datenum(2020,06,10);
		datesEasingContactRestrictions(8,1) = datenum(2020,06,10);
		datesEasingContactRestrictions(9,1) = datenum(2020,06,17);
		datesEasingContactRestrictions(10,1) = datenum(2020,06,01);
		datesEasingContactRestrictions(11,1) = datenum(2020,05,30);%nur 5 Personen
		datesEasingContactRestrictions(12,1) = datenum(2020,05,25);
		datesEasingContactRestrictions(13,1) = datenum(2020,06,05);
		datesEasingContactRestrictions(14,1) = datenum(2020,06,06);
		datesEasingContactRestrictions(15,1) = datenum(2020,05,28);
		datesEasingContactRestrictions(16,1) = datenum(2020,06,13);
	endfunction
  function datesLiftingTravelRestictionsEU = getDatesLiftingTravelRestictionsEU ()
		% dates where travel restictions for EU countreis were lifted
		datesLiftingTravelRestictionsEU = [datenum(2020,06,15) * ones(17,1) ,inf(17,1)];
	endfunction
  function datesEndOfVacation= getdatesEndOfVacation ()
    % Dates on which the vacations started
##    datesStartOfVacation = inf(17,2);
##		datesStartOfVacation(1,1) = datenum(2020,06,29);
##		datesStartOfVacation(2,1) = datenum(2020,06,25);
##		datesStartOfVacation(3,1) = datenum(2020,07,16);
##		datesStartOfVacation(4,1) = datenum(2020,07,16);
##		datesStartOfVacation(5,1) = datenum(2020,06,29);
##		datesStartOfVacation(6,1) = datenum(2020,06,07);
##		datesStartOfVacation(7,1) = datenum(2020,06,07);
##		datesStartOfVacation(8,1) = datenum(2020,07,30);
##		datesStartOfVacation(9,1) = datenum(2020,07,27);
##		datesStartOfVacation(10,1) = datenum(2020,06,07);
##		datesStartOfVacation(11,1) = datenum(2020,06,25);
##		datesStartOfVacation(12,1) = datenum(2020,06,25);
##		datesStartOfVacation(13,1) = datenum(2020,06,22);
##		datesStartOfVacation(14,1) = datenum(2020,07,20);
##		datesStartOfVacation(15,1) = datenum(2020,07,16);
##		datesStartOfVacation(16,1) = datenum(2020,07,20);

    % Dates on which the vacations ended
    datesEndOfVacation = inf(17,2);
		datesEndOfVacation(1,1) = datenum(2020,08,08);
		datesEndOfVacation(2,1) = datenum(2020,08,05);
		datesEndOfVacation(3,1) = datenum(2020,08,26);
		datesEndOfVacation(4,1) = datenum(2020,08,26);
		datesEndOfVacation(5,1) = datenum(2020,08,11);
		datesEndOfVacation(6,1) = datenum(2020,08,14);
		datesEndOfVacation(7,1) = datenum(2020,08,14);
		datesEndOfVacation(8,1) = datenum(2020,09,12);
		datesEndOfVacation(9,1) = datenum(2020,09,07);
		datesEndOfVacation(10,1) = datenum(2020,08,14);
		datesEndOfVacation(11,1) = datenum(2020,08,07);
		datesEndOfVacation(12,1) = datenum(2020,08,08);
		datesEndOfVacation(13,1) = datenum(2020,08,01);
		datesEndOfVacation(14,1) = datenum(2020,08,28);
		datesEndOfVacation(15,1) = datenum(2020,08,26);
		datesEndOfVacation(16,1) = datenum(2020,08,29);
	endfunction
  function datesMasks = getDatesMasks ()
    %Day on which the recommendation to wear a mask was given
    %only required for SIR model
    datesMasks = [datenum(2020,04,16) * ones(17,1) ,inf(17,1)];
  endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	##  startDatesExitRestrictions(1) = 0-calculationStartDate;
	function datesMajorEventsRestrictions = getDatesMajorEventsRestrictions ()
		% values equal to zero mean that there have been no exit restrictions
		datesMajorEventsRestrictions = [datenum(2020,03,8) * ones(17,1) ,inf(17,1)];
		##  datesMajorEventsRestrictions(10,1) = inf;
	endfunction
  function datesLockdownLight = getDatesLockdownLight ()
		% values equal to zero mean that there have been no exit restrictions
		datesLockdownLight = [datenum(2020,11,02) * ones(17,1) ,inf(17,1)];
		##  datesMajorEventsRestrictions(10,1) = inf;
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	##  startDatesExitRestrictions(1) = 0-calculationStartDate;
	function datesExitRestrictions = getDatesExitRestrictions ()
		% values equal to zero mean that there have been no exit restrictions
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
  
  % Idee: Maskenpflicht als Restriktion einführen?
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function factorsBeta = getFactorsBeta (parameter, t)
		if !exist("datesSchoolClosing", "var")
			persistent datesSchoolClosing = getDatesSchoolClosing (parameter.startDate);
			persistent datesContactRestrictions = getDatesContactRestrictions (parameter.startDate);
			persistent datesExitRestrictions = getDatesExitRestrictions (parameter.startDate);
			persistent datesMajorEventsRestrictions = getDatesMajorEventsRestrictions (parameter.startDate);
      persistent datesEasingContactRestrictions = getDatesEasingContactRestrictions (parameter.startDate);
      persistent datesEndOfVacation = getdatesEndOfVacation (parameter.startDate);
      persistent datesLiftingTravelRestictionsEU = getDatesLiftingTravelRestictionsEU (parameter.startDate);
      persistent datesLockdownLight = getDatesLockdownLight (parameter.startDate);
      persistent datesMasks = getDatesMasks (parameter.startDate);
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
    
    %check if parameter file contains easings of restrictions or not
    if isfield(parameter,"easingContactRestrictions")%~length(strfind(parameter{},"easingContactRestrictions")) == 0 
      factorsBeta(and(datesEasingContactRestrictions(:,1) < t,...
      t < datesEasingContactRestrictions(:,2))) *= parameter.easingContactRestrictions;
    end 
    if isfield(parameter,"endOfVacation")
      factorsBeta(and(datesEndOfVacation(:,1) < t,...
      t < datesEndOfVacation(:,2))) *= parameter.endOfVacation; 
    end
    if isfield(parameter,"liftingTravelRestictionsEU")
      factorsBeta(and(datesLiftingTravelRestictionsEU(:,1) < t,...
      t < datesLiftingTravelRestictionsEU(:,2))) *= parameter.liftingTravelRestictionsEU; 
    end
    if isfield(parameter,"lockdownLight")
      factorsBeta(and(datesLockdownLight(:,1) < t,...
      t < datesLockdownLight(:,2))) *= parameter.lockdownLight; 
    end
    if isfield(parameter,"masks")
      factorsBeta(and(datesMasks(:,1) < t,...
      t < datesMasks(:,2))) *= parameter.masks; 
    end
	endfunction
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function showFactorsBeta (parameter)
		datesSchoolClosing = getDatesSchoolClosing (parameter);
		datesContactRestrictions = getDatesContactRestrictions (parameter);
		datesExitRestrictions = getDatesExitRestrictions (parameter);
		datesMajorEventsRestrictions = getDatesMajorEventsRestrictions (parameter);
		datesEasingContactRestrictions = getDatesEasingContactRestrictions (parameter);
    datesEndOfVacation = getdatesEndOfVacation (parameter);
    datesLiftingTravelRestictionsEU = getDatesLiftingTravelRestictionsEU (parameter);
    datesLockdownLight = getDatesLockdownLight (parameter); 
    datesMasks = getDatesMasks (parameter);
    
		t= -10+min([datesExitRestrictions(!isinf(datesExitRestrictions));...
		datesContactRestrictions(!isinf(datesContactRestrictions));...
		datesSchoolClosing(!isinf(datesSchoolClosing));...
		datesMajorEventsRestrictions(!isinf(datesMajorEventsRestrictions));...
    datesEasingContactRestrictions(!isinf(datesEasingContactRestrictions));...
    datesEndOfVacation(!isinf(datesEndOfVacation));...
    datesLiftingTravelRestictionsEU(!isinf(datesLiftingTravelRestictionsEU));...
    datesLockdownLight(!isinf(datesLockdownLight));...
    datesMasks(!isinf(datesMasks))]):...  
		max([datesExitRestrictions(!isinf(datesExitRestrictions));...
		datesContactRestrictions(!isinf(datesContactRestrictions));...
		datesSchoolClosing(!isinf(datesSchoolClosing));...
		datesMajorEventsRestrictions(!isinf(datesMajorEventsRestrictions));...
    datesEasingContactRestrictions(!isinf(datesEasingContactRestrictions));...
    datesEndOfVacation(!isinf(datesEndOfVacation));...
    datesLiftingTravelRestictionsEU(!isinf(datesLiftingTravelRestictionsEU));...
    datesLockdownLight(!isinf(datesLockdownLight));...
    datesMasks(!isinf(datesMasks))])+10;
		
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
    
 		factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesEasingContactRestrictions(:,1) < t,...
		t < datesEasingContactRestrictions(:,2))) *= parameter.easingContactRestrictions; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
    
    factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesEndOfVacation(:,1) < t,...
		t < datesEndOfVacation(:,2))) *= parameter.startOfVacation; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
		
    factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesLiftingTravelRestictionsEU(:,1) < t,...
		t < datesLiftingTravelRestictionsEU(:,2))) *= parameter.liftingTravelRestictionsEU; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
    
    factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesLockdownLight(:,1) < t,...
		t < datesLockdownLight(:,2))) *= parameter.lockdownLight; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
    
    factorsBeta = ones(size(datesSchoolClosing, 1), length(t));
		factorsBeta(and(datesMasks(:,1) < t,...
		t < datesMasks(:,2))) *= parameter.masks; 
		factorsBeta = population * factorsBeta(1:end-1,:) / sum(population);
		plot(t + parameter.startDate, factorsBeta);
    
		legend("Schools closing", "Contact restrictions", "Exit restrictions",...
		"Major events restrictions", "Easing contact restrictions", "End of Vacation",...
    "Lifting Travel Restictions EU", "Lockdown Light", "location", "northeastoutside");
		datetick ("x", "dd.mm.yyyy", "keeplimits");
	endfunction
	
	function darkFigureCityWise = getDarkFigures(cities, parameter, statewise)			
		% anhand vergleich tote (25.4) zu infizierten 3 Wochen vorher, checken todo
    %alte Faktoren berechnet aus infectionnumbers und deadnumbers
		factorsDfCorrection = [ 1.06263   0.92143   1.13461   1.39972   1.01612   1.18285   0.68441   1.15415   1.13741   1.23239   0.57678 1.11515   0.64394   0.91610   0.71415   1.10816 1];
		%neu berechnete Faktoren für neue RKI Daten aus parseRKI
    %factorsDfCorrection = [1.035017091	1.008622058	1.132102404	1.135272761	0.983703052	1.219837459	0.729016264	1.072816675	1.150455454	1.249227021	0.605305814	1.242551543	0.541320639	0.826518817	0.731429415	1.336803532 1];
    %factorsDfCorrection = [1.07864   1.12911   1.09597   0.77945   0.89199   1.02703   0.75310   1.17178   1.21760   1.47352   0.58356   1.14254   0.54313   0.97425   0.77320   1.36513   1];
    %factorsDfCorrection = [1.04291   0.97987   1.13102   1.12368   0.98082   1.22953   0.72034   1.07941   1.14565   1.28686   0.59521   1.29059   0.51938    0.82505   0.71892   1.33076 1];
    if isfield(parameter,"factorsDF")%~length(strfind(parameter{},"easingContactRestrictions")) == 0 
      factorsDfCorrection = [parameter.factorsDF, 1];
    end 
    factorsDfCorrection /= mean(factorsDfCorrection);
		
		if nargin == 3
			if strcmp(statewise, "statewise")
				darkFigureCityWise = parameter.darkFigure * factorsDfCorrection;
				return;
			end
		end
		
		darkFigureCityWise = ones(length(cities), 1) * parameter.darkFigure;	
		
		##	factors = [1.012487,0.886786,1.074653,1.036395,0.972342,...
		##	1.078477,0.660660,1.080391,1.070779,1.207535,0.564052,...
		##	0.999975,0.666807,0.885511,0.649587,0.946921, 1];
		##	
		for i=1:length(cities)
			darkFigureCityWise(i) *= factorsDfCorrection(floor(cities{i}.ags/1000));
		endfor  	
	endfunction
	
	function deltaCityWise = getDelta(cities, parameter, statewise)			

		factorsDfCorrection = [ 1.06263   0.92143   1.13461   1.39972   1.01612   1.18285   0.68441   1.15415   1.13741   1.23239   0.57678 1.11515   0.64394   0.91610   0.71415   1.10816 1];

    if isfield(parameter,"factorsDelta")%~length(strfind(parameter{},"easingContactRestrictions")) == 0 
      factorsDfCorrection = [parameter.factorsDelta, 1];
    end 
    factorsDfCorrection /= mean(factorsDfCorrection);
		
		if nargin == 3
			if strcmp(statewise, "statewise")
				deltaCityWise = parameter.darkFigure * factorsDfCorrection;
				return;
			end
		end
		
		deltaCityWise = ones(length(cities), 1) * parameter.delta;	
		
		for i=1:length(cities)
			deltaCityWise(i) *= factorsDfCorrection(floor(cities{i}.ags/1000));
		endfor  	
	endfunction
	
  
function mortalityCityWise = getMortalityCourse (parameter, cities, t)
  t+=parameter.startDate;
  
  t_start = datenum(2020,03,02);
  t_end = datenum(2021,03,02);
  runTime = t_end - t_start;
  time = linspace(t_start, t_end, runTime);
  
  %Data from Excel  "Altersverteilung_Neuinfektionen" with shift by two weeks
  %Factors are calculated for calendar weeks, values in between are interpolated
  mortalityCW = [0.563860495,0.563860495,0.563860495, 0.59470227,0.720584392,1.019749573,1.35690959,1.522336374,...
	1.555816934, 1.51985618,	1.366810599,	1.230078643,	1.15428665,	1.042598338,	...
  0.912176096,	0.789443471,	0.727658453,	0.501364919,	0.519042372,	0.530875906,...
	0.517781852,	0.500192626,	0.473300394,	0.516665509,	0.418170276,	0.356420909,...
	0.296906246,	0.315596437,	0.33182108,	0.468621497,	0.518510454,	0.54216471,	...
  0.590009594,	0.600851564,	0.655646015,	0.712544848,	0.741409155, 0.77440109, ...
	0.839890139];
  timeCW = datenum(2020,03,08):7:datenum(2021,03,02);
  splineY = [mortalityCW, ones(1,round(length(time)/7)-length(mortalityCW))*mortalityCW(end)];
  splineDates = timeCW;
  splineMort = spline(splineDates, splineY);
  splineFinal = ppval(splineMort, time);
  
  mortalityCityWise = ones(length(cities), 1);
  %cosider age structure of county, too
  for i = 1:length(cities)  
    if parameter.considerAverageAgeCounty
      mortalityCityWise(i) = (4*ppval(splineMort, t)+(6)*cities{i}.devAgeAverage)/10;
    else
      mortalityCityWise(i) = (4*ppval(splineMort, t)+6)/10;
    end
  end
endfunction
	
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
		
		
		##		if length(beta_cross_county) == 1			
		##			beta_cross_county = ones(size(beta)) * beta_cross_county;
		##			if !parameter.reduceToStates		
		##				factorsBetaCC = [1.493720 0.431189 1.445615 0.516137 0.344211 0.824671 1.175521 0.537480 0.506872 0.472437 0.515953 1.664755 1.191301 0.681262 0.948501 1.114234];
		##				%factorsBetaCC /= mean(factorsBetaCC);
		##				for i=1:length(cities)
		##					beta_cross_county(i) *= factorsBetaCC(floor(cities{i}.ags/1000));
		##				end
		##			end		
		##		elseif length(beta_cross_county) == 16
		##			beta_cross_county(end+1) = 1;
		##			tmp = zeros(length(cities), 1);
		##			for i=1:length(cities)
		##				tmp(i) = beta_cross_county(floor(cities{i}.ags/1000));
		##			end
		##			beta_cross_county = tmp;
		##		else
		##			error("Wrong length of betacc");
		##		endif
		
	endfunction
	
	
	
	
	
