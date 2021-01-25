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
## @deftypefn {} {@var{retval} =} expandSIRH (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-26

function xFull = expandSIRH(x, t, cities, parameter, initialHistory) 
  
  courseOfDisease = sir_eqn_spatial("totalCourseOfDisease", parameter);
  lags = sir_eqn_spatial("lags");
  S = x;
  
  S_history = initialHistory';
  
  exposed = zeros(size(S));
  effectiveInfective = zeros(size(S));
  symptoms = zeros(size(S));
  hospital = zeros(size(S));
  icu = zeros(size(S));
  
  % concatenate time and solution with initialHistory of the solver
  PPS = spline(...
  [linspace(-max(lags), min(t), size(S_history, 1))'; t(2:end)],...
  [S_history; S(2:end,:)]');
		
  %%%%%%%%%%%%%%%  
  %Preparations for the consideration of age-related mortality
  %Attention! not yet finished
  if isfield(parameter,"considereAgeOfNewlyInfected")
      if parameter.considereAgeOfNewlyInfected
       courseOfDisease2 = zeros(length(t)*length(cities), length(courseOfDisease.dying));
       mortalityInitial = parameter.mortality;
       parameter.courseNew = true;
       for i=1:length(t)
        factorMortality = sir_eqn_spatial("mortalityFactor", parameter, cities, i);
        parameter.mortality = mortalityInitial*factorMortality;
        course = sir_eqn_spatial("totalCourseOfDisease", parameter);
        courseOfDisease2(i:i+length(cities)-1,:) = course.dying;
       end
      end 
  end
  %%%%%%%%%%%%%%%  
	
  for j=1:length(t)
    Sdiff = -diff(ppval(PPS, t(j) - lags), 1, 2);
    exposed(j,:) = Sdiff*cumsum(courseOfDisease.recovering + courseOfDisease.dying)';
    effectiveInfective(j, :) = Sdiff*courseOfDisease.infective';
    symptoms(j, :) = Sdiff*courseOfDisease.symptoms';
    hospital(j, :) = Sdiff*courseOfDisease.hospital';
    icu(j, :) = Sdiff*courseOfDisease.icu';
  end
  
  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate>737881)
    % integrate death and discovered
    darkFigures = sir_eqn_spatial ("getDarkFigure", parameter, cities);
    opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
    get_dDeath_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags), 1, 2)*courseOfDisease.dying'];
    dead = zeros(1,length(cities));
    for j = 1:length(cities)
      dead(j) = cities{j}.deathStart/cities{j}.population;
    end

    if isfield(parameter,"considereAgeOfNewlyInfected")
      if parameter.considereAgeOfNewlyInfected
        death = zeros(length(t), length(cities));
        for i=1:length(cities)
          indCity=i;
          get_dDeath_dT2 = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags)(indCity,:), 1, 2)*courseOfDisease2((round(t_eq)+1)*indCity,:)'];
          [t1, death1] = ode23(get_dDeath_dT2, [min(t), max(t)], dead(i), opt);
          death(:,i) = spline(t1, death1, t);
        end
      else
        [t1, death] = ode23(get_dDeath_dT, [min(t), max(t)], dead, opt);
        death = spline(t1, death, t)';
      end
    else
      [t1, death] = ode23(get_dDeath_dT, [min(t), max(t)], dead, opt);
      death = spline(t1, death, t)';
    end
    
##    [t1, death] = ode23(get_dDeath_dT, [min(t), max(t)], dead, opt);
##    
    
    darkFigures = sir_eqn_spatial ("getDarkFigure", parameter, cities);
    get_dDiscovering_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags), 1, 2)*...
    courseOfDisease.discovering']./darkFigures;    
    [t2, discovered] = ode23(get_dDiscovering_dT, [min(t), max(t)], (1-S_history(end,:))./darkFigures', opt);
    discovered = spline(t2, discovered, t)';
  else   
    % integrate death and discovered
    opt = odeset('NormControl', 'on', 'MaxStep', parameter.maxTimeStep);
    get_dDeath_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags), 1, 2)*courseOfDisease.dying'];
    [t1, death] = ode23(get_dDeath_dT, [-max(lags), max(t)], zeros(1, size(S,2)), opt);
    death = spline(t1, death, t)';
    
    darkFigures = sir_eqn_spatial ("getDarkFigure", parameter, cities);
    get_dDiscovering_dT = @(t_eq, x0) [-diff(ppval(PPS, t_eq - lags), 1, 2)*...
    courseOfDisease.discovering']./darkFigures;
    [t2, discovered] = ode23(get_dDiscovering_dT, [-max(lags), max(t)], zeros(1, size(S,2)), opt);
    discovered = spline(t2, discovered, t)';
  end 
  
  
  ##    names = {"Susceptible", "Exposed", "Removed", "Infectious",...
  ##    "has_symptoms", "in_hospital", "needs_intensice_care"};
  xFull = reshape([S; exposed; 1-S-exposed; effectiveInfective;...
  symptoms; hospital; icu; death; discovered],...
  size(S, 1), 9 * size(S,2));
  %SIRH = [S, X, X, effectiveInfective, symptoms, hospital, icu, dead, discovered];
endfunction












