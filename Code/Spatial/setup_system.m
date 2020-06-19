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
## @deftypefn {} {@var{retval} =} setup_system (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-14


%% possible spatial cases are
%% bayern|germany
function cities = setup_system (parameter)
  if nargin == 0
    parameter.spatial = "bayern";
    parameter.initial = "GitData";
    parameter.model = "SIREDmod";
  endif
  
  % read cities opportunities
  t1 = tic;
  spatialDiscretization = strsplit(parameter.spatial,"&");
  citiesTmp = {};
  for i=1:length(spatialDiscretization)
    switch spatialDiscretization{i}      
      case "bayern"
        filename = '../../Daten/population_bayern.txt';
      case "germany"
        filename = '../../Daten/population.txt';
    endswitch
    citiesTmp = [citiesTmp, read_population(filename, parameter)];
  end

  names = zeros(1, length(citiesTmp));
  for i=1:length(citiesTmp)
    names(i) = citiesTmp{i}.area * citiesTmp{i}.ags *...
    citiesTmp{i}.x * citiesTmp{i}.y;
  end
  [~, I, ~] = unique(names);
  I = sort(I);
  cities = cell(1, length(I));
  for i=1:length(cities)
    cities{i} = citiesTmp{I(i)};
  end
  if parameter.fullConsoleOut
    fprintf("Read in of population took %fs\n", toc(t1));
  end
  
  % count citizens
  N = 0;
  for i=1:length(cities)
    N+=cities{i}.population;
  end  
  t2 = tic;
  
  % set all cities susceptible as staring point
  for i=1:length(cities)
    cities{i}.SIR = [1,0,0,0,0];
  end      
  % adjust SIRED values to those asked for in the parameter file
  scenarios = strsplit(parameter.initial,"&");
  for i=1:length(scenarios)  
    switch scenarios{i}
      case "1000Uniform" % uniform
        % set 1000 infected people
        I_0 = 1000;
        R_0 = 1; 
        E_0 = 0;
        if (isfield(parameter,"rho") == 1)
          E_0 = parameter.rho * I_0;
        end
        D_0 = 0;
        S_0 = N - I_0 - R_0 - E_0 - D_0;
        % set S_0, I_0 and R_0 in all cities
        for i=1:length(cities)
          cities{i}.SIR = [S_0/N, I_0/N, R_0/N, E_0/N, D_0/N];
        end           
      case "100Burghausen" % start in Burghausen
        % set 100 infected in Burghausen
        for i=1:length(cities)
          if length(strfind(cities{i}.name, "Burghausen"))
            N_city = cities{i}.population;
            I_0 = 100;
            R_0 = 1;
            E_0 = 0;
            if (isfield(parameter,"rho") == 1)
              E_0 = parameter.rho * I_0;
            end
            D_0 = 0;
            S_0 = N_city - I_0 - R_0 - E_0 - D_0;
            cities{i}.SIR = [S_0/N_city, I_0/N_city, R_0/N_city, E_0/N_city, D_0/N_city];
          end
        end
      case "1000Munich" % start with 1000 infected in Munich
        for i=1:length(cities)
          if or(length(strfind(cities{i}.name, "München, Landeshauptstadt")), ...
            length(strfind(cities{i}.name, "M�nchen, Landeshauptstadt")))
            N_city = cities{i}.population;
            I_0 = 1000;
            R_0 = 1;
            E_0 = 0;
            if (isfield(parameter,"rho") == 1)
              E_0 = parameter.rho * I_0;
            end
            D_0 = 0;
            S_0 = N_city - I_0 - R_0 - E_0 - D_0;
            cities{i}.SIR = [S_0/N_city, I_0/N_city, R_0/N_city, E_0/N_city, D_0/N_city];
          end
        end
        
      case "allMunich" % everybody in munich is infected
        for i=1:length(cities)
          if or(length(strfind(cities{i}.name, "München, Landeshauptstadt")), ...
            length(strfind(cities{i}.name, "M�nchen, Landeshauptstadt")))
            cities{i}.SIR = [0, 1, 0,0,0];
          end
        end
      case "RKIcountyWise" % data according to rki
        agsstruc = read_data_rki();
        agsVector = zeros(length(agsstruc), 1);
        for i = 1:size(agsstruc,2)
          agsVector(i) = agsstruc{1,i}.ags;
        end  
        for j=1:length(cities)
          index = find(agsVector == cities{j}.ags);
          if index
            r0 = 0;
            i0 = agsstruc{index}.Ifrac;
            d0 = agsstruc{index}.Dfrac;
            e0 = 0;
            if (isfield(parameter,"rho") == 1)
              e0 = parameter.rho * i0;
            end
            s0 = 1-i0 -d0 -e0-r0;
            cities{j}.SIR = [s0, i0, r0, e0, d0];
          else
            fprintf("ags not found: %s\n", cities{j}.ags);
            error("AGS not found.");
          end
        endfor
      case "Ischgl" % Ischgl gets added to counties and everybody there is infected
        cities(length(cities)+1) = getIschgl();   
      case "Heinsberg" % start in Heinsberg 
        % maybe 100% infected in Heinsberg is too much, set to 10%
        for i=1:length(cities)
          if length(strfind(cities{i}.name, "Heinsberg, Stadt"))
            cities{i}.SIR = [0.9, 0, 0, 0.1, 0];
          end
        end
      case "Heiligendamm" % possible summer breakout starting in heiligendamm
        % that is located in Bad Doberan    
        for i=1:length(cities)
          if length(strfind(cities{i}.name, "Bad Doberan, Stadt"))
            cities{i}.SIR = [0.5, 0.5, 0, 0, 0];
          end
        end
			case "GitData" % data obtained from git repo https://github.com/jgehrcke/covid-19-germany-gae
				agsstruc = read_history_data(parameter.initalDistributionDate);
        agsstrucSimStart = read_history_data(parameter.startDate);
				shiftDays = 4;
        agsstrucSimStartP1 = read_history_data(parameter.startDate);
        agsstrucSimStartP2 = read_history_data(parameter.startDate + shiftDays);
				sumIsimStart = 0;
				sumIsimStartP1 = 0;
				sumIsimStartP2 = 0;
				sumIDistribution = 0;
        agsVector = zeros(length(agsstruc), 1);
        for i = 1:size(agsstruc,2)
          agsVector(i) = agsstruc{1,i}.ags;
					sumIsimStart += agsstrucSimStart{1,i}.Infected;
					sumIsimStartP1 += agsstrucSimStartP1{1,i}.Infected;
					sumIsimStartP2 += agsstrucSimStartP2{1,i}.Infected;
					sumIDistribution += agsstruc{1,i}.Infected;
				end  
				agsVectorCities = zeros(length(cities), 1);
				population = zeros(length(cities), 1);				
				for i=1:length(cities)
					agsVectorCities(i) = cities{i}.ags;
					population(i) = cities{i}.population;
				end
				##				beta = sir_eqn_spatial ("getBeta", parameter, 0, cities);
				darkFigure = sir_eqn_spatial ("getDarkFigure", parameter, cities);
        for j=1:length(cities)
          index = find(agsVector == cities{j}.ags);
					if  max(population ( find(agsVectorCities == cities{j}.ags))) >...
						cities{j}.population
						% Only the biggest city in each county gets the infected people
						continue;
					end
          if index
						reduction_factor =  sumIsimStart / sumIDistribution / 207*274; 						
						% e funktions fit durch beide Punkte -> basis 
						% Q = Q0 * base_Q^t
						% dQ_dt = Q0 * ln(base_Q) * base_Q^(t)
						% Q(2.3) = Q0 * base_Q^0 -> Q0 = Q(2.3)
						% Q(6.3) = Q0 * base_Q^4 -> base_Q = (Q(6.3) / W(2.3)) ^(1/4)
						base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
						n_city = cities{j}.population;
            r0 = 0;
            d0 = 0;
            i0 = agsstruc{index}.Infected/n_city*reduction_factor;
            e0 = i0*log(base_Q)*((darkFigure(j) - 1) / parameter.gamma1)*exp(1);			
            s0 = 1-i0 -d0 -e0-r0;
             cities{j}.SIR = [s0, i0, r0, e0, d0];
          else
            fprintf("ags not found: %s\n", cities{j}.ags);
            error("AGS not found.");
          end
        endfor
    otherwise
      error('unknown testcase');
  endswitch
end

if length(strfind(parameter.initial,"GitData")) == 0
  for i=1:length(cities)
    % change inital values of E and I due to the behaviour of the SIREDmod model
    % where only the E group is infects further
    if(cities{i}.SIR(4) == 0)
    cities{i}.SIR(4) =  cities{i}.SIR(2);
    cities{i}.SIR(2) = 0;    
  else
    %error('SIREDmod exposed numbers unequal 0 before changing of infected and exposed initial values');  
  end
end  
else 
for i=1:length(cities)
  if cities{i}.ags >=17000
    cities{i}.SIR(4) =  cities{i}.SIR(2);
    cities{i}.SIR(2) = 0;                
  endif
endfor
endif


if parameter.fullConsoleOut
fprintf("Setting of initial values took %fs\n", toc(t2));
end
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function city = getIschgl()
populationIncreaseFactor = 30;

city.name = "Ischgl";
city.area = 103.01;
city.population = 1617*populationIncreaseFactor;
city.biggestCityInCounty = 1617*populationIncreaseFactor;
city.density = 3;
city.x = 780.77;
city.y = 5214.6;
city.ags = 17000;
city.SIR = [0, 0, 0, 1, 0];
city.fracOld = 0.0;

endfunction
