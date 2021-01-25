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
%% 3cities|small|bayern|germany
function cities = setup_system (parameter)
  if nargin == 0
    parameter.spatial = "small";
    parameter.initial = "1000Munich";
    parameter.model = "SIR";
  endif
  
  % read cities opportunities
  % 3cities|small|bayern| (otherwise it is germany)
  t1 = tic;
  % print status
  if parameter.fullConsoleOut
    fprintf("Read in of population ...");
  end
  spatialDiscretization = strsplit(parameter.spatial,"&");
  citiesTmp = {};
  for i=1:length(spatialDiscretization)
    switch spatialDiscretization{i}      
      case "3cities"
        filename = '../../Daten/3cities.txt';
      case "small"
        filename = '../../Daten/population_bayern_small.txt';
      case "bayern"
        filename = '../../Daten/population_bayern.txt';
      case "germany"
        filename = '../../Daten/population.txt';
      case "europeanCapitals"
        filename = '../../Daten/europeanCapitals.txt';
    endswitch
    citiesTmp = [citiesTmp, read_population(filename, parameter)];
  end
  ##  % unique cities according to ags
  ##  names = cell(1, length(citiesTmp));
  ##  for i=1:length(citiesTmp)
  ##    names{i} = [num2str(citiesTmp{i}.area), num2str(citiesTmp{i}.ags),...
  ##    num2str(citiesTmp{i}.x), num2str(citiesTmp{i}.y)];
  ##  end
  % unique cities according to ags
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
  % update status
  if parameter.fullConsoleOut
    fprintf("\rRead in of population took %.2fs\n", toc(t1));
  end
  
  % count citizens
  N = 0;
  for i=1:length(cities)
    N+=cities{i}.population;
  end  
  t2 = tic;
  
  % set all cities susceptible as staring point
  if and(strcmp(parameter.initial,"RKIfiles"), parameter.startDate  > 737908)
    for i=1:length(cities)
      cities{i}.SIR = [0,0,0,0,0];
    end
  else
    for i=1:length(cities)
      cities{i}.SIR = [1,0,0,0,0];
    end
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
        % maybe 100% infected in Heinsberg is too much, set to 0.1
        for i=1:length(cities)
          if length(strfind(cities{i}.name, "Heinsberg, Stadt"))
            %if cities{i}.ags == 5370
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
        cities = getGitData(cities, parameter, 4);
      case "RKIfiles"
        cities = getRKIData(cities, parameter, 4);
      case "RKIdownload"
        error("Usage of RKIdownload (initial-)option. This option is depricated and was never really tested.");
        output_filename = "../../Results/RKI-data.csv";
        % output_filename = "../../Daten/cases-rki-by-ags_current.csv";  % This is the file updated by the parseRKIData.py Skript running on the Kerberos Server
        % only download and parse new data if it wasn't done already
        if ~isfile(output_filename)
          % download new data from RKI (same as for the RKI-dashboard) and parse it
          % from: https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
          t = tic;
          % print status
          if parameter.fullConsoleOut
            fprintf("downloading new data from RKI...");
          endif
          
          % download csv
          % this seems to be incomplete on 30.07. while it was still ok on 28.07. 
          % csv can also be downloaded here: https://www.arcgis.com/home/item.html?id=f10774f1c63e40168479a1feb6c7ca74
          % -> https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data
          % raw_filename = '../../Results/RKI-data-raw.csv';
          % urlwrite("https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv",raw_filename);
          % s = fileread(raw_filename);
          s = urlread("https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv");
          
          % data cols
          % 1:ObjectId 2:IdBundesland 3:Bundesland 4:Landkreis 5:Altersgruppe 6:Geschlecht 7:AnzahlFall
          % 8:AnzahlTodesfall 9:Meldedatum 10:IdLandkreis 11:Datenstand 12:NeuerFall 13:NeuerTodesfall
          % 14:Refdatum 15:NeuGenesen 16:AnzahlGenesen 17:IstErkrankungsbeginn 18:Altersgruppe2
          data = strsplit(s, '\n');
          
          % update status
          if parameter.fullConsoleOut
            fprintf("\rdownloading new data took %.2fs (%i entries)\n", toc(t), length(data));
          end
          
          % parse the data
          [cases, lastDayFound] = parseRKIData(data, parameter.fullConsoleOut);
          
          % write data to file
          csvwrite(output_filename, cases);
          % dlmwrite(output_filename, ";", cases);
          
        else
          % an existing data file was found; delete it if new data should be downloaded
          if parameter.fullConsoleOut
            fprintf('existing data file found (%s)\n', output_filename);
          end
          
          % get the last day from the found data
          lastDayFound = csvread(output_filename)(1,end);
        end
        
        % use the algorithm for GitData for further processing, as this new data has the same format now
        shiftDays = 7;  # the numbers of days after which a which a second measurement is done to fit the exponentials (If I understand correctly)
        if parameter.fullConsoleOut
          fprintf("setting startDate and initalDistributionDate to %s (totalRuntime=%d)\n", datestr(lastDayFound-shiftDays), parameter.totalRuntime)
        end
        parameter.startDate = lastDayFound - shiftDays;
        parameter.initalDistributionDate = lastDayFound - shiftDays;
        parameter.endDateOpt = parameter.startDate + parameter.totalRuntime;
        cities = getGitData(cities, parameter, shiftDays);
        
    otherwise
      error('unknown testcase');
  endswitch
end

switch parameter.model
  case {"SIR", "SIRH"}
    for i=1:length(cities)
      cities{i}.SIR(5) =  [];
      cities{i}.SIR(4) =  [];
    end
  case {"SIREDLiterature"}
    %change values for E group, since only Individuals from I group can infect others
    for i=1:length(cities)
      %cities{i}.SIR(4) =  0;
    end
  case "SIREDmod"
    if length(strfind(parameter.initial,"GitData")) == 0
      for i=1:length(cities)
        % change inital values of E and I due to the behaviour of the SIREDmod model
        % where only the E group is infects further
        if (cities{i}.SIR(4) == 0)
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
  case "SIREDYO"
    ageStruc = read_agestructure();
    ageVector = zeros(length(ageStruc), 1);
    for i = 1:size(ageStruc,2)
      ageVector(i) = ageStruc{1,i}.ags;
    end
    for j=1:length(cities)
      index = find(ageVector == cities{j}.ags);
      if index
        fracOld = ageStruc{index}.fracOld;
        cities{j}.fracOld = fracOld;
        Syng = cities{j}.SIR(1) * (1-fracOld);
        Iyng = cities{j}.SIR(2) * (1-fracOld);
        Ryng = cities{j}.SIR(3) * (1-fracOld);
        Eyng = cities{j}.SIR(4) * (1-fracOld);
        Dyng = cities{j}.SIR(5) * (1-fracOld);          
        Sold = cities{j}.SIR(1) * fracOld;
        Iold = cities{j}.SIR(2) * fracOld;
        Rold = cities{j}.SIR(3) * fracOld;
        Eold = cities{j}.SIR(4) * fracOld;
        Dold = cities{j}.SIR(5) * fracOld;
        
        cities{j}.SIR = [Syng, Iyng, Ryng, Eyng, Dyng, Sold, Iold, Rold, Eold, Dold];
      else
        fprintf("ags not found: %s\n", cities{j}.name);
        error("AGS not found.");
      endif
    end
endswitch

if parameter.fullConsoleOut
  fprintf("Setting of initial values took %.2fs\n", toc(t2));
end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function city = getIschgl ()
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cities = getGitData(cities, parameter, shiftDays)
% read the data for the relevant dates
%GitData is not used anymore
error("Git data is no longer used! Please use RKIfiles as start condition");
agsstruc = read_history_data(parameter.initalDistributionDate, 0, nan, parameter);
agsstrucSimStart = read_history_data(parameter.startDate, 0, nan, parameter);
agsstrucSimStartP1 = read_history_data(parameter.startDate, 0, nan, parameter);
agsstrucSimStartP2 = read_history_data(parameter.startDate + shiftDays, 0, nan, parameter);

% initialize variables
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
##        beta = sir_eqn_spatial ("getBeta", parameter, 0, cities);
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
    switch parameter.model
      case {"SIRED", "SIREDmod", "SIREDYO","SIREDLiterature"}
        cities{j}.SIR = [s0, i0, r0, e0, d0];
      case {"SIR", "SIRH"}
        % factor = 65*3.3;
        factor = 1;
        cities{j}.SIR = [1-i0*factor, i0*factor, r0, 0, 0];
        cities{j}.baseQ = base_Q;
    endswitch
  else
    fprintf("ags not found: %i (%s)\n", cities{j}.ags, cities{j}.name);
    error("AGS not found.");
  end
endfor
endfunction
function cities = getRKIData(cities, parameter, shiftDays)
% read RKIdata
filename_infected = '../../Daten/cases-rki-by-ags_current.csv';
filename_dead = '../../Daten/deaths-rki-by-ags_current.csv';
filename_ageStructure = '../../Daten/Deviation-from-average-age_AGS.csv';
delimiterIn = ',';
infected = importdata(filename_infected,delimiterIn);
dead = importdata(filename_dead,delimiterIn);
ageStruc = importdata(filename_ageStructure,delimiterIn);
lags = sir_eqn_spatial("lags");

%save data status
parameter.dataStatus = datestr(infected(1,end));
cities{1}.dataStatus = datestr(infected(1,end));

datavec = [];
datavec.ags = zeros(size(infected,1)-1,1);
for i = 2:size(infected,1)
  datavec.ags(i-1) = infected(i,1);
end
%find dates from parameter file
indexInitialDD = find(infected(1,:) == parameter.initalDistributionDate);
indexStart = find(infected(1,:) == parameter.startDate);
indexStartP1 = find(infected(1,:) == parameter.startDate);
indexStartP2 = find(infected(1,:) == parameter.startDate + shiftDays);

datavec.infInitialDD = zeros(size(infected,1)-1,1);
datavec.infStart = zeros(size(infected,1)-1,1);
datavec.infStartP1 = zeros(size(infected,1)-1,1);
datavec.infStartP2 = zeros(size(infected,1)-1,1);
datavec.recovered = zeros(size(infected,1)-1,1);

for i = 2:size(infected,1)
  datavec.infInitialDD(i-1) = infected(i,indexInitialDD);
  datavec.infExposedDD(i-1) = infected(i,indexInitialDD+3);
  %All simulations that start before 16-Mar-2020 take the initial distribution date into account
  if parameter.startDate > 737866
    %calculate people that are currently infected (recovery takes 14 days) and recovered
    %On average, it takes six days until an infected person is detected, after about 
    %2 days the person is infectious, i.e. the person was infectious four days before
    %he or she was detected
    %Idea: take new infections of the last four days, all others should be in quarantine by now
    %datavec.infStart(i-1) = infected(i,indexStart) - infected(i, indexStart-14);
    datavec.infStart(i-1) = infected(i,indexStart) - infected(i, indexStart-14);
    datavec.infExposedStart(i-1) = infected(i, indexStart+3) - infected(i, indexStart);;
    %we need all infected from the last 14 days so diff = 0
    datavec.diff(i-1) = infected(i,indexStart-14) - infected(i, indexStart-14);
    datavec.recovered(i-1) = infected(i, indexStart-14);
    datavec.dead(i-1) = dead(i,indexStart-14);
    datavec.deadStart(i-1) = dead(i,indexStart);
    %All simulations that start after 31-Mar-2020 can use the real history as start condition
    if parameter.startDate > 737881
      %datavec.Infhistory(i-1,1:29) = infected(i, indexStart-28:indexStart);
      datavec.Infhistory(i-1,4:26) = movmean(infected(i, indexStart-28:indexStart), [3 3] ,'Endpoints', 'discard');
      datavec.Infhistory(i-1,3) = movmean(infected(i, indexStart-28:indexStart-24), [2 2] ,'Endpoints', 'discard');
      datavec.Infhistory(i-1,27) = movmean(infected(i, indexStart-4:indexStart), [2 2] ,'Endpoints', 'discard');
      datavec.Infhistory(i-1,2) = movmean(infected(i, indexStart-28:indexStart-26), [1 1] ,'Endpoints', 'discard');
      datavec.Infhistory(i-1,28) = movmean(infected(i, indexStart-2:indexStart), [1 1] ,'Endpoints', 'discard');
      datavec.Infhistory(i-1,1) = infected(i, indexStart-28);
      datavec.Infhistory(i-1,29) = infected(i, indexStart);
    end
  else
    datavec.infStart(i-1) = infected(i,indexStart);
  end
  datavec.infStartP1(i-1) = infected(i,indexStartP1);
  datavec.infStartP2(i-1) = infected(i,indexStartP2);
end
sumIsimStartState = zeros(16,1);
sumIDistributionState = zeros(16,1);
sumIsimStart = 0;
sumIsimStartP1 = 0;
sumIsimStartP2 = 0;
sumIDistribution = 0;

ags = 1:16;
agsVector = zeros(size(datavec.ags,1), 1);
agsVector(:) = datavec.ags(:,1);
reduction = cell(1, length(ags));
for i=1:length(ags)
  reduction{i} = find(floor(agsVector/1e3) == ags(i));
end 
for i = 1:size(datavec.ags,1)
  sumIsimStart += datavec.infStart(i,1);
  sumIsimStartP1 += datavec.infStartP1(i,1);
  sumIsimStartP2 += datavec.infStartP2(i,1);
  sumIDistribution += datavec.infInitialDD(i,1);
end  
%Define reduction factor by state
%With a Germany-wide definition, there would be major deviations in the starting conditions
for i=1:length(ags)
  sumIsimStartState(i) = sum(datavec.infStart(reduction{i}));
  %Avoid division by zero
  if sumIsimStartState(i)==0
    sumIDistributionState(i) = 1;
  end
  sumIDistributionState(i) = sum(datavec.infInitialDD(reduction{i}));
end

agsVectorCities = zeros(length(cities), 1);
population = zeros(length(cities), 1);        
for i=1:length(cities)
  agsVectorCities(i) = cities{i}.ags;
  population(i) = cities{i}.population;
end
darkFigure = sir_eqn_spatial ("getDarkFigure", parameter, cities);
base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
index2 = 0;
for j=1:length(cities)
  index = find(agsVector == cities{j}.ags);
  stateInd = find(ags == floor(cities{j}.ags/1e3)); 
  cities{j}.history = zeros(40,1);
  cities{j}.deathStart = 0;
  cities{j}.devAgeAverage = ageStruc(index,2);
##  base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
  if index ~=index2
    n_city = 0;
    indags = find(agsVectorCities == cities{j}.ags);
    for i = indags(1):indags(end)
      n_city += cities{i}.population;
    end  
  end
  %index2 = index;
  %All simulations that start before 16-Mar-2020 take the initial distribution date into account
  if (parameter.startDate<737866)
    reduction_factor =  sumIsimStartState(stateInd) / sumIDistributionState(stateInd);            
    r0 = 0;
    d0 = 0;
    i0 = datavec.infInitialDD(index)/n_city*reduction_factor;
    e0l = datavec.infExposedDD(index)/n_city*reduction_factor;
    e0 = i0*log(base_Q)*((darkFigure(j) - 1) / parameter.gamma1)*exp(1);  
    s0 = 1-i0 -d0 -e0-r0;
    switch parameter.model
      case {"SIREDLiterature"}
        cities{j}.SIR = [1-i0-e0l, i0, r0, e0l, d0];
      case {"SIRED", "SIREDmod", "SIREDYO"}
        cities{j}.SIR = [s0, i0, r0, e0, d0];
      case {"SIR", "SIRH"}
        % factor = 65*3.3;
        factor = 1;
        cities{j}.SIR = [1-i0*factor, i0*factor, r0, 0, 0];
        cities{j}.baseQ = base_Q;
      %Todo SIR Startbedingungen
      %Todo SIREDLiterature
    endswitch
  else 
    %Avoid that values are unnecessarily recalculated
    if index ~=index2
      diff = (datavec.diff(index)*darkFigure(j))/n_city;
      r0 = (datavec.recovered(index)*darkFigure(j)-datavec.dead(index))/n_city;
      d0 = datavec.deadStart(index)/n_city;
      i0 = datavec.infStart(index)/n_city;%*reduction_factor;  
      e0l = datavec.infExposedStart(index)/n_city;
      e0 = i0*darkFigure(j);  
      s0 = 1-i0-d0-e0-r0-diff-(datavec.dead(index)/n_city);
    end
    switch parameter.model
    case {"SIREDLiterature"}
        cities{j}.SIR = [1-i0-e0l-r0-d0-(datavec.dead(index)/n_city), i0, r0, e0l, d0];
    case {"SIRED", "SIREDmod", "SIREDYO"}
      cities{j}.SIR = [s0, i0, r0, e0, d0];
    case {"SIR"}
      cities{j}.SIR = [1-i0-r0, i0, r0, 0, 0];%
    case {"SIRH"}
      %cities{j}.SIR = [1-i0*parameter.darkFigure, i0, r0, 0, 0];
      %[1-i0-d0-r0, i0, r0, 0, 0];
      %cities{j}.baseQ = base_Q;
      cities{j}.SIR = [s0, i0, r0, 0, 0];
      cities{j}.deathStart = d0;
      if (parameter.startDate>737881)
        %sum up population of all cities in county otherwise history would
        %be negative in case of small cities in big counties     
        if index ~=index2       
          i0h = datavec.Infhistory(index,:)/n_city;
          history = i0h*darkFigure(j);%*factor;%*factor;%+i0h*parameter.darkFigure; % +r0h;%+d0h
          splineH = pchip(1:29, history, 29 - linspace(lags(1), lags(end), 40));
        end
        cities{j}.history = splineH;
        %cities{j}.SIR = [1-i0*darkFigure(j), i0, r0, 0, 0];
        %cities{j}.SIR = [s0, i0, r0, 0, 0];
        %cities{j}.deathStart = d0;
      end
  endswitch
  
endif
##else
##fprintf("ags not found: %i (%s)\n", cities{j}.ags, cities{j}.name);
##error("AGS not found.");
##end
index2 = index;
endfor

##for j=1:length(cities)
##  index = find(agsVector == cities{j}.ags);
##  cities{j}.history = zeros(40,1);
##  cities{j}.deathStart = 0;
##  cities{j}.devAgeAverage = ageStruc(index,2);
##  if  max(population ( find(agsVectorCities == cities{j}.ags))) >...
##    cities{j}.population
##    % Only the biggest city in each county gets the infected people
##    continue;
##  end
##  if index
##    %if startDate is in early March and there are not enough infected in some 
##    %couties to do realistic simulation with inital values    
##    % e funktions fit durch beide Punkte -> basis 
##    % Q = Q0 * base_Q^t
##    % dQ_dt = Q0 * ln(base_Q) * base_Q^(t)
##    % Q(2.3) = Q0 * base_Q^0 -> Q0 = Q(2.3)
##    % Q(6.3) = Q0 * base_Q^4 -> base_Q = (Q(6.3) / W(2.3)) ^(1/4)
##    %cities{j}.devAgeAverage = ageStruc(index,2);
##    base_Q = (sumIsimStartP2/sumIsimStartP1)^(1/shiftDays);
##    n_city = 0;
##    indags = find(agsVectorCities == cities{j}.ags);
##    for i = indags(1):indags(end)
##      n_city += cities{i}.population;
##    end  
##    if (parameter.startDate<737866)
##      reduction_factor =  sumIsimStart / sumIDistribution;            
##      r0 = 0;
##      d0 = 0;
##      i0 = datavec.infInitialDD(index)/n_city*reduction_factor;
##      e0 = i0*log(base_Q)*((darkFigure(j) - 1) / parameter.gamma1)*exp(1);  
##      s0 = 1-i0 -d0 -e0-r0;
##      switch parameter.model
##        case {"SIRED", "SIREDmod", "SIREDYO","SIREDLiterature"}
##          cities{j}.SIR = [s0, i0, r0, e0, d0];
##        case {"SIR", "SIRH"}
##          % factor = 65*3.3;
##          factor = 1;
##          cities{j}.SIR = [1-i0*factor, i0*factor, r0, 0, 0];
##          cities{j}.baseQ = base_Q;
##      endswitch
##    else      
##      r0 = (datavec.recovered(index)*darkFigure(j)-datavec.dead(index))/n_city;
##      d0 = datavec.deadStart(index)/n_city;
##      i0 = datavec.infStart(index)/n_city;%*reduction_factor;  
##      e0 = i0*darkFigure(j);  
##      s0 = 1-i0 -d0-e0-r0;
##      switch parameter.model
##        %Todo SIREDLiterature
##      case {"SIRED", "SIREDmod", "SIREDYO","SIREDLiterature"}
##        cities{j}.SIR = [s0, i0, r0, e0, d0];
##      case {"SIR"}
##        cities{j}.SIR = [1-i0-r0, i0, r0, 0, 0];%
##      case {"SIRH"}
##        %cities{j}.SIR = [1-i0*parameter.darkFigure, i0, r0, 0, 0];
##        %[1-i0-d0-r0, i0, r0, 0, 0];
##        cities{j}.baseQ = base_Q;
##        if (parameter.startDate>737908)
##          %sum up population of all cities in county otherwise history would
##          %be negative in case of small cities in big counties            
##          i0h = datavec.Infhistory(index,:)/n_city;
##          history = i0h*darkFigure(j);%*factor;%*factor;%+i0h*parameter.darkFigure; % +r0h;%+d0h
##          splineH = pchip(1:29, history, 29 - linspace(lags(1), lags(end), 40));
##          cities{j}.history = splineH;
##          %cities{j}.SIR = [1-i0*parameter.darkFigure, i0, r0, 0, 0];
##          cities{j}.SIR = [1-i0*darkFigure(j), i0, r0, 0, 0];
##          cities{j}.deathStart = d0;
##        else
##          cities{j}.SIR = [1-i0*darkFigure(j), i0, r0, 0, 0];
##          cities{j}.deathStart = d0;
##        end
##    endswitch
##  endif
##else
##  fprintf("ags not found: %i (%s)\n", cities{j}.ags, cities{j}.name);
##  error("AGS not found.");
##end
##endfor
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cases, lastDayFound] = parseRKIData(data, fullConsoleOut)
% parse data
t = tic;
if fullConsoleOut
fprintf("parsing new data...");
endif

% init date variables
day0 = datenum(2020,01,27);  % the day before it all started in germany
firstDayFound = 1000000;
lastDayFound = -1;

% init cases matrix
cases = zeros(402,floor(now())-day0+1);  # +1 just in case
cases(1,2:end) = (day0):(day0+size(cases)(2)-2);  % set dates in the first line
currLandkreis = 0;  % start at 0 because data is sorted by Landkreis-id (from low to high)
currIndex = 1;

% go through all data-lines
for i=2:length(data)-1  % ignore the first (names) and the last (empty) line
% parse line
entries = strsplit(data{i}, ",");
if length(entries) < 19
  % error: incomplete line
  fprintf("\nincomplete line %i\n", i);
  fprintf("%i: `%s`\n", i-1, data{i-1});
  fprintf("%i: `%s`\n", i  , data{i  });
  error('incomplete line');
end
value = str2num(entries{7});
landkreis = str2num(entries{10});
day = datenum(entries{9}(1:10)) - day0;  % days since the first day
firstDayFound = min(firstDayFound, day);
lastDayFound = max(lastDayFound, day);
% switch to a new landkreis/city
if landkreis != currLandkreis
  % Berlin is splitted in RKI data but not in population data
  % so set all parts of Berlin equal to Berlin itself
  if landkreis >=11000 && landkreis <=11012  % `read_data_rki.m` uses 11011 as upper limit, but `InfectedDataHistory.txt` seems to also include 11012
    # inside berlin
    landkreis = 11000;
    if currLandkreis != landkreis
      # first time in Berlin
      currIndex += 1;
      cases(currIndex,1) = landkreis;
      currLandkreis = landkreis;
    endif
  else
    % normal action for a new landkreis
    currIndex += 1;  % this only works because the data is sorted by Landkreis
    cases(currIndex,1) = landkreis;  % write the Landkreis-id to the first column
    currLandkreis = landkreis;
    % update status
    if fullConsoleOut
      fprintf("\rparsing new data (Lk %i/402 - %s) ...                    ", currIndex, entries{4});
    endif
  end
end
% add the current case
cases(currIndex,day+1) += value;
end

% check if the dates are ok
if firstDayFound < 0 || firstDayFound >= 1000000
fprintf('\ninvalid first date: %i\n', firstDayFound);
error('invalid first date');
end

% update status
if fullConsoleOut
fprintf("\rparsing new data (sum up values) ...                         ");
end

% remove calumns full of 0s at the end
cases = cases(:,1:lastDayFound-firstDayFound+2);

% sum values, so that all cases until date i are saved (to match `InfectedDataHistory.txt`)
cases(2:end,2:end) = cumsum(cases(2:end,2:end), 2);

# cut off 32 days at the beginning to have the same start-day as `InfectedDataHistory.txt`
cases_old = cases;
cases = zeros(402, size(cases_old)(2)-32);
cases(:,1) = cases_old(:,1);
cases(:,2:end) = cases_old(:,34:end);

% update status
if fullConsoleOut
fprintf("\rparsing new data took %.2fs                                  \n", toc(t));
fprintf("found data from %s until %s\n", datestr(firstDayFound+day0), datestr(lastDayFound+day0));
end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cases, lastDayFound] = parseRKIData2(data, fullConsoleOut)
% this is WIP
% try to avoid for-loops

% parse data
t = tic;
if fullConsoleOut
fprintf("parsing new data...");
endif

% init date variables
day0 = datenum(2020,01,27);  % the day before it all started in germany
firstDayFound = 1000000;
lastDayFound = -1;

% convert string data to integers
%1 2 3 4 5 6      7 8    9         10                                                                 0  1  2  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
% [_,_,_,_,_,_,values,_,days,landkreise,_,_,_,_,_,_,_,_] = textscan('../../Results/RKI-data-raw.csv', "%u,%u,%s,%s,%s,%s,%d,%d,%s,%u,%s,%d,%d,%s,%d,%d,%d,%s");
fprintf("\rparsing new data: strsplit     ");  % debug
data2 = cellfun (@(x) strsplit(x, ","), data, "uniformoutput", false);
fprintf("\rparsing new data: zeros     ");  % debug
data3 = zeros(3, length(data));  # [values, days, landkreise]
fprintf("\rparsing new data: values     ");  % debug
data3(1) = cellfun (@(x) str2num(x), data2(:,7), "uniformoutput", false);  # values
fprintf("\rparsing new data: days     ");  % debug
data3(2) = cellfun (@(x) datenum(x), data2(:,9,1:10), "uniformoutput", false);  # days
fprintf("\rparsing new data: landkreise     ");  % debug
data3(3) = cellfun (@(x) str2num(x), data2(:,10), "uniformoutput", false);  # landkreise

% init cases matrix
cases = zeros(402,floor(now())-day0+1);  # +1 just in case
cases(1,2:end) = (day0+1):(day0+length(cases(1,:))-1);  % set dates in the first line

% the main part of the algorithm is still missing
% TODO

endfunction