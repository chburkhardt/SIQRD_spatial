## Copyright (C) 2020 iwtm96_user
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
## @deftypefn {} {@var{retval} =} read_case_history (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-04-29

function data = read_case_history_RKIfiles()


  filename_infected = '../../Daten/cases-rki-by-ags_current.csv';
  filename_dead = '../../Daten/deaths-rki-by-ags_current.csv';
  
  % !!!!!!!!!!!!!!!
  % what should be noted here is that the excel sheets shouldn't have blanks
  % representing 0 values before being transverted to .txt files
  % !!!!!!!!!!!!!!!
  %% stateNames sorted in the way as in the input file
  stateNames= {"Baden-Wuerttemberg","Freistaat Bayern","Berlin"...
    "Brandenburg","Freie Hansestadt Bremen","Freie und Hansestadt Hamburg"...
    "Hessen","Mecklenburg-Vorpommern","Niedersachsen",  "Nordrhein-Westfalen"...
    "Rheinland-Pfalz","Saarland","Freistaat Sachsen","Sachsen-Anhalt"...
    "Schleswig-Holstein","Freistaat Thueringen"};
    
  %% read infection and death numbers from txt or csv file  
  delimiterIn = ',';
  infected = importdata(filename_infected,delimiterIn);
  death = importdata(filename_dead,delimiterIn);
  
  infectionsvec = {};
  deathsvec = {};
  timesvec = {};
  dateInNumber = infected(1,5:end);
  for i = 1:length(dateInNumber)
    timesvec(length(timesvec)+1) = dateInNumber(i);
  end
  runningIndinit = 2;
  runningInd = 2;
  infections = {};
  deaths = {};
  for i = 1:16
    county = infected(runningInd,1);
    d = 1000;
    countyIndex = mod(county, d);
    y = (county - countyIndex) / d;
    runningInd = runningIndinit;
    while(i==y)%
       if runningInd > 402
         break
       end
       county = infected(runningInd,1);
       d = 1000;
       countyIndex = mod(county, d);
       y = (county - countyIndex) / d;
       if y ~= i
         break
       end
       runningInd = runningInd +1;
    end
    %cut first days, since simulations start from 02.03.20 and data starts
    %on 28.02.20
    for j = 5:size(infected,2)
      infSum = 0;
      deathSum = 0;
      for k = runningIndinit:runningInd-1
        infSum += infected(k,j);
        deathSum += death(k,j);
      end
      infections(j-4) = infSum;
      deaths(j-4) = deathSum;
    end
##    if i==1
##      timesvec(length(timesvec)+1) = dateInNumber;
##    end
    infectionsvec{length(infectionsvec)+1} = infections;
    deathsvec{length(deathsvec)+1} = deaths;
    runningIndinit = runningInd;
  end
  
  %daysDieing = 12;
  daysDieing = 15;
  dfVec = {};  
  darkFigures = {};
  for j = 1:16
    %for i=1:length(deathsvec{j})-12
    for i=1:length(deathsvec{j})-15
      if i > daysDieing
        darkFigures(i) = cell2mat(deathsvec{j}(i+daysDieing))/(cell2mat(infectionsvec{j}(i)) * 0.006); %0.006=mortality        
      end
    end
    for i=1:daysDieing
      darkFigures(i) = darkFigures(daysDieing+1);
    end
    %for i=length(deathsvec{j})-12:length(deathsvec{j})
    for i=length(deathsvec{j})-15+1:length(deathsvec{j})
      darkFigures(i) = darkFigures(length(deathsvec{j})-(daysDieing));
    end
    dfVec{length(dfVec)+1} = darkFigures;
  end

  % convert date read to a statewise format
  dataTmp = cell(16,1);
  for i=1:length(dataTmp)
  
    dataTmp{i}.time=timesvec;
    for j = 1:length(infectionsvec{1,i})
      dataTmp{i}.infected(j)=infectionsvec{i}(j);
    end
    for k = 1:length(deathsvec{1,i})
      dataTmp{i}.dead(k)=deathsvec{i}(k);
    end
    for l = 1:length(dfVec{1,i})
      dataTmp{i}.df(l)=dfVec{i}(l);
    end
  end
  
  data = cell(16,1);
  data{1} = dataTmp{8};
  data{2} = dataTmp{9};
  data{3} = dataTmp{11};
  data{4} = dataTmp{12};
  data{5} = dataTmp{4};
  data{6} = dataTmp{2};
  data{7} = dataTmp{6};
  data{8} = dataTmp{13};
  data{9} = dataTmp{3};
  data{10} = dataTmp{5};
  data{11} = dataTmp{7};
  data{12} = dataTmp{10};
  data{13} = dataTmp{14};
  data{14} = dataTmp{15};
  data{15} = dataTmp{1};
  data{16} = dataTmp{16};
  
  for i=1:16
    data{i}.name=stateNames{i};
  end  
  I = zeros(17,length(data{1}.infected));
  D = zeros(17,length(data{1}.infected));
  I(1,:) = cell2mat(data{1}.time(:));
  D(1,:) = cell2mat(data{1}.time(:));
  for i = 2:17
    I(i,:) = cell2mat(data{i-1}.infected(:));
    D(i,:) = cell2mat(data{i-1}.dead(:));
  end
  
  %write .txt file with cumulated infection numbers for each state
##  folder = ["../../Results/"];
##  fidI = fopen([folder, "/InfectedStatewise.txt"], "w");
##  fidD = fopen([folder, "/DeadStatewise.txt"], "w");
##  dlmwrite (fidI, I, "delimiter", ",", "newline", "\n");
##  dlmwrite (fidD, D, "delimiter", ",", "newline", "\n");
##  fclose(fidI);
##  fclose(fidD);
  
endfunction
