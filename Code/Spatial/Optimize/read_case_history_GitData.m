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

function data = read_case_history_test ()
  filename_infected = '../../Daten/cases-rki-by-states_19082020.txt';
  filename_dead = '../../Daten/deathNumbers-rki-by-states_19082020.txt';
  
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
  
  %% read infection numbers from txt file
  fid = fopen(filename_infected);
  infectionsvec = {};
  timesvec = {};
  %skip first 5 lines as this is the header
  for i=1:2
    line = fgetl(fid);
  end
  while(ischar(line))
    entries = strsplit(line, ";");
    if (length(entries)==18)
      numbers = strsplit(entries{1},"-");
      % convert the date to a number that can be calculated with
      dateInNumber =datenum(str2num(numbers{3}),str2num(numbers{2}),str2num(numbers{1}));
      infections = {};
      for i = 1:16
        infections(i) = str2num(entries{i+1});
      end
      infectionsvec{length(infectionsvec)+1} = infections;
      timesvec(length(timesvec)+1) = dateInNumber;
    end
    line = fgetl(fid);
  end
  fclose(fid);

  
  %% read death numbers from txt file
  fid_2 = fopen(filename_dead);
  deathsvec = {};
  timesvecdead = {};
  %skip first 5 lines as this is the header
  for i=1:2
    line = fgetl(fid_2);
  end
  while(ischar(line))
    entries = strsplit(line, ";");
    if (length(entries)==18)
      %% convert the date to a number that can be calculated with
      numbers = strsplit(entries{1},"-");
      dateInNumber = datenum(str2num(numbers{3}),str2num(numbers{2}),str2num(numbers{1}));
      dead = {};
      for i = 1:16
        dead(i) = str2num(entries{i+1});
      end
      deathsvec{length(deathsvec)+1} = dead;
      timesvecdead{length(timesvecdead)+1} = dateInNumber;
    end
    line = fgetl(fid_2);
  end
  fclose(fid_2);
  
  % convert date read to a statewise format
  dataTmp = cell(16,1);
  for i=1:length(dataTmp)
  
    dataTmp{i}.time=timesvec;
    for j = 1:min(length(infectionsvec))
      dataTmp{i}.infected(j)=infectionsvec{j}(i);
    end
    for k = 1:min(length(deathsvec), length(infectionsvec))
      dataTmp{i}.dead(k)=deathsvec{k}(i);
    end
  end
  
  data = cell(16,1);
  data{1} = dataTmp{8};
  data{2} = dataTmp{9};
  data{3} = dataTmp{16};
  data{4} = dataTmp{11};
  data{5} = dataTmp{4};
  data{6} = dataTmp{2};
  data{7} = dataTmp{6};
  data{8} = dataTmp{12};
  data{9} = dataTmp{3};
  data{10} = dataTmp{5};
  data{11} = dataTmp{7};
  data{12} = dataTmp{10};
  data{13} = dataTmp{13};
  data{14} = dataTmp{14};
  data{15} = dataTmp{1};
  data{16} = dataTmp{15};
  
  for i=1:16
    data{i}.name=stateNames{i};
  end
  
endfunction
