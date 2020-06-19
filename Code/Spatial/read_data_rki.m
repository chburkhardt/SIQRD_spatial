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
## @deftypefn {} {@var{retval} =} read_population_rki (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-04-11

function datavec = read_data_rki (specialDay)
  if nargin == 0
    filename = '../../Daten/RKI_Corona_Landkreise.txt';
    
    
    datavec = {};
    casesBerlin = 0;
    populationBerlin = 0;
    deathBerlin = 0;
    
    fid = fopen(filename);
    % skip first line, as this is the legend
    line = fgetl(fid);
    line = fgetl(fid);
    while(ischar(line))
    
    entries = strsplit(line, ",");
    if size(entries,2) == 40
      %data.name = entries{8};
      data.ags = str2num(entries{6});
      data.Ifrac = str2num(entries{30}) / str2num(entries{24});
      data.Dfrac = str2num(entries{31}) / str2num(entries{24});
      datavec(length(datavec)+1) = data;     
    else
      % Berlin is splitted in RKI data but not in population data
      data.ags = str2num(entries{2});
      if data.ags >=11000 && data.ags <=11011
        casesBerlin += str2num(entries{9});
        deathBerlin += str2num(entries{10});
        populationBerlin += str2num(entries{5});
      endif  
    endif
    line = fgetl(fid);
  endwhile
  fclose(fid);
  
  % add Berlin
  data.ags = 11000;
  data.Ifrac = casesBerlin / populationBerlin;
  data.Dfrac = deathBerlin / populationBerlin;
  
  datavec(length(datavec)+1) = data;
elseif
    filename = '../../Daten/InfectedOverAGS_RKI_28032020_11042020_25042020.txt';
    
    if (specialDay == datenum(2020,03,28))
      posInfected = 2;
    elseif (specialDay == datenum(2020,04,11))
      posInfected = 3;
    elseif (specialDay == datenum(2020,04,25))
      posInfected = 4;
    else   
      error("Day to read agswise RKI infection numbers in, not in database");
    endif  
    
    datavec = {};   
    fid = fopen(filename);
    % skip first line, as this is the legend
    line = fgetl(fid);
        line = fgetl(fid);
    line = fgetl(fid);
    while(ischar(line))
    entries = strsplit(line, ";");
    if size(entries,2) >= 3
      data.ags = str2num(entries{1});
      data.Infected = str2num(entries{posInfected});
      datavec(length(datavec)+1) = data;     
    endif
    line = fgetl(fid);
  endwhile
  fclose(fid);
endif

endfunction
