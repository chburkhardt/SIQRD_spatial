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
## @deftypefn {} {@var{retval} =} read_population (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-03

function cities = read_population (filename, parameter)
  % read the population data of all cities in the given file

  % check if the result is already cached
  if exist("citiesTmp", "var")
    if strcmp(oldFilename, filename)
      cities = citiesTmp;
      return;
    end 
  endif

  % init
  persistent oldFilename;
  oldFilename = filename;
  s0=1;
  i0=0;
  r0=0;
  persistent citiesTmp;
  citiesTmp = {};
  
  % read file
  fid = fopen(filename);
  % go through file line by line
  line = fgetl(fid);
  while (ischar(line))
    entries = strsplit(line, "\t");
    % parse only valid lines
    if size(entries, 2) == 21
      city.name = entries{8};
      city.area = str2num(strrep(entries{9}, ",", "."));
      city.population = str2num(strrep(entries{11}, " ", "")); 
      if strcmp(entries{18}, "Ausland")      
        city.biggestCityInCounty = str2num(strrep(entries{12}, " ", "")); 
      else
        city.biggestCityInCounty = str2num(strrep(entries{11}, " ", "")); 
      end
      if city.population == 0
        line = fgetl(fid);
        continue;
      end
      %      city.longitude = str2num(strrep(entries{16}, ",", "."));
      %      city.latitude = str2num(strrep(entries{17}, ",", "."));
      city.density = str2num(entries{20});
      %AGS_temp = [entries{3} entries{4} entries{5}];
      %city.ags = substr(AGS_temp,2);
      
      
      [x, y] = latLon2XY(str2num(strrep(entries{17}, ",", ".")),...
      str2num(strrep(entries{16}, ",", ".")));
      city.x = x;
      city.y = y;    
      city.ags = str2num([entries{3} entries{4} entries{5}]); 
      city.verbandschluessel = str2num(entries{6});
      city.gemeindekennzahl = str2num(entries{7});
      
      city.SIR = [s0, i0, r0];
      city.fracOld = 0.0;
      
      citiesTmp(length(citiesTmp)+1) = city;
    endif    
    line = fgetl(fid);
  endwhile
  fclose(fid);
  ##if nargin == 0
  ##  show_map(citiesTmp);
  ##endif
  cities = citiesTmp;
endfunction
%
%
%
function [x, y] = latLon2XY(lat, lon)
  circumferenceAtLat = cos(lat*0.01745329251)*111.305;
  x = lon*circumferenceAtLat;
  y = lat*110.919;
endfunction
