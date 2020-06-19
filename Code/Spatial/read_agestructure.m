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
## @deftypefn {} {@var{retval} =} read_agestructure() (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-04-13

function agevec = read_agestructure
  filename = '../../Daten/Bevoelkerung_ueber65_ueber80.txt';
  
  agevec = {};
  
  fid = fopen(filename);
  line = fgetl(fid);
  line = fgetl(fid);
  while(ischar(line))
    entries = strsplit(line, ",");
    if or(size(entries,2) == 4,size(entries,2)==3)
      agestructure.ags = str2num(entries{1});
      agestructure.fracOld = str2num(entries{size(entries,2)})/100;
      agevec(length(agevec)+1) = agestructure;
    endif
    line = fgetl(fid);
    
  endwhile
fclose(fid);

endfunction
