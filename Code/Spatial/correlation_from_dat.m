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
## @deftypefn {} {@var{retval} =} correlation_from_dat (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-06-29

function r_squared = correlation_from_dat (filename_sim_data, filename_rki_data, InfectedOrDead)
if nargin == 0
  filename_sim_data = '../../sirhCountywise.dat';
  filename_rki_data = '../../rkiDataTable.dat';
  InfectedOrDead = 'Dead'
end
if nargin == 2
    InfectedOrDead = 'Infected'
end

sim_data = dlmread(filename_sim_data,'');
rki_data = dlmread(filename_rki_data,'');
%% get vector with infected(discovered) or dead from sim (57 days x 16 Bundesländer)
%% get vector with infected or dead from rki (57 days x 16 Bundesländer)
sim_state_daywise = zeros(16,size(sim_data,1));
rki_state_daywise = zeros(16,size(sim_data,1));
for i=1:16
  switch InfectedOrDead
    case 'Infected'
      sim_state_daywise(i,:) = sim_data(:,1+i*9);
      rki_state_daywise(i,:) = rki_data(:,i*2);
    case 'Dead'
      sim_state_daywise(i,:) = sim_data(:,i*9);
      rki_state_daywise(i,:) = rki_data(:,1+i*2);
    end
end
%% kick entries of first colum out, as this are just 0 values
sim_state_daywise(:,[1]) = [];
rki_state_daywise(:,[1]) = [];

%% calc correlation factor for each federal-state
r_squared = zeros(16,1);
p_values = zeros(16,1);
for i=1:16
  [r1, p1] = corrcoef(sim_state_daywise(i,:)',rki_state_daywise(i,:)');
  r_squared(i) = r1(1,2)*r1(1,2);
  p_values(i) = p1(1,1);
end
    
  


endfunction
