## Copyright (C) 2020 andreaskerg
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
## @deftypefn {} {@var{retval} =} visualizeOptResult (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andreaskerg <andreaskerg@andreaskerg-VirtualBox>
## Created: 2020-05-12

function retval = visualizeOptResult (input1, input2)
  fid = fopen("../../Daten/protokoll1500runs.txt");
  data = fscanf(fid, ["Res: %f beta: %f | gamma1: %f | gamma2: %f |",...
  " mortality: %f | darkFigure: %f | startDate: %f | majorEvents: %f |",...
  " schoolClosing: %f | contactRestrictions: %f | \n"],inf);
  fclose(fid);
  data = reshape(data, 10, length(data)/10)';
  paranames = {"Res", "beta", "gamma1", "gamma2", " mortality", "darkFigure", ...
  "startDate", "majorEvents", "schoolClosing", "contactRestrictions"};
  lb=min(data, [], 1);
  ub=max(data, [], 1);
  
  resolution = 10;  
  
  close all;
  for k = 1:8
    for l = (k):8
      para1 = k+1;
      para2 = l+2;
      halfaxis = (ub([para1, para2]) - lb([para1, para2]))/10;      
      [X,Y] = meshgrid(linspace(lb(para1), ub(para1), resolution),...
      linspace(lb(para2),ub(para2), resolution));
      Z = nan(size(X));
      for i=1:size(X,1)*size(X,2)
        idxRel = ((data(:,para1) - X(i))/halfaxis(1)).^2+...
        ((data(:,para2) - Y(i)) / halfaxis(2)).^2 < 1;
        res = data(idxRel, 1);
        if and(numel(res), min(res) < 0.2)
          Z(i) = min(res);
        else
          Z(i) = nan;
        endif
      endfor
      
      if l<=4
        figure(1);
        row = k;
        col = l;
      elseif k>=5
        figure(3);
        row = k-4;
        col = l-4;
      else
        figure(2);
        row = k;
        col=l-4;       
      endif
      subplot(4,4,(row-1)*4+col), surf(X,Y,Z);
      hold on;
      plot3(data(:,para1), data(:,para2), data(:,1), '*');
      xlabel(paranames(para1));
      ylabel(paranames(para2));
      zlabel("Residuum");
      caxis([0,0.05]);
      zlim([0, 0.1]);
      colormap summer;
      shading interp;
      view(2);
    endfor
  endfor
  
endfunction













