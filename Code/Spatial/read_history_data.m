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
## @deftypefn {} {@var{retval} =} read_history_rki_data (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: iwtm96_user <iwtm96_user@VMWARE-IWTM96>
## Created: 2020-05-07

function out = read_history_data (dateForInitialDistribution, triggerAllHistory, max_T, parameter)
  % defauls: datenum(2020,03,06), 0, nan, ...
  % initialize missing arguments
  % if nargin == 0
  %   dateForInitialDistribution = datenum(2020,03,06);
  %   triggerAllHistory = 0;
  % elseif nargin == 1
  %   triggerAllHistory = 0;
  % end
  % if nargin < 3
  %   max_T = nan;
  % end
  if nargin < 4
    error("read_history_data arguments missing");
  end

  % fprintf('read_history_data(%d, %d, %d, %s)\n', dateForInitialDistribution, triggerAllHistory, max_T, ...);  % for debugging

  % initialize settings based on the data source
  switch parameter.initial
    case 'RKIdownload'
      filename = '../../Results/RKI-data.csv';
      delimiter = ',';
      startdateOfData = dateForInitialDistribution;  % TODO: seet this to the actual start date of the data (2020-01-27?)
      enddateOfData = startdateOfData+1;  % TODO: set this to the actual last date if necessary
    otherwise
      % set default filename and delimiter
      filename = '../../Daten/InfectedDataHistory.txt';
      delimiter = ';';
      startdateOfData = datenum(2020,03,02);
      enddateOfData = datenum(2020,05,03);
  end

  % make sure that the dates make sense
  if dateForInitialDistribution < startdateOfData || enddateOfData < dateForInitialDistribution
    fprintf('Uncommon date order: start %s, init: %s, end: %s\n', datestr(startdateOfData), datestr(dateForInitialDistribution), datestr(enddateOfData))
  endif

  if exist("datavec", "var")
    if and(dateForInitialDistribution == dateForInitialDistributionOld,...
      triggerAllHistory == triggerAllHistoryOld,max_T==max_TOld)
      out = datavec;
      return
    end 
  end
  
  persistent dateForInitialDistributionOld;
  persistent triggerAllHistoryOld;
  persistent max_TOld;
  
  dateForInitialDistributionOld = dateForInitialDistribution;
  triggerAllHistoryOld=triggerAllHistory;
  max_TOld=max_T;

  persistent datavec;
  if triggerAllHistory == 0
    datavec = {};
    position = dateForInitialDistribution - startdateOfData + 1;
    fid = fopen(filename);
    % skip first line, as this is the legend
    line = fgetl(fid);
    while !feof(fid)   
      line = fgetl(fid);
      entries = strsplit(line, delimiter);
      if size(entries,2) > 2      
        data.ags = str2num(entries{1});
        data.Infected = str2num(entries{1 + position});
        datavec(end+1) = data;  
      end
    endwhile
    fclose(fid);
  else  
    globaldata = {};
    positionStart = startdateOfData - dateForInitialDistribution + 1;
    % how does the above make sense? (its the original); replaced with the following line:
    % positionStart = dateForInitialDistribution - startdateOfData + 1;
    positionEnd = enddateOfData - startdateOfData + 1;
    fid = fopen(filename);
    % skip first line, as this is the legend
    line = fgetl(fid);
    while !feof(fid)   
      line = fgetl(fid);
      entries = strsplit(line, delimiter);
      if size(entries,2) > 2      
        data.ags = str2num(entries{1});
        for i=positionStart:1:positionEnd
          data.Infected(i) = str2num(entries{1 + i});
          data.t(i) = i-positionStart;
        end 
        globaldata(end+1) = data;  
      end
    endwhile
    fclose(fid);  
    
    datavec = [];
    numberTimestepsNeeded = min(length(globaldata{1}.t),max_T);
    datavec.t = zeros(numberTimestepsNeeded,1);
    datavec.Infected = zeros(numberTimestepsNeeded,length(globaldata));
    datavec.ags = zeros(length(globaldata),1);
    counter = 1;
    for k=1:length(globaldata{1}.t)
      if k <= max_T+1
        for j=1:length(globaldata)
          if j == 1
            datavec.t(k) = globaldata{j}.t(k);
          end
          if k == 1
            datavec.ags(j) = globaldata{j}.ags(k); 
          end
          datavec.Infected(k,j) = globaldata{j}.Infected(k);
        end
      end
    end
  end
  out = datavec;
endfunction
