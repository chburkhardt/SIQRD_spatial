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

function out = read_history_data (dateForInitialDistribution, triggerAllHistory, max_T)
  if nargin == 0
    dateForInitialDistribution =datenum(2020,03,06);
    triggerAllHistory = 0;
  elseif nargin == 1
    triggerAllHistory = 0;
  end	
	if nargin < 3
		max_T = nan;
	end
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
	filename = '../../Daten/InfectedDataHistory.txt';
	if triggerAllHistory == 0
		datavec = {};
		startdateOfData = datenum(2020,03,02);
		position = dateForInitialDistribution - startdateOfData + 1;
		fid = fopen(filename);
		% skip first line, as this is the legend
		line = fgetl(fid);
		while !feof(fid)   
			line = fgetl(fid);
			entries = strsplit(line, ";");
			if size(entries,2) > 2      
				data.ags = str2num(entries{1});
				data.Infected = str2num(entries{1 + position});
				datavec(end+1) = data;  
			end
		endwhile
		fclose(fid);
	else  
		globaldata = {};
		startdateOfData = datenum(2020,03,02);
		enddateOfData = datenum(2020,05,03);
		positionStart = startdateOfData - dateForInitialDistribution + 1;
		positionEnd = enddateOfData - startdateOfData + 1;
		fid = fopen(filename);
		% skip first line, as this is the legend
		line = fgetl(fid);
		while !feof(fid)   
			line = fgetl(fid);
			entries = strsplit(line, ";");
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
