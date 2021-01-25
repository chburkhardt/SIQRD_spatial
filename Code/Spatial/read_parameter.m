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
## @deftypefn {} {@var{retval} =} read_parameter (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-04-20

function parameter = read_parameter (filename)
  fid = fopen(filename);
  line = fgetl(fid);
  while(ischar(line))
  % ignore comments
  entries = strsplit(line(1:(strfind(line,"%") - 1)));
  
  %ignore lines not starting with "set"
  if ~strcmp(entries{1}, "set") 
    line = fgetl(fid);
    continue;
  end
  
  switch entries{2}
    case "model"
      parameter.model = getValue(entries{4}, "string");
    case "analytically"
      parameter.analytically = getValue(entries{4}, "boolean");
    case "reduce"
      parameter.reduce = getValue(entries{4}, "boolean");
    case "reduceToStates"
      parameter.reduceToStates = getValue(entries{4}, "boolean");
    case "keepAGS"
      parameter.keepAGS = getValue(entries{4}, "numericalArray");
    case "beta"
      parameter.beta = getValue(entries{4}, "numerical");
    case "betaSWscaling"
      parameter.betaSWscaling = getValue(entries{4}, "numerical");
    case "betaStatewise"
      parameter.betaStatewise = getValue(entries{4}, "numericalArray");
    case "gamma1"
      parameter.gamma1 = getValue(entries{4}, "numerical");
    case "gamma2"
      parameter.gamma2 = getValue(entries{4}, "numerical");
    case "darkFigure"
      parameter.darkFigure = getValue(entries{4}, "numerical");
    case "mortality"
      parameter.mortality = getValue(entries{4}, "numerical");
    case "schoolClosing"
      parameter.schoolClosing = getValue(entries{4}, "numerical");
    case "contactRestrictions"
      parameter.contactRestrictions = getValue(entries{4}, "numerical");
    case "exitRestrictions"
      parameter.exitRestrictions = getValue(entries{4}, "numerical");
    case "majorEvents"
      parameter.majorEvents = getValue(entries{4}, "numerical");
    case "masks"
      parameter.masks = getValue(entries{4}, "numerical");
##    case "betaStateWise"
##      parameter.betaStateWise = getValue(entries{4}, "boolean");
    case "c_s"
      parameter.c_s = getValue(entries{4}, "numerical");
    case "c_i"
      parameter.c_i = getValue(entries{4}, "numerical");
    case "beta_cross_county"
      parameter.beta_cross_county = getValue(entries{4}, "numerical");
    case "beta_cc_statewise"
      parameter.beta_cc_statewise = getValue(entries{4}, "numericalArray");
    case "n_traffic"
      parameter.n_traffic = getValue(entries{4}, "numerical");
    case "spatial"
      parameter.spatial = getValue(entries{4}, "string");
    case "initial"
      parameter.initial = getValue(entries{4}, "string");
    case "totalRuntime"
      parameter.totalRuntime = getValue(entries{4}, "numerical");
    case "startDate"
      parameter.startDate = getValue(entries{4}, "date");
    case "endDateOpt"
      parameter.endDateOpt = getValue(entries{4}, "date");
    case "initalDistributionDate"
      parameter.initalDistributionDate = getValue(entries{4}, "date");
    case "maxTimeStep"
      parameter.maxTimeStep = getValue(entries{4}, "numerical");
    case "nRuns"
      parameter.nRuns = getValue(entries{4}, "numerical");
    case "showDiagrams"
      parameter.showDiagrams = getValue(entries{4}, "boolean"); 
    case "saveDiagrams"
      parameter.saveDiagrams = getValue(entries{4}, "boolean");
    case "saveVTK"
      parameter.saveVTK = getValue(entries{4}, "boolean");
    case "saveConnections"
      parameter.saveConnections = getValue(entries{4}, "boolean");
    case "resolutionVTK"
      parameter.resolutionVTK = getValue(entries{4}, "numerical");
    case "blowUp"
      parameter.blowUp = getValue(entries{4}, "boolean");
    case "folderName"
      parameter.folderName = getValue(entries{4}, "string");
    case "fullConsoleOut"
      parameter.fullConsoleOut = getValue(entries{4}, "boolean");
    case "MobNetworkModel"
      parameter.MobNetworkModel = getValue(entries{4}, "boolean");
    case "OutputComparisonRKI"
      parameter.OutputComparisonRKI = getValue(entries{4},"boolean");
    case "discoveredInVTK"
      parameter.discoveredInVTK = getValue(entries{4},"boolean");
    case "vtkStyleMap"
      parameter.vtkStyleMap = getValue(entries{4},"boolean");
    case "seed"
      parameter.seed = getValue(entries{4},"boolean");      
    case "seedAGS"
      parameter.seedAGS = getValue(entries{4},"numericalArray");
    case "seedDate"
      parameter.seedDate = getValue(entries{4},"dateArray");
    case "seedSize"
      parameter.seedSize = getValue(entries{4},"numericalArray");
    case "easingContactRestrictions"
      parameter.easingContactRestrictions = getValue(entries{4},"numerical");
    case "endOfVacation"
      parameter.endOfVacation = getValue(entries{4},"numerical");
    case "liftingTravelRestictionsEU"
      parameter.liftingTravelRestictionsEU = getValue(entries{4},"numerical"); 
    case "lockdownLight"
      parameter.lockdownLight = getValue(entries{4},"numerical"); 
    case "considerAverageAgeCounty"
      parameter.considerAverageAgeCounty = getValue(entries{4},"boolean");
    case "considereAgeOfNewlyInfected"
      parameter.considereAgeOfNewlyInfected = getValue(entries{4},"boolean");
    case "factorsDF"
      parameter.factorsDF = getValue(entries{4}, "numericalArray");
  endswitch
  line = fgetl(fid);
end
fclose(fid);
parameter = expandParameters(parameter);
sensecheck(parameter);
endfunction

function parameterArray = expandParameters(parameter)
names = fieldnames(parameter);
parameterArray = {};

% add one parameterstruct to parameterArray for all manifestations of the
% first parameter
for i=1:length(parameter.(names(1){1}))
  if iscell(parameter.(names(1){1})(i))
    tmpParameter.(names(1){1}) = cell2mat(parameter.(names(1){1})(i));
  else
    tmpParameter.(names(1){1}) = parameter.(names(1){1})(i);
  end
  parameterArray{end + 1} = tmpParameter;      
end  

for i=2:length(names) % loop over all remaining parameters
  % duplicate all parametersets as often as the ith parameter has manifestations
  oldLength = length(parameterArray);
  parameterArray = repmat(parameterArray, 1, length(parameter.(names(i){1})));
  for j=1:length(parameter.(names(i){1}))
    % add the jth manifestation of parameter(name(i)) to all entries of the array
    for k=1:oldLength      
      if iscell(parameter.(names(i){1})(j))
        parameterArray{(j-1)*oldLength + k}.(names(i){1}) = cell2mat(parameter.(names(i){1})(j));     
			else
        parameterArray{(j-1)*oldLength + k}.(names(i){1}) = parameter.(names(i){1})(j);        
      end      
    end  
  end  
end
end


function retval = getValue(entry, type)
switch type
  case "boolean"
    if or(numel(strfind(entry, ":")), numel(strfind(entry, ";")))
      retval = {true, false};
    elseif strcmp(entry, "true")
      retval = true;
    else
      retval = false;
    end
    
    
  case "string"
    if numel(strfind(entry, ";"))
      retval = strsplit(entry, ";");
    else
      retval = {entry};
    end
    
    
  case "date"
    if numel(strfind(entry, ";"))
      retval = datenum(strsplit(entry, ";")); 
    else
      retval = datenum(strrep(entry, "_", " "));
    end   
    
   case "dateArray"
    if numel(strfind(entry, ";"))
      retval = datenum(strsplit(entry, ";")); 
    else
      retval = {datenum(strrep((strsplit(entry, ",")), "_", " "))};
    end   

  case "numerical"
    if numel(strfind(entry, ";"))
      retval = str2num(char(strsplit(entry, ";")));
    elseif numel(strfind(entry, ":"))
      tmp = str2num(char(strsplit(entry, ":")));
      if length(tmp) != 3
        fprintf ("Entry was: %s\n", entry);
        error("Wrong format for range");
      end
      retval = tmp(1):tmp(2):tmp(3);
    else
      retval = str2num(entry);
    end    
    
  case "numericalArray"
    if numel(strfind(entry, ";"))
      retval = str2num(char(strsplit(entry, ";")));
    elseif numel(strfind(entry, ":"))
      tmp = str2num(char(strsplit(entry, ":")));
      if length(tmp) != 3
        fprintf ("Entry was: %s\n", entry);
        error("Wrong format for range");
      end
      retval = tmp(1):tmp(2):tmp(3);
    else
      retval = {str2double(strsplit(entry, ","))};
    end
endswitch
if !exist("retval", "var")
  error("unknown type");
end
endfunction

function retval = sensecheck(parameter)
for i=1:length(parameter)
	if and(parameter{i}.blowUp, parameter{i}.OutputComparisonRKI)
		error("Blow up cannot be used with OutputComparisonRKI");
	endif
	
##	if parameter{i}.vtkStyleMap	
##		if or(!parameter{i}.reduce, parameter{i}.reduceToStates)
##			error("vtkStyleMap only works for countywise calculations");
##		endif
##	endif	
endfor
endfunction










