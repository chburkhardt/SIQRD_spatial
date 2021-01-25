function factors = read_in_deaths()
  filename_infected_Age = '../../Daten/cases-rki-by-age-ags_current.csv';
  filename_infected = '../../Daten/cases-rki-by-ags_current.csv';
  delimiterIn = ',';
  infectedByAge = importdata(filename_infected_Age,delimiterIn);
  infected = importdata(filename_infected,delimiterIn);
  
  lengthAgeGroup = 402;
  
  mortalityAgeGroup = [0 0 0.0001 0.0017 0.034 0.1382];

  normalMortailty = 0.016957926;
  
  infectedAgeGroups = cell(6,1);
  
  c = 1;
  for i = 1:lengthAgeGroup:size(infectedByAge,1)-402
    infectedAgeGroups{c}.infectedByAge = infectedByAge(i:c*lengthAgeGroup,:);
    c += 1;
  end
  
  mortality = zeros(402, size(infectedByAge,2));
  for i = 1:length(infectedAgeGroups)
    infectedAgeGroups{i}.nomalMortInAgeGroup = zeros(402, size(infectedByAge,2));
    infectedAgeGroups{i}.nomalMortInAgeGroup(1,:) = infectedAgeGroups{i}.infectedByAge(1,:);
    infectedAgeGroups{i}.nomalMortInAgeGroup(:,1) = infectedAgeGroups{i}.infectedByAge(:,1);
    infectedAgeGroups{i}.nomalMortInAgeGroup(2:end,2:end) =...
    infectedAgeGroups{i}.infectedByAge(2:end,2:end)./infected(2:end,2:end) * mortalityAgeGroup(i);
    infectedAgeGroups{i}.nomalMortInAgeGroup(isnan(infectedAgeGroups{i}.nomalMortInAgeGroup))=0;
  end
  
  mortality = zeros(402, size(infectedByAge,2));
  mortality(1,:) = infectedAgeGroups{1}.infectedByAge(1,:);
  mortality(:,1) = infectedAgeGroups{1}.infectedByAge(:,1);
  for i = 1:length(infectedAgeGroups)
    mortality(2:end,2:end) += infectedAgeGroups{i}.nomalMortInAgeGroup(2:end,2:end);
  end
  
  factorsMortality = mortality/normalMortailty;
  factorsMortality(1,:) = mortality(1,:);
  factorsMortality(:,1) = mortality(:,1);
  
  factors = factorsMortality;
  
  reduceToStates = true;
  
  if reduceToStates 
    ags = 1:16;
    reduction = cell(1, length(ags));
    factorsMortalitySW = zeros(length(ags)+1, size(infectedByAge,2));
    factorsMortalitySW(1,:) = mortality(1,:);
    factorsMortalitySW(2:end,1) = ags;
    for i = 1:length(reduction)
      reduction{i} = find(floor(factorsMortality(2:end,1)/1e3) == ags(i));
      factorsMortalitySW(i+1,2:end) = sum(factorsMortality(reduction{i}+1,2:end),1)/length(reduction{i});
    end
    factors = factorsMortalitySW;
  end
  
endfunction
