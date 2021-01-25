function ageDistributionNewInfections = getAgeDistributionNewInfections ()
  filename_infected_Age = '../../Daten/cases-rki-by-age-ags_current.csv';
  filename_infected = '../../Daten/cases-rki-by-ags_current.csv';
  delimiterIn = ',';
  infectedByAge = importdata(filename_infected_Age,delimiterIn);
  infected = importdata(filename_infected,delimiterIn);
  
  newInfectionsSum = diff(sum(infected(2:end,2:end),1));
  
  lengthAgeGroup = 402;
  
  infectedAgeGroups = cell(6,1);
  
  c = 1;
  for i = 1:lengthAgeGroup:size(infectedByAge,1)-402
    infectedAgeGroups{c}.infectedByAge = infectedByAge(i:c*lengthAgeGroup,:);
    c += 1;
  end
  
##  mortality = zeros(402, size(infectedByAge,2));
##  for i = 1:length(infectedAgeGroups)
##    infectedAgeGroups{i}.nomalMortInAgeGroup = zeros(402, size(infectedByAge,2));
##    infectedAgeGroups{i}.nomalMortInAgeGroup(1,:) = infectedAgeGroups{i}.infectedByAge(1,:);
##    infectedAgeGroups{i}.nomalMortInAgeGroup(:,1) = infectedAgeGroups{i}.infectedByAge(:,1);
##    infectedAgeGroups{i}.nomalMortInAgeGroup(2:end,2:end) =...
##    infectedAgeGroups{i}.infectedByAge(2:end,2:end)./infected(2:end,2:end) * mortalityAgeGroup(i);
##    infectedAgeGroups{i}.nomalMortInAgeGroup(isnan(infectedAgeGroups{i}.nomalMortInAgeGroup))=0;
##  end
  
  for i = 1:length(infectedAgeGroups)
    infectedAgeGroups{i}.sumAgeGroup = sum(infectedAgeGroups{i}.infectedByAge(2:end,2:end),1);
    infectedAgeGroups{i}.newInfections = diff(infectedAgeGroups{i}.sumAgeGroup);
    infectedAgeGroups{i}.proportionNewInfs = infectedAgeGroups{i}.newInfections./newInfectionsSum;
  end
  
  
  filename = "../../Results/newInfectionsAge.dat";
	%filename = "../../Results/rkiDataTableNew.dat";
	titles = {"0-4J", "5-14J", "15-34J", "35-59J", "60-79J", "80+J"};
	nVals = 2;
	% time x (time, vals per county x counties)

  time = infectedAgeGroups{1}.infectedByAge(1,2:end);

  out = zeros(length(time), length(infectedAgeGroups) + 1);
  out(:,1) = time - datenum([2020 3 2]);
  for i=1:length(infectedAgeGroups)
    out(2:end, i+1) = infectedAgeGroups{i}.proportionNewInfs;
  endfor
  
  fid = fopen(filename, 'w'); 
  fprintf(fid, "%s ", "time");
  for j=1:length(titles)
      fprintf(fid, "%s ", [strrep(titles{j}," ","_")]);
  endfor
  
  fprintf(fid, "\n");
  for i=1:size(out, 1)
    fprintf(fid, "%f ", out(i,:));
    fprintf(fid, "\n");
  endfor  
  fclose(fid);
 
   %clf;
   figure;
   Contours=0:0.1:1;
   colormap ("viridis");
   labelNames = {"0-4 Y", "5-14 Y", "15-34 Y", "35-59 Y", "60-79 Y", "80+ Y"}; 
   img = out(2:end,2:end)';
   x = 0:289;
   y = 1:6;
   h = imagesc (x, y, img); 
   set(gca, 'YTick', 1:6); 
   set(gca,'YTickLabel',labelNames);
   set(gca,'YDir','normal');
   set(gca,'Position',[0.1 0.1 0.9 0.4])
   colorbar;   
   axis fill
  
endfunction
