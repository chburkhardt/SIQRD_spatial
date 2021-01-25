function out = writeTexTableDFs(input1, input2)
  RKIread = read_case_history_RKIfiles();
  
  %rewrite the cell, since states are arranged differently in RKIread and statewiseSIR
  RKIdata = cell(16,1);
  RKIdata{1} = RKIread{15};
  RKIdata{2} = RKIread{6};
  RKIdata{3} = RKIread{9};
  RKIdata{4} = RKIread{5};
  RKIdata{5} = RKIread{10};
  RKIdata{6} = RKIread{7};
  RKIdata{7} = RKIread{11};
  RKIdata{8} = RKIread{1};
  RKIdata{9} = RKIread{2};
  RKIdata{10} = RKIread{12};
  RKIdata{11} = RKIread{3};
  RKIdata{12} = RKIread{4};
  RKIdata{13} = RKIread{8};
  RKIdata{14} = RKIread{13};
  RKIdata{15} = RKIread{14};
  RKIdata{16} = RKIread{16};

  for i=1:length(RKIdata)
    RKIdata{i}.df = cell2mat(RKIdata{i}.df);
  end
  
	filename = "../../Results/DataTableDFs.dat";
	titles = {"DF"};
	nVals = 1;
  
  out = zeros(length(RKIdata{1}.time), length(RKIdata) * nVals + 1);
  out(:,1) = cell2mat(RKIdata{1}.time) - datenum([2020 3 2]);
  for i=1:length(RKIdata)
    out(:, (1 + nVals*(i-1)+1):(1 + nVals*i)) = [RKIdata{i}.df'];
  endfor
  
  fid = fopen(filename, 'w'); 
  fprintf(fid, "%s ", "time");
  for i=1:length(RKIdata)
    for j=1:length(titles)
      fprintf(fid, "%s ", [strrep(titles{j}," ","_"), "_",...
      strrep(RKIdata{i}.name," ","_")]);
    endfor
  endfor
  fprintf(fid, "\n");
  for i=1:size(out, 1)
    fprintf(fid, "%f ", out(i,:));
    fprintf(fid, "\n");
  endfor  
  fclose(fid);
	
	
endfunction
