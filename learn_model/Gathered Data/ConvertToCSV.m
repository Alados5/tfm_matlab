function ConvertToCSV(ExpDate)
%
% ConvertToCSV(ExpDate)
%
% Converts the TXT files in the 'Converted TXT' subfolder of the input
% experiment session into .CSV files. If the final folder for the CSVs
% does not exist, it is also created beforehand.
%
% ExpDate must be a string indicating a date in 'yyyy_mm_dd' format.

dirname = ['.\Session_',ExpDate, '\Converted CSV'];

if ~isfolder(['.\Session_',ExpDate])
    error("There is no data folder for the indicated session (date)");
elseif ~isfolder(['.\Session_',ExpDate,'\Converted TXT'])
    error("There is no 'Converted TXT' folder inside this session folder");
elseif ~isfolder(dirname)
    copyfile(['.\Session_',ExpDate,'\Converted TXT'], dirname);    
end

fileList = dir([dirname, '\*.txt']);

for i = 1:numel(fileList)
    file = fullfile(dirname, fileList(i).name);
    [tempDir, tempFile] = fileparts(file); 
    status = copyfile(file, fullfile(tempDir, [tempFile, '.csv']));
    if status == 1
        delete(file)  
    else
        error("File could not be converted");
    end 
end

end