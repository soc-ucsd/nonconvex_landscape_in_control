% Specify the path to your text file
% file_path = ["data/1dim.txt","data/2dim.txt","data/3dim.txt",...
%     "data/AC8.txt","data/HE1.txt","data/REA2.txt"];
% 
% save_path = ["data/1dim.mat","data/2dim.mat","data/3dim.mat",...
%     "data/AC8.mat","data/HE1.mat","data/REA2.mat"];
file_paths = ["data/dim1.txt","data/dim2.txt","data/dim3.txt"];

save_paths = ["data/dim1.mat","data/dim2.mat","data/dim3.mat"];

save_names = ["dim1","dim2","dim3"];

for i=1:length(file_paths)

    % Open the file for reading
    fid = fopen(file_paths(i), 'r');

    % Check if the file was opened successfully
    if fid == -1
        error('Error opening the file');
    end

    % Initialize an array to store the extracted values
    values = [];

    % Read each line of the file
    while ~feof(fid)
        % Read a line from the file
        line = fgetl(fid);

        % Check if the line contains "f="
        if contains(line, 'f = ')
            % Extract the value after "f="
            value_str = regexp(line, 'f = ([\d.]+)', 'match', 'once');

            value_str = regexprep(value_str, '[^0-9.]', '');

            % Convert the string to a number and add it to the array
            value = str2double(value_str);

            % Check if the conversion was successful
            if ~isnan(value)
                values = [values, value];
            end
        end
    end

    % Close the file
    fclose(fid);

    % % Display the extracted values
    % disp('Extracted values:');
    % disp(values);

    % Save the values to a file (optional)
    save_path = save_paths(i);
    save_name = save_names(i);
    save(save_path, 'values');
end

