function road_profile = read_profile_2c(file_name)
%READ_PROFILE_2c reads road profile from two-column text file 
%
%   function road_profile = read_profile_2c(file_name)
%
%   file_name ... file name of a text file with 2 columns: stationing, height
%               and N rows
%   road_profile ... outputs two column table Nx2 with values: stationing, height

fi = fopen(file_name);
D = textscan(fi,'%f %f');
fclose(fi);
road_profile = [D{1},D{2}];
