%% this m file is using one single DLC file to detect whisker base of both Right whiskers and Left whiskser (in Mirror R format)
data_folder = '\\dk-server.dk.ucsd.edu\afassihizakeri\topview data\2023_02_22_ 163923';
tracker = WhiskerTracker(data_folder);
tracker.main()