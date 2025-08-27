    function run_magia_gui()

    % Main GUI Window
    f = figure('Name','Tool Run Config', 'Position', [500 200 600 580], 'Resize', 'on');

    % --- PET toolbox Logo Text ---
    uicontrol(f, 'Style', 'text', ...
        'String', 'PET Toolbox v0.1', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'ForegroundColor', [0.1 0.3 0.6], ...
        'Position', [180 530 240 30], ...
        'HorizontalAlignment', 'center');

    %% Inputs and Outputs Panel
    panelIO = uipanel(f, 'Title', 'Input/Output Directories', ...
        'FontSize', 12, 'Position', [0.05 0.52 0.9 0.35]);

    % PET Directory
    uicontrol(panelIO, 'Style', 'text', 'Position', [10 140 120 20], ...
        'String', 'PET Directory:', 'HorizontalAlignment', 'left');
    dirEdit = uicontrol(panelIO, 'Style', 'edit', 'Position', [140 140 280 25], ...
        'Enable', 'inactive', 'String', '');
    uicontrol(panelIO, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [430 140 70 25], ...
        'Callback', @(~,~)browse_data_parent_directory());

    % Subject(s)
    uicontrol(panelIO, 'Style', 'text', 'Position', [10 100 120 20], ...
        'String', 'Subject(s):', 'HorizontalAlignment', 'left');
    sbEdit = uicontrol(panelIO, 'Style', 'edit', 'Position', [140 100 280 25], ...
        'Enable', 'inactive', 'String', '');
    uicontrol(panelIO, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [430 100 70 25], ...
        'Callback', @(~,~)browse_subjects());


    % Sessions
    uicontrol(panelIO, 'Style', 'text', 'Position', [10 60 120 20], ...
        'String', 'Sessions:', 'HorizontalAlignment', 'left');
    session1Check = uicontrol(panelIO, 'Style', 'checkbox', ...
        'Position', [140 60 120 25], 'String', 'Use session 1', 'Value', 1);
    session2Check = uicontrol(panelIO, 'Style', 'checkbox', ...
        'Position', [280 60 120 25], 'String', 'Use session 2', 'Value', 0);

    % Output Folder
    uicontrol(panelIO, 'Style', 'text', 'Position', [10 20 120 20], ...
        'String', 'Output Folder:', 'HorizontalAlignment', 'left');
    workfolderEdit = uicontrol(panelIO, 'Style', 'edit', ...
        'Position', [140 20 280 25], 'String', '/vf/users/yanw4/GUI_MAGIA');
    uicontrol(panelIO, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [430 20 70 25], ...
        'Callback', @(~,~)browse_output_folder());

    %% Pipeline Settings Panel
    panelPipe = uipanel(f, 'Title', 'Pipeline Settings', ...
        'FontSize', 12, 'Position', [0.05 0.15 0.9 0.35]);

    % Tracer
    uicontrol(panelPipe, 'Style', 'text', 'Position', [10 150 120 20], ...
        'String', 'Tracer (pet_tr):', 'HorizontalAlignment', 'left');
    tracerPopup = uicontrol(panelPipe, 'Style', 'popupmenu', ...
        'String', {'[11c]nnc','[11c]raclopride'}, 'Value', 1, ...
        'Position', [140 150 120 25]);

    % modeling_options.txt
    uicontrol(panelPipe, 'Style', 'text', 'Position', [10 110 120 20], ...
        'String', 'Options:', 'HorizontalAlignment', 'left');
    optionsEdit = uicontrol(panelPipe, 'Style', 'edit', ...
        'Position', [140 110 280 25], 'String', './demo_modeling_options.txt');
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [430 110 70 25], ...
        'Callback', @(~,~)browse_options());
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Preview', ...
        'Position', [140 80 100 25], ...
        'Callback', @(~,~)preview_modeling_options());
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Create New', ...
        'Position', [260 80 100 25], ...
        'Callback', @(~,~)setup_modeling_options_gui());

    % specs.txt
    uicontrol(panelPipe, 'Style', 'text', 'Position', [10 40 120 20], ...
        'String', 'Specs:', 'HorizontalAlignment', 'left');
    specsEdit = uicontrol(panelPipe, 'Style', 'edit', ...
        'Position', [140 40 280 25], 'String', './demo_specs.txt');
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Browse...', ...
        'Position', [430 40 70 25], ...
        'Callback', @(~,~)browse_specs());
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Preview', ...
        'Position', [140 10 100 25], ...
        'Callback', @(~,~)preview_specs_file());
    uicontrol(panelPipe, 'Style', 'pushbutton', 'String', 'Create New', ...
        'Position', [260 10 100 25], ...
        'Callback', @(~,~)setup_specs_gui());

    %% Action Buttons

    % Run button (moved down from Y=120 to Y=40)
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Run MAGIA', ...
        'Position', [120 40 160 45], ...
        'FontSize', 12, 'BackgroundColor', [0.8 1 0.8], ...
        'Callback', @(~,~)run_all());

    % SPM Visualization (moved down from Y=60 to Y=-20)
    uicontrol(f, 'Style', 'pushbutton', ...
        'String', 'Visualize NIfTI Files (SPM)', ...
        'Position', [300 40 200 45], ...
        'FontSize', 10, 'BackgroundColor', [0.9 0.9 1], ...
        'Callback', @(~,~)visualize_nii_files());


    %% Callbacks

    function browse_data_parent_directory()
        selectedDir = uigetdir(pwd, 'Select Data Parent Directory');
        if selectedDir ~= 0
            dirEdit.String = selectedDir;
            setappdata(f, 'DataParentDir', selectedDir);
        end
    end

    % Browse for subfolders inside selected parent directory
    function browse_subjects()
        parentDir = getappdata(f, 'DataParentDir');
        if isempty(parentDir)
            errordlg('Please select the PET directory first.','Missing Directory');
            return;
        end
    
        % Get list of subfolders in the selected directory
        allItems = dir(parentDir);
        subfolders = allItems([allItems.isdir] & ~startsWith({allItems.name}, '.'));
    
        [selectionIdx, ok] = listdlg( ...
            'PromptString', 'Select Subject Folders:', ...
            'ListString', {subfolders.name}, ...
            'SelectionMode', 'multiple', ...
            'ListSize', [300, 200]);
    
        if ok
            subjects = {subfolders(selectionIdx).name};
            sbEdit.String = strjoin(subjects, ', ');
            assignin('base', 'subjects', subjects);
            assignin('base', 'directory', parentDir);
        end
    end


    function browse_output_folder()
        folder = uigetdir('', 'Select Output Folder');
        if folder ~= 0 
            set(workfolderEdit, 'String', folder);
        end
    end


    function browse_options()
        [file, path] = uigetfile('*.txt', 'Select modeling_options.txt');
        if isequal(file, 0)
            return;
        else
            set(optionsEdit, 'String', fullfile(path, file));
        end
    end

    function preview_modeling_options()
        filepath = get(optionsEdit, 'String');
        if ~isfile(filepath)
            errordlg('The specified modeling_options.txt file does not exist.', 'File Not Found');
            return;
        end

        try
            % Read text file line by line
            fid = fopen(filepath, 'r');
            rawLines = {};
            line = fgetl(fid);
            while ischar(line)
                rawLines{end+1} = line; %#ok<AGROW>
                line = fgetl(fid);
            end
            fclose(fid);

            % Split lines into key-value pairs
            data = {};
            for i = 1:length(rawLines)
                tokens = strsplit(rawLines{i}, ':');
                if numel(tokens) >= 2
                    key = strtrim(tokens{1});
                    value = strtrim(strjoin(tokens(2:end), ':'));
                    data{end+1, 1} = key; %#ok<AGROW>
                    data{end, 2} = value;
                end
            end

            % Create a new figure to display the table
            previewFig = figure('Name', 'Modeling Options Preview', ...
                                'Position', [600 400 500 300], ...
                                'Resize', 'off', ...
                                'NumberTitle', 'off', ...
                                'MenuBar', 'none');

            uitable(previewFig, 'Data', data, ...
                    'ColumnName', {'Field', 'Value'}, ...
                    'ColumnWidth', {150, 320}, ...
                    'Position', [10 10 480 280], ...
                    'RowName', []);

        catch ME
            errordlg(['Failed to load .txt file: ' ME.message], 'Error');
        end
    end


    function browse_specs()
        [file, path] = uigetfile('*.txt', 'Select specs');
        if isequal(file, 0)
            return;
        else
            set(specsEdit, 'String', fullfile(path, file));
        end
    end

    function preview_specs_file()
        filepath = get(specsEdit, 'String');
        if ~isfile(filepath)
            errordlg('The specified specs_srtm.txt file does not exist.', 'File Not Found');
            return;
        end
    
        try
            fileContent = fileread(filepath);
            lines = strsplit(fileContent, newline);
            lines = lines(~cellfun('isempty', lines)); % Remove empty lines
    
            data = {};
            for i = 1:numel(lines)
                line = strtrim(lines{i});
                parts = strsplit(line, ':');
                if numel(parts) >= 2
                    key = strtrim(parts{1});
                    value = strtrim(strjoin(parts(2:end), ':')); % handle colons inside value
                    data(end+1, :) = {key, value};
                end
            end
    
            % Create a figure window
            previewFig = figure('Name', 'Specs File Preview', ...
                                'Position', [620 420 450 300], ...
                                'Resize', 'off', ...
                                'NumberTitle', 'off', ...
                                'MenuBar', 'none');
    
            uitable(previewFig, 'Data', data, ...
                    'ColumnName', {'Field', 'Value'}, ...
                    'ColumnWidth', {150, 270}, ...
                    'Position', [10 10 430 280], ...
                    'RowName', []);
        catch ME
            errordlg(['Failed to read specs file: ' ME.message], 'Error');
        end
    end

    function visualize_nii_files()
        [files, path] = uigetfile('*.nii', 'Select NIfTI file(s) to view', 'MultiSelect', 'on');
        if isequal(files, 0)
            return;
        end
        if ischar(files)
            files = {files}; % ensure it's a cell array
        end
        filepaths = fullfile(path, files);
        
        spm('defaults', 'PET');
        
        spm_image('init', char(filepaths)); % Launch spm_image with selected files
    end


    function run_all()

        data_parent_directory = get(dirEdit,'String');
        sb_raw = get(sbEdit, 'String');
        tracerList = get(tracerPopup, 'String');
        pet_tr = tracerList{get(tracerPopup, 'Value')};
        workfolder = get(workfolderEdit, 'String');
        modeling_options_file = get(optionsEdit, 'String');
        specs_file = get(specsEdit, 'String');
        sessions = [];
        if get(session1Check, 'Value'), sessions = [sessions 1]; end
        if get(session2Check, 'Value'), sessions = [sessions 2]; end

        if isempty(sessions)
            errordlg('Please select at least one session.', 'Input Error');
            return;
        end

        % Clean tracer string
        if strcmp(pet_tr,'[11c]nnc'), pet_tr = 'nnc'; end
        if strcmp(pet_tr,'[11c]raclopride'), pet_tr = 'rac'; end

        % Parse subject list
        subjects = strtrim(strsplit(sb_raw, ','));


        try
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(subjects) == 1
                % === Run directly in MATLAB ===
                sb = subjects{1};
                run_magia_with_params(data_parent_directory, sb, pet_tr, workfolder, sessions, modeling_options_file, specs_file);
                msgbox('MAGIA run completed successfully.', 'Success');
                
            else
                % === Submit multiple sbatch jobs ===
                for i = 1:length(subjects)
                    sb = subjects{i};


                    cmd = sprintf('sbatch --output="%s/log_error_files/%s_%s.out" --error="%s/log_error_files/%s_%s.err" run_magia_parallel.sh "%s" "%s" "%s" "%s" "%s" "%s" "%s"', ...
                        workfolder, sb, pet_tr, workfolder, sb, pet_tr, ...
                        data_parent_directory, sb, pet_tr, workfolder, mat2str(sessions), modeling_options_file, specs_file);

                    % cmd = sprintf('sbatch run_magia_parallel.sh "%s" "%s" "%s" "%s" "%s" "%s" "%s"', data_parent_directory, sb, pet_tr, workfolder, mat2str(sessions), modeling_options_file, specs_file);
                    [status, result] = system(cmd);
                    if status == 0
                        fprintf('Submitted job for %s: %s\n', sb, strtrim(result));
                    else
                        warning('Failed to submit job for %s: %s\n', sb, result);
                    end
                end
                msgbox('All MAGIA jobs submitted via sbatch.', 'Submitted');
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        catch ME
            errordlg(['Error running MAGIA: ', ME.message], 'Execution Failed');
        end
    end




end
