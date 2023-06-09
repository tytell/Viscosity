function export_all_swim_kinematics

filenames_fromtable = ["lamprey5-1", "lamprey5-1", "lamprey5-2", "lamprey5-3", ...
    "lamprey5-4", "lamprey5-5", "lamprey5-5", "lamprey5-6", "lamprey5-7", ...
    "lamprey5-8", "lamprey5-9", "lamprey5-10", "lamprey5-11", ...
    "lamprey5-12", "lamprey5-13", "lamprey5-14", "lamprey5-15", "lamprey5-16", ...
    "lamprey5-17", "lamprey5-18", "lamprey5-19", "lamprey5-20", ...
    "lamprey5-21", "lamprey5-22", "lamprey5-23", "lamprey5-24", "lamprey5-25", ...
    "lamprey5-27", "lamprey5-28", "lamprey5-29", "lamprey5-30", "lamprey5-31", ...
    "lamprey5-31", "lamprey5-32", "lamprey5-32", "lamprey5-33", ...
    "lamprey5-33", "lamprey5-34", "lamprey5-34", "lamprey5-35", "lamprey5-35",...
    "lamprey5-35", "lamprey5-36", "lamprey5-36", "lamprey5-37", ...
    "lamprey5-37", "lamprey5-38", "lamprey5-39", "lamprey5-40", "lamprey5-41", ...
    "lamprey6-2", "lamprey6-3", "lamprey6-4", "lamprey6-5", ...
    "lamprey6-6", "lamprey6-7", "lamprey6-8", "lamprey6-9", "lamprey6-10", ...
    "lamprey6-11", "lamprey6-12", "lamprey6-13", "lamprey6-14", ...
    "lamprey6-15", "lamprey6-16", "lamprey6-17", "lamprey6-18", "lamprey6-19", ...
    "lamprey6-20", "lamprey6-21", "lamprey6-22", "lamprey6-23", ...
    "lamprey6-24", "lamprey6-25", "lamprey6-26", "lamprey6-27", "lamprey6-28", ...
    "lamprey6-28", "lamprey6-29", "lamprey6-29", "lamprey6-29", ...
    "lamprey6-30", "lamprey6-31", "lamprey6-32", "lamprey6-33", "lamprey6-33", ...
    "lamprey6-34", "lamprey6-35", "lamprey6-36", "lamprey6-37", ...
    "lamprey6-37", "lamprey6-37", "lamprey6-38", "lamprey6-38", "lamprey6-38", ...
    "lamprey6-39", "lamprey6-40", "lamprey6-41", "lamprey6-42", ...
    "lamprey6-43", "lamprey7-1", "lamprey7-2", "lamprey7-3", "lamprey7-4", ...
    "lamprey7-5", "lamprey7-6", "lamprey7-7", "lamprey7-8", ...
    "lamprey7-9", "lamprey7-10", "lamprey7-11", "lamprey7-12", "lamprey7-13", ...
    "lamprey7-14", "lamprey7-15", "lamprey7-16", "lamprey7-17", ...
    "lamprey7-18", "lamprey7-19", "lamprey7-20", "lamprey7-21", "lamprey7-22", ...
    "lamprey7-23", "lamprey7-24", "lamprey7-25", "lamprey7-26", ...
    "lamprey7-27", "lamprey7-28", "lamprey7-29", "lamprey7-30", "lamprey7-31", ...
    "lamprey7-32", "lamprey7-32", "lamprey7-33", "lamprey7-34", ...
    "lamprey7-35", "lamprey7-36", "lamprey7-37", "lamprey7-38", "lamprey7-39", ...
    "lamprey7-40", "lamprey7-41", "lamprey7-42", "lamprey7-42", ...
    "lamprey7-43", "lamprey8-1", "lamprey8-2", "lamprey8-3", "lamprey8-4", ...
    "lamprey8-5", "lamprey8-6", "lamprey8-7", "lamprey8-8", ...
    "lamprey8-9", "lamprey8-10", "lamprey8-11", "lamprey8-12", "lamprey8-13", ...
    "lamprey8-14", "lamprey8-15", "lamprey8-16", "lamprey8-17", ...
    "lamprey8-18", "lamprey8-19", "lamprey8-20", "lamprey8-21", "lamprey8-22", ...
    "lamprey8-23", "lamprey8-24", "lamprey8-25", "lamprey8-26", ...
    "lamprey8-27", "lamprey8-28", "lamprey8-29", "lamprey8-30", "lamprey8-30", ...
    "lamprey8-31", "lamprey8-32", "lamprey8-32", "lamprey8-33", ...
    "lamprey8-33", "lamprey8-34", "lamprey8-35", "lamprey8-36", "lamprey8-37", ...
    "lamprey8-38", "lamprey8-38", "lamprey8-39", "lamprey8-40", ...
    "lamprey8-41", "lamprey8-42", "lamprey8-43", "lamprey8-43", "lamprey9-1", ...
    "lamprey9-2", "lamprey9-4", "lamprey9-5", "lamprey9-6", ...
    "lamprey9-7", "lamprey9-7", "lamprey9-8", "lamprey9-9", "lamprey9-10", ...
    "lamprey9-11", "lamprey9-12", "lamprey9-13", "lamprey9-14", ...
    "lamprey9-15", "lamprey9-16", "lamprey9-17", "lamprey9-18", "lamprey9-19", ...
    "lamprey9-20", "lamprey9-21", "lamprey9-22", "lamprey9-23", ...
    "lamprey9-24", "lamprey9-25", "lamprey9-26", "lamprey9-27", "lamprey9-28", ...
    "lamprey9-29", "lamprey9-30", "lamprey9-31", "lamprey9-32", ...
    "lamprey9-33", "lamprey9-34", "lamprey9-35", "lamprey9-36", "lamprey9-37", ...
    "lamprey9-38", "lamprey9-39", "lamprey9-40", "lamprey9-41", ...
    "lamprey9-42", "lamprey9-43"];

datadir = '/Users/etytel01/Documents/2022/Viscosity/rawdata2';

filepathsall = getfilenames(datadir, 'recursive', 'include',{'.*\.mat'});

tok = regexp(filepathsall, '(lamprey\d-\d+\w*)(-?corr)?', 'once','tokens');
good = true(length(filepathsall), 1);

basename = cell(size(filepathsall));
corrtxt = cell(size(filepathsall));

for i = 1:length(tok)
    basename{i} = '';
    corrtxt{i} = '';

    if isempty(tok{i})
        good(i) = false;
        continue
    end
    
    if ~contains(tok{i}{1}, filenames_fromtable)
        good(i) = false;
        continue;
    end

    basename{i} = tok{i}{1};
    if ~isempty(tok{i}{2})
        corrtxt{i} = tok{i}{2};
    end
end

[basenames, ~, ind] = unique(basename);
for i = 1:length(basenames)
    isind = ind == i;
    if sum(isind) > 1
        good(isind) = corrtxt(isind) == "-corr";
    end
end

filenames = filepathsall(good);

progress(0, length(filenames), 'Processing!');
for i = 1:length(filenames)
    [pn, fn, ext] = fileparts(filenames{i});
    outfilename = fullfile(pn, fn + "-midline.csv");

    try
        export_swim_kinematics(filenames{i}, outfilename);
    catch
        warning('Could not process file %s', fn);
    end

    progress(i);
end

