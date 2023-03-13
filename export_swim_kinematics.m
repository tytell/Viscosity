function export_swim_kinematics(infile, outfile)

if nargin == 0
    [fn,pn] = uigetfile('*.mat', 'Choose file');
    if isempty(fn)
        return
    end
    infile = fullfile(pn,fn);
end

load(infile, 't','fr', 'mxmm', 'mymm');

npt = size(mxmm, 1);

pt = repmat((1:npt)', [1 size(mxmm,2)]);
t = repmat(t, [npt 1]);
fr = repmat(fr, [npt 1]);

tab = table(t(:), fr(:), pt(:), mxmm(:), mymm(:));

tab.Properties.VariableNames = {'t', 'frame', 'point', 'mxmm', 'mymm'};

if nargin == 0
    [~, fn, ext] = fileparts(fn);
    outfilename = fullfile(pn, fn + ".csv");
    
    [fn,pn] = uiputfile(outfilename, 'Choose output file');
    if isempty(fn)
        fprintf('Cancelled');
    end
    outfile = fullfile(pn,fn);
end

writetable(tab, outfile);
