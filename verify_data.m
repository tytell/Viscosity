function verify_data(files)

% files = getfilenames('rawdata/*','recursive');

good = true(length(files),1);
i = 2;
while i < length(files)
    [pn,fn,ext] = fileparts(files{i});
    if strcmp(files{i-1},fullfile(pn,[fn '-corr.mat']))
        good(i) = false;
        i = i+2;
    else
        i = i+1;
    end
end

h(1) = subplot(2,1,1);
h(2) = subplot(2,4,5);
h(3) = subplot(2,4,6);
h(4) = subplot(2,4,7);
h(5) = subplot(2,4,8);

files = files(good);

good = true(size(files));
for i = 1:length(files)
    fprintf('** %s\n', files{i});
    load(files{i},'indpeak', 'per','wavelen','wavevel','amp','fishlenmm');
    
    plot(h(1), 1:20,indpeak);
    
    plot(h(2), per(end,:), 'o');
    ylabel(h(2),'Period');
    
    plot(h(3), amp([1 end],:)'/fishlenmm, 'o');
    ylabel(h(3),'Amplitude');
    
    plot(h(4), wavelen(end,:)/fishlenmm, 'o');
    ylabel(h(4),'Wavelen');
    
    plot(h(5), wavevel/fishlenmm, 'o');
    ylabel(h(5),'Wavevel');
    
    good(i) = inputyn('Good? ','default',true);
end
    