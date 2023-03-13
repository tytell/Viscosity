function correct_viscosity_kinematics

if (~getvar('filenames','corrected'))
    filenames = getfilenames('rawdata','recursive','include',{'*.mat'});
    corrected = 1;
end

for i = corrected:length(filenames)
    correct_kinematics(filenames{i});
    corrected = i;
    putvar filenames corrected;
end
