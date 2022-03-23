function show_viscosity_kinematics

rho = 1000;     %kg/m^3
height = 8;     %mm

load viscosity_kinematics;

[~,~,viscind] = unique(visc);

comspeed(comspeed > 3.5) = NaN;

comspeedmn = nanmean(comspeed);
freqmn = nanmean(freq);
tailampmn = nanmean(tailamp);
wavelenmn = nanmean(wavelen);
Stmn = nanmean(freq .* tailamp ./ comspeed);

tailperpspeedmn = NaN(size(comspeedmn));
for i = 1:size(mid.t,3)
    t1 = mid.t(:,:,i);
    tx1 = mid.mx(20,:,i) .* len(i);
    ty1 = mid.my(20,:,i) .* len(i);
    dx1 = mean(diff(mid.mx(17:20,:,i)));
    dy1 = mean(diff(mid.my(17:20,:,i)));
    
    tailvelx1 = deriv(t1,tx1);
    tailvely1 = deriv(t1,ty1);
    
    tailspeed1 = sqrt(tailvelx1.^2 + tailvely1.^2);
    tailperpspeed1 = (-tailvelx1.*dy1 + tailvely1.*dx1) ./ sqrt(dx1.^2 + dy1.^2);
    tailperpspeedmn(i) = nanmean(abs(tailperpspeed1));
    tailspeedmn(i) = nanmean(abs(tailspeed1));
    tailspeedsqmn(i) = nanmean(tailvelx1.^2 + tailvely1.^2);
end

Remn = rho * (tailspeedmn/1000).*(len/1000)./(visc/100);
Cdmn = 1 + sqrt(8./Remn);
dragmn = 0.5 * Cdmn .* rho .* ((len/1000) .* height/1000) .* (tailperpspeedmn/1000).^2;

if (inputyn('Save CSV file?','default',false))
    cols = {'indiv','trial','viscosity','rep','len',...
        'comspeed','tailamp','freq','wavelen','tailperpspeed','St','Re','Cd','drag'};
    
    fid = fopen('viscosity_kinematics.csv','w');
    fprintf(fid,'%s, ',cols{:});
    fprintf(fid,'\n');
    
    fprintf(fid,'%d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
        [indiv; trial; visc; rep; len; ...
        comspeedmn; tailampmn; freqmn; wavelenmn; tailperpspeedmn; ...
        Stmn; Remn; Cdmn; dragmn]);
    fclose(fid);
end

comspeedmnvisc = accumarray(viscind',comspeedmn',[],@nanmean);
freqmnvisc = accumarray(viscind',freqmn',[],@nanmean);
wavelenmnvisc = accumarray(viscind',wavelenmn',[],@nanmean);
tailampmnvisc = accumarray(viscind',tailampmn',[],@nanmean);
Stmnvisc = accumarray(viscind',Stmn',[],@nanmean);
Remnvisc = accumarray(viscind',Remn',[],@nanmean);
Cdmnvisc = accumarray(viscind',Cdmn',[],@nanmean);
dragmnvisc = accumarray(viscind',dragmn',[],@nanmean);

figureseries('Basic kinematics');
clf;

h(1) = subplot(2,2,1);
bar(comspeedmnvisc,'k');
hold on;
plotgroups(viscind,comspeedmn,{indiv},{'mf'},'xoffset',0.1, 'MarkerSize',3);
hold off;
xtick({'1x','6x','20x'});
ylabel({'Swimming speed','(L/s)'});


h(2) = subplot(2,2,2);
bar(tailampmnvisc,'k');
hold on;
plotgroups(viscind,tailampmn,{indiv},{'mf'},'xoffset',0.1, 'MarkerSize',3);
hold off;
xtick({'1x','6x','20x'});
ylabel({'Tail amplitude','(L)'});


h(3) = subplot(2,2,3);
bar(freqmnvisc,'k');
hold on;
plotgroups(viscind,freqmn,{indiv},{'mf'},'xoffset',0.1, 'MarkerSize',3);
hold off;
xtick({'1x','6x','20x'});
ylabel({'Frequency','(Hz)'});


h(4) = subplot(2,2,4);
bar(wavelenmnvisc,'k');
hold on;
plotgroups(viscind,wavelenmn,{indiv},{'mf'},'xoffset',0.1, 'MarkerSize',3);
hold off;
xtick({'1x','6x','20x'});
ylabel({'Wave length','(L)'});

set(h,'Box','off','TickDir','out');

figureseries('Forces');
clf;

hf(1) = subplot(2,2,1);
superboxplot(viscind,log10(Remn),'notch',false,'col','k');
ytick(log10([20 50 100 200 500 1000 2000 5000 10000]),...
     {'20','50','100','200','500','1000','2000','5000','10000'});
ylim(log10([50 15000]));
xtick([1 2 3],{'1x','6x','20x'});
ylabel({'Re'});


hf(2) = subplot(2,2,2);
superboxplot(viscind,Cdmn,'notch',false,'col','k');
xtick([1 2 3],{'1x','6x','20x'});
ylabel({'C_d'});


hf(3) = subplot(2,2,3);
superboxplot(viscind,(tailperpspeedmn),'notch',false,'col','k');
xtick([1 2 3],{'1x','6x','20x'});
ylabel({'Tail speed','(mm/s)'});


hf(4) = subplot(2,2,4);
superboxplot(viscind,dragmn*1000,'notch',false,'col','k');
xtick([1 2 3],{'1x','6x','20x'});
ylabel({'Drag force','(mN)'});
