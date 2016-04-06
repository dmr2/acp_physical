function [TAS,sTAS,years,modellist]=readTAS(scen,smoothwin,baseyears,years,subdir);

%  [TAS,sTAS,years,modellist]=readTAS(scen,smoothwin,baseyears,subdir)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Feb 12 00:34:46 EST 2014

defval('scen','rcp85');
defval('smoothwin',19);
defval('subdir',fullfile('IFILES/atm/global_tas'));
defval('years',1860:2099);
defval('baseyears',[1860 1900]);

% read TAS

clear TAS;
clear sTAS;

pd=pwd;
cd(subdir);
cd(scen);
files=dir('*');
modellist={};
jj=1;
for ii=1:length(files)
	if (files(ii).isdir) && (files(ii).name(1)~='.')
		if exist([files(ii).name '/global_mean_tas'],'dir')
			disp(files(ii).name);
			cd(files(ii).name);
			cd global_mean_tas;
			files2=dir('*.txt');
			if length(files2)>0
				dat=importdata(files2(1).name); dat=dat.data;
				modellist{jj}=lower(files(ii).name);
				TAS(:,jj)=interp1(dat(:,1),dat(:,2),years);
				jj=jj+1;
			end
			cd ..
		end
		cd ..
	end
end
cd(pd);
sub=find((years<=baseyears(end)).*(years>=baseyears(1)));

for jj=1:size(TAS,2);
	TAS(:,jj)=TAS(:,jj)-mean(TAS(sub,jj));
end

TAS(end+1,:) = TAS(end,:) + (TAS(end,:)-TAS(end-1,:));
years(end+1) = years(end)+(years(end)-years(end-1));

sTAS = NaN*TAS;

for jj=1:size(TAS,2);
	sub1 = find(~isnan(TAS(:,jj)));
	sTAS(sub1,jj) = smooth(TAS(sub1,jj),smoothwin);	
end
