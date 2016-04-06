function [Tproj,samppattern,sampresidual,weatherproj,targetquantcenters,sTAS,TASyrs]=WeightModelsByMAGICCMCPR(MAGICCyears,Mproj,scen,doplot,TASbaseyears,psyears,Nsamps,subdir)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Feb 9 12:01:19 EST 2014

defval('doplot',1);
defval('scen','rcp85');
defval('Nsamps',100);
defval('subdir','');

defval('psyears',2000:2100);
defval('TASbaseyears',[1981 2010]);
if length(TASbaseyears)==1
	TASbaseyears=[1 1]*TASbaseyears;
end
defval('MAGICCbaseyears',TASbaseyears);
defval('smoothwin',19);
defval('weatheryears',1981:2010);

sub=find((MAGICCyears>=MAGICCbaseyears(1)).*(MAGICCyears<=MAGICCbaseyears(2)));
reftemp=mean(Mproj(sub,:),1);
MAGICCproj=bsxfun(@minus,Mproj,reftemp);

if smoothwin>1
	for i=1:size(MAGICCproj,2)
		MAGICCproj(:,i) = smooth(MAGICCproj(:,i),smoothwin);
	end
end

[TAS,sTAS,TASyrs,TASmodellist]=readTAS(scen,smoothwin,TASbaseyears,[MAGICCyears(1):psyears(end)],subdir);

%%%

targetquantbounds = linspace(0,1,Nsamps+1);
targetquantcenters=.5*(targetquantbounds(1:end-1)+targetquantbounds(2:end));
[jk,ia,ib]=intersect(MAGICCyears,psyears);
Tproj =  quantile(MAGICCproj(ia,:),targetquantcenters,2);
psyears=psyears(ib);

% patterns and residuals
modelpool = repmat(TASmodellist,1,ceil(Nsamps/length(TASmodellist)));
samppattern=modelpool(randperm(length(modelpool),Nsamps));
sampresidual=modelpool(randperm(length(modelpool),Nsamps));

% weather years
weatherproj=nan(size(Tproj));
[jk,ia,ib]=intersect(weatheryears,psyears);
weatherproj(ib,:)=repmat(jk(:),1,size(Tproj,2));
subneedsweather=find(psyears>weatheryears(end));
Nforweather=length(subneedsweather);
weatheryearpool=repmat(weatheryears,1,ceil(Nforweather/length(weatheryears)));

for ii=1:size(Tproj,2)
    weatherproj(subneedsweather,ii)=weatheryearpool(randperm(length(weatheryearpool),Nforweather));
end

%%%

if doplot

	clf;
	subplot(2,2,1);
	plot(psyears,Tproj); hold on;
	xlim(psyears([1 end]));
	title([scen]);
	ylabel(['C above ' num2str(MAGICCbaseyears(1)) '-' num2str(MAGICCbaseyears(2))]);
	pdfwrite([scen '_MAGICCproj_MCPR']);

	fid=fopen([scen '_MCPR.tsv'],'w');
    fprintf(fid,'quantile\tpattern\tresidual');
    fprintf(fid,'\t%0.0f',psyears);
    fprintf(fid,'\tweather_%0.0f',psyears(subneedsweather));
    fprintf(fid,'\n');
    for ii=1:Nsamps
        fprintf(fid,'%0.3f',targetquantcenters(ii));
        fprintf(fid,['\t' samppattern{ii} '\t' sampresidual{ii}]);
        fprintf(fid,'\t%0.2f',Tproj(:,ii));
        fprintf(fid,'\t%0.0f',weatherproj(subneedsweather,ii));
        fprintf(fid,'\n');
    end
    fclose(fid);

end

%%

