function [patternwts,Tctr,Tedges,targetquantcenters,sTAS,TASyrs]=WeightModelsByMAGICCSMME(MAGICCyears,Mproj,scen,doplot,TASbaseyears,calibyear,psyears,targetquantedges,subdir)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Feb 9 12:02:18 EST 2014

defval('doplot',1);
defval('scen','rcp85');
defval('subdir','');

defval('targetquantedges',[.0001 .08 .12 .20 .40 .60 .80 .88 .92 .98 .9999]);
defval('psyears',2000:2100);
defval('calibyear',2090);
defval('TASbaseyears',[1981 2010]);
if length(TASbaseyears)==1
	TASbaseyears=[1 1]*TASbaseyears;
end
defval('MAGICCbaseyears',TASbaseyears);
defval('smoothwin',19);
defval('npatternperbin',2);

sub=find((MAGICCyears>=MAGICCbaseyears(1)).*(MAGICCyears<=MAGICCbaseyears(2)));
reftemp=mean(Mproj(sub,:),1);
MAGICCproj=bsxfun(@minus,Mproj,reftemp);

if smoothwin>1
	for i=1:size(MAGICCproj,2)
		MAGICCproj(:,i) = smooth(MAGICCproj(:,i),smoothwin);
	end
end


calibrow=find(MAGICCyears==calibyear);

[TAS,sTAS,TASyrs,TASmodellist]=readTAS(scen,smoothwin,TASbaseyears,[MAGICCyears(1):psyears(end)],subdir);
modelT = sTAS(find(TASyrs==calibyear),:);
sub=find(~isnan(modelT)); modelT=modelT(sub); sTAS=sTAS(:,sub); TAS=TAS(:,sub); TASmodellist=TASmodellist(sub);

%%%

testtargquants = [.001 .05:.05:.95 .999];
y = quantile(MAGICCproj(calibrow,:),testtargquants);

targetquantcenters = (targetquantedges(1:end-1)+targetquantedges(2:end))/2;
targetquantwts = diff(targetquantedges);
Tedges =  quantile(MAGICCproj(calibrow,:),targetquantedges);
Tctr =  quantile(MAGICCproj(calibrow,:),targetquantcenters);; 
[nmodels,modelbin]=histc(modelT,Tedges); nmodels=nmodels(1:end-1);
npatterns = (npatternperbin-nmodels) .* (nmodels<npatternperbin);

fullqlevs=[.001 .01:.01:.99 .999];
fullqvals = quantile(MAGICCproj(calibrow,:),fullqlevs);
mdlquants=interp1(fullqvals,fullqlevs,modelT);

patternquant=[];
patternT=[];
patternmodel=[];

for i=1:length(npatterns)
	patternquant=[patternquant repmat(targetquantcenters(i),1,npatterns(i))];
	patternT=[patternT repmat(Tctr(i),1,npatterns(i))];
	patternmodel=[patternmodel repmat(0,1,npatterns(i))];
	sub=find(modelbin==i);
	for jj=1:length(sub)
		patternmodel=[patternmodel sub(jj)];
		patternquant=[patternquant mdlquants(sub(jj))];
		patternT = [patternT modelT(sub(jj))];
	end
end

[s,si]=sort(patternT);
patternT=patternT(si); patternquant=patternquant(si); patternmodel=patternmodel(si);

[neff,patternbin]=histc(patternT,Tedges); neff=neff(1:end-1);
patternbinwts=targetquantwts./neff*sum(neff);
patternbinwts(find(~isfinite(patternbinwts)))=0;
patternwts = patternbinwts(patternbin);

clear patternTpathway patternscaled;
patlist={};
for i=1:length(patternquant)
	if patternmodel(i)==0
		patternTpathway(i,:) = quantile(MAGICCproj,patternquant(i),2);
		patlist={patlist{:},['pattern' num2str(i) '_' sprintf('%0.3f',100*patternquant(i)) '*']};
		patternscaled(i)=1;
	else
		mdlpattern=interp1(TASyrs,sTAS(:,patternmodel(i)),MAGICCyears);
		sub2=find(isnan(mdlpattern));
		[jk,jki]=intersect(MAGICCyears(sub2),psyears);
		if length(jk)>0
			patternscaled(i)=1;
		else
			patternscaled(i)=0;
		end
		mdlpattern(sub2)=quantile(MAGICCproj(sub2,:),patternquant(i),2);
		patternTpathway(i,:)=mdlpattern;
		if patternscaled(i)
			patlist={patlist{:},[TASmodellist{patternmodel(i)} '*']};
		else
			patlist={patlist{:},[TASmodellist{patternmodel(i)}]};
		end
	end	
end	

if doplot
	clf;
	subplot(2,1,1);
	plot(y,testtargquants,'r','linew',3); hold on;
	plot(Tedges,[0 cumsum(nmodels)/sum(nmodels)],'linew',2); hold on;
	%plot([0 Tctr],[0 cumsum(targetquantwts)],'r','linew',2);
	%plot(Tedges,[0 cumsum(nmodels.*modelbinwts)/sum(modelwts)],'g','linew',2); hold on;
	plot(Tedges,[0 cumsum(neff.*patternbinwts)/(sum(patternwts))],'c','linew',2); 
	plot(modelT,1.05*ones(size(modelT)),'bo','MarkerFaceColor','b');
	sub=find(patternmodel==0);
	plot(patternT(sub),1.05*ones(size(sub)),'kd','MarkerFaceColor','k');
	ylim([0 1.1]); %xlim([2 11]);
	xlabel(['Temperature in ' num2str(calibyear) ' (C)']); ylabel('CDF');
	title(scen);
	legend('MAGICC','CMIP5','CMIP5 weighted','models','patterns',	'Location','Southeast');

%	dTedges=diff(Tedges);
%	subplot(2,1,2);
%	plot((y(1:end-1)+y(2:end))/2,diff(testtargquants)./diff(y),'r','linew',3); hold on;
%	plot([Tctr],diff([0 cumsum(nmodels)/sum(nmodels)])./dTedges,'linew',2); hold on;
%	%plot([Tctr],diff([0 cumsum(targetquantwts)])./dTedges,'r','linew',2);
%%	plot([Tctr],diff([0 cumsum(nmodels.*modelbinwts)/sum(modelwts)])./dTedges,'g','linew',2); hold on;
%	plot([Tctr],diff([0 cumsum(neff.*patternbinwts)/(sum(patternwts))])./dTedges,'c','linew',2); hold on;
%	%xlim([2 11]);
%	xlabel(['Temperature in ' num2str(calibyear) ' (C)']); ylabel('PDF');

	pdfwrite([scen '_' num2str(calibyear) '_SMME_modelweighting']);

	[jk,ia,ib]=intersect(psyears,MAGICCyears);
	fid=fopen([scen '_' num2str(calibyear) '_SMME.tsv'],'w');
	fprintf(fid,['quantile\tmodel\tweight']);
	fprintf(fid,'\t%0.0f',jk);
	fprintf(fid,'\n');
	for i=1:length(patlist)
		fprintf(fid,'%0.2f',patternquant(i));
		fprintf(fid,['\t' patlist{i}]);
		fprintf(fid,'\t%0.4f',patternwts(i)/sum(patternwts));
		fprintf(fid,'\t%0.2f',patternTpathway(i,ib));
		fprintf(fid,'\n');
	end
	fclose(fid);
	

	clf;
	subplot(2,2,1);
	plot(MAGICCyears,patternTpathway,'Color',[.6 .6 .6]); hold on;
	plot(TASyrs,sTAS,'r'); hold on;
	plot(MAGICCyears,quantile(MAGICCproj,.5,2),'b','linew',3); hold on;
	plot(MAGICCyears,quantile(MAGICCproj,[.167 .833],2),'b'); hold on;
	plot(MAGICCyears,quantile(MAGICCproj,[.05 .95],2),'b--'); hold on;
	plot(MAGICCyears,quantile(MAGICCproj,[.01 .99],2),'b:'); hold on;
	xlim(psyears([1 end]));
	title([scen]);
	ylabel(['C above ' num2str(MAGICCbaseyears(1)) '-' num2str(MAGICCbaseyears(2))]);
	pdfwrite([scen '_' num2str(calibyear) '_SMME_MAGICCproj']);
end


end

%%

