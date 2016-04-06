% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Jun 03 23:08:58 EDT 2014

% configuration variables

mf=mfilename('fullpath'); mfsl=strfind(mf,'/'); mypath=mf(1:mfsl(end)-1);
addpath(mypath,fullfile(mypath,'../lib/MFILES'),fullfile(mypath,'MFILES'));

IFILES = fullfile(mypath, '../../IFILES/atm/');
subdir=fullfile(IFILES,'global_tas');
scens={'rcp85','rcp60','rcp45','rcp26'};
magiccfiles={'IPCCAR5climsens_rcp85_DAT_SURFACE_TEMP_BO_15Nov2013_185227.OUT','IPCCAR5climsens_rcp6_DAT_SURFACE_TEMP_BO_16Nov2013_064508.OUT','IPCCAR5climsens_rcp45_DAT_SURFACE_TEMP_BO_16Nov2013_070923.OUT','IPCCAR5climsens_rcp3pd_DAT_SURFACE_TEMP_BO_16Nov2013_085858.OUT'};

% Generate Weights

clear Myears Mproj;
for ii=1:length(scens)
	disp(scens{ii});
	mfile=fullfile(IFILES,magiccfiles{ii});

	dat=importdata(mfile,' ',25);
	Myears{ii}=dat.data(:,1); Mproj{ii}=dat.data(:,2:end);

	[patternwts,Tctr,Tedges,targetquantcenters]=WeightModelsByMAGICCSMME(Myears{ii},Mproj{ii},scens{ii},[],[1981 2010],[],1950:2200,[],subdir);

    [Tproj,samppattern,sampresidual,weatherproj,targetquantcenters2]=WeightModelsByMAGICCMCPR(Myears{ii},Mproj{ii},scens{ii},[],[1981 2010],[1950:2200],[],subdir);

end

subBASE=find((Myears{1}>=1981).*(Myears{1}<=2010));
sub30=find((Myears{1}>=2020).*(Myears{1}<=2039));
sub50=find((Myears{1}>=2040).*(Myears{1}<=2059));
sub90=find((Myears{1}>=2080).*(Myears{1}<=2099));
sub100=find(Myears{1}==2100);
sub150=find((Myears{1}>=2140).*(Myears{1}<=2159));
sub200=find(Myears{1}==2200);

for ii=1:length(scens)
    delta{ii}=bsxfun(@minus,Mproj{ii},mean(Mproj{ii}(subBASE,:),1));
end

clear quants50 quants90 quants100 quants150 quants200 quantsall;
quantlevs=[.01 .05 .167 .5 .833 .95 .99];
for ii=1:length(scens)
    quants0(ii,:)=quantile(mean(delta{ii}(subBASE,:),1),quantlevs);
    quants30(ii,:)=quantile(mean(delta{ii}(sub30,:),1),quantlevs);
    quants50(ii,:)=quantile(mean(delta{ii}(sub50,:),1),quantlevs);
    quants90(ii,:)=quantile(mean(delta{ii}(sub90,:),1),quantlevs);
    quants100(ii,:)=quantile(delta{ii}(sub100,:),quantlevs);
    quants150(ii,:)=quantile(mean(delta{ii}(sub150,:),1),quantlevs);
    quants200(ii,:)=quantile(delta{ii}(sub200,:),quantlevs);
    quantsall{ii}=quantile(delta{ii}',quantlevs);
end

fid=fopen('globalT.tsv','w');
fprintf(fid,'global mean T, degrees C above 1981-2010, based on distribution of MAGICC runs\n\n');   
for ppp=1:length(scens)
    fprintf(fid,[scens{ppp} ' (%0.0f runs)\n'],size(Mproj{ppp},2));
    fprintf(fid,'level\t1981-2010\t2020-2039\t2040-2059\t2080-2099\t2100\t2140-2159\t2200\t');
    fprintf(fid,'\t%0.0f',Myears{1});
    for nn=1:length(quantlevs)
        fprintf(fid,'\n');
        fprintf(fid,'%0.1f',100*quantlevs(nn));
        fprintf(fid,'\t%0.2f',quants0(ppp,nn));
        fprintf(fid,'\t%0.2f',quants30(ppp,nn));
        fprintf(fid,'\t%0.2f',quants50(ppp,nn));
        fprintf(fid,'\t%0.2f',quants90(ppp,nn));
        fprintf(fid,'\t%0.2f',quants100(ppp,nn));
        fprintf(fid,'\t%0.2f',quants150(ppp,nn));
        fprintf(fid,'\t%0.2f',quants200(ppp,nn));
        fprintf(fid,'\t');
        fprintf(fid,'\t%0.2f',quantsall{ppp}(nn,:));
    end
    
    fprintf(fid,'\n\n');
end
fclose(fid);


fid=fopen('globalT_F.tsv','w');
fprintf(fid,'global mean T, degrees F above 1981-2010, based on distribution of MAGICC runs\n\n');   
for ppp=1:length(scens)
    fprintf(fid,[scens{ppp} ' (%0.0f runs)\n'],size(Mproj{ppp},2));
    fprintf(fid,'level\t1981-2010\t2020-2039\t2040-2059\t2080-2099\t2100\t2140-2159\t2200\t');
    fprintf(fid,'\t%0.0f',Myears{1});
    for nn=1:length(quantlevs)
        fprintf(fid,'\n');
        fprintf(fid,'%0.1f',100*quantlevs(nn));
        fprintf(fid,'\t%0.2f',quants0(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants30(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants50(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants90(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants100(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants150(ppp,nn)*1.8);
        fprintf(fid,'\t%0.2f',quants200(ppp,nn)*1.8);
        fprintf(fid,'\t');
        fprintf(fid,'\t%0.2f',quantsall{ppp}(nn,:)*1.8);
    end
    
    fprintf(fid,'\n\n');
end
fclose(fid);
