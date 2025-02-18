%   script nam_narr_2swan.m
%
%   inputs-     3hour NAM 12km data grib files and/or
%               3hour NARR 32 km grib files and can  
%               combine and interpolate to a common grid.
%               All data is read thru THREDDs from
%               http://nomads.ncdc.noaa.gov/data.php?name=access
% 
%   user selects:   - to get NAM and/or NARR data.
%                     or read data from a netcdf file.
%                   - time interval [nam_start:nam_end]
%                   - spatial interval [roms grid or generic grid]
%                   - variables to be inlucded
%
%   output-     SWAN ASCII wind forcing file
%   needs-      native matlab to read opendap
%
% 29Sept2014 - jcwarner
%

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) Select which variables to include in this ascii forcing file.
%  put a '1' if you want to include it, '0' otherwise.
get_Wind=1;       % surface u- and v- winds (m/s)

%(2) Enter name of output SWAN forcing file
SWAN_NAMNARR_name='swan_namnarr_1020Sep2012.dat';

%(3) Enter start and end dates - needed if using get_NARR or get_NAM
namnarr_start = datenum('10-Sep-2012');
namnarr_end   = datenum('20-Sep-2012');

%(4) Select to interpolate to a roms grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_swan_grid = 0;
interpto_user_grid = 1;
%
if (interpto_swan_grid)
  model_grid='Sandy_roms_grid.coord';
  nx=84;
  ny=24;
elseif (interpto_user_grid)
  lon_rho=[255:0.1:310]-360;
  lat_rho=[ 10:0.1:50 ];  % Create a 0.1 degree lon-lat grid
else
  disp('pick a grid')
end

%5) Select which data to obtain: NAM, NARR, or both.
get_NARR=1;  %NARR-A grid 221 32km data
get_NAM=1;   %NAM grid 218 12km data
%
% -- or  --
%
% read in data from a netcdf file
get_netcdf_file=0;
if (get_netcdf_file)
  ncfile='roms_namnarr_30Sept10Oct2012.nc';
end
%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%get some grid info
if (interpto_swan_grid)
  grd=textread(model_grid);
  gridsize=length(grd)/2;
  zz=grd(1:gridsize);
  lon_rho=reshape(zz,nx,ny);
  zz=grd(gridsize+1:end);
  lat_rho=reshape(zz,nx,ny);
elseif (interpto_user_grid)
  lon_rho=repmat(lon_rho,length(lat_rho),1)';
  lat_rho=repmat(lat_rho',1,size(lon_rho,1))';
else
  disp('pick a grid')
end
[Lp,Mp]=size(lon_rho);
L=Lp-1;
M=Mp-1;

% now figure out what year they want
NAMNARR_time=[namnarr_start:3/24:namnarr_end];
ntimes=length(NAMNARR_time);
Time=NAMNARR_time-datenum(1858,11,17,0,0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% pre allocate some arrays
%
if (get_Wind)
  Uwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
  Vwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
%
if (get_NARR)
  disp('going to get NARR-A grid 221 32km data');
%
  for mm=1:ntimes
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NARR-A grid 221 32km data at ',dd]);
    url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/',dd(1:6),'/',dd(1:8),'/narr-a_221_',dd(1:8),'_',dd(10:11),'00_000.grb'];
    if (mm==1)
      x=ncread(url,'x');
      y=ncread(url,'y');
      clo=-107.0;   clat=50.0;
      earth_rad=6367.470;
      [X,Y]=meshgrid(x,y);
      m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
      [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
    end
%
    if (get_Wind)
      var=squeeze(ncread(url,'u_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Uwind(:,:,mm)=cff;
%
      var=squeeze(ncread(url,'v_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Vwind(:,:,mm)=cff;
%
      if (~get_NAM)
        if (interpto_swan_grid)    % do this below, dont need it twice
          cffx=Uwind(:,:,count).*cos(angle)+Vwind(:,:,count).*sin(angle);
          cffy=Vwind(:,:,count).*cos(angle)-Uwind(:,:,count).*sin(angle);
          Uwind(:,:,count)=cffx;
          Vwind(:,:,count)=cffy;
        end
      end
    end
  end
  save NARR_data.mat
end
%
if (get_NAM)
  disp('going to get NAM grid 218 12km data');
%
  for mm=1:ntimes
    nstp=mod(mm,8);
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NAM grid 218 12km data at ',dd]);
          %http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/201210/20121028/namanl_218_20121028_0000_000.grb
    %url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',dd(10:11),'00_000.grb'];
    if ismember(nstp,[1 2])
      first='0000';
    elseif ismember(nstp,[3 4])
      first='0600';
    elseif ismember(nstp,[5 6])
      first='1200';
    elseif ismember(nstp,[7 0])
      first='1800';
    end
    if ismember(nstp,[1 3 5 7])
      second='000';
    else
      second='003';
    end
    url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb'];
    try
      if (mm==1)
        x=ncread(url,'x');
        y=ncread(url,'y');
        clo=-95.0;   clat=25.0;
        earth_rad=6367.470;
        [X,Y]=meshgrid(x,y);
        m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
        [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
    %
    % find the indices of the lon_rho lat_rho grid that are inside the NAM
    % data. we will just use these points from NAM and take the rest from NARR.
    %
        disp('computing mask to merge NARR and NAM')
        mask=zeros(size(lon_rho));
        X=[nlon(:,1); nlon(end,:)' ;nlon(end:-1:1,end); nlon(1,end:-1:1)'];
        Y=[nlat(:,1); nlat(end,:)' ;nlat(end:-1:1,end); nlat(1,end:-1:1)'];
        zz=inpolygon(lon_rho,lat_rho,X, Y);
        mask(zz==1)=1;
      end
    %
      if (get_Wind)
        var=squeeze(ncread(url,'u_wind_height_above_ground'));
        var=squeeze(var(:,:,1));
        var=var.';
        F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
        zz=F(lon_rho,lat_rho);
        zz(isnan(zz))=0;
        cff=squeeze(Uwind(:,:,mm)).*(1-mask)+zz.*mask;
        Uwind(:,:,mm)=cff;
    %
        var=squeeze(ncread(url,'v_wind_height_above_ground'));
        var=squeeze(var(:,:,1));
        var=var.';
        F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
        zz=F(lon_rho,lat_rho);
        zz(isnan(zz))=0;
        cff=squeeze(Vwind(:,:,mm)).*(1-mask)+zz.*mask;
        Vwind(:,:,mm)=cff;
      end
    catch ME
      disp(['cldnt get that data at ', url])
    end
  end
end
%
% write data to ascii file
%
if(get_netcdf_file)
  netcdf_load(ncfile);
  Time=wind_time;
  get_Wind=1;
end
%
fid = fopen(SWAN_NAMNARR_name,'w');
for i=1:length(Time)
  if (get_Wind)
    disp(['Writing winds for SWAN at ',datestr(Time(i)+datenum(1858,11,17,0,0,0))])
    uswan=squeeze(Uwind(:,:,i)');
    vswan=squeeze(Vwind(:,:,i)');
    fprintf(fid,'%10.2f\n',uswan');
    fprintf(fid,'%10.2f\n',vswan');
  end
end
fclose(fid);

%
disp(['------------ wrote ',SWAN_NAMNARR_name,' ------------']);


