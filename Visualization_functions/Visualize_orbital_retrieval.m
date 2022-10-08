% Emission Model
% Divya Kumawat, 09/2022
%% Visualizing Orbital retrievals

function Visualize_orbital_retrieval(Data,SMAP_final_retrieved,mask_snowcover_moistsoil,mask_wet_snow_retrievable_area)

SMAP = Data.SMAP;
Mapped_SNODAS = Data.SNODAS;

numberOfAttempts = 5;
attempt = 0;
info = [];
serverURL = 'http://basemap.nationalmap.gov/ArcGIS/services/USGSImageryOnly/MapServer/WMSServer?';
while(isempty(info))
    try
        info = wmsinfo(serverURL);
        orthoLayer = info.Layer(1);
    catch e

        attempt = attempt + 1;
        if attempt > numberOfAttempts
            throw(e);
        else
            fprintf('Attempting to connect to server:\n"%s"\n', serverURL)
        end
    end
end
% You can change it to 2048 or 4096 for a better quality of
% topographic map.
imageLength = 4096;
fs = 14;

figure('color','w');
lw = 1.5;
ec = [115 114 113]./255;
set(gcf,'Position',[325 175 1125 625])

S = shaperead('./Visualization_functions/USA_shape_file/cb_2018_us_nation_5m.shp','UseGeoCoords', true);
ltlim = [25 49] + [-1 1];
lnlim = [-125 -52] + [-1 1];

%figure 1 (background)-------------------------

latlim = [25 50];
lonlim = [-124 -65];
[A,R] = wmsread(orthoLayer,'Latlim',latlim,'Lonlim',lonlim,...
    'ImageHeight',imageLength,'ImageWidth',imageLength);
US_states = shaperead('usastatelo','UseGeoCoords', true);


ax1 = usamap('conus');
set(ax1,'Position',[-0.0078125,0.13542,0.5,0.73177])
setm(ax1,'MapProjection','mercator')
geoimg1 = geoshow(A,R);
tightmap; plabel('off'); mlabel('off');
alpha 0.6

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

%figure 2 (soil moisture)-------------------

cmap1  = cmapgen([253 166 17; 254 252 16; 152 202 14; 46 142 22;17 255 255;0 0 230; 60 60 86],[30 30 20 15 20 30 30]);


ax2 = usamap('conus');
set(ax2,'Position',[0.5,0.13542,0.5,0.73177])
setm(ax2,'MapProjection','mercator')
geoimg4 = geoshow(SMAP.lat_preserved,SMAP.lon_preserved, SMAP_final_retrieved.sm_r,'DisplayType','texturemap');        % to show coarseness
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off
geoimg4.AlphaDataMapping = 'none';
geoimg4.FaceAlpha = 'texturemap';
alpha(geoimg4,double((~isnan(SMAP_final_retrieved.sm_r))));
caxis(gca,[0 0.62]);
colormap(gcf,(cmap1./255));
set(ax1,'Position',[0,-0.04,0.4,0.8])
set(ax2,'Position',[0,-0.04,0.4,0.8],'visible','off')
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);


%figure 1 (background)-------------------------

ax1 = usamap('conus');
set(ax1,'Position',[-0.0078125,0.13542,0.5,0.73177])
setm(ax1,'MapProjection','mercator')
geoimg1 = geoshow(A,R);
tightmap; plabel('off'); mlabel('off');
alpha 0.6

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', 'k', 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox);
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

%figure 2 (soil moisture)-------------------

ax2 = usamap('conus');
set(ax2,'Position',[0.5,0.13542,0.5,0.73177])
setm(ax2,'MapProjection','mercator')
geoimg4 = geoshow(SMAP.lat,SMAP.lon, SMAP.soil_moist,'DisplayType','texturemap');        % to show coarseness
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus);
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off
geoimg4.AlphaDataMapping = 'none';
geoimg4.FaceAlpha = 'texturemap';
alpha(geoimg4,double((~isnan(SMAP.soil_moist))));
caxis(gca,[0 0.62]);
colormap(ax2,(cmap1./255));

ax3 = usamap('conus');
set(ax3,'Position',[0.3,0.13542,0.5,0.73177])
setm(ax3,'MapProjection','mercator')
geoshow('worldlakes.shp','FaceColor','W')
geoimg3 = geoshow(Mapped_SNODAS.Lat,Mapped_SNODAS.Lon,double(mask_snowcover_moistsoil),'DisplayType','texturemap');
tightmap;plabel('off'); mlabel('off');
geoimg3.AlphaDataMapping = 'none';
geoimg3.FaceAlpha = 'texturemap';
alpha(geoimg3,double((mask_snowcover_moistsoil)));
colormap(ax3,[191 247 245]./255) % Light Blue

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

ax4 = usamap('conus');
set(ax4,'Position',[0.3,0.13542,0.5,0.73177])
setm(ax4,'MapProjection','mercator')
geoshow('worldlakes.shp','FaceColor','W')
geoimg5 = geoshow(Mapped_SNODAS.Lat,Mapped_SNODAS.Lon,double((mask_wet_snow_retrievable_area)),'DisplayType','texturemap');
tightmap;plabel('off'); mlabel('off');
geoimg5.AlphaDataMapping = 'none';
geoimg5.FaceAlpha = 'texturemap';
alpha(geoimg5,double((mask_wet_snow_retrievable_area)));
colormap(ax4,[107 174 214]./255) % Dark blue Color

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

h2 = colorbar(ax2);
h2.Location = 'southoutside';
h2.TickLabelInterpreter = 'Latex';
h2.Label.Interpreter = 'Latex';
h2.Label.String = '$\theta$ [$\rm m^{3}.m^{-3}$]';
h2.Box = 'On';
h2.Label.FontSize = fs;
h2.FontSize = fs;
h2.Position = [0.0132,0.16,0.373,0.0144];

set(ax1,'Position',[0,0.35,0.4,0.8])
set(ax2,'Position',[0,0.35,0.4,0.8],'visible','off')
set(ax3,'Position',[0,0.35,0.4,0.8],'visible','off')
set(ax4,'Position',[0,0.35,0.4,0.8],'visible','off')
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);

% FINAL TAU------------------
lw = 1.5;
ec = [115 114 113]./255;

S = shaperead('./Visualization_functions/USA_shape_file/cb_2018_us_nation_5m.shp','UseGeoCoords', true);
ltlim = [25 49] + [-1 1];
lnlim = [-125 -52] + [-1 1];

%figure 1 (background)-------------------------

latlim = [25 50];
lonlim = [-124 -65];
[A,R] = wmsread(orthoLayer,'Latlim',latlim,'Lonlim',lonlim,...
    'ImageHeight',imageLength,'ImageWidth',imageLength);
US_states = shaperead('usastatelo','UseGeoCoords', true);


ax1 = usamap('conus');
set(ax1,'Position',[-0.0078125,0.13542,0.5,0.73177])
setm(ax1,'MapProjection','mercator')
geoimg1 = geoshow(A,R);
tightmap; plabel('off'); mlabel('off');
alpha 0.6

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

%figure 2 (tau)-------------------

cmap2 = cmapgen([255 255 255;... % white
    155 110 90;...  % brown
    190 220 120;... % light green
    0 50 0],...     % dark green
    [20,40,180,120]);


ax2 = usamap('conus');
set(ax2,'Position',[0.5,0.13542,0.5,0.73177])
setm(ax2,'MapProjection','mercator')
geoimg4 = geoshow(SMAP.lat,SMAP.lon, SMAP_final_retrieved.tau_r,'DisplayType','texturemap');
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off
geoimg4.AlphaDataMapping = 'none';
geoimg4.FaceAlpha = 'texturemap';
alpha(geoimg4,double((~isnan(SMAP_final_retrieved.tau_r))));
caxis(gca,[0 1.5]);
colormap(ax2,(cmap2./255));

set(ax1,'Position',[0.38,-0.04,0.4,0.8])
set(ax2,'Position',[0.38,-0.04,0.4,0.8],'visible','off')
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);

%figure 1 (background)-------------------------

latlim = [25 50];
lonlim = [-124 -65];
[A,R] = wmsread(orthoLayer,'Latlim',latlim,'Lonlim',lonlim,...
    'ImageHeight',imageLength,'ImageWidth',imageLength);
US_states = shaperead('usastatelo','UseGeoCoords', true);


ax1 = usamap('conus');
set(ax1,'Position',[-0.0078125,0.13542,0.5,0.73177])
setm(ax1,'MapProjection','mercator')
geoimg1 = geoshow(A,R);
tightmap; plabel('off'); mlabel('off');
alpha 0.6

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', 'k', 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

%figure 2 (soil moisture)-------------------
ax2 = usamap('conus');
set(ax2,'Position',[0.5,0.13542,0.5,0.73177])
setm(ax2,'MapProjection','mercator')
geoimg4 = geoshow(SMAP.lat,SMAP.lon, SMAP.tau,'DisplayType','texturemap');        
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off
geoimg4.AlphaDataMapping = 'none';
geoimg4.FaceAlpha = 'texturemap';
alpha(geoimg4,double((~isnan(SMAP.soil_moist))));
caxis(gca,[0 1.5]);
colormap(ax2,(cmap2./255));

ax3 = usamap('conus');
set(ax3,'Position',[0.3,0.13542,0.5,0.73177])
setm(ax3,'MapProjection','mercator')
geoshow('worldlakes.shp','FaceColor','W')
geoimg3 = geoshow(Mapped_SNODAS.Lat,Mapped_SNODAS.Lon,double(mask_snowcover_moistsoil),'DisplayType','texturemap');
tightmap;plabel('off'); mlabel('off');
geoimg3.AlphaDataMapping = 'none';
geoimg3.FaceAlpha = 'texturemap';
alpha(geoimg3,double((mask_snowcover_moistsoil)));
colormap(ax3,[191 247 245]./255) % Light Blue

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

ax4 = usamap('conus');
set(ax4,'Position',[0.3,0.13542,0.5,0.73177])
setm(ax4,'MapProjection','mercator')
geoshow('worldlakes.shp','FaceColor','W')

geoimg5 = geoshow(Mapped_SNODAS.Lat,Mapped_SNODAS.Lon,double((mask_wet_snow_retrievable_area)),'DisplayType','texturemap');
tightmap;plabel('off'); mlabel('off');
geoimg5.AlphaDataMapping = 'none';
geoimg5.FaceAlpha = 'texturemap';
alpha(geoimg5,double((mask_wet_snow_retrievable_area)));
colormap(ax4,[107 174 214]./255)

g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);
geoshow(US_states,'edgecolor', ec, 'facecolor', 'none');
ltbox = ltlim([1 2 2 1 1]);
lnbox = lnlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 1);
[xbox, ybox] = mfwdtran(ltbox, lnbox); 
[xaus, yaus] = mfwdtran(S.Lat, S.Lon); 
[xmask, ymask] = polybool('-', xbox, ybox, xaus, yaus); 
[f,v] = poly2fv(xmask, ymask);
hmask = patch(gca,'faces', f, 'vertices', v, ...
    'facecolor', 'w', 'edgecolor', 'none');
tightmap; plabel('off'); mlabel('off');
framem off
gridm off

h2 = colorbar(ax2);
h2.Location = 'southoutside';
h2.TickLabelInterpreter = 'Latex';
h2.Label.Interpreter = 'Latex';
h2.Label.String = '$\tau$ [-]';
h2.Box = 'On';
h2.Label.FontSize = fs;
h2.FontSize = fs;
h2.Position = [0.393,0.16,0.3751,0.0144];
set(ax1,'Position',[0.38,0.35,0.4,0.8])
set(ax2,'Position',[0.38,0.35,0.4,0.8],'visible','off')
set(ax3,'Position',[0.38,0.35,0.4,0.8],'visible','off')
set(ax4,'Position',[0.38,0.35,0.4,0.8],'visible','off')
g = geoshow(S, 'edgecolor', 'k', 'facecolor', 'none','LineWidth',lw);

annotation(gcf,'rectangle',...
    [0.0143333333333333 0.5712 0.373222222222222 0.3728],'LineWidth',1.5);
annotation(gcf,'rectangle',...
    [0.393666666666667 0.571200000000001 0.373222222222222 0.3728],...
    'LineWidth',1.5);
annotation(gcf,'rectangle',...
    [0.0143333333333333 0.1756 0.373222222222222 0.3728],'LineWidth',1.5);
annotation(gcf,'rectangle',...
    [0.394555555555556 0.1756 0.373222222222222 0.3728],'LineWidth',1.5,...
    'FaceAlpha',1.5);

annotation(gcf,'textbox',...
    [0.346777777777778 0.571520001068116 0.053926063043096 0.062079998931885],...
    'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.724777777777779 0.568320001068115 0.0551112521489419 0.062079998931885],...
    'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.351222222222222 0.171520001068116 0.052740841724252 0.062079998931885],...
    'String',{'(c)'},...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.726555555555556 0.169920001068116 0.055111252148942 0.062079998931885],...
    'String',{'(d)'},...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.395000000000001 0.167829200272883 0.167704056139335 0.102570799727118],...
    'String',['bias = 0.012',sprintf('\n'),'RMSE = 0.122'],...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.393222222222223 0.555029200272883 0.167704056139335 0.102570799727118],...
    'String',['bias = 0.019',sprintf('\n'),'RMSE = 0.154'],...
    'Interpreter','latex',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'EdgeColor','none');
end