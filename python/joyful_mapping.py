import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
#from gprm.utils.raster import to_anchor_plate
import pygmt
import xarray as xr
import rioxarray as rio
from gprm.utils.raster import to_anchor_plate
import joyful_geochemistry as joy



def reconstruct(df, reconstruction_model, reconstruction_time, anchor_plate_id=201, valid_time_filter=True):
    
    df = reconstruction_model.assign_plate_ids(df, copy_valid_times=valid_time_filter)
    if valid_time_filter:
        df = df.query('age<=FROMAGE')
    df = df.replace(np.nan,-99999.9)
    df_r = reconstruction_model.reconstruct(df, reconstruction_time=reconstruction_time, anchor_plate_id=anchor_plate_id)
    if df_r is not None:
        df_r = df_r.replace(-99999.9,np.nan)
    
    return df_r


def select_orogens(df, gdf=None, orogen_names=None, continent_names=None, region=None):
    if gdf is None:
        from gprm.datasets import Geology
        gdf = Geology.fetch_GlobalTectonicMap()
    
    if continent_names is not None:
        if isinstance(continent_names, list):
            continent_names = '|'.join(continent_names)
        gdf = gdf[gdf['continent'].str.contains(continent_names, na=False)]
         
    if orogen_names is not None:
        if isinstance(orogen_names, list):
            orogen_names = '|'.join(orogen_names)
        gdf = gdf[gdf['lastorogen'].str.contains(orogen_names, na=False)]

    if region is not None:
        gdf = gdf.cx[region[0]:region[1], region[2]:region[3]]

    return df.clip(gdf)


def load_pilger_volcanics(excelfile=None):

    if excelfile is None:
        excelfile = '/Users/simon/OneDrive/Andes_works//datafiles/2022-2-NGG0Q7_Pilger_Andean-Igneous-Radiometric-Dates.xlsx'
    
    pilger_Igneous = pd.read_excel(excelfile, sheet_name='Igneous-Master')

    pilger_Igneous['Longitude'] = pilger_Igneous['Longitude'].astype(str)
    pilger_Igneous['Longitude'] = pilger_Igneous['Longitude'].str.replace("\)",'', regex=True)
    pilger_Igneous['Longitude'] = pilger_Igneous['Longitude'].astype(np.float64)
    pilger_Igneous['Latitude'] = pilger_Igneous['Latitude'].astype(np.float64) #.str.replace('[^\w]','', regex=True)
    pilger_Igneous['Age'] = pd.to_numeric(pilger_Igneous['Age'], errors='coerce')
    pilger_Igneous = pilger_Igneous[['Longitude', 'Latitude', 'Age']]
    
    return gpd.GeoDataFrame(pilger_Igneous, geometry=gpd.points_from_xy(pilger_Igneous.Longitude, pilger_Igneous.Latitude), crs=4326)


def plot_volcanics(fig, volcanics, time_min, time_max, 
                   reconstruction_time=None,
                   reconstruction_model=None, anchor_plate_id=0,
                   color='black', style='c0.2c', transparency=50,
                   perspective=None, region=None, projection=None):

    tmp = volcanics[(volcanics.Age>time_min) & (volcanics.Age<time_max)]
    #tmpM = pilger_Metamorphic[(pilger_Metamorphic.Age>time_min) & (pilger_Metamorphic.Age<time_max)]
    
    if reconstruction_model is not None:
        if reconstruction_time is None:
            reconstruction_time = (time_min+time_max)/2
        tmp = reconstruction_model.assign_plate_ids(tmp)
        tmp_r = reconstruction_model.reconstruct(tmp, reconstruction_time=reconstruction_time, 
                                                 anchor_plate_id=anchor_plate_id)
        print(reconstruction_time)
    else:
        tmp_r = tmp

    fig.plot(x=tmp_r.geometry.x, y=tmp_r.geometry.y, style=style, fill=color, transparency=transparency, 
             region=region, projection=projection, perspective=perspective)     



def plot_elevation_basemap(fig, grid=None, cmap='hot', 
                           region='d', projection='M10c', perspective=[240, 35], coastlines=True):
    

    fig.basemap(region=region, projection=projection, perspective=perspective, frame='afg')

    if grid is not None:
        pygmt.makecpt(cmap=cmap, series=[-1000,5000,500], background='o', reverse=True)
        fig.grdimage(grid=grid, cmap=True, projection=projection, perspective=perspective, transparency=30, nan_transparent=True)
    
    if coastlines:
        fig.coast(shorelines='1p,gray', transparency=60, projection=projection, perspective=perspective)
    

        
def add_labels(fig, reconstruction_time=None,
               label_string=None,
               x=2.5, y=1.15, font='48p',
               add_colorbar=False, 
               colorbar_position='JBL+jBL+o0.5c/1.6c+w12c/0.8c+h',
               colorbar_title='Elevation [m]'):
    
    if reconstruction_time is not None:
        fig.text(x=x ,y=y, text='{:0.0f} Ma'.format(reconstruction_time),
                 region='0/1/0/1', projection='x10c', font=font, no_clip=True)
    elif label_string is not None:
        fig.text(x=x ,y=y, text=label_string,
                 region='0/1/0/1', projection='x10c', font=font, no_clip=True)
    if add_colorbar:
        with pygmt.config(FONT_ANNOT_PRIMARY='16p', FONT_LABEL='24p'):
            fig.colorbar(position=colorbar_position, frame=['x+l{:s}'.format(colorbar_title)])

    return


def plot_elevations_as_columns(fig, binned_df, cmap='hot',
                               column_delta_height=None, column_marker_size='0.3c', column_pen='0.3p,black',
                               region='d', projection_3d='Z3c', perspective=[240, 35]):
    
    fig.plot3d(x=binned_df.x, 
               y=binned_df.y, 
               z=binned_df.median_elevation, 
               style='O{:s}'.format(column_marker_size), 
               region=region+[0,6000], projection=projection_3d, perspective=perspective, pen=column_pen)

    pygmt.makecpt(cmap=cmap, series=[-1000,5000,500], background='o', reverse=True)
    if column_delta_height is None:
        fig.plot3d(data = pd.DataFrame(data={'x':binned_df.x, 
                                             'y':binned_df.y, 
                                             'z1':binned_df.q75,  # zvalue for the top of the column
                                             'z2':binned_df.median_elevation, # value assigning color
                                             'z3':binned_df.q25}), # zvalue for the base of the column,
                style='o{:s}+b'.format(column_marker_size), 
                cmap=True,
                region=region+[0,6000], projection=projection_3d, perspective=perspective, pen=column_pen)
    else:
        fig.plot3d(data = pd.DataFrame(data={'x':binned_df.x, 
                                             'y':binned_df.y, 
                                             'z1':binned_df.median_elevation+column_delta_height,  # zvalue for the top of the column
                                             'z2':binned_df.median_elevation, # value assigning color
                                             'z3':binned_df.median_elevation-column_delta_height}), # zvalue for the base of the column,
                style='o{:s}+b'.format(column_marker_size), 
                cmap=True,
                region=region+[0,6000], projection=projection_3d, perspective=perspective, pen=column_pen)
    
    return


def geochem_timeslice(df, reconstruction_time, time_bin_size, 
                      calibration='luffi',
                      reconstruction_model=None, anchor_plate_id=0):
    
    df_filt = joy.filter_the_database(df, filter_method=calibration, 
                                      age_min=reconstruction_time-time_bin_size/2., 
                                      age_max=reconstruction_time+time_bin_size/2.)
    
    if df_filt.empty:
        return None
    
    if reconstruction_model is not None:
        df_filt_r = reconstruct(df_filt, reconstruction_model, reconstruction_time=reconstruction_time, 
                                anchor_plate_id=anchor_plate_id)
    else:
        print('No reconstruction model specified, results will assume present-day coordinates...')
        df_filt_r = df_filt

    return df_filt_r
    
'''
def binned_elevation_estimates(df_filt_r, mohometer_selection, gc_interpolator_dict, space_bin_size):
    # This should be df_filt_r???
    elevation_estimates = joy.get_elevations(df_filt_r, gc_interpolator_dict=gc_interpolator_dict)
    
    return joy.bin_elevations(df_filt_r.geometry, elevation_estimates, space_bin_size)
'''


def timeslice_plot(df, reconstruction_time,
                   time_bin_size, space_bin_size, 
                   fig, reconstruction_model=None, raster_sequence=None,  
                   anchor_plate_id=0, raster_anchor_plate_id=None,
                   calibration='luffi', mohometer_selection=None, gc_interpolator_dict=None,
                   residuals=False, volcanics=False, return_type=False,
                   column_marker_size='0.3c', plot_basemap=True, coastlines=True,
                   region='d', projection='M15c', 
                   perspective=None, projection_3d='Z3c'):

    if raster_anchor_plate_id is None:
        raster_anchor_plate_id = anchor_plate_id
    
    df_filt_r = geochem_timeslice(df, reconstruction_time, time_bin_size, calibration=calibration,
                                  reconstruction_model=reconstruction_model, anchor_plate_id=anchor_plate_id)
    
    if df_filt_r is not None:
        elevation_estimates = joy.get_elevations(df_filt_r, 
                                                gc_interpolator_dict=gc_interpolator_dict,
                                                calibration=calibration,
                                                mohometer_selection=mohometer_selection)
        
        binned_df = joy.bin_elevations(df_filt_r.geometry, elevation_estimates, space_bin_size)
        #binned_df = binned_elevation_estimates(df_filt_r, gc_interpolator_dict, space_bin_size)
    else:
        binned_df = []
        elevation_estimates = []
        residuals_bdf = []
        elevations_residuals = []
    
    if raster_sequence is not None:
        if reconstruction_model is None:
            filtered_topo = raster_sequence
        else:
            filtered_topo = pygmt.grdfilter(to_anchor_plate(raster_sequence[reconstruction_time], 
                                                    reconstruction_model, 
                                                    reconstruction_time, 
                                                    anchor_plate_id, 
                                                    old_anchor_plate_id=raster_anchor_plate_id,
                                                    region=region, spacing=0.2),
                                    distance='4',
                                    filter='m250+g0.5',
                                    spacing=0.25,
                                    coltypes='g')
    else:
        filtered_topo = None

    if plot_basemap:
        plot_elevation_basemap(fig, filtered_topo, cmap='cubhelix',
                               region=region, projection=projection, perspective=perspective, coastlines=coastlines)
    
    if volcanics is not None:
        plot_volcanics(fig, 
                       volcanics,
                       reconstruction_time-(time_bin_size/2), 
                       reconstruction_time+(time_bin_size/2), 
                       reconstruction_time=reconstruction_time,
                       reconstruction_model=reconstruction_model, 
                       anchor_plate_id=anchor_plate_id,
                       perspective=perspective)
    
    if len(binned_df)>0:
        if residuals:

            DEM = pygmt.grdtrack(points=pd.DataFrame(data={'x':df_filt_r.geometry.x,'y':df_filt_r.geometry.y}),
                                 grid=filtered_topo, 
                                 newcolname='z', 
                                 no_skip=True)['z']
            elevations_residuals = elevation_estimates.sub(DEM, axis=0)

            residuals_bdf = joy.bin_elevations(df_filt_r.geometry, elevations_residuals, bin_size_degrees=space_bin_size)
            #resd = elevations_residuals.drop(columns=['bin_latitude', 'bin_longitude', 'geometry'])

            #elevations_residuals = pygmt.grdtrack(points=binned_df,#pd.DataFrame(data={'x':bdf.x,'y':bdf.y}),
            #                                      grid=filtered_topo, 
            #                                      newcolname='z')
            
            pygmt.makecpt(cmap='polar', series=[-3000,3000,500], background='o')
            residual_marker_size = '{:f}{:s}'.format(float(column_marker_size[:-1])*1.5, column_marker_size[-1])
            fig.plot(x=residuals_bdf.x, y=residuals_bdf.y, 
                     fill=residuals_bdf.median_elevation, # np.array(elevations_residuals.z)-np.array(elevations_residuals.median_elevation), 
                     style='s{:s}'.format(residual_marker_size), pen='0.2,black', perspective=perspective, cmap=True)

        else:
            plot_elevations_as_columns(fig, binned_df, cmap='cubhelix',
                                       column_marker_size=column_marker_size, column_pen='0.2p,black',
                                       column_delta_height=500, #300,
                                       region=region, perspective=perspective, projection_3d=projection_3d)

    if return_type is not None:
        if return_type=='raw_residuals':
            return elevations_residuals
        elif return_type=='binned_residuals':
            return residuals_bdf
        elif return_type=='raw_elevations':
            return elevation_estimates
        elif return_type=='binned_elevations':
            return binned_df



def get_raster_stats_in_polygon_time_series(raster_list, reconstruction_model, polygon_gdf,
                                            plate_id_to_rotate_to=None, old_anchor_plate_id=0,
                                            region=[-180,180,-90,90], spacing=0.2, return_type='averages'):
    
    my_time_series_stats = []
    my_time_series_histogram = []
    
    for reconstruction_time in raster_list:

        #print(raster_list[reconstruction_time])
        if plate_id_to_rotate_to is not None:
            rotated_grid = to_anchor_plate(raster_list[reconstruction_time], 
                                        reconstruction_model, 
                                        reconstruction_time, 
                                        plate_id_to_rotate_to,
                                        old_anchor_plate_id=old_anchor_plate_id,
                                        region=region, spacing=spacing)
        else:
            rotated_grid = xr.load_dataarray(raster_list[reconstruction_time])

        # To use the clipping functionality of rioxarray, we need to set some metadato to the grid object 
        # specifically the coordinate projection system and the names of the coordinate fields
        rotated_grid.rio.set_crs("epsg:4326")
        rotated_grid.rio.set_spatial_dims('lon', 'lat')

        # Now apply the polygon mask to the rotated grid
        # Note here we also specify the crs of the geometries though this is not strictly necessary where the 
        # crs is the same as the grid
        masked_rotated_grid = rotated_grid.rio.clip([polygon_gdf.geometry], crs="epsg:4326")
            #index = rotated_grid.data>=200.
        masked_rotated_grid.data[masked_rotated_grid.data<-500.] = np.nan
        
        my_time_series_stats.append((reconstruction_time, 
                                     np.nanmin(masked_rotated_grid.data),
                                     np.nanmax(masked_rotated_grid.data),
                                     np.nanmedian(masked_rotated_grid.data),
                                     np.nanmean(masked_rotated_grid.data),
                                     np.nanquantile(masked_rotated_grid.data, q=0.25),
                                     np.nanquantile(masked_rotated_grid.data, q=0.75)))
        
        h,bin_edges = np.histogram(masked_rotated_grid.data[~np.isnan(masked_rotated_grid.data)], 
                                   bins=np.arange(-500.,8000,100.), density=True)
        my_time_series_histogram.append(100*h)
        
    my_time_series_histogram = np.flipud(np.cumsum(np.flipud(np.vstack(my_time_series_histogram).T), axis=0))
    my_time_series_histogram[my_time_series_histogram==0] = np.nan
    
    #return pd.DataFrame(data=zip(my_time_series_time, my_time_series_min_elev, my_time_series_max_elev, my_time_series_median_elev, my_time_series_mean_elev),
    #                    columns=['age', 'min_elevation', 'max_elevation', 'median_elevation', 'mean_elevation']), my_time_series_histogram
    return pd.DataFrame(data=my_time_series_stats,
                        columns=['age', 'min_elevation', 'max_elevation', 'median_elevation', 'mean_elevation', 'q25_elevation', 'q75_elevation']), my_time_series_histogram


def add_boxplots(df_binned, ax, 
                 color='black',
                 colnames = ['la_yb_elevation', 'sr_y_elevation', 'gd_yb_elevation'],
                 style='box',
                 widths=5):
    
    midpoints = [key.mid for key in list(df_binned.groups.keys())]

    median_elevations = []
    elevation_violins = []
    for group in df_binned:
        if not group[1].empty:
            median_elevations.append(group[1][colnames].stack().median())
            elevation_violins.append(group[1][colnames].stack())
            group[0]
        else:
            median_elevations.append(np.nan)
            elevation_violins.append([np.nan,np.nan])

    if style=='box':
        parts = ax.boxplot(elevation_violins, positions=midpoints, widths=widths, 
                               patch_artist=True, manage_ticks=False)
        for pc in parts['boxes']:
            pc.set_facecolor(color)
            pc.set_alpha(0.5)
            pc.set_edgecolor(color)
    
    if style=='violin':
        parts = ax.violinplot(elevation_violins, positions=midpoints, widths=widths, bw_method=.25)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.5)
            pc.set_edgecolor(color)

        

def residuals_crossplot(binned_df, residuals_bdf, fname):

    #resd = elevations_residuals.drop(columns=['bin_latitude', 'bin_longitude', 'geometry'])
    
    fig,ax = plt.subplots(figsize=(6,4.5))
    ax.plot([-1000,8000],[-1000,8000], 'k-')
    if len(binned_df)>0:
        ax.scatter(binned_df.median_elevation-residuals_bdf.median_elevation, binned_df.median_elevation,
                   c=residuals_bdf.median_elevation, vmin=-3000, vmax=3000, cmap='seismic', edgecolor='k')
    #plt.bar(bdf.median_elevation,height=residuals_bdf.median_elevation, width=20,color='k')
    ax.grid()
    ax.set_xlabel('DEM Elevation [m]')
    ax.set_ylabel('Elevation Proxy from geochemistry [m]')
    ax.set_xlim(-1000,6000)
    ax.set_ylim(-1000,6000)
    if len(binned_df)>0:
        ax.text(4200,600,'RMS = {:0.1f}m'.format(np.sqrt(np.mean(residuals_bdf.median_elevation**2))))
        ax.text(4200,200,'STD = {:0.1f}m'.format(residuals_bdf.median_elevation.std()))
    #ax.set_title('{:s}, {:s}'.format(calibration,mohometer_description_string))
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

