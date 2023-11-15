import numpy as np
import pandas as pd
import geopandas as gpd
import pygmt
from collections import OrderedDict
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import median_abs_deviation


def geochem_from_csv(csvfile,
                     longitude_field_name=None,
                     latitude_field_name=None):

    df = pd.read_csv(csvfile)
 
    #if (longitude_field_name is not None) | (latitude_field_name is not None):
    #    df.rename(columns={longitude_field_name: 'Longitude',
    #                       latitude_field_name: 'Latitude'})

    df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[longitude_field_name], 
                                                          df[latitude_field_name]), crs=4326)
    return compute_ratios(df, long_list=True)

'''
def fix_fe(df):
    
    # Based on the method outlined here: 
    # https://serc.carleton.edu/research_education/cyberinfrastructure/Harkers/part_4.html
    
    fe_fields = df.copy()[['feo', 'feo_tot', 'fe2o3', 'fe2o3_tot']]

    feo_tot = fe_fields['feo_tot'].to_numpy()

    # We will compute the result for all rows where the feo_tot is nan
    ind = np.isnan(feo_tot)

    # TODO check whether divide by 
    corrected_value = (fe_fields.fe2o3*0.8998 + fe_fields.fe2o3_tot*0.8998 + fe_fields.feo)/2.

    feo_tot[ind] = corrected_value[ind]
    
    return feo_tot
'''
def fix_fe(df):
    """
    If FeOT was reported, and then FeOT = FeOT. 
    If FeO was reported and Fe2O3 was not, FeOT = FeO. 
    If only Fe2O3T was reported, FeOT = 0.9 * Fe2O3T. 
    If both FeO and Fe2O3 were reported, FeOT = FeO + 0.9 * Fe2O3. 
    If neither FeO nor Fe2O3T were reported and Fe2O3 was reported, FeOT = 0.9 * Fe2O3
    """
    for column in ['feo', 'feo_tot', 'fe2o3', 'fe2o3_tot']:
        if column not in df.columns:
            df[column] = np.nan
            
    
    ind1 = df['feo_tot'].isna()
    
    ind2 = df['feo_tot'].isna() & ~df['feo'].isna() & df['fe2o3'].isna()
    
    ind3 = ~df['fe2o3_tot'].isna() & df['feo_tot'].isna() & df['feo'].isna() & df['fe2o3'].isna()
    
    ind4 = df['feo_tot'].isna() & ~df['feo'].isna() & ~df['fe2o3'].isna()
    
    ind5 = df['fe2o3_tot'].isna() & df['feo_tot'].isna() & df['feo'].isna() & ~df['fe2o3'].isna()
    
    fixed_feo_tot = np.array(df['feo_tot'])
    fixed_feo_tot[ind2] = df.loc[ind2,'feo']
    fixed_feo_tot[ind3] = df['fe2o3_tot'][ind3] * 0.9
    fixed_feo_tot[ind4] = df['feo'].astype(float)[ind4] + df['fe2o3'].astype(float)[ind4] * 0.9
    fixed_feo_tot[ind5] = df['fe2o3'].astype(float)[ind5] * 0.9
    
    #print(fixed_feo_tot.values)
    df['feo_tot'] = fixed_feo_tot
    
    return df



# TODO add a fix MnO, fix P2O5 function after Gale et al??




def normalize(df, columns, reference_composition='NMORB_SM89'):
    """
    Normalise columns in geochemical database (assuming column names after Gard et al, 2019)
    using one of the reference compositions available through the pyrolite function
    'pyrolite.geochem.norm'
    """
    from pyrolite.geochem.norm import get_reference_composition
    
    f = get_reference_composition(reference_composition)
    
    for column_name in columns:
        value = f.comp[column_name.capitalize()[:-4]].value
        
        df['{:s}_n'.format(column_name)] = df[column_name] / value
        
    return df

# For Chondrite SM89, reference values are:
# La 0.237
# Yb 0.17
    
def compute_ratios(df, long_list=False):
    
    df['la_yb'] = df['la_ppm']/df['yb_ppm']
    df['gd_yb'] = df['gd_ppm']/df['yb_ppm']
    df['sr_y'] = df['sr_ppm']/df['y_ppm']
    
    df['ln_la_yb'] = np.log(df['la_yb'])
    df['ln_gd_yb'] = np.log(df['gd_yb'])
    df['ln_sr_y'] = np.log(df['sr_y'])
    
    if long_list:
        # ordering based on Luffi and Ducea table #1
        df['ce_yb'] = df['ce_ppm']/df['yb_ppm']
        df['lu_hf'] = df['lu_ppm']/df['hf_ppm']
        df['nd_yb'] = df['nd_ppm']/df['yb_ppm']
        df['la_y'] = df['la_ppm']/df['y_ppm']
        df['nb_y'] = df['nb_ppm']/df['y_ppm']
        df['sm_yb'] = df['sm_ppm']/df['yb_ppm']
        df['ce_y'] = df['ce_ppm']/df['y_ppm']
        df['nd_y'] = df['nd_ppm']/df['y_ppm']
        df['th_yb'] = df['th_ppm']/df['yb_ppm']
        df['la_sm'] = df['la_ppm']/df['sm_ppm']
        df['ba_v'] = df['ba_ppm']/df['v_ppm']
        df['zr_y'] = df['zr_ppm']/df['y_ppm']
        df['ba_sc'] = df['ba_ppm']/df['sc_ppm']
        df['a'] = df['na2o'] + df['k2o']
        df['nb_yb'] = df['nb_ppm']/df['yb_ppm']
        df['dy_yb'] = df['dy_ppm']/df['yb_ppm']
        df['a_cao'] = df['a']/df['cao']
        df['th_y'] = df['th_ppm']/df['y_ppm']
        df['ni_sc'] = df['ni_ppm']/df['sc_ppm']
        df['cr_sc'] = df['cr_ppm']/df['sc_ppm']
        df['zr_ti'] = df['zr_ppm']/df['ti_ppm']
        df['ni_v'] = df['ni_ppm']/df['v_ppm']
        df['cr_v'] = df['cr_ppm']/df['v_ppm']    
    
    return df

    
def filter_the_database(df, filter_method, nans_to_zeros=True, age_min=-1e9, age_max=1e9):
    
    # Regardless of filtering method, we will always:
    # 1. remove points without long/lat
    #df_filt = df.dropna(subset=['Longitude', 'Latitude'])
    df_filt = df.copy(deep=True)
    # 2. Assume that if there is no age, then it is a very young rock
    # then, remove samples with age greater than 1 Myr
    if nans_to_zeros:
        df_filt['age'] = df_filt['age'].replace(np.nan, 0.)
    
    # Age Filter
    df_filt = df_filt[(df_filt['age']<=age_max) & (df_filt['age']>=age_min)]

    print('Number of samples after basic filtering {:d}'.format(len(df_filt)))

    if filter_method in ['Chapman','Hu']:

        # Chapman (2015), Hu et al (2017)
        df_filt = df_filt[(df_filt['sio2']>=55)&(df_filt['sio2']<=70)]
        print('Number of samples with 55<=sio2<=70 = {:d}'.format(len(df_filt)))
        df_filt = df_filt[(df_filt['mgo']>=1)&(df_filt['mgo']<=4)]
        print('Number of these samples with 1<=mgo<=4 = {:d}'.format(len(df_filt)))
        df_filt = df_filt.dropna(subset=['rb_ppm','sr_ppm'])
        df_filt['rb/sr'] = df_filt['rb_ppm']/df_filt['sr_ppm']
        df_filt = df_filt[(df_filt['rb/sr']>=0.05)&(df_filt['rb/sr']<=0.25)]
        print('Number of these samples with 0.05<=rb/sr<=0.25 = {:d}'.format(len(df_filt)))

    if filter_method == 'FarnerLee':
        
        # Farner and Lee (2017)
        df_filt = df_filt.dropna(subset=['sio2'])
        print('Number of these samples with a valid sio2 = {:d}'.format(len(df_filt)))


        mee_list = ['sio2', 'mgo', 'feo_tot', 'cao', 'k2o', 'tio2', 'al2o3', 'mno', 'na2o', 'p2o5']
        for column in mee_list:
            df_filt[column] = df_filt[column].replace(np.nan, 0.)
        
        major_element_oxide_sum = df_filt['sio2']+df_filt['mgo']+df_filt['feo_tot']+df_filt['cao']+df_filt['k2o']+df_filt['tio2']+df_filt['al2o3']+df_filt['mno']+df_filt['na2o']+df_filt['p2o5']

        #print(df_filt.columns.values)
        df_filt=df_filt[major_element_oxide_sum>=98.]
        print('Number of these samples with major element sum > 98%= {:d}'.format(len(df_filt)))
        
        #df_filt = df_filt[(df_filt['mgo']>=4)]
        #print('Number of these samples with mgo >= 4 {:d}'.format(len(df_filt)))
        
    if filter_method in ['Luffi','luffi']:
        #print('TODO implement a specific alkaline/subalkaline boundary')
        df_filt = df_filt.query('`sio2`>=45 & `sio2`<=80')
        
        #'''
        # Ref Rickwood
        # Maximum boundary...
        #alkaline_boundary = pd.DataFrame(data={'sio2': [40, 45, 50, 55, 60, 65, 70, 100],
        #                                       'A': [0.45, 2.8, 4.75, 6.5, 8.0, 9.6, 11.1, 19]})

        # Irvine and Baragar 1971 (as defined in Rickwood, 1989)
        alkaline_boundary  = pd.DataFrame(
            data={'sio2': [39, 41.56, 43.28, 45.47, 48.18, 51.02, 53.72, 56.58, 60.47, 66.82, 77.15, 100],
            'A': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12.5]})
        
        buffer_percentage = 1.5
        A_interp = np.interp(df_filt['sio2'], 
                             alkaline_boundary['sio2'], 
                             alkaline_boundary['A']+buffer_percentage)

        #A_test = A_interp-df_filt['a']

        df_filt = df_filt.query('`a` < @A_interp')
        #'''

        #pass


    print('Final number of samples passed = {:d}'.format(len(df_filt)))
    
    df_filt = df_filt.reset_index(drop=True)
    return df_filt


def get_elevations_for_samples(df):
    
    df = df.dropna(subset=['Longitude', 'Latitude'])

    grid = pygmt.datasets.load_earth_relief(region='d', resolution='02m', registration='gridline')

    res = pygmt.grdtrack(grid=grid, points = df[['Longitude','Latitude']], newcolname='elevation')

    df = df.reset_index(drop=True)
    df['elevation'] = res['elevation']
    
    return df



def get_data_within_polygon(gc_gdf, polygon_gdf):
    
    return gpd.clip(gc_gdf, polygon_gdf)
    
    

def bin_gc_data(df_input, upper_value, lower_value=None, 
                elev_binsize=100, spatial_binsize=0.1, age_binsize = 20., polygon_file=None):
    
    # TODO make the binning field(s) optional
    
    df_filt = df_input
    
    if polygon_file is not None:
        gc_gdf = gpd.GeoDataFrame(df_filt, geometry=gpd.points_from_xy(df_filt.longitude, df_filt.latitude), crs=4326)
        polygon_gdf = gpd.read_file(polygon_file)
        df_filt = gpd.clip(gc_gdf, polygon_gdf)

    df_filt = df_filt.reset_index(drop=True)

    # If we have defined both an upper and lower value, we compute the ratio
    # If no lower value is defined, just take the upper_value directly for the next 
    # calculations
    if lower_value is not None:
        subset = df_filt.dropna(subset=['elevation',upper_value,lower_value])
        subset = subset[subset[upper_value]>=0]
        subset = subset[subset[lower_value]>=0]
        subset['ratio'] = subset[upper_value]/subset[lower_value]
        
    else:
        subset = df_filt.dropna(subset=['elevation',upper_value])
        subset = subset[subset[upper_value]>0]
        subset['ratio'] = subset[upper_value]

    # Because it is possible to have issues with 'divide by zero' when computing
    # a ratio, we check for infinity values and remove them
    subset['ratio'] = subset['ratio'].replace([np.inf,-np.inf], np.nan)
    subset = subset.dropna(subset=['ratio'])
    subset = subset.reset_index(drop=True)
    
    num_samples = len(subset)
    
    #gc_clip_blockmedian_elev=subset[['Longitude','Latitude','elevation']]
    #res_elev = pygmt.blockmedian(gc_clip_blockmedian_elev,region='d',spacing='10k',C=True)
    #gc_clip_blockmedian_geochem = subset[['Longitude','Latitude','ratio']]
    #res_geochem = pygmt.blockmedian(gc_clip_blockmedian_geochem,region='d',spacing='10k',C=True)

    # TODO
    # Replace blockmedian with a step to round longs and lats to a specified factor
    subset['Longitude_bin'] = np.round(subset['Longitude']/spatial_binsize)*spatial_binsize
    subset['Latitude_bin'] = np.round(subset['Latitude']/spatial_binsize)*spatial_binsize
    
    subset['elevation_bin'] = np.round(subset['elevation']/elev_binsize)*elev_binsize
    
    subset['age_bin'] = np.round(subset['age']/age_binsize)*age_binsize

        
    return num_samples, subset


def get_elevation(values, ratio, return_error_bounds=False):
    #print(values)
    values = values.mask(values<=0.)
    values = values.mask(np.isinf(values))
    #values[values<=0.] = np.nan
    #values[np.isinf(values)] = np.nan
    if ratio=='gd_yb':
        return np.log(values / 1.8) / 1.157e-4  # From Farner and Lee 2017
    elif ratio=='la_yb':
        #values = (values/0.237)*0.17 # normalize to chondrite
        return np.log(values / 6.91) / 2.651e-4  # From Farner and Lee 2017
    elif ratio=='la_yb_hu':
        values = (values/0.237)*0.17 # normalize to chondrite
        return np.log(values / 2.61) / 0.41e-3  # From Hu++2020
    elif ratio=='sr_y' or ratio=='sr_y_hu':
        return (values-4.71) / 0.0105  # From Hu++2020
    

def get_elevations(df, gc_interpolator_dict=None, calibration='luffi', mohometer_selection=None, 
                   elevation_min=None, elevation_max=None):
    '''
    given a dataframe with various columns of geochemical data, return
    the estimated elevations
    '''
    print('TODO implement min/max elevation cutoffs')
    if calibration in ['Hu', 'FarnerLee']:
        elevation_estimates = get_two_elevations(df, calibration=calibration, elevation_min=elevation_min, elevation_max=elevation_max)
    elif calibration == 'luffi':
        if gc_interpolator_dict is None:
            ValueError('Interpolator dictionary must be supplied for luffi calibration')
        elevation_estimates = get_luffi_elevations(df, gc_interpolator_dict, 
                                                   elevation_min=elevation_min, elevation_max=elevation_max)
    else:
        ValueError('Unknown calibration')

    if mohometer_selection is None:
        return elevation_estimates
    elif (isinstance(mohometer_selection, str) or (isinstance(mohometer_selection, list))):
        return elevation_estimates.loc[:,mohometer_selection]
    else: # if mohometer selection is an integer
        return elevation_estimates.iloc[:,:mohometer_selection]
           


def get_two_elevations(df, calibration='Hu', elevation_min=-2000, elevation_max=6000):
    # assuming we are basing our elevations directly on a previous study, this
    # will be either 
    # 1. Hu++2020, with calibrations for (La/Yb)n and Sr/Y
    # 2. Farner+Lee2017, with calibrations for La/Yb and Gd/Yb
    
    if calibration=='Hu':
        df['la_yb_elevation'] = get_elevation(df['la_yb'], 'la_yb_hu')
        df['sr_y_elevation'] = get_elevation(df['sr_y'], 'sr_y')
        df = df[['la_yb_elevation', 'sr_y_elevation']]
    
    elif calibration=='FarnerLee':
        df['la_yb_elevation'] = get_elevation(df['la_yb'], 'la_yb')
        df['gd_yb_elevation'] = get_elevation(df['gd_yb'], 'gd_yb')
        df = df[['la_yb_elevation', 'gd_yb_elevation']]
    
    #df.loc[df['la_yb_elevation']>6000., 'la_yb_elevation'] = np.nan
    #df.loc[df['sr_y_elevation']>6000., 'sr_y_elevation'] = np.nan
    #df.loc[df['gd_yb_elevation']>6000., 'gd_yb_elevation'] = np.nan
        
    #df.loc[df['la_yb_elevation']<-2000., 'la_yb_elevation'] = np.nan
    #df.loc[df['sr_y_elevation']<-2000., 'sr_y_elevation'] = np.nan
    #df.loc[df['gd_yb_elevation']<-2000., 'gd_yb_elevation'] = np.nan
    
    #return df[['la_yb_elevation', 'sr_y_elevation', 'gd_yb_elevation']]
    if elevation_min is not None:
        df[df<elevation_min] = np.nan
    if elevation_max is not None:
        df[df>elevation_max] = np.nan
    return df


def time_binned_elevations(elevations, bin_size):
    
    median_ages=np.arange(-bin_size/2,500,bin_size)
    df_binned = elevations.groupby([pd.cut(elevations['age'], bins=median_ages)])
    
    return df_binned


############################################################
# Luffi and Ducea functions

def create_gc_interpolator(fname):
    
    dat = np.loadtxt(fname, delimiter=',')

    dat = dat.reshape((np.unique(dat[:,0]).shape[0], np.unique(dat[:,1]).shape[0], 3))

    points = (dat[:,0,0],dat[0,:,1])
    f = RegularGridInterpolator(points,
                                dat[:,:,2], method='linear',
                                bounds_error=False)
    
    return f


def make_gc_interpolator_dict(model_dir):
    '''
    Create grid interpolator objects for each of the geochemical 'sensors'
    covered by Luffi and Ducea
    '''
    
    table_dict = OrderedDict({
        'la_yb': 'LaYb_model.csv',
        'ce_yb': 'CeYb_model.csv',
        'lu_hf': 'LuHf_CS_model.csv',
        'nd_yb': 'NdYb_model.csv',
        'la_y': 'LaY_CS_model.csv',
        'nb_y': 'NbY_CS_model.csv',
        'sm_yb': 'SmYb_CS_model.csv',
        'ce_y': 'CeY_CS_model.csv',
        'nd_y': 'NdY_CS_model.csv',
        'th_yb': 'ThYb_model.csv',
        'hf_ppm': 'Hf_CS_model.csv',
        'la_sm': 'LaSm_CS_model.csv',
        'gd_yb': 'GdYb_CS_model.csv',
        'ba_v': 'BaV_model.csv',
        'mno': 'MnO_CS_model.csv',
        'k2o': 'K2O_CS_model.csv',
        'zr_y': 'ZrY_model.csv',
        'ba_sc': 'BaSc_CS_model.csv',
        'a': 'A_CS_model.csv',
        'nb_yb': 'NbYb_model.csv',
        'dy_yb': 'DyYb_CS_model.csv',
        'a_cao': 'ACaO_model.csv',
        'th_y': 'ThY_CS_model.csv',
        'cao': 'CaO_CS_model.csv',
        'rb_ppm': 'Rb_CS_model.csv',
        'ni_sc': 'NiSc_model.csv',
        'sc_ppm': 'Sc_model.csv',
        'ga_ppm': 'Ga_model.csv',
        'pb_ppm': 'Pb_model.csv',
        'ba_ppm': 'Ba_CS_model.csv',
        'u_ppm': 'U_CS_model.csv',
        'feo_tot': 'FeOt_CS_model.csv',
        'cr_sc': 'CrSc_model.csv',
        'ni_ppm': 'Ni_model.csv',
        'sr_y': 'SrYx_model.csv',
        'sr_ppm': 'Sr_CS_model.csv',
        'zr_ti': 'ZrTi_CS_model.csv',
        'ni_v': 'NiV_model.csv',
        'cr_v': 'CrV_CS_model.csv',
    })
 
    gc_interpolator_dict = OrderedDict()

    for gc in table_dict:
        #print(table_dict[gc])

        gc_interpolator_dict[gc] = create_gc_interpolator('{:s}/{:s}'.format(model_dir,table_dict[gc]))

    return gc_interpolator_dict


def get_luffi_elevations(df, gc_interpolator_dict, elevation_min=None, elevation_max=None):

    results = pd.DataFrame()
    #for col in df.columns:
    #    if col not in list(gc_interpolator_dict.keys()):
    for col in gc_interpolator_dict.keys():
        if col not in df.columns:
            pass
        else:
            interpolator = gc_interpolator_dict[col]
            result = interpolator((df['mgo'],df[col]))
            
            results['{:s}_elevation'.format(col)] = result*1000.
    
    return results


def bin_elevations(geometry, elevations_df, bin_size_degrees):
    
    #ppdat = elevations_df.join(coords_df) #df_filt[['Longitude', 'Latitude']])
    ppdat = gpd.GeoDataFrame(elevations_df, geometry=geometry, crs=4326)
    
    ppdat['bin_latitude'] = np.round(ppdat.geometry.y/bin_size_degrees) * bin_size_degrees
    ppdat['bin_longitude'] = np.round(ppdat.geometry.x/bin_size_degrees) * bin_size_degrees
    
    #ppdat = ppdat.drop(columns=['Longitude', 'Latitude'])

    p_groups = ppdat.groupby(by=['bin_longitude','bin_latitude'])

    binned_list = []
    for g in p_groups:
        data = g[1].drop(columns=['bin_longitude', 'bin_latitude', 'geometry']).stack()
        binned_list.append([g[0][0], 
                            g[0][1], 
                            data.median(),
                            median_abs_deviation(data),
                            data.quantile(q=0.25),
                            data.quantile(q=0.75)])

    binned_df = pd.DataFrame(data=binned_list, 
                             columns=['x', 'y', 'median_elevation', 'elevation_mad', 'q25', 'q75'])#.dropna(subset=['median_elevation'])

    return binned_df


def assign_bin_ages(ages, time_bin_size):

    return (np.round((ages-time_bin_size/2)/time_bin_size)*time_bin_size) + (time_bin_size/2)



def thickness2elevation(thickness, rho_crust=2.75, rho_mantle=3.3, thickness0=32000):
    
    elevation = (1 - rho_crust / rho_mantle) * (thickness-thickness0)
    
    return elevation


def elevation2thickness(elevation, rho_crust=2.75, rho_mantle=3.3, thickness0=32000):

    thickness = elevation / (1 - rho_crust / rho_mantle)
    
    return thickness + thickness0
    