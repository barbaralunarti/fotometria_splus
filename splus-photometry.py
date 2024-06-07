# %%

%pip install pyerfa
%pip install git+https://github.com/astropy/pyregion.git
%pip install splusdata==4.0
%pip install aplpy==2.1.0
%pip install astropy==5.3.4
%pip install pandas
%pip install numpy==1.26.2
%pip install matplotlib
%pip install setuptools==69.0.3

# %%

import splusdata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
import aplpy
import pickle
import astropy.constants as c
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.wcs import WCS
from IPython.display import display

# %%

# Opens the list of galaxies

csv = open('galaxies_splus.csv', errors='ignore')
df = pd.read_csv(csv)
df

# %%

conn = splusdata.connect('username','password')

# %%

counttrue = 0
countfalse = 0
found_ind = []
for i in range(0,len(df)):
    amostra = df.iloc[i]
    RA = amostra['RA']
    DEC = amostra['DEC']
    coord = SkyCoord(RA, DEC, unit = (u.hourangle, u.deg))
    res = conn.checkcoords(coord.ra.degree, coord.dec.degree)
    if(res['distance'] > 0.7):
        countfalse += 1
    else:
        counttrue += 1
        found_ind.append(i)
print('Amostras encontradas: ', counttrue)
print('Amostras não encontradas: ', countfalse)

# %%

gal_sample = df.iloc[found_ind]

c = SkyCoord(ra=gal_sample['RA'], dec=gal_sample['DEC'], unit=(u.hourangle, u.deg))
type(c.ra.value.tolist())

gal_sample['RA_ICRS'] = c.ra.value.tolist()
gal_sample['DEC_ICRS'] = c.dec.value.tolist()
gal_sample

# %%

gal_sample.to_pickle('gal_sample.pkl')

# %%

loaded_gal_sample = pd.read_pickle('gal_sample.pkl')
loaded_gal_sample

# %%

# Query Splus using the sample data

querytable = []

for i in range(0,len(gal_sample)):
    inp_ra = gal_sample['RA_ICRS'].iloc[i]
    inp_dec = gal_sample['DEC_ICRS'].iloc[i]
    radius = (60.0 * u.arcsec).to(u.degree).value

    query = f"SELECT * FROM idr4_dual.idr4_detection_image \
    WHERE 1=CONTAINS( POINT('ICRS',{inp_ra}, {inp_dec}), CIRCLE('ICRS', ra, dec, {radius}) )"

    querytable.append(conn.query(query))
    
# %%

for i in range(0,len(querytable)):
    with open(f'{gal_sample["Object"].iloc[i]}.pkl', mode='wb+') as file:
        pickle.dump(querytable[i],file)

# %%

loaded_list=[]
for i in range(0,len(loaded_gal_sample)-1):
    with open(f'{loaded_gal_sample["Object"].iloc[i]}.pkl', mode='rb') as file:
        loaded_list.append(pickle.load(file))

# %%

frequencies = ['FWHM_g','FWHM_i','FWHM_J0378','FWHM_J0395','FWHM_J0410','FWHM_J0430','FWHM_J0515','FWHM_J0660','FWHM_J0861','FWHM_r','FWHM_u','FWHM_z']
new_frequencies = ['G','I','F378','F395','F410','F430','F515','F660','F861','R','U','Z']

# %%

select_table = loaded_list[0]
select_table

df_data = select_table.to_pandas()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
display(df_data)

# %%

# ISO aperture

gal_index = 8
select_table = loaded_list[gal_index] # select the table
index_select = 29 # select the index of the table
hdu = conn.get_cut(select_table[index_select]['RA'], select_table[index_select]['DEC'], 200, 'R')
w = WCS(hdu[1].header)
gc = aplpy.FITSFigure(hdu,hdu=1)
gc.show_grayscale(invert=True)

iso = select_table['ISOarea'][index_select]
a = select_table['A'][index_select]
b = select_table['B'][index_select]
s = iso / (np.pi * a * b) # scaling parameter
print("Scaling parameter:", s)
w = (s*a)
h = (s*b)
gc.show_ellipses(select_table['RA'][index_select], select_table['DEC'][index_select],
                width=(w),
                height=(h),
                angle=-Angle(select_table['THETA'][index_select], unit=u.degree),
                edgecolor='green')

print(df['Object'][gal_index])
# %%

list_to_convert = [16, 97, 16, 18, 4, 23, 10, 37, 29, 10, 27, 11, 9, 9, 15, 18, 12, 4, 18, 15, 24, 8, 6, 16, 31, 22, 1, 30, 1, 28, 5, 6, 17, 19, 21, 14, 7, 6, 1, 19, 3, 26, 5, 22, 34, 9, 24, 24, 20, 15, 9, 18, 17, 23, 10, 19, 16, 11, 31, 22, 21, 24, 43, 32, 1, 6, 14, 8, 5, 34, 10, 18, 20, 24, 21, 27, 13, 9, 1, 21, 23, 16, 3, 18, 22, 17, 13, 17, 5, 4, 5, 39, 5, 5, 16, 17, 19, 12, 14, 0, 17, 25, 11, 19, 12, 6, 32, 20, 13, 25, 23, 21, 26, 46, 17, 29, 13, 25]
converted_list_ids = []
for i,obj_i in enumerate(list_to_convert):
    if i == len(list_to_convert)-1: ## missing last object
        break
    print(i)
    select_table = loaded_list[i]
    converted_list_ids.append(select_table[obj_i]['ID'])

# %%

# ISO aperture for all the filters

obj_index = 138
select_table = loaded_list[obj_index] # select the table of the galaxy (object)
table_index = 25 # select the index to center the image
index_ellipses = [25, 16] # draw the ellipses of the index
nrows = 4
ncols = 3

apertures_directory="apertures"
if not os.path.exists(apertures_directory):
    os.makedirs(apertures_directory)

obj_directory = os.path.join(apertures_directory, f"{loaded_gal_sample.iloc[obj_index]['Object']}")
if not os.path.exists(obj_directory):
    os.makedirs(obj_directory)

fig, axs = plt.subplots(nrows, ncols, figsize=(24, 36))
for ax_row in axs:
    for ax in ax_row:
        ax.axis('off')
    
fig.suptitle(f"{loaded_gal_sample.iloc[obj_index]['Object']}", fontsize=22)

for j in range(0,len(frequencies)):
    row = j // ncols
    col = j % ncols

    try:
        hdu = conn.get_cut(select_table[table_index]['RA'], select_table[table_index]['DEC'], 200, new_frequencies[j])
        w = WCS(hdu[1].header)
        gc = aplpy.FITSFigure(hdu,hdu=1, figure=fig, subplot=(nrows, ncols, j + 1))
        gc.show_grayscale(invert=True)
    except Exception as e:
        print("Exception:", str(e))
        continue
    
    for index in index_ellipses:
        iso = select_table['ISOarea'][index]
        a = select_table['A'][index]
        b = select_table['B'][index]
        s = iso / (np.pi * a * b) # scaling parameter
        print("Scaling parameter:", s)
        w = (s*a)
        h = (s*b)
            
        gc.show_ellipses(select_table['RA'][index], select_table['DEC'][index],
                        width=(w),
                        height=(h),
                        angle=-Angle(select_table['THETA'][index], unit=u.degree),
                        edgecolor='green')
        
        plt.title(f"{new_frequencies[j]}", fontsize=18)    
        aperture_filename = os.path.join(obj_directory, f"ISO.png")
        plt.tight_layout()
plt.close()
fig.savefig(aperture_filename)

# %%

# AUTO aperture (Kron radius)

select_table = loaded_list[0] # select the table
index_select=10 # select the index of the table
hdu = conn.get_cut(select_table[index_select]['RA'], select_table[index_select]['DEC'], 200, 'R')
w = WCS(hdu[1].header)
gc = aplpy.FITSFigure(hdu,hdu=1)
gc.show_grayscale(invert=True)

a = select_table['A'][index_select]
b = select_table['B'][index_select]
k = select_table['KRON_RADIUS'][index_select]

gc.show_ellipses(select_table['RA'][index_select], select_table['DEC'][index_select],
                width=(k * a), # major axis
                height=(k * b), # minor axis
                angle=-Angle(select_table['THETA'][index_select], unit=u.degree),
                edgecolor='red')

# %%

# AUTO aperture (Kron radius) for all the filters

obj_index = 138
select_table = loaded_list[obj_index] # select the table of the galaxy (object)
table_index = 25 # select the index to center the image
index_ellipses = [25, 16]  # draw the ellipses of the index
nrows = 4
ncols = 3

apertures_directory="apertures"
if not os.path.exists(apertures_directory):
    os.makedirs(apertures_directory)

obj_directory = os.path.join(apertures_directory, f"{loaded_gal_sample.iloc[obj_index]['Object']}")
if not os.path.exists(obj_directory):
    os.makedirs(obj_directory)

fig, axs = plt.subplots(nrows, ncols, figsize=(24, 36))
for ax_row in axs:
    for ax in ax_row:
        ax.axis('off')
    
fig.suptitle(f"{loaded_gal_sample.iloc[obj_index]['Object']}", fontsize=22)

for j in range(0,len(frequencies)):
    row = j // ncols
    col = j % ncols

    try:
        hdu = conn.get_cut(select_table[table_index]['RA'], select_table[table_index]['DEC'], 200, new_frequencies[j])
        w = WCS(hdu[1].header)
        gc = aplpy.FITSFigure(hdu,hdu=1, figure=fig, subplot=(nrows, ncols, j + 1))
        gc.show_grayscale(invert=True)
    except Exception as e:
        print("Exception:", str(e))
        continue
    
    for index in index_ellipses:    
        a = select_table['A'][index]
        b = select_table['B'][index]
        k = select_table['KRON_RADIUS'][index]

        gc.show_ellipses(select_table['RA'][index], select_table['DEC'][index],
                        width=(k * a), # major axis
                        height=(k * b), # minor axis
                        angle=-Angle(select_table['THETA'][index], unit=u.degree),
                        edgecolor='red')
        
        plt.title(f"{new_frequencies[j]}", fontsize=18)    
        aperture_filename = os.path.join(obj_directory, f"AUTO.png")
        plt.tight_layout()
plt.close()
fig.savefig(aperture_filename)

# %%

# PETRO aperture (Petro radius)

select_table = loaded_list[0] # select the table
index_select=10 # select the index of the table
hdu = conn.get_cut(select_table[index_select]['RA'], select_table[index_select]['DEC'], 200, 'R')
w = WCS(hdu[1].header)
gc = aplpy.FITSFigure(hdu,hdu=1)
gc.show_grayscale(invert=True)

a = select_table['A'][index_select]
b = select_table['B'][index_select]
p = select_table['PETRO_RADIUS'][index_select]

gc.show_ellipses(select_table['RA'][index_select], select_table['DEC'][index_select],
                width=(p * a), # major axis
                height=(p * b), # minor axis
                angle=-Angle(select_table['THETA'][index_select], unit=u.degree),
                edgecolor='blue')

# %%

# PETRO aperture (Petro radius) for all the filters

obj_index = 138
select_table = loaded_list[obj_index] # select the table of the galaxy (object)
table_index = 25 # select the index to center the image
index_ellipses = [25, 16]  # draw the ellipses of the index
nrows = 4
ncols = 3

apertures_directory="apertures"
if not os.path.exists(apertures_directory):
    os.makedirs(apertures_directory)

obj_directory = os.path.join(apertures_directory, f"{loaded_gal_sample.iloc[obj_index]['Object']}")
if not os.path.exists(obj_directory):
    os.makedirs(obj_directory)

fig, axs = plt.subplots(nrows, ncols, figsize=(24, 36))
for ax_row in axs:
    for ax in ax_row:
        ax.axis('off')
    
fig.suptitle(f"{loaded_gal_sample.iloc[obj_index]['Object']}", fontsize=22)

for j in range(0,len(frequencies)):
    row = j // ncols
    col = j % ncols

    try:
        hdu = conn.get_cut(select_table[table_index]['RA'], select_table[table_index]['DEC'], 200, new_frequencies[j])
        w = WCS(hdu[1].header)
        gc = aplpy.FITSFigure(hdu,hdu=1, figure=fig, subplot=(nrows, ncols, j + 1))
        gc.show_grayscale(invert=True)
    except Exception as e:
        print("Exception:", str(e))
        continue
    
    for index in index_ellipses:    
        a = select_table['A'][index]
        b = select_table['B'][index]
        p = select_table['PETRO_RADIUS'][index]

        gc.show_ellipses(select_table['RA'][index], select_table['DEC'][index],
                        width=(p * a), # major axis
                        height=(p * b), # minor axis
                        angle=-Angle(select_table['THETA'][index], unit=u.degree),
                        edgecolor='blue')
        
        plt.title(f"{new_frequencies[j]}", fontsize=18)    
        aperture_filename = os.path.join(obj_directory, f"PETRO.png")
        plt.tight_layout()
plt.close()
fig.savefig(aperture_filename)

# %%

hdu.info()

# %%

hdu[1].header