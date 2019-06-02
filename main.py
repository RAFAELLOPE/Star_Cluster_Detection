#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 11:23:54 2019

@author: rafa
"""
import sys
import numpy as np
import pandas as pd
import GaiaData as gd
import DiasCatalog as dc
import datetime
# Third-party dependencies
from astropy import units as u
from astropy.coordinates import Angle, Longitude, Latitude, Distance
from astropy.coordinates import SkyCoord
from astropy.table import Table
# Set up matplotlib and use a nicer set of plot parameters
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.visualization import astropy_mpl_style
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
#Sklear algorithm
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_score
#To suppress warnings
import warnings
warnings.filterwarnings("ignore")

# DATA EXTRACTION #############################################################
def extract_data(r, err_pos, min_g_mag, cluster=None, coordinates=None):
    #Define variables needed to make the query 
    data = None
    dias_catalog = dc.DiasCatalog()
    if (coordinates != None) or (cluster != None):
        if coordinates == None:
            point = dias_catalog.get_ra_dec(cluster)
        else:
            point = coordinates
        #Create the GAIA query
        attributes = ['source_id', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmdec','phot_g_mean_mag']
        gaia_query = gd.GaiaData(attributes)
        data = gaia_query.astrometric_query_circle(point, r, err_pos, min_g_mag)
    return(data, point)

# DATA PREPROCESSING ##########################################################
def preprocessing(data, sample_factor = None):
    #Preprocessing
    data_metrics = data
    data_metrics['parallax'] = data_metrics['parallax'].to(u.parsec, equivalencies=u.parallax())
    #Adapt data metrics to a numpy array
    np_data_metrics = np.transpose(np.array([data_metrics['ra'], data_metrics['dec'], data_metrics['parallax'],
                                             data_metrics['pmra'], data_metrics['pmdec']]))
    #Sample data
    if(sample_factor != None):
        np.random.seed(0)
        idx = np.random.choice(np_data_metrics.shape[0], size=int(sample_factor*np_data_metrics.shape[0]),
                               replace= False)
        np_data_metrics = np_data_metrics[idx,:]
    #Change coordinates from Spherical to Cartesian coordinate system
    ra = Longitude(np_data_metrics[:,0], unit=u.deg)
    ra.wrap_angle = 180 * u.deg
    dec = Latitude(np_data_metrics[:,1], unit=u.deg)
    dist = Distance(np_data_metrics[:,2], unit=u.parsec)
    sphe_coordinates = SkyCoord(ra, dec, distance = dist, frame='icrs', representation_type='spherical')
    cart_coordinates = sphe_coordinates.cartesian
    #Adapt data to normalize it correctly
    data_sphe_adapted = np.transpose(np.array([sphe_coordinates.ra, sphe_coordinates.dec, sphe_coordinates.distance]))
    data_cart_adapted = np.transpose(np.array([cart_coordinates.x, cart_coordinates.y, cart_coordinates.z]))
    data_pm_adapted = np_data_metrics[:,3:5]   
    data_all_adapted = np.append(data_cart_adapted, data_pm_adapted, axis=1)
    return(data_sphe_adapted, data_cart_adapted, data_all_adapted)


def get_distances_for_ML(X, Y):
    distance = 0
    ra1 = X[0]*u.deg
    ra2 = Y[0]*u.deg
    dec1 = X[1]*u.deg
    dec2 = Y[1]*u.deg
    if(len(X) == 3):
        dist1 = X[2]*u.parsec
        dist2 = Y[2]*u.parsec
        point1 = SkyCoord(ra1, dec1, dist1)
        point2 = SkyCoord(ra2, dec2, dist2)
        distance = point1.separation_3d(point2)
    else:
        point1 = SkyCoord(ra1, dec1)
        point2 = SkyCoord(ra2, dec2)
        distance = point1.separation(point2)
    return(distance.value)


def get_distance_matrix(data, scale=False, metric='euclidean'):
    if scale:
        scaler = MinMaxScaler()
        data_scaled = scaler.fit_transform(data)
    else:
        data_scaled = data
    if metric != 'euclidean':
        dist_matrix = pairwise_distances(data_scaled, metric=get_distances_for_ML, n_jobs=-1)
    else:
        dist_matrix = pairwise_distances(data_scaled, metric='euclidean', n_jobs=-1)
    return(dist_matrix)


################# DATA EVALUATON ####################################################

def DBSCAN_result(param_scores, dist_matrix, data, center, radius):
    dias_catalog = dc.DiasCatalog()
    cluster_centers = []
    sorted_param_scores = param_scores.sort_values(by=['local_score'], ascending=False)
    opt_epsilon = sorted_param_scores.iloc[0,0]
    opt_min_pts = sorted_param_scores.iloc[0,1]
    db = DBSCAN(eps=opt_epsilon, min_samples=opt_min_pts, metric='precomputed', n_jobs=-1).fit(dist_matrix)
    labels = db.labels_
    for i in range(len(set(labels))-1):
        cluster = data[np.where(labels == i)]
        cluster_center = cluster.mean(axis=0)
        cluster_centers.append(cluster_center)
    matches = dias_catalog.get_match(center, radius, cluster_centers)
    for m in matches:
        print('Cluster found: %s'%(m[0]))
    plot_clusters(data, labels)
    plot_score(param_scores)


def DBSCAN_eval(data, center, r, sample_factor, ftype, dim3, eps_range, min_samples_range, scale=False, metric = 'euclidean', pm=False):
    data_sphe, data_cart, data_all = preprocessing(data, sample_factor)
    data_search = data_sphe[:,:2]
    if ftype == 'cart':
        if dim3:
            data = data_cart
        else:
            data = data_cart[:,:2]
    elif ftype == 'sphe':
        if dim3:
            data = data_sphe
        else:
            data = data_sphe[:,:2]
    elif ftype == 'pm':
        data = data_all
    else:
        print('Specify tpye of dataset: cart, sphe, pm')
        sys.exit()
    
    dist_matrix = get_distance_matrix(data, scale, metric)
    
    if pm:
        data = data[:,:3]
    
    dias_catalog = dc.DiasCatalog()
    num_cum = len(dias_catalog.get_clusters(center, r))
    tmp_sscores = []
    matches = []
    sscores = pd.DataFrame(columns=['epsilon', 'minpts', 'local_score'])
    cluster_centers = []
    for eps in eps_range:
        for min_samples in min_samples_range:
            db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed", n_jobs=-1).fit(dist_matrix)
            labels = db.labels_
            for i in range(len(set(labels))-1):
                cluster = data_search[np.where(labels == i)]
                cluster_center = cluster.mean(axis=0)
                cluster_centers.append(cluster_center)
            matches = dias_catalog.get_match(center, r, cluster_centers)
            tmp_sscores.append(eps)
            tmp_sscores.append(min_samples)
            if(len(matches) > 0 and len(set(labels)) > 1):
                tmp_sscores.append(len(matches)/(num_cum + (len(cluster_centers)-len(matches))))
            else:
                tmp_sscores.append(0.0)
            sscores.loc[len(sscores)] = tmp_sscores
            tmp_sscores = []
            cluster_centers = []
    return(sscores, dist_matrix, data_search)

# DATA PLOTTING ###############################################################
def plot_clusters(X, labels):
    # Black removed and is used for noise instead.
    size = 6.0
    f = plt.figure(figsize=(20,20))
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 0.5, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
            size = 3.5

        class_member_mask = (labels == k)

        xy = X[class_member_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=size)
        size=6.0

    plt.title('Estimated number of clusters: %5d' % int(len(set(labels))-1))
    plt.xlabel("Right Ascension (deg)")
    plt.ylabel("Declination (deg)")
    f.savefig('plot_cluster_execution_%s.png'%(datetime.datetime.now()))

    
    
def plot_score(param_scores):
    np_param_scores = param_scores.values
    eps = np.sort(np.array(list(set(param_scores['epsilon']))))
    Nmin = np.sort(np.array(list(set(param_scores['minpts']))))
    Z = np.empty((len(Nmin), len(eps)))
    fig, ax = plt.subplots(constrained_layout = True)
    X, Y = np.meshgrid(eps, Nmin)
    for i, n in enumerate(Nmin):
        for j, e in enumerate(eps):
            Z[i,j] = np_param_scores[np.where((param_scores['epsilon'] == e) & (param_scores['minpts'] == n)),2]
    extend = "neither"
    cmap = plt.cm.get_cmap('hot')
    CS = ax.contourf(X,Y,Z, cmap=cmap, extend=extend)
    fig.colorbar(CS)
    ax.set_xlabel('Epsilon')
    ax.set_ylabel('Nmin')
    ax.set_title('DBSCAN matchin M')
    fig.savefig('plot_score_execution_%s.png'%(datetime.date.today()))


def check_errors(f):
    def check(*args, **kwargs):
        try:
            f(*args, **kwargs)
        except:
            print("An error has ocurred")
    return(check)

@check_errors
def plot_bar(*args):
    if args[1] == 'eps':
        df = args[0][['epsilon','local_score']]
        df['bucket'] = pd.cut(df.epsilon, args[2])
        title = 'Máximo M en función de epsilon'
        x_label = 'epsilon'
    else:
        df = args[0][['minpts','local_score']]
        df['bucket'] = pd.cut(df.minpts, args[2])
        title = 'Máximo M en función de Nmin'
        x_label = 'Nmin'
    newdf = df[['bucket','local_score']].groupby('bucket').max()   
    ax = newdf.plot(kind='bar', title=title, colormap = 'jet', fontsize=7, legend=False)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Max(M)")
    ax.grid()    



### MAIN ######################################################################
def process_line(ftype, params):
    ra = params['ra']
    dec = params['dec']
    radius = params['r']
    err_pos = params['err_pos']
    min_g_mag = params['min_mag_g']
    sample_factor = params['sample']
    scale = params['norm'] == 'X'
    metric = params['distance']
    dim3 = params['dim3'] == 'X'
    eps_min = params['eps_min']
    eps_max = params['eps_max']
    eps_num = params['eps_num']
    min_pts_min = params['min_pts_min']
    min_pts_max = params['min_pts_max']
    min_pts_num = params['min_pts_num']
    eps_range = np.linspace(eps_min, eps_max, eps_num)
    min_pts_range = np.linspace(min_pts_min, min_pts_max, min_pts_num)
    
    center = [ra, dec]
    data, center = extract_data(radius, err_pos, min_g_mag, coordinates = center)
    param_scores, dist_matrix, data_search = DBSCAN_eval(data, center, radius, sample_factor, ftype, dim3, eps_range, min_pts_range, scale, metric)

    return(param_scores, dist_matrix, data_search, center, radius)
    
    

def main():
    if len(sys.argv) == 3:
        file = sys.argv[1]
        ftype = sys.argv[2]
        try:
            params = pd.read_csv(file)
        except:
            print('The file does not exist')
            sys.exit(0)
        for i in range(len(params)):
            param_scores, dist_matrix, data_search, center, radius = process_line(ftype, params.iloc[i,:])
            param_scores.to_csv('execution%s_%s.csv'%(str(i), datetime.date.today()))
            DBSCAN_result(param_scores, dist_matrix, data_search, center, radius)



if __name__ == "__main__":
    main()