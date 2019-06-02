#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 14:57:32 2018

@author: rafa
"""

from astroquery.vizier import Vizier
from astropy.coordinates import Angle, Longitude, Latitude
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt


class DiasCatalog(object):
    def __init__(self, *attributes):
        self.attributes = attributes
        Vizier.ROW_LIMIT = -1
        dias_catalog = Vizier.get_catalogs('/B/ocl/clusters')
        if attributes == ():
            self.table = dias_catalog['B/ocl/clusters']
        else:
            self.table = dias_catalog['B/ocl/clusters'][attributes]
        
    
    def get_ra_dec(self, cluster):
        """Method to get right ascension and declination from a cluster
            
        Parameters
        --------
        cluster: Name of the cluster to get its position
        Returns
        --------
        SkyCoord object containing ra and dec in degrees
        """
        try:
            row = self.table[np.where(self.table['Cluster'] == cluster)]
            ra = Angle(row['RAJ2000'], unit=u.hourangle)
            dec = Angle(row['DEJ2000'], unit=u.deg)
            return([float(ra.deg), float(dec.deg)])
        except:
            print('Cluster does not exist')

    def get_cluster_names(self):
        """Method to get the name of the clusters in the catalog
        Returns
        --------
        print out name of the clusters in the catalog
        """
        print(self.table['Cluster'])
    
    def get_diam(self, cluster):
        """Method to get right ascension and declination from a cluster
            
        Parameters
        --------
        cluster: Name of the cluster to get its diameter
        Returns
        --------
        SkyCoord object containing ra and dec in degrees
        """
        try:
            row = self.table[np.where(self.table['Cluster'] == cluster)]
            diam = float(row['Diam'])*u.arcmin
            return(diam.to('deg').value)
        except:
            print('Cluster does not exist')
    
    
    def get_cluster_coordinates(self):
        ra = Angle(self.table['RAJ2000'], unit=u.hourangle)
        ra = ra.deg
        ra = Longitude(ra, unit=u.deg)
        dec = Latitude(self.table['DEJ2000'], unit=u.deg )
        return(SkyCoord(ra, dec))
    
    def plot_all(self):
        coordinates = self.get_cluster_coordinates()
        ra = coordinates.ra
        dec = coordinates.dec
        fig, ax = plt.subplots()
        ax.scatter(ra, dec, s=0.5)
        #for i, txt in enumerate(np.array(self.table['Cluster'], dtype='S')):
        #    ax.annotate(txt, (ra[i].value, dec[i].value))
        plt.show()
    
    def get_clusters(self, center, r):
        clusters = []
        ra_center = Longitude(center[0], unit=u.deg)
        dec_center = Latitude(center[1], unit=u.deg)
        center = SkyCoord(ra_center, dec_center)
        for i in range(self.table['Cluster'].shape[0]):
            ra = Angle(self.table[i]['RAJ2000'], unit=u.hourangle)
            ra = ra.deg
            ra = Longitude(ra, unit=u.deg)
            dec = Latitude(self.table[i]['DEJ2000'], unit=u.deg)
            cluster = SkyCoord(ra, dec)
            diam = float(self.table[i]['Diam'])*u.arcmin
            if(np.isnan(diam)):
                diam = 0.05
            else:
                diam = diam.to('deg').value
            if(cluster.separation(center).value <= r):
                clusters.append((self.table['Cluster'][i], cluster, diam))
        return(clusters)
    
    def get_rate(self, center, r, clusters):
        '''
        Method to get a metric between 0 and 1 to measure the success of our 
        algorithm.
        
        Parameters
        ---------
        center : List
                 Center of the sky
        r : float
            Radius of the sky
        
        cluster: List
                 Center point of a cluster detected
                 
        Returns
        -------
        rate: float
              Metric between 0 and 1 that measures how good our algorithm worked
        '''
        num_cum = int(0)
        matches = int(0)
        rate = 0.0
        ra = Longitude(center[0], unit=u.deg)
        dec = Latitude(center[1], unit=u.deg)
        center = SkyCoord(ra, dec)
#        ra_cluster = Longitude(cluster[0], unit=u.deg)
#        dec_cluster = Latitude(cluster[1], unit=u.deg)
#        cluster = SkyCoord(ra_cluster, dec_cluster)
        for i in range(self.table['Cluster'].shape[0]):
            ra_dias = Angle(self.table[i]['RAJ2000'], unit=u.hourangle)
            ra_dias = ra_dias.deg
            ra_dias = Longitude(ra_dias, unit=u.deg)
            dec_dias = Latitude(self.table[i]['DEJ2000'], unit=u.deg)
            cluster_dias = SkyCoord(ra_dias, dec_dias)
            diam = float(self.table[i]['Diam'])*u.arcmin
            if(np.isnan(diam)):
                diam = 0.05
            else:
                diam = diam.to('deg').value
            if(cluster_dias.separation(center).value <= r):
                num_cum += 1
                for cluster in clusters:
                    ra_cluster = Longitude(cluster[0], unit=u.deg)
                    dec_cluster = Latitude(cluster[1], unit=u.deg)
                    cluster = SkyCoord(ra_cluster, dec_cluster)
                    if(cluster_dias.separation(cluster).value <= diam/2.0):
                        matches += 1
                        break
        try:
            rate = matches/num_cum
            print('Number of cluster in the sky %5d' % num_cum)
            print('Number of matches: %5d' % matches)
        except:
            print('No star cluster detected in Dias catalog in the part of\n the sky specified')
        return(rate)
        
        
        
        
    def get_match(self, center, r, clusters_detected):
        existing_clusters = self.get_clusters(center, r)
        tot_matches = []
        for cluster in clusters_detected:
            ra_cluster_detected = Longitude(cluster[0], unit=u.deg)
            dec_cluster_detected = Latitude(cluster[1], unit=u.deg)
            cluster_detected = SkyCoord(ra_cluster_detected, dec_cluster_detected)
            l_tmp = [match for match in existing_clusters if match[1].separation(cluster_detected).value <= match[2]/2.0]
            tot_matches = list(set(tot_matches + l_tmp))
        return tot_matches
            
        
        

         
#dias_catalog = DiasCatalog()
#center = [93.6084, 14.6927]
#r = 5.0
#cumulos_detectados = [[89.5583, 16.55944], [92.1, 13.965]]
#matches = dias_catalog.get_match(center, r, cumulos_detectados)
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
        
        
        