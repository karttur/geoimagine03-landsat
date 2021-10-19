'''
Created on 7 Oct 2018
Last Updated 12 Feb 2021

@author: thomasgumbricht
'''

__version__ = '0.9.0'

VERSION = tuple( int(x) for x in __version__.split('.') )

metadataD = { 'name':'landsat', 'author':'Thomas Gumbricht', 'author_email':'thomas.gumbricht@gmail.com',
             'title':'landsat', 'label':'Landsat specific processing.',
             'image':'avg-trmm-3b43v7-precip_3B43_trmm_2001-2016_A'}