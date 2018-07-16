# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 09:53:31 2018

@author: Guillermo
"""


def remove_none_elements_from_list(list):
    return [e for e in list if e != '']

f1=open('GMRX0157.607','r')
texto=f1.readlines()
traz=open('trazabilidad.txt','a')

for i in range(20,len(texto)):
    REFSV = remove_none_elements_from_list(texto[i].split(' '))
    traz.write(REFSV[7]+", ")
    
traz.close()