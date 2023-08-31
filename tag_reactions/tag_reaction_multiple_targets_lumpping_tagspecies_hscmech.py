#!/usr/bin/env python
"""
To tag reactions leading to loss and production of Ox and HOx individually. 

Usage:
python tag_reaction_multiple_targets.py

Revision History:
2023-01-30 - hansen.cao@york.ac.uk - initial version
"""

import sys
import os
import re
import datetime
import numpy as np


def extractReactions(familynames,familynametags,
                     dict_reactant1,dict_reactant2,dict_product1,dict_product2,
                     dict_altproduct2,dict_othaltproduct2,dict_product3,
                     dict_exlproduct1,dict_exlreactant1,dict_exlreactant2,dict_exlreactant3):
    """
    Extract reactions leading to Ox loss in the gckpp_Monitor.F90 module.
    
    Args:
    -----
    familynames (str) : Names of the families you are intested in
    """

    os.system('rm reactions_need_to_be_tagged_'+familynametags[0]+'.csv')
    # Filenames
    file_reactions = 'gckpp_Monitor.F90'
    file_species   = 'gckpp_Parameters.F90'
    file_updated   = 'reactions_need_to_be_tagged_'+familynametags[0]+'.csv'

    # Rread file with reactions
    with open(file_reactions, 'r') as f:
        lines = f.readlines()

    # Read files with species
    with open(file_species, 'r') as f:
        params = f.readlines()

    # Open gckpp_Util.F90
    fo = open(file_updated,"a")

    
    nn   = 0
    ntot = 0

    # get all lines with a reaction
    rxt = [i for i in lines if '-->' in i]

    #print(rxt)
    #exit()
    # walk through all reactions
    irct = 0
    for irxt in rxt:
        irct += 1
        irxt = irxt.replace("\n", "")
        irxt = irxt.replace("'", "")
        #---------------------------------------------------
        # Also remove extra text for safety's sake
        #  -- Bob Yantosca (03 May 2022)
        irxt = irxt.replace("/)", "")
        irxt = irxt.replace(",", "")
        irxt = irxt.replace("&", "")
        irxt = irxt.replace("!", "")
        irxt = irxt.replace("up to ", "")
        #---------------------------------------------------
        spl = irxt.split('-->')
        #print(spl[0])
        #print(spl[1])
        #exit()
        # get reaction number and reactant if OH is on left-hand side of reaction
        #if ' LOx ' in spl[1]:
        for name in familynames:
          if ' '+name+' ' in spl[1]:
            # write to file
            nn   += 1
            ntot += 1
            #fo.write(istr)
            #fo.write(irxt.lstrip() + '\n')
            find=0
            tmpstr=''
            
            for key in dict_reactant1:
                 if (find == 1):
                    break
                 #print(key)
                 #print(dict_reactant1[key])
                 #print(type(dict_reactant1[key]))
                 if (dict_reactant1[key][0] in spl[0] and \
                    (not dict_exlreactant1[key][0] in spl[0]) and \
                    (not dict_exlreactant2[key][0] in spl[0]) and \
                    (not dict_exlreactant3[key][0] in spl[0]) and \
                     dict_product1[key][0] in spl[1] and \
                     dict_product3[key][0] in spl[1] and \
                    (not dict_exlproduct1[key][0] in spl[1]) and \
                    (dict_product2[key][0] in spl[1] or dict_altproduct2[key][0] in spl[1] or dict_othaltproduct2[key][0] in spl[1])):
                    for val in dict_reactant2[key]:
                      if val in spl[0]: 
                
                        tmpstr = key
                        find = 1
                        print('find ',key,'for EQ ',irxt)
                        break
              

            if (find == 0):
                tmpstr = familynametags[0]+'OTHPTW'
                print('not find key for EQ ',irxt)
            fo.write(irxt +' LUMPSPC ' +tmpstr+ '\n')
            break


    # all done
    print('Reactivity consists of '+str(ntot)+' reactions')
    print('Written to '+file_updated)
    fo.close()


def getPathways(familynametag,mechName):

    if (familynametag in 'TOxP'):
      dict_reactant1={familynametag+'FromNORO2':[' NO '],familynametag+'FromNOMO2':[' NO '],
                      familynametag+'FromNOHO2':[' NO '],
                      familynametag+'FromHOBrSS':[' '],familynametag+'FromHOClSS':[' ']}
  
      dict_reactant2={familynametag+'FromNORO2':[' MCO3 ',' ETO2 ',' OTHRO2 ',' A3O2 ',
                                      ' PO2 ', ' R4O2 ',' R4N1 ', ' ATO2 ',
                                      ' KO2 ', ' B3O2 ',' PRN1 ', ' RCO3 ',
                                    ' CH2OO ',' CH3CHOO ', ' PIO2 ',' LIMO2 ',
                                    ' OLNN ', ' OLND ', ' IHOO1 ',' IHOO4 ',
                                    ' IHPOO1 ',' IHPOO2 ', ' IHPOO3 ', ' IEPOXAOO ',
                                    ' IEPOXBOO ',' ICHOO ',' HPALD1OO ',' HPALD2OO ',
                                    ' ISOPNOO1 ',' ISOPNOO2 ',' IDHNDOO1 ',' IDHNDOO2 ',
                                    ' IDHNBOO ', ' INO2B ', ' INO2D ',' IHPNBOO ',
                                    ' IHPNDOO ', ' ICNOO ',' IDNOO ',' C4HVP1 ',
                                    ' C4HVP1 ', ' C4HVP2 ',' MVKOHOO ',' MCROHOO ',
                                    ' MACR1OO ', ' MACRNO2 ',' ETOO ',' ETHN ',
                                    ' AROMRO2 ',' BZCO3 ', ' BENZO2 ' ],
                      familynametag+'FromNOMO2':[' MO2 '],
                      familynametag+'FromNOHO2':[' HO2 '],
                      familynametag+'FromHOBrSS':[' HOBrCl ',' HOBrSSC '],
                      familynametag+'FromHOClSS':[' HOClSSA ',' HOClSSC ']}
  
      dict_exlreactant1={familynametag+'FromNORO2':[' None '],familynametag+'FromNOMO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],
                      familynametag+'FromHOBrSS':[' None '],
                      familynametag+'FromHOClSS':[' None ']}
  

      dict_exlreactant2={familynametag+'FromNORO2':[' None '],familynametag+'FromNOMO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],
                      familynametag+'FromHOBrSS':[' None '],
                      familynametag+'FromHOClSS':[' None ']}
  

      dict_exlreactant3={familynametag+'FromNORO2':[' None '],familynametag+'FromNOMO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],
                      familynametag+'FromHOBrSS':[' None '],
                      familynametag+'FromHOClSS':[' None ']}
  

      dict_exlproduct1={familynametag+'FromNORO2':[' None '],familynametag+'FromNOMO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],
                      familynametag+'FromHOBrSS':[' None '],
                      familynametag+'FromHOClSS':[' None ']}
  

      dict_product1={familynametag+'FromNORO2':[' POx '],familynametag+'FromNOMO2':[' POx '],
                      familynametag+'FromNOHO2':[' POx '],
                      familynametag+'FromHOBrSS':[' POx '],
                      familynametag+'FromHOClSS':[' POx ']}
  
  
      dict_product3={familynametag+'FromNORO2':[' POx '],familynametag+'FromNOMO2':[' POx '],
                      familynametag+'FromNOHO2':[' POx '],
                      familynametag+'FromHOBrSS':[' POx '],
                      familynametag+'FromHOClSS':[' POx ']}
  
  
      dict_product2={familynametag+'FromNORO2':[' NO2 '],familynametag+'FromNOMO2':[' NO2 '],
                      familynametag+'FromNOHO2':[' NO2 '],
                      familynametag+'FromHOBrSS':[' HOBr '],
                      familynametag+'FromHOClSS':[' HOCl ']}
  
  
      dict_altproduct2={familynametag+'FromNORO2':[' etc.'],familynametag+'FromNOMO2':[' NO2 '],
                      familynametag+'FromNOHO2':[' NO2 '],
                      familynametag+'FromHOBrSS':[' HOBr '],
                      familynametag+'FromHOClSS':[' HOCl ']}
  
  
      dict_othaltproduct2={familynametag+'FromNORO2':[' etc.'],familynametag+'FromNOMO2':[' NO2 '],
                      familynametag+'FromNOHO2':[' NO2 '],
                      familynametag+'FromHOBrSS':[' HOBr '],
                      familynametag+'FromHOClSS':[' HOCl ']}
  
  

    elif (familynametag in 'TOxL'):
       dict_reactant1={familynametag+'FromO3H2O':[' O1D '],familynametag+'FromO3HO2':[' O3 '],
                       familynametag+'FromO3OH':[' O3 '],
                       familynametag+'FromHOBrhv':[' HOBr '],
                       familynametag+'FromHOBrHX':[' HOBr '],familynametag+'FromBrOXO':[' BrO '],
                       familynametag+'FromHOClHX':[' HOCl '],
                       familynametag+'FromBrOOH':[' BrO '],
                       familynametag+'FromClxOyOH':[' OH '],
                       familynametag+'FromHOIhv':[' HOI '], familynametag+'FromOIOhv':[' OIO '],
                       familynametag+'FromIOXO':[' IO '],
                       familynametag+'FromHOClhv':[' HOCl '],
                       familynametag+'FromCl2O2hv':[' Cl2O2 '],
                       familynametag+'FromClOMO2':[' ClO '],
                       familynametag+'FromClNOyHX':[' '],
                       familynametag+'FromHOClToSS':[' HOCl '], familynametag+'FromHOBrToSS':[' HOBr '],
                       familynametag+'FromBrNO3ToSS':[' BrNO3 '],
                       familynametag+'FromBrNO3HCl':[' BrNO3 '],
                       familynametag+'FromOwHOBrBrO':[' O '],
                       familynametag+'FromIxOyNzHOIToSS':[' '],
                       familynametag+'FromIxOyNzHOITodSS':[' '],
                       familynametag+'FromClOClO':[' 2 ClO '],
                       familynametag+'FromO3wH':[' O3 '],
                       familynametag+'FromOHwO':[' OH '],
                       familynametag+'FromNO2wO':[' NO2 '],
                       familynametag+'FromHO2wO':[' HO2 ']}

   
   
       dict_reactant2={familynametag+'FromO3H2O':[' H2O '],
                       familynametag+'FromO3HO2':[' HO2 '],
                       familynametag+'FromO3OH':[' OH '],
                       familynametag+'FromHOBrhv':[' HOBr '],
                       familynametag+'FromHOBrHX':[' HBr ',' HCl '],familynametag+'FromBrOXO':[' BrO ',' ClO ',' IO '],
                       familynametag+'FromHOClHX':[' HBr ',' HCl '],
                       familynametag+'FromBrOOH':[' OH '],
                       familynametag+'FromClxOyOH':[' ClO ',' OClO ',' Cl2O2 '],
                       familynametag+'FromHOIhv':[' HOI '],familynametag+'FromOIOhv':[' OIO '],
                       familynametag+'FromIOXO':[' ClO '],
                       familynametag+'FromHOClhv':[' HOCl '], familynametag+'FromClOMO2':[' MO2 '],
                       familynametag+'FromCl2O2hv':[' Cl2O2 '],
                       familynametag+'FromClNOyHX':[' ClNO3 ',' ClNO2 '],
                       familynametag+'FromHOClToSS':[' HOCl ',' SALACL ', ' SALCCL ',' SO2 ',' BrSALA ',' BrSALC '], 
                       familynametag+'FromHOBrToSS':[' HOBr ',' SALACL ', ' SALCCL ',' SO2 ',' BrSALA ',' BrSALC '],
                       familynametag+'FromBrNO3ToSS':[' BrNO3 '],
                       familynametag+'FromBrNO3HCl':[' HCl '],
                       familynametag+'FromOwHOBrBrO':[' HOBr ', ' BrO '],
                       familynametag+'FromIxOyNzHOITodSS':[' HOI ',' I2O2 ',' I2O3 ',' I2O4 ',' IONO ', ' IONO2 '],
                       familynametag+'FromIxOyNzHOIToSS':[' HOI ',' IONO ', ' IONO2 '],
                       familynametag+'FromClOClO':[' 2 ClO '],
                       familynametag+'FromO3wH':[' H '],
                       familynametag+'FromOHwO':[' O '],
                       familynametag+'FromNO2wO':[' O '],
                       familynametag+'FromHO2wO':[' O ']}
   
       dict_exlreactant1={familynametag+'FromO3H2O':[' None '],familynametag+'FromO3HO2':[' None '],
                       familynametag+'FromO3OH':[' None '],
                       familynametag+'FromHOBrhv':[' None '],
                       familynametag+'FromHOBrHX':[' None '],familynametag+'FromBrOXO':[' O '],
                       familynametag+'FromHOClHX':[' None '],
                       familynametag+'FromBrOOH':[' None '],
                       familynametag+'FromClxOyOH':[' None '],
                       familynametag+'FromHOIhv':[' None '], familynametag+'FromOIOhv':[' None '],
                       familynametag+'FromIOXO':[' None '],
                       familynametag+'FromHOClhv':[' None '], familynametag+'FromClOMO2':[' None '],
                       familynametag+'FromClNOyHX':[' None '],
                       familynametag+'FromHOClToSS':[' HCl '], familynametag+'FromHOBrToSS':[' HCl '],
                       familynametag+'FromCl2O2hv':[' None '],
                       familynametag+'FromBrNO3ToSS':[' HCl '],
                       familynametag+'FromBrNO3HCl':[' None '],
                       familynametag+'FromOwHOBrBrO':[' None '],
                       familynametag+'FromIxOyNzHOIToSS':[' None '],
                       familynametag+'FromIxOyNzHOITodSS':[' None '],
                       familynametag+'FromClOClO':[' None '],
                       familynametag+'FromO3wH':[' None '],
                       familynametag+'FromOHwO':[' None '],
		       familynametag+'FromNO2wO':[' None '],
                       familynametag+'FromHO2wO':[' None ']}
  
       dict_exlreactant2={familynametag+'FromO3H2O':[' None '],familynametag+'FromO3HO2':[' None '],
                       familynametag+'FromO3OH':[' None '],
                       familynametag+'FromHOBrhv':[' None '],
                       familynametag+'FromHOBrHX':[' None '],familynametag+'FromBrOXO':[' O '],
                       familynametag+'FromHOClHX':[' None '],
                       familynametag+'FromBrOOH':[' None '],
                       familynametag+'FromClxOyOH':[' None '],
                       familynametag+'FromHOIhv':[' None '], familynametag+'FromOIOhv':[' None '],
                       familynametag+'FromIOXO':[' None '],
                       familynametag+'FromHOClhv':[' None '], familynametag+'FromClOMO2':[' None '],
                       familynametag+'FromCl2O2hv':[' None '],
                       familynametag+'FromClNOyHX':[' None '],
                       familynametag+'FromHOClToSS':[' HBr '], familynametag+'FromHOBrToSS':[' HBr '],
                       familynametag+'FromBrNO3ToSS':[' HCl '],
                       familynametag+'FromBrNO3HCl':[' None '],
                       familynametag+'FromOwHOBrBrO':[' None '],
                       familynametag+'FromIxOyNzHOIToSS':[' None '],
                       familynametag+'FromIxOyNzHOITodSS':[' None '],
                       familynametag+'FromClOClO':[' None '],
                       familynametag+'FromO3wH':[' None '],
                       familynametag+'FromOHwO':[' None '],
		       familynametag+'FromNO2wO':[' None '],
                       familynametag+'FromHO2wO':[' None ']}
  
       dict_exlreactant3={familynametag+'FromO3H2O':[' None '],familynametag+'FromO3HO2':[' None '],
                       familynametag+'FromO3OH':[' None '],
                       familynametag+'FromHOBrhv':[' None '],
                       familynametag+'FromHOBrHX':[' None '],familynametag+'FromBrOXO':[' O '],
                       familynametag+'FromHOClHX':[' None '],
                       familynametag+'FromBrOOH':[' None '],
                       familynametag+'FromClxOyOH':[' None '],
                       familynametag+'FromHOIhv':[' None '], familynametag+'FromOIOhv':[' None '],
                       familynametag+'FromIOXO':[' None '],
                       familynametag+'FromHOClhv':[' None '], familynametag+'FromClOMO2':[' None '],
                       familynametag+'FromCl2O2hv':[' None '],
                       familynametag+'FromClNOyHX':[' None '],
                       familynametag+'FromHOClToSS':[' O '], familynametag+'FromHOBrToSS':[' O '],
                       familynametag+'FromBrNO3ToSS':[' HCl '],
                       familynametag+'FromBrNO3HCl':[' None '],
                       familynametag+'FromOwHOBrBrO':[' None '],
                       familynametag+'FromIxOyNzHOIToSS':[' None '],
                       familynametag+'FromIxOyNzHOITodSS':[' None '],
                       familynametag+'FromClOClO':[' None '],
                       familynametag+'FromO3wH':[' None '],
                       familynametag+'FromOHwO':[' None '],
		       familynametag+'FromNO2wO':[' None '],
                       familynametag+'FromHO2wO':[' None ']}
  
   
   
       dict_exlproduct1={familynametag+'FromO3H2O':[' None '],familynametag+'FromO3HO2':[' None '],
                       familynametag+'FromO3OH':[' None '],
                       familynametag+'FromHOBrhv':[' None '],
                       familynametag+'FromHOBrHX':[' None '],familynametag+'FromBrOXO':[' HO2 '],
                       familynametag+'FromHOClHX':[' None '],
                       familynametag+'FromBrOOH':[' None '],
                       familynametag+'FromClxOyOH':[' None '],
                       familynametag+'FromHOIhv':[' None '], familynametag+'FromOIOhv':[' None '],
                       familynametag+'FromIOXO':[' None '],
                       familynametag+'FromHOClhv':[' None '], familynametag+'FromClOMO2':[' None '],
                       familynametag+'FromCl2O2hv':[' None '],
                       familynametag+'FromClNOyHX':[' None '],
                       familynametag+'FromHOClToSS':[' None '], familynametag+'FromHOBrToSS':[' None '],
                       familynametag+'FromBrNO3ToSS':[' None '],
                       familynametag+'FromBrNO3HCl':[' None '],
                       familynametag+'FromOwHOBrBrO':[' None '],
                       familynametag+'FromIxOyNzHOIToSS':[' None '],
                       familynametag+'FromIxOyNzHOITodSS':[' None '],
                       familynametag+'FromClOClO':[' None '],
                       familynametag+'FromO3wH':[' None '],
                       familynametag+'FromOHwO':[' None '],
		       familynametag+'FromNO2wO':[' None '],
                       familynametag+'FromHO2wO':[' None ']}
   
   
       dict_product1={familynametag+'FromO3H2O':[' LOx '],familynametag+'FromO3HO2':[' LOx '],
                       familynametag+'FromO3OH':[' LOx '],
                       familynametag+'FromHOBrhv':[' LOx '],
                       familynametag+'FromHOBrHX':[' LOx '],familynametag+'FromBrOXO':[' LOx '],
                       familynametag+'FromHOClHX':[' LOx '],
                       familynametag+'FromBrOOH':[' LOx '],
                       familynametag+'FromClxOyOH':[' LOx '],
                       familynametag+'FromHOIhv':[' LOx '],familynametag+'FromOIOhv':[' LOx '],
                       familynametag+'FromIOXO':[' LOx '],
                       familynametag+'FromHOClhv':[' LOx '], familynametag+'FromClOMO2':[' LOx '],
                       familynametag+'FromCl2O2hv':[' LOx '],
                       familynametag+'FromClNOyHX':[' LOx '],
                       familynametag+'FromHOClToSS':[' LOx '], familynametag+'FromHOBrToSS':[' LOx '],
                       familynametag+'FromBrNO3ToSS':[' LOx '],
                       familynametag+'FromBrNO3HCl':[' LOx '],
                       familynametag+'FromOwHOBrBrO':[' LOx '],
                       familynametag+'FromIxOyNzHOIToSS':[' LOx '],
                       familynametag+'FromIxOyNzHOITodSS':[' LOx '],
                       familynametag+'FromClOClO':[' LOx '],
                       familynametag+'FromO3wH':[' LOx '],
                       familynametag+'FromOHwO':[' LOx '],
		       familynametag+'FromNO2wO':[' LOx '],
                       familynametag+'FromHO2wO':[' LOx ']}
   
       dict_product3={familynametag+'FromO3H2O':[' LOx '],familynametag+'FromO3HO2':[' LOx '],
                       familynametag+'FromO3OH':[' LOx '],
                       familynametag+'FromHOBrhv':[' Br '],
                       familynametag+'FromHOBrHX':[' LOx '],familynametag+'FromBrOXO':[' LOx '],
                       familynametag+'FromHOClHX':[' LOx '],
                       familynametag+'FromBrOOH':[' Br '],
                       familynametag+'FromClxOyOH':[' LOx '],
                       familynametag+'FromHOIhv':[' I '],familynametag+'FromOIOhv':[' O2 '],
                       familynametag+'FromIOXO':[' O2 '],
                       familynametag+'FromHOClhv':[' Cl '], familynametag+'FromClOMO2':[' ClOO '],
                       familynametag+'FromCl2O2hv':[' Cl '],
                       familynametag+'FromClNOyHX':[' LOx '],
                       familynametag+'FromHOClToSS':[' LOx '], familynametag+'FromHOBrToSS':[' LOx '],
                       familynametag+'FromBrNO3ToSS':[' LOx '],
                       familynametag+'FromBrNO3HCl':[' LOx '],
                       familynametag+'FromOwHOBrBrO':[' LOx '],
                       familynametag+'FromIxOyNzHOIToSS':[' LOx '],
                       familynametag+'FromIxOyNzHOITodSS':[' LOx '],
                       familynametag+'FromClOClO':[' LOx '],
                       familynametag+'FromO3wH':[' OH '],
                       familynametag+'FromOHwO':[' H '],
		       familynametag+'FromNO2wO':[' NO '],
                       familynametag+'FromHO2wO':[' OH ']}
    
  
       dict_product2={familynametag+'FromO3H2O':[' OH '],familynametag+'FromO3HO2':[' OH '],
                       familynametag+'FromO3OH':[' HO2 '],
                       familynametag+'FromHOBrhv':[' OH '],
                       familynametag+'FromHOBrHX':[' Br2 '],familynametag+'FromBrOXO':[' Br '],
                       familynametag+'FromHOClHX':[' Cl2 '],
                       familynametag+'FromBrOOH':[' HO2 '],
                       familynametag+'FromClxOyOH':[' LOx '],
                       familynametag+'FromHOIhv':[' OH '],familynametag+'FromOIOhv':[' I '],
                       familynametag+'FromIOXO':[' I '],
                       familynametag+'FromHOClhv':[' OH '], familynametag+'FromClOMO2':[' ClOO '],
                       familynametag+'FromCl2O2hv':[' ClOO '],
                       familynametag+'FromClNOyHX':[' LOx '],
                       familynametag+'FromHOClToSS':[' HOClSSA '], familynametag+'FromHOBrToSS':[' HOBrCl '],
                       familynametag+'FromBrNO3ToSS':[' HOBrCl '],
                       familynametag+'FromBrNO3HCl':[' BrCl '],
                       familynametag+'FromOwHOBrBrO':[' LOx '],
                       familynametag+'FromIxOyNzHOITodSS':[' ISALA '],
                       familynametag+'FromIxOyNzHOIToSS':[' IBr '],
                       familynametag+'FromClOClO':[' LOx '],
                       familynametag+'FromO3wH':[' O2 '],
                       familynametag+'FromOHwO':[' O2 '],
		       familynametag+'FromNO2wO':[' O2 '],
                       familynametag+'FromHO2wO':[' O2 ']}
    
    
       dict_altproduct2={familynametag+'FromO3H2O':[' OH'],familynametag+'FromO3HO2':[' OH '],
                       familynametag+'FromO3OH':[' HO2 '],
                       familynametag+'FromHOBrhv':[' OH '],
                       familynametag+'FromHOBrHX':[' BrCl '],familynametag+'FromBrOXO':[' Br2 '],
                       familynametag+'FromHOClHX':[' BrCl '],
                       familynametag+'FromBrOOH':[' HO2 '],
                       familynametag+'FromClxOyOH':[' LOx '],
                       familynametag+'FromHOIhv':[' OH '],familynametag+'FromOIOhv':[' I '],
                       familynametag+'FromIOXO':[' Cl '],
                       familynametag+'FromHOClhv':[' OH '], familynametag+'FromClOMO2':[' ClOO '],
                       familynametag+'FromCl2O2hv':[' ClOO '],
                       familynametag+'FromClNOyHX':[' LOx '],
                       familynametag+'FromHOClToSS':[' HOClSSC '], familynametag+'FromHOBrToSS':[' HOBrSSC '],
                       familynametag+'FromBrNO3ToSS':[' HOBrSSC '],
                       familynametag+'FromBrNO3HCl':[' BrCl '],
                       familynametag+'FromOwHOBrBrO':[' LOx '],
                       familynametag+'FromIxOyNzHOITodSS':[' ISALC '],
                       familynametag+'FromIxOyNzHOIToSS':[' ICl '],
                       familynametag+'FromClOClO':[' LOx '],
                       familynametag+'FromO3wH':[' OH '],
                       familynametag+'FromOHwO':[' O2 '],
		       familynametag+'FromNO2wO':[' O2 '],
                       familynametag+'FromHO2wO':[' O2 ']}
    
    
       dict_othaltproduct2={familynametag+'FromO3H2O':[' OH'],familynametag+'FromO3HO2':[' OH '],
                       familynametag+'FromO3OH':[' HO2 '],
                       familynametag+'FromHOBrhv':[' OH '],
                       familynametag+'FromHOBrHX':[' Br2 '],familynametag+'FromBrOXO':[' BrCl '],
                       familynametag+'FromHOClHX':[' Cl2 '],
                       familynametag+'FromBrOOH':[' HO2 '],
                       familynametag+'FromClxOyOH':[' LOx '],
                       familynametag+'FromHOIhv':[' OH '],familynametag+'FromOIOhv':[' I '],
                       familynametag+'FromIOXO':[' ICl '],
                       familynametag+'FromHOClhv':[' OH '], familynametag+'FromClOMO2':[' ClOO '],
                       familynametag+'FromCl2O2hv':[' ClOO '],
                       familynametag+'FromClNOyHX':[' LOx '],
                       familynametag+'FromHOClToSS':[' LOx '], familynametag+'FromHOBrToSS':[' LOx '],
                       familynametag+'FromBrNO3ToSS':[' HOBrCl '],
                       familynametag+'FromBrNO3HCl':[' BrCl '],
                       familynametag+'FromOwHOBrBrO':[' LOx '],
                       familynametag+'FromIxOyNzHOITodSS':[' AERI '],
                       familynametag+'FromIxOyNzHOIToSS':[' ICl '],
                       familynametag+'FromClOClO':[' LOx '],
                       familynametag+'FromO3wH':[' OH '],
                       familynametag+'FromOHwO':[' O2 '],
		       familynametag+'FromNO2wO':[' O2 '],
                       familynametag+'FromHO2wO':[' O2 ']}

    elif (familynametag in 'THOxP'):

      dict_reactant1={familynametag+'FromO1DH2O':[' O1D '],familynametag+'FromNOHO2':[' HO2 '],
                      familynametag+'FromO3HO2':[' HO2 '],familynametag+'FromH2O2hv':[' H2O2 '],
                      familynametag+'FromOVOCsROOHhv':[' '],
                      familynametag+'FromHOBrhv':[' HOBr '],familynametag+'FromHBrwO':[' HBr '],
                      familynametag+'FromHOClhv':[' HOCl '],familynametag+'FromHOIhv':[' HOI ']}
  
      dict_reactant2={familynametag+'FromOVOCsROOHhv':[' MP ',' GLYC ',' PRPN ',' ETP ',
                                      ' RA3P ', ' RB3P ',' R4P ', ' PP ',
                                      ' RP ', ' MAP ',' ATOOH ', ' PIP ',
                                    ' HMHP ',' HPETHNL ', ' MVKHP ',' MVKPC ',
                                    ' MCRENOL ', ' MCRHP ', ' MACR1OOH ',' MVKN ',
                                    ' MCRHNB ',' RIPA ', ' RIPB ', ' RIPC ',
                                    ' RIPD ',' HPALD1 ',' HPALD2 ',' HPALD3 ',
                                    ' HPALD4 ',' IHN1 ',' IHN4 ',' INPB ',
                                    ' INPD ', ' ICN ', ' ICPDH ',' IDHDP ',
                                    ' IDHPE ', ' IDCHP ',' ITHN ',' ITCN ',
                                    ' ETHP ', ' BZCO3H '],
                      familynametag+'FromO1DH2O':[' H2O '],
                      familynametag+'FromNOHO2':[' NO '],
                      familynametag+'FromO3HO2':[' O3 '],
                      familynametag+'FromH2O2hv':[' H2O2 '],
                      familynametag+'FromHOBrhv':[' HOBr '],familynametag+'FromHBrwO':[' O ',' O1D '],
                      familynametag+'FromHOClhv':[' HOCl '],familynametag+'FromHOIhv':[' HOI ']}
  
  
      dict_exlreactant1={familynametag+'FromH2O2hv':[' O '],familynametag+'FromO3HO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],familynametag+'FromO1DH2O':[' None '],
                      familynametag+'FromOVOCsROOHhv':[' OH '],
                      familynametag+'FromHOBrhv':[' None '],familynametag+'FromHBrwO':[' None '],
                      familynametag+'FromHOClhv':[' None '],familynametag+'FromHOIhv':[' None ']}
  
      dict_exlreactant2={familynametag+'FromH2O2hv':[' None '],familynametag+'FromO3HO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],familynametag+'FromO1DH2O':[' None '],
                      familynametag+'FromOVOCsROOHhv':[' None '],
                      familynametag+'FromHOBrhv':[' None '],familynametag+'FromHBrwO':[' None '],
                      familynametag+'FromHOClhv':[' None '],familynametag+'FromHOIhv':[' None ']}
  
      dict_exlreactant3={familynametag+'FromH2O2hv':[' None '],familynametag+'FromO3HO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],familynametag+'FromO1DH2O':[' None '],
                      familynametag+'FromOVOCsROOHhv':[' None '],
                      familynametag+'FromHOBrhv':[' None '],familynametag+'FromHBrwO':[' None '],
                      familynametag+'FromHOClhv':[' None '],familynametag+'FromHOIhv':[' None ']}
  

      dict_exlproduct1={familynametag+'FromH2O2hv':[' None '],familynametag+'FromO3HO2':[' None '],
                      familynametag+'FromNOHO2':[' None '],familynametag+'FromO1DH2O':[' None '],
                      familynametag+'FromOVOCsROOHhv':[' None '],
                      familynametag+'FromHOBrhv':[' None '],familynametag+'FromHBrwO':[' None '],
                      familynametag+'FromHOClhv':[' None '],familynametag+'FromHOIhv':[' None ']}
  

      dict_product1={familynametag+'FromH2O2hv':[' PHOx '],familynametag+'FromO3HO2':[' PHOx '],
                      familynametag+'FromNOHO2':[' PHOx '],familynametag+'FromO1DH2O':[' PHOx '],
                      familynametag+'FromOVOCsROOHhv':[' PHOx '],
                      familynametag+'FromHOBrhv':[' PHOx '],familynametag+'FromHBrwO':[' PHOx '],
                      familynametag+'FromHOClhv':[' PHOx '],familynametag+'FromHOIhv':[' PHOx ']}
  

  
      dict_product3={familynametag+'FromH2O2hv':[' PHOx '],familynametag+'FromO3HO2':[' PHOx '],
                      familynametag+'FromNOHO2':[' PHOx '],familynametag+'FromO1DH2O':[' PHOx '],
                      familynametag+'FromOVOCsROOHhv':[' PHOx '],
                      familynametag+'FromHOBrhv':[' PHOx '],familynametag+'FromHBrwO':[' PHOx '],
                      familynametag+'FromHOClhv':[' PHOx '],familynametag+'FromHOIhv':[' PHOx ']}
  

  
      dict_product2={familynametag+'FromH2O2hv':[' OH '],familynametag+'FromO3HO2':[' OH '],
                      familynametag+'FromNOHO2':[' OH '],familynametag+'FromO1DH2O':[' OH '],
                      familynametag+'FromOVOCsROOHhv':[' OH '],
                      familynametag+'FromHOBrhv':[' OH '],familynametag+'FromHBrwO':[' OH '],
                      familynametag+'FromHOClhv':[' OH '],familynametag+'FromHOIhv':[' OH ']}
  

  
      dict_altproduct2={familynametag+'FromH2O2hv':[' OH '],familynametag+'FromO3HO2':[' OH '],
                      familynametag+'FromNOHO2':[' OH '],familynametag+'FromO1DH2O':[' OH '],
                      familynametag+'FromOVOCsROOHhv':[' etc. '],
                      familynametag+'FromHOBrhv':[' OH '],familynametag+'FromHBrwO':[' OH '],
                      familynametag+'FromHOClhv':[' OH '],familynametag+'FromHOIhv':[' OH ']}
  
  
      dict_othaltproduct2={familynametag+'FromH2O2hv':[' OH '],familynametag+'FromO3HO2':[' OH '],
                      familynametag+'FromNOHO2':[' OH '],familynametag+'FromO1DH2O':[' OH '],
                      familynametag+'FromOVOCsROOHhv':[' OH '],
                      familynametag+'FromHOBrhv':[' OH '],familynametag+'FromHBrwO':[' OH '],
                      familynametag+'FromHOClhv':[' OH '],familynametag+'FromHOIhv':[' OH ']}
  

    elif (familynametag in 'THOxL'):

      dict_reactant1={familynametag+'FromOHHOy':[' OH '],familynametag+'FromOHNOy':[' OH '],
                      familynametag+'FromOHCH4':[' OH '],familynametag+'FromOHCO':[' OH '],
                      familynametag+'FromOHC1VOC':[' OH '],familynametag+'FromOHC2VOC':[' OH '],
                      familynametag+'FromOHHBr':[' OH '],familynametag+'FromOHrBrx':[' OH '],
                      familynametag+'FromOHHI':[' OH '],familynametag+'FromOHrIx':[' OH '],
                      familynametag+'FromOHHCl':[' OH '],familynametag+'FromOHrClx':[' OH '],
                      familynametag+'FromOHSSClm':[' OH ']}

      dict_reactant2={familynametag+'FromOHHOy':[' H2 ',' O3 ',' H2O2 ',' 2 OH ', ' HO2 ', ' MO2 ', ' NO3 ', ' O '],
                      familynametag+'FromOHNOy':[' NO ',' NO2 ',' HNO2 ',' HNO3 ',' HNO4 ', ' R4N2 ', ' ETHLN ',
                                                 ' PROPNN ', ' MONITS ', ' MONITU ', ' HONIT ', ' MENO3 ', ' ETNO3 ',
                                                 ' IPRNO3 ', ' NPRNO3 ', ' IHN2 ', ' IHN3 ', ' IHN1 ', ' IHN4 ', ' INPB ',
                                                 ' INPD ', ' ICN ', ' IDN ', ' MVKN ', ' MCRHN ', ' MPAN ', ' ITCN ',
                                                 ' ITHN ', ' ETHN ', ' BZPAN ', ' NPHEN '],
                      familynametag+'FromOHCH4':[' CH4 '],
                      familynametag+'FromOHCO':[' CO '],
                      familynametag+'FromOHC1VOC':[' MP ', ' CH2O ', ' HCOOH ', ' MOH '],
                      familynametag+'FromOHC2VOC':[' ATOOH ',' ALD2 ', ' C2H6 ',' C3H8 ', ' ALK4 ',' ACTA ', ' RCHO ',
                                                   ' ACET ', ' MEK ',  ' EOH ', ' ROH ',  ' PRPE ', ' GLYC ', ' GLYX ',
                                                   ' MGLY ', ' HAC ',  ' PRPN ',' ETP ',  ' RA3P ', ' RB3P ', ' R4P ',
                                                   ' RP ',   ' PP ',   ' MAP ', ' MTPA ', ' MTPO ', ' LIMO ', ' PIP ',
                                                   ' ISOP ', ' IDC ',  ' RIPA ',' RIPB ', ' RIPC ', ' RIPD ', 
                                                   ' IEPOXD ', ' IEPOXA ', ' IEPOXB ',    ' MVK ',  ' MACR ', ' MCRENOL ',
                                                   ' MVKDH ', ' MVKHC ', ' MCRDH ', ' MACR1OOH ', ' HMML ', ' ICPDH ',
                                                   ' IDCHP ', ' PYAC ', ' HMHP ', ' C2H4 ', ' C2H2 ', ' ETHP ', ' BENZ ',
                                                   ' TOLU ', ' XYLE ', ' PHEN ', ' CSL ', ' MCT ', ' BALD ', ' BZCO3H ',
                                                   ' BENZP ',' AROMP4 ', ' AROMP5 ' ],
                      familynametag+'FromOHHBr':[' HBr '],familynametag+'FromOHrBrx':[' Br2 ', ' BrO '],
                      familynametag+'FromOHHI':[' HI '],familynametag+'FromOHrIx':[' I2 ', ' HOI '],
                      familynametag+'FromOHHCl':[' HCl '],familynametag+'FromOHrClx':[' Cl2 ', ' ClO ', ' OClO ', ' Cl2O2 ', ' HOCl ', ' ClNO2 ',' ClNO3 '],
                      familynametag+'FromOHSSClm':[' SALACL ', ' SALCCL ']}
                      
           
      dict_exlreactant1={familynametag+'FromOHHOy':[' None '],familynametag+'FromOHNOy':[' None '],
                      familynametag+'FromOHCH4':[' None '],familynametag+'FromOHCO':[' None '],
                      familynametag+'FromOHC1VOC':[' None '],familynametag+'FromOHC2VOC':[' None '],
                      familynametag+'FromOHHBr':[' None '],familynametag+'FromOHrBrx':[ ' None '],
                      familynametag+'FromOHHI':[' None '],familynametag+'FromOHrIx':[' None '],
                      familynametag+'FromOHHCl':[' None '],familynametag+'FromOHrClx':[ ' None '],
                      familynametag+'FromOHSSClm':[' None ']}

      dict_exlreactant2={familynametag+'FromOHHOy':[' None '],familynametag+'FromOHNOy':[' None '],
                      familynametag+'FromOHCH4':[' None '],familynametag+'FromOHCO':[' None '],
                      familynametag+'FromOHC1VOC':[' None '],familynametag+'FromOHC2VOC':[' None '],
                      familynametag+'FromOHHBr':[' None '],familynametag+'FromOHrBrx':[ ' None '],
                      familynametag+'FromOHHI':[' None '],familynametag+'FromOHrIx':[' None '],
                      familynametag+'FromOHHCl':[' None '],familynametag+'FromOHrClx':[ ' None '],
                      familynametag+'FromOHSSClm':[' None ']}

      dict_exlreactant3={familynametag+'FromOHHOy':[' None '],familynametag+'FromOHNOy':[' None '],
                      familynametag+'FromOHCH4':[' None '],familynametag+'FromOHCO':[' None '],
                      familynametag+'FromOHC1VOC':[' None '],familynametag+'FromOHC2VOC':[' None '],
                      familynametag+'FromOHHBr':[' None '],familynametag+'FromOHrBrx':[ ' None '],
                      familynametag+'FromOHHI':[' None '],familynametag+'FromOHrIx':[' None '],
                      familynametag+'FromOHHCl':[' None '],familynametag+'FromOHrClx':[ ' None '],
                      familynametag+'FromOHSSClm':[' None ']}

      dict_exlproduct1={familynametag+'FromOHHOy':[' None '],familynametag+'FromOHNOy':[' None '],
                      familynametag+'FromOHCH4':[' None '],familynametag+'FromOHCO':[' None '],
                      familynametag+'FromOHC1VOC':[' None '],familynametag+'FromOHC2VOC':[' None '],
                      familynametag+'FromOHHBr':[' None '],familynametag+'FromOHrBrx':[ ' None '],
                      familynametag+'FromOHHI':[' None '],familynametag+'FromOHrIx':[' None '],
                      familynametag+'FromOHHCl':[' None '],familynametag+'FromOHrClx':[ ' None '],
                      familynametag+'FromOHSSClm':[' None ']}


      dict_product1={familynametag+'FromOHHOy':[' LHOx '],familynametag+'FromOHNOy':[' LHOx '],
                      familynametag+'FromOHCH4':[' LHOx '],familynametag+'FromOHCO':[' LHOx '],
                      familynametag+'FromOHC1VOC':[' LHOx '],familynametag+'FromOHC2VOC':[' LHOx '],
                      familynametag+'FromOHHBr':[' LHOx '],familynametag+'FromOHrBrx':[ ' LHOx '],
                      familynametag+'FromOHHI':[' LHOx '],familynametag+'FromOHrIx':[' LHOx '],
                      familynametag+'FromOHHCl':[' LHOx '],familynametag+'FromOHrClx':[ ' LHOx '],
                      familynametag+'FromOHSSClm':[' LHOx ']}

      dict_product3={familynametag+'FromOHHOy':[' LHOx '],familynametag+'FromOHNOy':[' LHOx '],
                      familynametag+'FromOHCH4':[' LHOx '],familynametag+'FromOHCO':[' LHOx '],
                      familynametag+'FromOHC1VOC':[' LHOx '],familynametag+'FromOHC2VOC':[' LHOx '],
                      familynametag+'FromOHHBr':[' LHOx '],familynametag+'FromOHrBrx':[ ' LHOx '],
                      familynametag+'FromOHHI':[' LHOx '],familynametag+'FromOHrIx':[' LHOx '],
                      familynametag+'FromOHHCl':[' LHOx '],familynametag+'FromOHrClx':[ ' LHOx '],
                      familynametag+'FromOHSSClm':[' LHOx ']}

      dict_product2={familynametag+'FromOHHOy':[' LHOx '],familynametag+'FromOHNOy':[' LHOx '],
                      familynametag+'FromOHCH4':[' LHOx '],familynametag+'FromOHCO':[' LHOx '],
                      familynametag+'FromOHC1VOC':[' LHOx '],familynametag+'FromOHC2VOC':[' LHOx '],
                      familynametag+'FromOHHBr':[' LHOx '],familynametag+'FromOHrBrx':[ ' LHOx '],
                      familynametag+'FromOHHI':[' LHOx '],familynametag+'FromOHrIx':[' LHOx '],
                      familynametag+'FromOHHCl':[' LHOx '],familynametag+'FromOHrClx':[ ' LHOx '],
                      familynametag+'FromOHSSClm':[' LHOx ']}

  
      dict_altproduct2={familynametag+'FromOHHOy':[' LHOx '],familynametag+'FromOHNOy':[' LHOx '],
                      familynametag+'FromOHCH4':[' LHOx '],familynametag+'FromOHCO':[' LHOx '],
                      familynametag+'FromOHC1VOC':[' LHOx '],familynametag+'FromOHC2VOC':[' LHOx '],
                      familynametag+'FromOHHBr':[' LHOx '],familynametag+'FromOHrBrx':[ ' LHOx '],
                      familynametag+'FromOHHI':[' LHOx '],familynametag+'FromOHrIx':[' LHOx '],
                      familynametag+'FromOHHCl':[' LHOx '],familynametag+'FromOHrClx':[ ' LHOx '],
                      familynametag+'FromOHSSClm':[' LHOx ']}


      dict_othaltproduct2={familynametag+'FromOHHOy':[' LHOx '],familynametag+'FromOHNOy':[' LHOx '],
                      familynametag+'FromOHCH4':[' LHOx '],familynametag+'FromOHCO':[' LHOx '],
                      familynametag+'FromOHC1VOC':[' LHOx '],familynametag+'FromOHC2VOC':[' LHOx '],
                      familynametag+'FromOHHBr':[' LHOx '],familynametag+'FromOHrBrx':[ ' LHOx '],
                      familynametag+'FromOHHI':[' LHOx '],familynametag+'FromOHrIx':[' LHOx '],
                      familynametag+'FromOHHCl':[' LHOx '],familynametag+'FromOHrClx':[ ' LHOx '],
                      familynametag+'FromOHSSClm':[' LHOx ']}

  
  

    else:
        print('No database for the tag species yet :',familynametag)
        exit()
    return dict_reactant1,dict_reactant2,dict_exlreactant1,dict_exlreactant2,dict_exlreactant3,dict_exlproduct1,dict_product1,dict_product3,dict_product2,dict_altproduct2,dict_othaltproduct2


def tagReactions(familyname,familynametag,MolarWeight):
    """
    tagging reactions in the custom_tagged.eqn.
    
    Args:
    -----
    familyname (str) : Name of family species
    familynametag (str) : tag Name (loss or production) of family species
    MolarWeight (float) : molar weight of family species [unit: g]
    """


    # Filenames
    file_reactions = 'reactions_need_to_be_tagged_'+familynametag+'.csv'
    file_species   = 'gckpp_Parameters.F90'
    file_org   = 'custom.eqn'
    file_updated   = 'custom_tagged.eqn'

    kppfile_org   = 'custom.kpp'
    kppfile_updated   = 'custom_tagged.kpp'

    spcfile_org   = 'species_database.yml'
    spcfile_updated   = 'species_database_tagged.yml'
    restartfile   = 'GEOSChem.Restart.20190701_0000z.nc4'

    strlist    = ['Br2FromSS',   'BrClFromSS', 'Cl2FromSS',
                  'PBr2FromSS',  'PBrClFromSS','PCl2FromSS',
                  'Br2ToSS',   'BrClToSS', 'Cl2ToSS',
                  'PBr2ToSS',  'PBrClToSS','PCl2ToSS',
                  'HBrFromSS','HBrToSS',
                  'PHBrFromSS','PHBrToSS']


    bihalogenmw_g= [159.80,      115.45,       70.90,
                    159.80,      115.45,       70.90,
                    159.80,      115.45,       70.90,
                    159.80,      115.45,       70.90,
                    80.91,       80.91,
                    80.91,       80.91]


    bihalogen=dict(zip(strlist,bihalogenmw_g))

    # Rread file with reactions
    with open(file_reactions, 'r') as f:
        lines = f.readlines()

    # Read files with species
    with open(file_species, 'r') as f:
        params = f.readlines()


    
    nn   = 0
    ntot = 0

    # get all lines with a reaction
    #rxt = [i for i in lines if '-->' in i]


    rctcefarr=np.zeros((2000))
    equationstrs = ["" for x in range(2000)]
    eqntags = ["" for x in range(2000)]
    #print(rctcefarr)
    #exit()
    # walk through all reactions
    irct = 0
    #----------------
    # read reactions which leads to Ox loss
    for irxt in lines:
        #---------------------------------------------------
        #spl = irxt.split('index')
        spl = re.split('index| LUMPSPC ',irxt)
        print(spl[0])
        print(spl[1])
        print(spl[2])
        rn = int(spl[1])
        #print(rn)
        # get reaction number and reactant if OH is on left-hand side of reaction
        #if ' LOx ' in spl[0]:
        if ' '+familyname.strip()+' ' in spl[0]:
            substr=re.split('-->|\+|-',spl[0])
            #print(substr)
            #substr_LOx = [isubstr for isubstr in substr if 'LOx' in isubstr]
            substr_LOx = [isubstr for isubstr in substr if familyname.strip() in isubstr]
            #print(substr_LOx)
            aa=substr_LOx[0].split(' ')
            aa[:] = [x for x in aa if x]
            #print(aa)
            rctcef='1.00'
            #if (aa[0] not in 'LOx'):
            if (aa[0] not in familyname.strip()):
               rctcef=aa[0]
            #print(rctcef)
            rctcefarr[rn-1]=rctcef
            equationstrs[rn-1]=spl[0].strip()
            eqntags[rn-1]=spl[2].strip()
    #--------------------
    # Save # that line to write out later.
    # update .eqn file
    with open(file_org, 'r') as f:
        utillines = f.readlines()
    nequation=0
    inequation=0
    with open(file_updated, 'w') as f:
        for line in utillines:
            if (inequation == 1) and (not line.startswith('//') and (' : ' in line)):
                nequation = nequation + 1
                sublines=line.split(' : ')
                #print(sublines)
                #exit()
                if (rctcefarr[nequation-1] != 0.0):
                   #f.write(sublines[0]+'+ '+"{:05.3f}".format(rctcefarr[nequation-1])+'TOxL'+"{:04d}".format(nequation)+' : '+sublines[1].rstrip()+'\n')
                   #f.write(sublines[0]+'+ '+"{:05.3f}".format(rctcefarr[nequation-1])+familynametag.strip()+"{:04d}".format(nequation)+' : '+sublines[1].rstrip()+'\n')
                   f.write(sublines[0]+' + '+"{:05.3f}".format(rctcefarr[nequation-1])+eqntags[nequation-1]+' : '+sublines[1].rstrip()+'\n')
                   #print(sublines[0]+'+ '+"{:05.3f}".format(rctcefarr[nequation-1])+'TOxL'+"{:04d}".format(nequation)+' : '+sublines[1].rstrip()+'\n')
                   #print(sublines[0]+'+ '+"{:05.3f}".format(rctcefarr[nequation-1])+familynametag.strip()+"{:04d}".format(nequation)+' : '+sublines[1].rstrip()+'\n')
                   print(sublines[0]+' + '+"{:05.3f}".format(rctcefarr[nequation-1])+eqntags[nequation-1]+' : '+sublines[1].rstrip()+'\n')
                else:
                   f.write(line.rstrip()+'\n')

            elif '#EQUATIONS' in line:
                inequation = 1
                f.write(line.rstrip()+'\n')

            elif '#DEFVAR' in line:
                f.write(line.rstrip()+'\n')

                for ieqn in list(filter(None, list(set(eqntags)))):
                       f.write(ieqn.strip()+'   = IGNORE; {Dummy species to track reactions leading to '+familyname.strip()+ '}'+'\n')
            else:
                f.write(line.rstrip()+'\n')


    print('Written to '+file_updated)
    f.close()

    #-----------
    # update .kpp file
    with open(kppfile_org, 'r') as f:
        utillines = f.readlines()
    with open(kppfile_updated, 'w') as f:
        for line in utillines:

            if '#FAMILIES' in line:
                f.write(line.rstrip()+'\n')

                for ieqn in list(filter(None, list(set(eqntags)))):
                       f.write('P'+ieqn.strip()+' : '+ieqn.strip()+';'+'\n')
            else:
                f.write(line.rstrip()+'\n')

    print('Written to '+kppfile_updated)
    f.close()

    #-----------
    # update species_database file
    with open(spcfile_org, 'r') as f:
        utillines = f.readlines()
    with open(spcfile_updated, 'w') as f:
        for line in utillines:
            if '# NOTE:' in line:
                f.write(line.rstrip()+'\n')

                for ieqn in list(filter(None, list(set(eqntags)))):
                       f.write(ieqn.strip()+':'+'\n')
                       f.write('  FullName: Dummy species to track '+familyname.strip()+' via reaction '+ieqn.strip()+'\n')
                       f.write('  Is_Gas: true'+'\n')
                       f.write('  MW_g: '+"{:05.2f}".format(MolarWeight)+'\n')

                       f.write('P'+ieqn.strip()+':'+'\n')
                       f.write('  FullName: Dummy species to track '+familyname.strip()+' via reaction '+ieqn.strip()+'\n')
                       f.write('  Is_Gas: true'+'\n')
                       f.write('  MW_g: '+"{:05.2f}".format(MolarWeight)+'\n')

                if familyname.strip() in ['LHOx','PHOx']:
                       f.write(familyname.strip()+':'+'\n')
                       f.write('  FullName: Dummy species to track '+familyname.strip()+' via reaction '+familyname.strip()+'\n')
                       f.write('  Is_Gas: true'+'\n')
                       f.write('  MW_g: '+"{:05.2f}".format(MolarWeight)+'\n')

                if familyname.strip() in ['PHOx']:
                   for key,val in bihalogen.items():
                       f.write(key+':'+'\n')
                       f.write('  FullName: Dummy species to track '+key+' via reaction '+key+'\n')
                       f.write('  Is_Gas: true'+'\n')
                       f.write('  MW_g: '+"{:05.2f}".format(val)+'\n')


            else:
                f.write(line.rstrip()+'\n')

    print('Written to '+spcfile_updated)
    f.close()
   
    #--------------
    #add new species to restart file


    for ieqn in list(filter(None, list(set(eqntags)))):
          newspc=ieqn
          os.system('cdo selvar,SpeciesRst_HOBrSSC '+restartfile+' '+newspc+'.nc4')
          os.system('ncrename -h -v SpeciesRst_HOBrSSC,SpeciesRst_'+newspc+' '+newspc+'.nc4')
          os.system('ncks -h -A -M '+newspc+'.nc4 '+restartfile)
          os.system('rm '+newspc+'.nc4')

          newspc='P'+ieqn
          os.system('cdo selvar,SpeciesRst_HOBrSSC '+restartfile+' '+newspc+'.nc4')
          os.system('ncrename -h -v SpeciesRst_HOBrSSC,SpeciesRst_'+newspc+' '+newspc+'.nc4')
          os.system('ncks -h -A -M '+newspc+'.nc4 '+restartfile)
          os.system('rm '+newspc+'.nc4')

          #exit()
    
    if familyname.strip() in ['LHOx','PHOx']:

          newspc=familyname.strip()
          os.system('cdo selvar,SpeciesRst_HOBrSSC '+restartfile+' '+newspc+'.nc4')
          os.system('ncrename -h -v SpeciesRst_HOBrSSC,SpeciesRst_'+newspc+' '+newspc+'.nc4')
          os.system('ncks -h -A -M '+newspc+'.nc4 '+restartfile)
          os.system('rm '+newspc+'.nc4')

    if familyname.strip() in ['PHOx']:
  
        for key,val in bihalogen.items():
            
            os.system('cdo selvar,SpeciesRst_HOBrSSC '+restartfile+' '+key+'.nc4')
            os.system('ncrename -h -v SpeciesRst_HOBrSSC,SpeciesRst_'+key+' '+key+'.nc4')
            os.system('ncks -h -A -M '+key+'.nc4 '+restartfile)
            os.system('rm '+key+'.nc4')



    os.system('cp custom_tagged.kpp custom.kpp')
    os.system('cp custom_tagged.eqn custom.eqn')
    os.system('cp species_database_tagged.yml species_database.yml')

def main():
    """
    Main program.  Gets number of arguments and calls writeReactivity
    """

    # Get the name of the mechansim
    mechName = 'hsc'
    #mechName = 'standard'

    os.system('rm custom_tagged.kpp custom_tagged.eqn species_database_tagged.yml')
    os.system('rm GEOSChem.Restart.20190701_0000z.nc4')
    os.system('cp custom_org.kpp custom.kpp')
    os.system('cp custom_org.eqn custom.eqn')
    os.system('cp species_database_org.yml species_database.yml')
    os.system('cp /users/mjr583/scratch/GC/14.1.0/rundirs/new_base_4x5_tracers/Restarts/GEOSChem.Restart.20190701_0000z.nc4 .')

    familynames=['POx','LOx','LHOx','PHOx']
    familynametags=['TOxP','TOxL','THOxL','THOxP']
    MWg=[48.00,48.00,17.00,17.00]

    #familynames=['LHOx','PHOx']
    #familynametags=['THOxL','THOxP']
    #MWg=[17.00,17.00]

    #familynames=['PHOx']
    #familynametags=['THOxP']
    #MWg=[17.00]


    for ifamily in range(0,len(familynames)):
        dict_reactant1,dict_reactant2,dict_exlreactant1,dict_exlreactant2,dict_exlreactant3,dict_exlproduct1,dict_product1,dict_product3,dict_product2,dict_altproduct2,dict_othaltproduct2=getPathways(familynametags[ifamily],mechName)
        extractReactions([familynames[ifamily]],[familynametags[ifamily]],
                         dict_reactant1,dict_reactant2,dict_product1,dict_product2,
                         dict_altproduct2,dict_othaltproduct2, dict_product3,
                         dict_exlproduct1,dict_exlreactant1,dict_exlreactant2,dict_exlreactant3)
    #exit()

    for ifamily in range(0,len(familynames)):
        tagReactions(familynames[ifamily],familynametags[ifamily],MWg[ifamily])
    

if __name__ == '__main__':
    main()
