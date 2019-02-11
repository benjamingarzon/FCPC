#!/home/benjamin.garzon/Software/anaconda2/bin/python
import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import imsave
from PIL import Image, ImageFont, ImageDraw
import os

#widthBW = 3430 # 8.7cm x 600 dpi
width = 3425 # 8.7cm x 1000 dpi
width_h = 7007 # 17.8cm x 1000 dpi
width0 = 3070

normalfont = ImageFont.truetype("/usr/share/fonts/liberation/LiberationSerif-Bold.ttf", size = 20)
smallfont = ImageFont.truetype("/usr/share/fonts/liberation/LiberationSerif-Bold.ttf", size = 30)

letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"]
def combine_figs(fig_list, nrows, ncols, outputname, width = 0, letters = letters, with_letter = True, letter_color = ['black' for l in letters], fig_order = "col", letter_pos = (20, 10), font = normalfont):

    rows = []
    for row in range(nrows):
        cols = []
        for col in range(ncols):

            if fig_order == "col":
                myletter = letters[row*ncols + col]
                myfig = fig_list[row*ncols + col]
                myletter_color = letter_color[row + col*nrows]
            else:
                myletter = letters[row + col*nrows]
                myfig = fig_list[row + col*nrows]
                myletter_color = letter_color[row + col*nrows]

            newim = Image.open(myfig).convert('RGB') 
            draw = ImageDraw.Draw(newim) 
	    
            if with_letter:
                   
                 if myletter_color == 'white':
                     myfill=(255,255,255)
                 else:
                     myfill=(0,0,0)

                 draw.text(letter_pos, myletter, font = font, fill = myfill)
            del draw
            print(np.array(newim).shape)
            cols.append(np.array(newim)[:, :, :3])      

        rows.append(np.concatenate(tuple(cols), axis=1))

    im = np.concatenate(tuple(rows), axis=0)
    im = Image.fromarray(im)

    if width !=0:
        wpercent = (width/float(im.size[0]))
        height = int((float(im.size[1])*float(wpercent)))
        im = im.resize((width, height), Image.ANTIALIAS)
        print("Resizing image")
        print(im.size)

    imsave(outputname, im)

# definitions here

ext='.tif'
ext='.jpg'

outputname = ['/home/benjamin.garzon/Data/NSPN/module_stats/Figure%d%s'%(x+1, ext) for x in range(8) ]
suppname = ['/home/benjamin.garzon/Data/NSPN/module_stats/FigureS%d%s'%(x+1, ext) for x in range(2) ]


# decision acuity accuracy
fig_list_1 = [
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_decAc/Accuracy.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQmatrix/Accuracy.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQvocab/Accuracy.png' 
]

# modules
fig_list_2 = ['/home/benjamin.garzon/Data/NSPN/module_stats/modules/modules.sum/000%d.png'%(x) for x in range(7) ]

# connections
fig_list_3 = [
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_decAc/decAc_connections.png'
]

# degree
fig_list_4 = [
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_decAc/decAc_regions.pos.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_decAc/decAc_regions.neg.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQmatrix/IQmatrix_regions.pos.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQmatrix/IQmatrix_regions.neg.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQvocab/IQvocab_regions.pos.png',
  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_IQvocab/IQvocab_regions.neg.png'
]

#fig_list_S1 = [
#  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_IQmatrix/Accuracy.png',
#  '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_IQvocab/Accuracy.png'
#]

modules = [" SMT ", "IMOFC", "BGTMP", "DLPFC", "DMPFC", " VIS ", " PDMN"]

module_ordering = [2,3,4,1,6,0,5]
modules = [modules[i] for i in module_ordering]
fig_list_2 = [fig_list_2[i] for i in module_ordering]

#  
combine_figs(fig_list_1, 1, 3, outputname[0], with_letter = False, width = width)

combine_figs(fig_list_2, 1, 7, 
outputname[1], width = width_h, letters = modules, 
letter_color = ['white' for l in modules], letter_pos = (75, 1800), font = smallfont)

combine_figs(fig_list_3, 1, 1, outputname[2], with_letter = False, width = width_h)

combine_figs(fig_list_4, 1, 6, 
outputname[3], width = width_h, 
letter_color = 6*['white'], font = smallfont)

#combine_figs(fig_list_S1, 1, 2, suppname[0], width = width_h)

