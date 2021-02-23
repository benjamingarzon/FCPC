#!/usr/bin/env python

# a flexible version of FSL's overlay
import sys
from nilearn import plotting, datasets, image
import numpy as np
import matplotlib.pyplot as plt


if True:
    bg_img = sys.argv[1] 
    img = sys.argv[2]
    sl = int(sys.argv[3])
    outpng = sys.argv[4]
    lthresh = float(sys.argv[5])
    useatlas = sys.argv[6]
    annotate = sys.argv[7]
    alpha = sys.argv[8]
    alpha = 0.9 if alpha == '' else float(alpha)

else:
    import os
    os.chdir('/home/benjamin.garzon/Data/NSPN/module_stats/modules/modules_optimal_trans.sum')
    bg_img = "BG.nii.gz" 
    img = "vol0000.nii.gz"
    sl = int(100)
    outpng = "test"
    lthresh = float(4)
    useatlas = "empty.nii.gz"
    annotate = "0"
    alpha = 0.4
    
alpha = 0.3
    #bgrange = sys.argv[2]
#ulthresh = sys.argv[5]
img = image.load_img(img)
useatlas = '1'
img.get_fdata()[ img.get_fdata() < 0] = 0
(x, y, z) = image.coord_transform(0, 0, sl, img.affine)
print(sl, x, y, z)

smith = datasets.fetch_atlas_smith_2009()
if useatlas == '0':	
    atlas = image.load_img(smith.bm10)
    atlas = image.index_img(atlas, (7, 8, 9))

    print(useatlas)
    display = plotting.plot_prob_atlas(atlas, 
        bg_img=bg_img, 
        #colorbar=True, 
        black_bg=True, 
        display_mode='z', cut_coords = (z, ), # display_mode="z", cut_coords = 3, 
        annotate=annotate == '1', 
        view_type='filled_contours', draw_cross = False, alpha = 0.5)
    display.add_overlay(img, threshold = lthresh, alpha = alpha)


else:
    if useatlas == '1':	
        fig = plt.figure()
        display = plotting.plot_anat(bg_img, 
	    black_bg=True, 
	    display_mode='z', cut_coords = (z, ),
	    annotate= annotate == '1', draw_cross = False, figure = fig)
        display.add_contours(img, 
                             filled = True, 
                             alpha = alpha, 
                             cmap = 'nipy_spectral')

    else:

        fig = plt.figure()
        display = plotting.plot_stat_map( img, 
	    bg_img=bg_img, 
	    colorbar=False,
	    black_bg=True,
	    display_mode='z', 
        cut_coords = (z, ), 
	    annotate = annotate == '1', 
        draw_cross = False,  
        threshold = lthresh, 
        figure = fig)
        #display.add_contours(atlas, levels=[.5], colors='b', alpha = alpha, 
	    #filled = True)
#    display.add_edges(atlas)
#    display.add_overlay(atlas)
#        display.add_overlay(img, cmap= plotting.cm.black_red, figure = fig)
#        plotting.show()


#display = plotting.plot_stat_map(img, bg_img = bg_img, display_mode="ortho", 
#                  colorbar=False, 
#                  threshold = lthresh, annotate = False, black_bg = True,
#                  draw_cross = False)
#display.add_overlay(atlas,
#                    cmap=plotting.cm.black_blue)

display.savefig(outpng)
display.close()

from imageio import imwrite
from PIL import Image, ImageFont, ImageDraw

im = Image.open(outpng + '.png')
h = im.height
w = im.width
rh = 0.1
rw = 0.15
im = im.crop((int(w*rw), int(h*rh), int(w*(1-rw)), h))
imwrite(outpng + '.png', im)

#cut_coords= (0, 0, sl),

#plotting.show()
#plot_prob_atlas


