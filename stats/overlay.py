#!/home/benjamin.garzon/Software/anaconda/bin/python

# a flexible version of FSL's overlay
import sys
from nilearn import plotting, datasets, image
import numpy as np

bg_img = sys.argv[1] 
img = sys.argv[2]
sl = int(sys.argv[3])
outpng = sys.argv[4]
lthresh = float(sys.argv[5])
atlas = sys.argv[6]
annotate = sys.argv[7]
alpha = sys.argv[8]
alpha = 0.9 if alpha == '' else float(alpha)
#bgrange = sys.argv[2]
#ulthresh = sys.argv[5]
img = image.load_img(img)

img.get_data()[ img.get_data() < 0] = 0
(x, y, z) = image.coord_transform(0, 0, sl, img.affine)

smith = datasets.fetch_atlas_smith_2009()
if atlas == '0':	
    atlas = image.load_img(smith.bm10)
    atlas = image.index_img(atlas, (7, 8, 9))

    display = plotting.plot_prob_atlas(atlas, 
        bg_img=bg_img, 
        #colorbar=True, 
        black_bg=True, 
        display_mode='z', cut_coords = (z, ), # display_mode="z", cut_coords = 3, 
        annotate=annotate == '1', 
        view_type='filled_contours', draw_cross = False, alpha = 0.5)
    display.add_overlay(img, threshold = lthresh, alpha = alpha)


else:
    if atlas == '1':	
        display = plotting.plot_anat(bg_img, 
	    black_bg=True, 
	    display_mode='z', cut_coords = (z, ),
	    annotate= annotate == '1', draw_cross = False)
        display.add_overlay(img, threshold = lthresh, alpha = alpha, cmap = 'hot')

    else:
        display = plotting.plot_stat_map(img, 
	    bg_img=bg_img, 
	    colorbar=False, 
	    black_bg=True, 
	    display_mode='z', cut_coords = (z, ),
	    annotate= annotate == '1', draw_cross = False,  threshold = lthresh)
        display.add_contours(atlas, levels=[.5], colors='b', alpha = alpha, 
	    filled = True)
#    display.add_edges(atlas)
#    display.add_overlay(atlas)
#    display.add_overlay(img)



#display = plotting.plot_stat_map(img, bg_img = bg_img, display_mode="ortho", 
#                  colorbar=False, 
#                  threshold = lthresh, annotate = False, black_bg = True,
#                  draw_cross = False)
#display.add_overlay(atlas,
#                    cmap=plotting.cm.black_blue)

display.savefig(outpng)
display.close()
#cut_coords= (0, 0, sl),

#plotting.show()
#plot_prob_atlas


