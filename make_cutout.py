import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_adaptive
import copy 

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Run make_cutout for overlapping images.'
                                     'Usage: python make_cutout.py im1 im2 --outdir --outsuffix',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required args
    parser.add_argument('im1',type=str,help='First image')
    parser.add_argument('im2',type=str,help='Second image')
    

    # optional args
    parser.add_argument('--outdir',type=str,default='./',help='output directory')
    parser.add_argument('--outsuffix',type=str,default='overlapped',help='suffix added behind')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true', help='Overwrite existing files?')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):
	# ## Try resample to common wcs
	# im1 = fits.open('sapphires/Obs41_F200W_sci.fits')
	# im2 = fits.open('~/var_agn/jades_goodss/goodss/goodss_F200W_sci_bksub.fits')

	# im1_wcs = WCS(im1[0].header)
	# im2_wcs = WCS(im2[0].header)

	# im_hdus = [im[0] for im in [im1,im2]]
	# wcs_out, shape_out = find_optimal_celestial_wcs(im_hdus,resolution=0.03*u.arcsec)

	# array1, footprint1 = reproject_adaptive(im1[0],wcs_out, shape_out=shape_out,
	# 									conserve_flux=True)
	# new_header = copy.deepcopy(im1[0].header)
	# temp_header = wcs_out.to_header()
	# for i in temp_header:
	# 	new_header[i] = temp_header[i]
	# new_header['PC1_1'] = 1
	# new_header['PC1_2'] = 0
	# new_header['PC2_1'] = 0
	# new_header['PC2_2'] = 1

	# hdulist = fits.HDUList([fits.PrimaryHDU(data=array1,header=new_header)])
	# hdulist.writeto('Obs41_F200W_sci_reproject.fits',overwrite=True)

	# array2, footprint2 = reproject_adaptive(im2[0],wcs_out, shape_out=shape_out,
	# 									conserve_flux=True)
	# new_header = copy.deepcopy(im2[0].header)
	# temp_header = wcs_out.to_header()
	# for i in temp_header:
	# 	new_header[i] = temp_header[i]
	# new_header['PC1_1'] = 1
	# new_header['PC1_2'] = 0
	# new_header['PC2_1'] = 0
	# new_header['PC2_2'] = 1

	# hdulist = fits.HDUList([fits.PrimaryHDU(data=array2,header=new_header)])
	# hdulist.writeto('goodss_F200W_sci_bksub_reproject.fits',overwrite=True)

	## Try resample to smaller array
	im1 = fits.open(args.im1)
	im2 = fits.open(args.im2)

	im1_wcs = WCS(im1[0].header)
	im2_wcs = WCS(im2[0].header)

	im1_total_shape, im2_total_shape = np.count_nonzero(np.isfinite(im1[0].data) & (im1[0].data!=0.)), np.count_nonzero(np.isfinite(im2[0].data) & (im2[0].data!=0.))
	if im1_total_shape>im2_total_shape:
		wcs_out, shape_out = find_optimal_celestial_wcs(im2[0],resolution=0.03*u.arcsec)
		# wcs_out, shape_out = im2_wcs, (im2[0].data.shape[0],im2[0].data.shape[1])
	else:
		wcs_out, shape_out = find_optimal_celestial_wcs(im1[0],resolution=0.03*u.arcsec)
		# wcs_out, shape_out = im1_wcs, (im1[0].data.shape[0],im1[0].data.shape[1])

	array1, footprint1 = reproject_adaptive(im1[0],wcs_out, shape_out=shape_out,
										conserve_flux=True)

	array2, footprint2 = reproject_adaptive(im2[0],wcs_out, shape_out=shape_out,
										conserve_flux=True)

	# mask1 = np.isnan(array1)
	# mask2 = array2==0.
	# mask = mask1 | mask2
	# array1[mask] = np.nan
	# array2[mask] = np.nan 



	new_header = copy.deepcopy(im1[0].header)
	temp_header = wcs_out.to_header()
	for i in temp_header:
		new_header[i] = temp_header[i]
	new_header['PC1_1'] = 1
	new_header['PC1_2'] = 0
	new_header['PC2_1'] = 0
	new_header['PC2_2'] = 1

	hdulist = fits.HDUList([fits.PrimaryHDU(data=array1,header=new_header)])
	name1 = args.im1.split('/')[-1].split('.fits')[0]
	hdulist.writeto(f'{args.outdir}/{name1}_{args.outsuffix}.fits',overwrite=True)


	new_header = copy.deepcopy(im2[0].header)
	temp_header = wcs_out.to_header()
	for i in temp_header:
		new_header[i] = temp_header[i]
	new_header['PC1_1'] = 1
	new_header['PC1_2'] = 0
	new_header['PC2_1'] = 0
	new_header['PC2_2'] = 1

	hdulist = fits.HDUList([fits.PrimaryHDU(data=array2,header=new_header)])
	name2 = args.im2.split('/')[-1].split('.fits')[0]
	hdulist.writeto(f'{args.outdir}/{name2}_{args.outsuffix}.fits',overwrite=True)


if __name__ == '__main__':
    main(parse_args())
