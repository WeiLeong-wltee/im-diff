import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_adaptive
import copy 
from astropy.nddata.utils import Cutout2D


def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Run make_cutout for overlapping images.'
                                     'Usage: python make_cutout.py im1 im2',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required args
    parser.add_argument('im1',type=str,help='First image')
    parser.add_argument('im2',type=str,help='Second image')
    

    # optional args
    parser.add_argument('--px_scale',type=float,default=0.03,help='pixel scale')
    parser.add_argument('--xlo',type=int,default=None)
    parser.add_argument('--xup',type=int,default=None)
    parser.add_argument('--ylo',type=int,default=None)
    parser.add_argument('--yup',type=int,default=None)
    parser.add_argument('-mo','--mutual_overlapped',default=False,action='store_true',help='Consider overlapped region only')
    parser.add_argument('--outdir',type=str,default='./',help='output directory')
    parser.add_argument('--outsuffix',type=str,default='overlapped',help='suffix added behind')
    parser.add_argument('-nan','--nan',default=False,action='store_true',help='Use nan for empty region. Default is 0.')
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
		wcs_out, shape_out = find_optimal_celestial_wcs(im2[0],resolution=args.px_scale*u.arcsec)
	else:
		wcs_out, shape_out = find_optimal_celestial_wcs(im1[0],resolution=args.px_scale*u.arcsec)

	array1, footprint1 = reproject_adaptive(im1[0],wcs_out, shape_out=shape_out,
										conserve_flux=True)

	array2, footprint2 = reproject_adaptive(im2[0],wcs_out, shape_out=shape_out,
										conserve_flux=True)

	if args.mutual_overlapped:
		mask = np.isnan(array1) | np.isnan(array2) | (array1==0) | (array2==0)
		if args.nan:
			array1 = np.where(mask,np.nan,array1)
			array2 = np.where(mask,np.nan,array2)
		else:
			array1 = np.where(mask,0,array1)
			array2 = np.where(mask,0,array2)

	new_header = copy.deepcopy(im1[0].header)
	temp_header = wcs_out.to_header()
	for i in temp_header:
		new_header[i] = temp_header[i]
	new_header['PC1_1'] = 1
	new_header['PC1_2'] = 0
	new_header['PC2_1'] = 0
	new_header['PC2_2'] = 1
	new_header1 = new_header
	new_header = copy.deepcopy(im2[0].header)
	temp_header = wcs_out.to_header()
	for i in temp_header:
		new_header[i] = temp_header[i]
	new_header['PC1_1'] = 1
	new_header['PC1_2'] = 0
	new_header['PC2_1'] = 0
	new_header['PC2_2'] = 1
	new_header2 = new_header


	if isinstance(args.xlo,int) and isinstance(args.xup,int) and isinstance(args.ylo,int) and isinstance(args.yup,int):
		print ('Run cutout2D.')
		xcenter, ycenter = (args.xlo+args.xup)/2, (args.ylo+args.yup)/2
		position = (xcenter,ycenter)
		size = (int(np.ceil(args.xup-args.xlo)),int(np.ceil(args.yup-args.ylo)))
		temp1 = Cutout2D(array1,position,size,wcs=wcs_out)
		temp2 = Cutout2D(array2,position,size,wcs=wcs_out)
		array1 = temp1.data; array2= temp2.data; 
		new_header = copy.deepcopy(im1[0].header)
		temp_header = temp1.wcs.to_header()
		for i in temp_header:
			new_header[i] = temp_header[i]
		new_header['PC1_1'] = 1
		new_header['PC1_2'] = 0
		new_header['PC2_1'] = 0
		new_header['PC2_2'] = 1
		new_header1 = new_header

		new_header = copy.deepcopy(im2[0].header)
		temp_header = temp2.wcs.to_header()
		for i in temp_header:
			new_header[i] = temp_header[i]
		new_header['PC1_1'] = 1
		new_header['PC1_2'] = 0
		new_header['PC2_1'] = 0
		new_header['PC2_2'] = 1
		new_header2 = new_header

	hdulist = fits.HDUList([fits.PrimaryHDU(data=array1,header=new_header1)])
	name1 = args.im1.split('/')[-1].split('.fits')[0]
	hdulist.writeto(f'{args.outdir}/{name1}_{args.outsuffix}.fits',overwrite=True)

	hdulist = fits.HDUList([fits.PrimaryHDU(data=array2,header=new_header2)])
	name2 = args.im2.split('/')[-1].split('.fits')[0]
	hdulist.writeto(f'{args.outdir}/{name2}_{args.outsuffix}.fits',overwrite=True)


if __name__ == '__main__':
    main(parse_args())
