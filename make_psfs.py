import numpy as np
import glob
from photutils.aperture import CircularAnnulus,aperture_photometry
from photutils.centroids import centroid_1dg
from astropy.io import fits
from astropy.visualization import SqrtStretch,LogStretch,LinearStretch,MinMaxInterval, SqrtStretch,ImageNormalize

def profile_2D(psf):
    x_marg = [np.sum(psf[:, i]) for i in range(psf.shape[0])]
    y_marg = [np.sum(psf[j, :]) for j in range(psf.shape[1])]
    return np.array(x_marg), np.array(y_marg)

def bg_subtract(image, rms_clean_level):
    im_shape = image.shape
    m1, m2 = int((image.shape[0] / 2) - 5), int((image.shape[0] / 2) + 5)

    mask_array = image.copy()
    mask_array[m1:m2, m1:m2] = np.nan

    center = image.shape[0] / 2

    # measure the background in an annulus around star
    aper_annulus = CircularAnnulus((center, center), r_in=27, r_out=30)
    phot = aperture_photometry(image, aper_annulus)
    bkg_mean = phot['aperture_sum'][0] / aper_annulus.area

    # subtract background
    im_new = image - bkg_mean

    # if there are remaining positive peaks, set them to local rms
    if rms_clean_level != 0:
        mask_new = im_new.copy()
        mask_new[m1:m2, m1:m2] = np.nan
        mask_new_flat = mask_new.flatten()
        mask_new_flat_fin = mask_new_flat[np.isnan(mask_new_flat) == False]
        rms_new = np.sqrt(np.average(mask_new_flat_fin ** 2))
        mask_new_flat[mask_new_flat > (rms_clean_level * rms_new)] = rms_new

        new_cutout_bkgsub = np.reshape(mask_new_flat, im_shape)
        new_cutout_bkgsub[m1:m2, m1:m2] = im_new[m1:m2, m1:m2]

    else:
        new_cutout_bkgsub = im_new

    return new_cutout_bkgsub

data_dir = 'psfs_jades_ceers'
default_input_field = [f'ceers_nircam{i}' for i in range(1,11)] + ['jades']
filter_list = ['F090W','F115W','F125W','F150W','F160W','F200W','F356W','F444W']

for fil in filter_list:
    goodstarlist = [np.loadtxt(f'{data_dir}/psfs_{field}/{fil}_good_star_list.dat',dtype=int).flatten().tolist() for field in default_input_field
     if glob.glob(f'{data_dir}/psfs_{field}/{fil}_good_star_list.dat')]
    goodfieldlist = [field for field in default_input_field if glob.glob(f'{data_dir}/psfs_{field}/{fil}_good_star_list.dat')]
    temp = []
    [temp.append(f'{data_dir}/star_library_{field}/{fil}_{_id}.fits') for index,field in enumerate(goodfieldlist) for _id in goodstarlist[index]]
    psf1 = fits.open(temp[0])[0]
    cutout_size = psf1.data.shape[0]
    goodstararr = np.zeros((len(temp),cutout_size,cutout_size))
    for i in range(len(temp)):
        goodstararr[i] += fits.getdata(temp[i])
    sum_psf = np.nansum(goodstararr,axis=0)
    med_psf = np.nanmedian(goodstararr,axis=0)
    #re-centroid the final product one last time - don't want an offset for psf-matching
    x_offset_sum, y_offset_sum = centroid_1dg(sum_psf)
    xdiff_sum, ydiff_sum = np.abs(cutout_size / 2 - x_offset_sum), np.abs(cutout_size / 2 - y_offset_sum)

    new_sum_psf = np.roll(sum_psf, int(xdiff_sum), axis=0)
    new_sum_psf = np.roll(new_sum_psf, int(ydiff_sum), axis=1)

    x_offset_med, y_offset_med = centroid_1dg(med_psf)
    xdiff_med, ydiff_med = np.abs(cutout_size / 2 - x_offset_med), np.abs(cutout_size / 2 - y_offset_med)

    new_med_psf = np.roll(med_psf, int(round(xdiff_med)), axis=0)
    new_med_psf = np.roll(new_med_psf, int(round(ydiff_med)), axis=1)

    sum_psf_bg_norm = bg_subtract(new_sum_psf, rms_clean_level=0)
    sum_psf_bg_norm /= np.sum(sum_psf_bg_norm)

    med_psf_bg_norm = bg_subtract(new_med_psf, rms_clean_level=0)
    med_psf_bg_norm /= np.sum(med_psf_bg_norm)

    sum_psf = sum_psf_bg_norm
    med_psf = med_psf_bg_norm

    sum_psf_profile2d = profile_2D(sum_psf)
    med_psf_profile2d = profile_2D(med_psf)

    hdr = fits.Header()
    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = sum_psf.shape[0]
    hdr['NAXIS2'] = sum_psf.shape[1]
    hdr['PIXSCALE'] = psf1.header['PIXSCALE']
    new_hdu = fits.PrimaryHDU(sum_psf, header=hdr)
    new_hdulist = fits.HDUList([new_hdu])
    new_psf = new_hdulist.writeto(f'{data_dir}/{fil}_sum_psf.fits', overwrite=True)

    hdr = fits.Header()
    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = med_psf.shape[0]
    hdr['NAXIS2'] = med_psf.shape[1]
    hdr['PIXSCALE'] = psf1.header['PIXSCALE']
    new_hdu = fits.PrimaryHDU(med_psf, header=hdr)
    new_hdulist = fits.HDUList([new_hdu])
    new_psf = new_hdulist.writeto(f'{data_dir}/{fil}_med_psf.fits', overwrite=True)

