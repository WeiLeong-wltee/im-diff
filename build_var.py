import os,glob,sys
import warnings

from unfold_jwst import msgs, utils

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Run imaging difference pipeline'
                                     'Usage: python build_var.py im1 im2 --psf1 psf1 --psf2 psf2 --outim_name im_diff.fits',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required args
    parser.add_argument('im1',type=str,help='First image')
    parser.add_argument('im2',type=str,help='Second image')
    # additional args
    parser.add_argument('--psf1',type=str,default=None,help='First image PSF. Provide filtername if using webbpsf.')
    parser.add_argument('--psf2',type=str,default=None,help='Second image PSF. Provide filtername if using webbpsf.')
    parser.add_argument('--outim_name',type=str,default='im_diff.fits',help='Output name for difference image')

    # optional args
    parser.add_argument('-h','--hotpants',default=False,action='store_true',help='im-diff with hotpants')
    parser.add_argument('-z','--pyzogy',default=True,action='store_true',help='im-diff with pyzogy')
    parser.add_argument('--mask1',type=str,default=None,help='First image mask')
    parser.add_argument('--mask2',type=str,default=None,help='Second image mask')
    parser.add_argument('--n_stamps',type=int,default=1,help='Number of stamps to use while fitting background')
    parser.add_argument('--sigma_cut',type=float,default=5,help='Threshold (in standard deviations) to extract a star from the image ("thresh" in "sep.extract")')
    parser.add_argument('--max_iterations',type=int,default=5,help='Maximum number of iterations to reconvolve the images for gain matching')
    parser.add_argument('--use_pixels',default=False,action='store_true',help='Use pixels instead of using sep extraction.')
    parser.add_argument('--percent',type=float,default=99,help='Remove pixels less than "percent" percentile above sky level to speed fitting.')
    parser.add_argument('--pixstack_limit',type=int,default=None,help='Number of active object pixels in sep, set with sep.set_extract_pixstack')
    parser.add_argument('--sep_deblend_nthresh',type=int,default=32,help='sep deblending threshold for neighbour connecting pixels, mainly affect the faint source extraction')
    parser.add_argument('--fit_noise_mode',type=str,default='iqr',help='Background noise std, can be quantile "iqr" or median "mad" or sep "sep"')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true', help='Overwrite existing files?')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):

    # PSF FWHM for nircam filters. Inherited from unfold_jwst/scripts/nircam_image.py
    filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M',
               'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
               'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
    psf_fwhm = [0.987, 1.103, 1.298, 1.553, 1.628, 1.770, 1.801, 1.494, 1.990,
                2.060, 2.141, 2.304, 2.341, 1.340, 1.444, 1.585, 1.547, 1.711,
                1.760, 1.830, 1.901, 2.165, 2.179, 2.300, 2.302, 2.459, 2.507, 2.535, 2.574]
    dict_utils = {filters[i]: {'FWHM': psf_fwhm[i]} for i in range(len(filters))}

    msgs.info(f'Running imaging difference pipeline for {args.im1.split("/")[-1]} and {args.im2.split("/")[-1]}.')
    ## Check running mode
    mode = 'hotpants' if args.hotpants else 'PyZOGY'
    msgs.info(f'Running with {mode}.')

    outname = args.outim_name.split('.fits')[0] + '.fits'
    if os.path.exists(outname) and not args.overwrite:
        msgs.info(f'{outname} exists. Skip image subtraction.')
    else:
        if mode == 'hotpants':
            msgs.info('Do not support hotpants yet. Terminated.')
            pass
            # import subprocess
            # inim = args.im1; tmplim = args.im2; outim = outname
            # command = f'hotpants -inim {inim} -tmplim {tmplim} -outim {outim}'
            # subprocess.run(command,shell=True)
        else:
            if glob.glob(args.psf1) and glob.glob(args.psf2):
                pass 
            else:
                import webbpsf
                ## Only include NIRCam now
                nc = webbpsf.NIRCam();
                nc.filter =  args.psf1; res1 = nc.calc_psf(oversample=4,fov_pixels=101);args.psf1 = res1['DET_DIST']
                nc.filter =  args.psf2; res2 = nc.calc_psf(oversample=4,fov_pixels=101);args.psf2 = res2['DET_DIST']
                res1, res2, nc = None,None,None

            from PyZOGY.subtract import run_subtraction
            run_subtraction(args.im1, args.im2, args.psf1, args.psf2,output=outname,
            science_mask=args.mask1,reference_mask=args.mask2,n_stamps=args.n_stamps,normalization="reference",
            corrected=False,matched_filter=False,sigma_cut=args.sigma_cut,
            max_iterations=args.max_iterations,pixstack_limit=args.pixstack_limit,sep_deblend_nthresh=args.sep_deblend_nthresh,
            fit_noise_mode=args.fit_noise_mode)
    
if __name__ == '__main__':
    main(parse_args())