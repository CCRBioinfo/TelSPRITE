#!/usr/bin/env python
#
# As a work of the United States Government, this project is in the public 
# domain within the United States.
#
# Additionally, we waive copyright and related rights in the work worldwide 
# through the CC0 1.0 Universal public domain dedication. 
#
#
#
# usage: segment_all_channels.py [-h] -i INDIR -o OUTDIR [--min MIN] 
#                                 [--max MAX] [--dapi DAPI] [--debug] [-t] [-b]
#
# Segment nuclei in TIF files
#
# options:
#   -h, --help           show this help message and exit
#   -i, --indir INDIR    Directory to process, containing TIF images with 
#                        fields of nuclei
#   -o, --outdir OUTDIR  Directory for segmented files
#   --min MIN            minimum nucleus area in pixels (default 10000)
#   --max MAX            maximum nucleus area in pixels (default 100000)
#   --dapi DAPI          channel for dapi (zero indexed)
#   --debug              Set debug to TRUE
#   -t, --test           Test the current settings
#   -b, --blank          set the background of the segment to zero

__author__ = "Sarah Clatterbuck Soper"
__license__ = "Public Domain"

import numpy as np
import os
import sys
import glob
import argparse
import pathlib
from skimage import (
    segmentation,
    morphology,
    filters,
    measure,
    feature,
    io
)
from PIL import Image
from skimage.segmentation import clear_border, expand_labels
from scipy import ndimage as ndi
from scipy.stats import scoreatpercentile

#
# Parse command line arguments
#
parser = argparse.ArgumentParser(description="Segment nuclei in TIF files")
parser.add_argument(
    "-i", "--indir", 
    required=True, 
    type=pathlib.Path, 
    help="Directory to process, containing TIF images with fields of nuclei"
)
parser.add_argument(
    "-o",
    "--outdir",
    required=True,
    type=pathlib.Path,
    help="Directory for segmented files"
)
parser.add_argument(
    "--min",
    type=int,
    default="10000",
    help="minimum nucleus area in pixels (default 10000)"
)
parser.add_argument(
    "--max",
    type=int,
    default="100000",
    help="maximum nucleus area in pixels (default 100000)"
)
parser.add_argument(
    "--dapi", 
    type=int, 
    default="3", 
    help="channel for dapi (zero indexed)"
)
parser.add_argument(
    "--debug", 
    action='store_true',
    help="Set debug to TRUE"
)
parser.add_argument(
    "-t", "--test", 
    action='store_true',
    help="Test the current settings"
)
parser.add_argument(
    "-b", "--blank", 
    action='store_true',
    help="set the background of the segment to zero"
)
args = parser.parse_args()

#
# function reads an image, processes it, saves segmented nuclei to 
# new directory specified in outdir
#

def open_file(fqn):
    try:
        img1 = io.imread(fqn)
    except OSError as err:
        print("OS error:", err)
    except ValueError as err:
        print("Got a value error.", err)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise
    return(img1)

# takes multi-channel image,
# scales and filters,
# returns corrected image
def scale_img(raw):
    new = np.zeros(raw.shape, dtype="uint8")
    for x in range(raw.shape[2]):
        v = raw[...,x]
        top = scoreatpercentile(v.flatten(), 99)
        v[v>top] = top
        v = (v - v.min()) / (v.max() - v.min())
        v *= 255
        v = v.astype(np.uint8)
        new[:,:,x] = v
    return(new)

def threshold_img(raw, blue):
    dapi = raw[:, :, blue]
    
    # These images are not well illuminated across the field
    # Using local threshold with a block of about a third of the
    # image width seems to look good.

    if (round(dapi.shape[0]*0.3))%2 == 0:
        block = dapi.shape[0]*0.3 + 1
    else:
        block = round(dapi.shape[0]*0.3)
    
    regions = dapi > filters.threshold_local(dapi, block, 'mean')

    return(regions)

def segment_img(raw, binary):
    #
    # operating on binary threshold
    #
    clear = clear_border(binary)  # Remove edge touching grains
    clear = morphology.remove_small_objects(clear, max_size=10000)
    clear = morphology.remove_small_holes(clear, max_size=128)
    
    labels = measure.label(clear)
    
    kernel = morphology.disk(15)
    distance = ndi.distance_transform_edt(labels)
    image_max = ndi.maximum_filter(distance, size=15, mode='constant')
    coords = feature.peak_local_max(image_max, footprint=kernel, labels=labels)
    
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    
    ws = segmentation.watershed(-distance, markers, mask=labels)
    ws = expand_labels(ws, distance=10)
    return(ws)

def process_file(raw, labels, fn, outdir, blank):   
    props = measure.regionprops(labels, intensity_image=raw)

    for prop in props:
        print('Label: {} Area: {}'.format(prop.label, prop.area))
        if (prop.area > args.min) and (prop.area < args.max):
            print("Saving image " + str(prop.label) + 
                  " from " + fn + " to " + str(outdir))
            qfn = os.path.join(outdir, str(prop.label) + "_" + fn)
            minr, minc, maxr, maxc = prop.bbox
            if (blank==True):
                bl = np.zeros_like(raw)
                coords= prop.coords
                idx = tuple(coords.T)  # (array of dim0 indices, array of dim1 indices, array of dim2 indices)
                bl[idx] = raw[idx]
                roi = bl[minr:maxr, minc:maxc]
            else:
                roi = raw[minr:maxr, minc:maxc]
            io.imsave(qfn, roi, check_contrast=False)
#
# Main loop.  
# Iterate through files in indir, processing each one 
#
def main() -> int:

    for filename in glob.iglob("*.tif", root_dir=args.indir, recursive=False):
        fqn = os.path.join(args.indir, filename)
        
        if (args.debug==True):
            print(args)
        
        img = open_file(fqn)
        normalized = scale_img(img)
        mask = threshold_img(normalized, args.dapi)
        nuclei = segment_img(normalized, mask)
        
        if (args.test==True):
            I  = Image.fromarray(nuclei)
            qfn = os.path.join(args.outdir, filename)
            I.save(qfn)
            continue
            
        # adding error handling just in case!
        try:
            process_file(img, nuclei, filename, args.outdir, args.blank)
        except:
            pass
    return 0

if __name__ == '__main__':
    sys.exit(main())  

