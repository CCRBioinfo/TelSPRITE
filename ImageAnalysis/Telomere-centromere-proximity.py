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
# usage: Telomere-centromere-proximity.py [-h] -i INDIR -o OUTDIR [--min [N]]
#                                        [--max [N]] [--tolerance [N]]
#                                        [--one [N]] [--two [N]] [--dapi [N]]
#                                         [--debug [N]]
# 
# Quantify proximity of spots in two channels 
#
# options:
#   -h, --help           show this help message and exit
#   -i, --indir INDIR    Directory to process, containing TIF images of single nuclei
#   -o, --outdir OUTDIR  Directory for output files
#   --min [N]            minimum focus area in pixels (default 9)
#   --max [N]            maximum focus area in pixels (default 10000)
#   --tolerance [N]      Neighbor gap tolerance in pixels (default 5)
#   --one [N]            first channel for spots, typically Cen (zero indexed)
#   --two [N]            second channel for spots, typically Tel (zero indexed)
#   --dapi [N]           channel for DAPI
#   --debug [N]          Debug Setting, set to > 0 for debugging messages

__author__ = "Sarah Clatterbuck Soper"
__license__ = "Public Domain"

import numpy as np
import os
import sys
import cv2
import glob
import argparse
import pathlib
from skimage import io
from scipy.stats import percentileofscore, mode
from skimage.filters import threshold_multiotsu, threshold_otsu
from skimage.measure import regionprops, regionprops_table
from skimage.morphology import dilation
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
from statistics import mean, median
from alive_progress import alive_bar
from contextlib import contextmanager
import threading
import _thread
import bigfish.stack as stack

#
# Parse command line arguments
#
parser = argparse.ArgumentParser(description="Quantify proximity of spots in two channels")
parser.add_argument(
    "-i", 
    "--indir", 
    required=True,
    type=pathlib.Path,
    help="Directory to process, containing TIF images of single nuclei"
)
parser.add_argument(
    "-o",
    "--outdir",
    required=True,
    type=pathlib.Path,
    help="Directory for output files",
)
parser.add_argument(
    "--min",
    metavar="N",
    type=int,
    nargs="?",
    default="9",
    help="minimum focus area in pixels (default 9)",
)
parser.add_argument(
    "--max",
    metavar="N",
    type=int,
    nargs="?",
    default="10000",
    help="maximum focus area in pixels (default 10000)",
)
parser.add_argument(
    "--tolerance",
    metavar="N",
    type=int,
    nargs="?",
    default="5",
    help="Neighbor gap tolerance in pixels (default 5)",
)
parser.add_argument(
    "--one", metavar="N", type=int, 
    nargs="?", default="1", 
    help="first channel for spots, typically Cen (zero indexed)"
)
parser.add_argument(
    "--two", metavar="N", type=int, 
    nargs="?", default="0", 
    help="second channel for spots, typically Tel (zero indexed)"
)
parser.add_argument(
    "--dapi", metavar="N", type=int, 
    nargs="?", default="2", help="channel for DAPI"
)
parser.add_argument(
    "--debug", metavar="N", type=int, 
    nargs="?", default="0", help="Debug Setting, set to > 0 for debugging messages"
)
args = parser.parse_args()


class TimeoutException(Exception):
    def __init__(self, msg=''):
        self.msg = msg

@contextmanager
def time_limit(seconds, msg=''):
    timer = threading.Timer(seconds, lambda: _thread.interrupt_main())
    timer.start()
    try:
        yield
    except KeyboardInterrupt:
        raise TimeoutException(f'Timed out processing {msg}')
    finally:
        # if the action ends in specified time, timer is canceled
        timer.cancel()

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def spot_channel_stats(ch_one, ch_two, debug):
    real_overlap = np.logical_and(ch_one, ch_two)
    real_overlap_count = np.count_nonzero(real_overlap)
    ch_one_count = np.count_nonzero(ch_one)
    ch_two_count = np.count_nonzero(ch_two)
    real_fraction_overlap = real_overlap_count / np.count_nonzero(ch_two)

    if debug > 0:
        print(f'real overlap is {str(real_overlap_count)}')
        print(f'{str(ch_one_count)} {str(ch_two_count)} \n')

        #
        # Get properties for the labeled regions of the two sets of foci
        #

    regions1 = regionprops(ch_one)
    kernel = np.ones((3, 3), np.uint8)

    spotnum_output = []
    dilateimagelist = []

    for props in regions1:
        # if debug > 0:
        #     print(f'Label is {str(props.label)}\t')
        spot1_test_image = np.zeros_like(ch_one)
        spot1_test_image[props.coords[:,0], props.coords[:,1]] = props.label
        spot1_test_image2 = dilation(spot1_test_image, kernel)
        dilateimagelist.append(spot1_test_image2)
        dil_overlap = np.logical_and(spot1_test_image2, ch_two)
        adjacent_spots = ch_two[dil_overlap]
        unique_spots = list(set(adjacent_spots))
        spotnum_output.append(len(unique_spots))

    real_mean = mean(spotnum_output)
    real_median = median(spotnum_output)

    return(real_overlap_count, real_fraction_overlap, real_mean, real_median)

def prep_channel(channel):
    try:
        (mvalue, count) = mode(channel[channel>0])
        channel[channel==0] = mvalue
        channel = stack.log_filter(channel, sigma=3)
    except OSError as err:
        print("OS error:", err)
    except ValueError as err:
        print("Got a value error.", err)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise
    #print("Prep complete")
    return(channel)

#
# function reads an image, processes it
#

def process_file(fname, indir, outdir, one, two, dapi, tolerance,
                 debug, FH):

    bn = os.path.splitext(os.path.basename(fname))[0]
    fqn = os.path.join(indir, fname)
    segfn = str(bn) + "_segments.png"
    segfqn = os.path.join(outdir, segfn)

    if is_non_zero_file(segfqn):
        print(f"File {bn} appears to have been processed already, skipping...\n")
        return 0

    if debug > 0:
        print("Processing " + fname)

    try:
        img1 = io.imread(fqn)
    except OSError as err:
        print("OS error:", err)
    except ValueError as err:
        print("Got a value error.", err)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise

    spot1 = img1[:, :, one]
    spot2 = img1[:, :, two]
    dapi = img1[:, :, dapi]

    if debug > 0:
        print("About to prep")

    spot1 = prep_channel(spot1)
    spot2 = prep_channel(spot2)

    #
    # Divide the image into three intensity classes for the spot channels.
    # Nucleus uses the basic two class otsu threshold
    #
    spot1_thresholds = threshold_multiotsu(spot1, classes=3)
    spot2_thresholds = threshold_multiotsu(spot2, classes=3)
    #
    # Now cleaning up spot segmentation, collapse two less intense
    # bins into one bin.  Spots are third, most intense bin
    #
    spot1_regions = np.digitize(spot1, bins=spot1_thresholds)
    spot2_regions = np.digitize(spot2, bins=spot2_thresholds)

    spot1_regions[spot1_regions<2] = 0
    spot2_regions[spot2_regions<2] = 0
    #
    # Segmenting and cleaning up the nucleus
    #
    nuc_thresholds = threshold_otsu(dapi)
    binary = dapi > (0.9 * nuc_thresholds)  # Create a binary image using the threshold
    kernel = np.ones((3, 3), np.uint8)      # create a 3x3 kernel for opening function
    filled = ndi.binary_fill_holes(binary)  # Fill small holes
    open_img = ndi.binary_opening(filled, structure=kernel)     # Remove small objects
    open_img = open_img.astype(np.uint8)
    dilate = cv2.dilate(open_img, kernel, iterations=1)         # Dilate/erode to smoothen
    nucleus = cv2.erode(dilate, kernel, iterations=1)

    #
    # Set up segmentation plot and save it
    #
    pfig, ax = plt.subplots(ncols=3, figsize=(10, 5))

    ax[0].imshow(spot1_regions)
    ax[0].set_title('Spots 1')
    ax[0].axis('off')
    ax[1].imshow(spot2_regions)
    ax[1].set_title('Spots 2')
    ax[1].axis('off')
    ax[2].imshow(nucleus)
    ax[2].set_title('nucleus')
    ax[2].axis('off')


    plt.savefig(segfqn, bbox_inches='tight')
    plt.close()

    if debug > 0:
        print("Wrote the segmentations")

    #
    # Now label the spots.
    # Spot markers are in markers1 and markers2
    #
    spot2_regions = np.uint8(spot2_regions)
    spot1_regions = np.uint8(spot1_regions)
    ret, markers2 = cv2.connectedComponents(spot2_regions)
    ret, markers1 = cv2.connectedComponents(spot1_regions)

    (opc, opf, osmn, osmd) = spot_channel_stats(markers1, markers2, debug)

    if debug > 0:
        print(f'real overlap is {str(opc)} {str(opf)} {str(osmn)} {str(osmd)}\n')

    # Erode the nucleus boundary a bit to keep randomized foci in the nucleus

    img_u8 = nucleus.astype(np.uint8)

    erode_nucleus = cv2.erode(img_u8, kernel, iterations=1)
    ret, nuc = cv2.connectedComponents(erode_nucleus)
    nuc_props = regionprops_table(nuc, properties=('centroid',
                                                   'orientation',
                                                   'area',
                                                   'coords'))

    range_output = []
    rng = np.random.default_rng()
    test_mean_list = []
    test_median_list = []
    regions1 = regionprops(markers1)
    regions2 = regionprops(markers2)

    if debug > 0:
        print("about to randomize")

    # We make 1000 test images
    for i in range(1, 1000):
        # create the test image
        # sys.stdout.write(".")
        # sys.stdout.flush()
        test_image = np.zeros_like(spot2)

        for props in regions2:
            # We try to place the spot randomly in the nucleus.
            # If it doesn't work the first time we just try again
            for attempt in range(10):
                try:
                    movetoindex = rng.choice(np.arange(0,len(nuc_props['coords'][0])))
                    new_centroid = nuc_props['coords'][0][movetoindex]
                    x0, y0 = props.centroid
                    x1, y1 = [round(x0), round(y0)] - new_centroid

                    new_region = props.coords - (x1, y1)

                    test_image[new_region[:,0], new_region[:,1]] = props.label
                except:
                    next
                else:
                    break
            else:
                pass

        (topc, topf, tosmn, tosmd) = spot_channel_stats(markers1, test_image, debug)

        range_output.append(topf)
        test_mean_list.append(tosmn)
        test_median_list.append(tosmd)

    #print("")
    ol_pctile = percentileofscore(range_output, opf)
    mn_pctile = percentileofscore(test_mean_list, osmn)

    md_count = median(range_output)
    md_mean = median(test_mean_list)

    FH.write(f'{str(bn)}\t{str(len(regions1))}\t{str(len(regions2))}\t')
    FH.write(f'{opf:.3f}\t{ol_pctile:.3f}\t{md_count:.3f}\t')
    FH.write(f'{osmn:.3f}\t{mn_pctile:.3f}\t{md_mean:.3f}\n')

    plt.clf()

    pfig, ax = plt.subplots(ncols=2, figsize=(10, 5))

    N, bins, patches = ax[0].hist(range_output, 20)
    ax[0].axvline(opf, color='r')
    ax[0].annotate(f'Percentile = {ol_pctile:.1f}',
                  xy=(0, N.max()),
                  xytext=(4,0),
                  textcoords='offset points')
    N, bins, patches = ax[1].hist(test_mean_list, 20)
    ax[1].axvline(osmn, color='r')
    ax[1].annotate(f'Percentile = {mn_pctile:.1f}',
                  xy=(0, N.max()),
                  xytext=(4,0),
                  textcoords='offset points')

    histfn = str(bn) + "_hist.png"
    histfqn = os.path.join(outdir, histfn)
    plt.savefig(histfqn, bbox_inches='tight')
    plt.close()

def exit_gracefully():
    print("Processing complete, exiting")
    sys.exit()

#
# Main loop.
# Iterate through files in indir, processing each one
#
def main() -> int:
    print("")

    outtxt1 = os.path.join(args.outdir, 'image_list.txt')

    try:
        FH = open(outtxt1, 'a', buffering=1)
        FH.write("Filename\t#Spots1\t#Spots2\t")
        FH.write("Overlap\tOverlap Pctile\tMedian of Test Range\t")
        FH.write("Mean Overlapped\tMean Percentile\tMean of test list\n")
    except OSError as err:
        print("OS error:", err)
    except ValueError as err:
        print("Got a value error:", err)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise

    total = len(glob.glob("*.tif", root_dir=args.indir, recursive=False))

    print(f"Processing {str(total)} images...\n")

    with alive_bar(total) as pbar:
        for filename in glob.iglob("*.tif", root_dir=args.indir, recursive=False):
            try:
                with time_limit(300, filename):
                    process_file(filename, args.indir, args.outdir, args.one,
                                 args.two, args.dapi, args.tolerance, args.debug,
                                 FH)
            except TimeoutException as e:
                print(f"Timed out! {e} \n")
            except:
                pass
            pbar()  # call `bar()` at the end
    return 0

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
    finally:
        exit_gracefully()

