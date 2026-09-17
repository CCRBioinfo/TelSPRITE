
#Telomere-Centromere proximity analysis

We wished to determine whether centromeres and telomeres are colocalized more often than expected by chance.  To make this assessment, we compared the observed localization with a distribution of 1000 randomizations of the position of the telomeres (the centromeres were held constant).

The script records overlap both in terms of pixels overlapped and in terms of foci overlapped, as explained below.

*Telomere-centromere-proximity.py* takes a directory of segmented nuclei (one nucleus per file)
and saves a list including the following information for each nucleus:

|Filename|The base name of the file (not fully qualified and no file extension)|
|\#Spots1|\# of spots in channel 1 (count not pixels)|
|#Spots2|# of spots in channel 2 (count not pixels)|
|Overlap|divide pixels of overlap by pixels in channel 2 (actual)|
|Overlap Pctile|percentile of actual “Overlap" relative to 1000x randomized overlaps |
|Median of Test Range|median of range of “Overlap" in 1000x randomized overlaps|
|Mean Overlapped|\#Spots1 (count) / \# of overlapping Spots 2 (count) = for each channel 1 spot, the mean # of channel 2 spots that contact|
|Mean Percentile|percentile of actual "Mean Overlapped” relative to 1000x randomized|
|Mean of test list|MEDIAN of “Mean Overlapped” for the 1000x randomized|

Additionally, for each nucleus the script saves a QC image with segmentations of each channel, as well as an image containing histograms of both pixel count and focus count overlap.

To segment images, we used *segment_all_channels.py*.  This script takes a folder of multichannel TIFF images and outputs a folder with one cropped nucleus image per file.

`usage: Telomere-centromere-proximity.py [-h] -i INDIR -o OUTDIR [--min [N]] [--max [N]] [--tolerance [N]] [--one [N]] [--two [N]] [--dapi [N]] [--debug [N]]

Quantify proximity of spots in two channels

options:
  -h, --help           show this help message and exit
  -i, --indir INDIR    Directory to process, containing TIF images of single nuclei
  -o, --outdir OUTDIR  Directory for output files
  --min [N]            minimum focus area in pixels (default 9)
  --max [N]            maximum focus area in pixels (default 10000)
  --tolerance [N]      Neighbor gap tolerance in pixels (default 5)
  --one [N]            first channel for spots, typically Cen (zero indexed)
  --two [N]            second channel for spots, typically Tel (zero indexed)
  --dapi [N]           channel for DAPI
  --debug [N]          Debug Setting, set to > 0 for debugging messages`

usage: segment_all_channels.py [-h] -i INDIR -o OUTDIR [--min MIN] [--max MAX] [--dapi DAPI] [--debug] [-t] [-b]

`Segment nuclei in TIF files

options:
  -h, --help           show this help message and exit
  -i, --indir INDIR    Directory to process, containing TIF images with fields of nuclei
  -o, --outdir OUTDIR  Directory for segmented files
  --min MIN            minimum nucleus area in pixels (default 10000)
  --max MAX            maximum nucleus area in pixels (default 100000)
  --dapi DAPI          channel for dapi (zero indexed)
  --debug              Set debug to TRUE
  -t, --test           Test the current settings
  -b, --blank          set the background of the segment to zero`

As a work of the United States Government, this project is in the public 
domain within the United States.

Additionally, we waive copyright and related rights in the work worldwide 
through the CC0 1.0 Universal public domain dedication. 