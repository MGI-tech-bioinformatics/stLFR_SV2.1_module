# stLFR_SV2.1

## Introuction:

This process is applicable to stLFR technology and similar linked read data. Currently running on stLFR data, theoretically it is also applicable to other linked read data. You can analyze and test the barcode by converting it to the "read_id # XXX_XXX_XXX" format on read ID.

This process uses barcode information to detect breakpoint signals of structural variation (SV), such as: equilibrium translocation, INV, partial missing repetitions, and more complex structural breakpoints, which can be combined with CNV results and phase results The accuracy of deriving the true structural variation of chromosomes is limited by the distribution density of linked reads on DNA molecules. For example, stLFR kits built with 1.5ng starting volume can produce data that can guarantee SV detection accuracy above 20K.

## Directory Structure:
'''
bin data example LFR-sv.pl Readme.docx run.sh tools
'''

Main program: LFR-sv

bin: related programs needed to run the process

lib: possibly missing support libraries

data: some pre-made database files, not required, other database files can be specified

example: Examples of related files, not required

tools: gadget directory, not required

Readme.docx: this document

run.sh: run script example

## Parameter Description：

**Name:**
        LFR-sv
        version 2.1
 
**Function:**
        Detect the SVs from stLFR WGS data
        
**Usage:**
        perl LFR-sv.pl -bam prefix.sorted.markdup.bam -o out.result.dir
        
**Options:**
        -bam <string>   original sorted and markduped bam file,if the index dose not exist, will be created.[necessary]
        
        -out <string>   output SV dir.[necessary](warning:if exists, the output dir will be cleaned first!!!!!)
        
        -ncpu <int>     thread number for running pipeline.[default 1]
        
        -bar_th <int> at least N read pairs in one barcode.[default 8]
        
        -seg_th <int> at least N read pairs in one segment.[default 4]
        
        -gap <int> define the gap size which should be considered as different segment[default 20000]
        
        -is <int> proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]
        
        -bin <int> bin size for cluster the segments.[default 2000]
        
        -merge1 <int> N continue outline bins could be considered as the same break point and will be merged into one evidence.[default 5]
        
        -merge2 <int> SVs nearby under N binsize will be considered as one event.[default 5]
        
        -mmax <int> the max SVs allowed in one event.[default 4]
        
        -low1 <float/int> single end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.95]
        
        -low2 <float/int> end to end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.9999]
        
        -ex1 <float> when low1 is a float of 0-1, exclude the bins which depth under ex1.[default 0.2]
        
        -ex2 <float> when low2 is a float of 0-1, exclude the bins which depth under ex2.[default 0.2]
        
        -phase <string> formatted phase result directory including phased barcode and region by chromosome. For details, see tools/gen_phase/phase.readme and the description of the phase section below. [default NULL]
        
        -bl <string> black list file(BED format).[default NULL]
        
        -cl <string> sorted control list file(BEDPE format).[default NULL](Be sure the chromosome and position are sorted in one line!!!)
        
        -sc <int> allow max sv counts for the same position in one direction.[default 1]
        
        -human <Y/N> for Homo sapiens,keep only [1234567890XYM] chromosome.[default N]
        
        -qc1 <float> valid read pair ratio for SV detection.[default 0.60]
        
        -qc2 <float> average read pair count for one barcode.[default 30]
        
        -qc3 <float> average segment end count for one bin.[default 8]
        
        -rlen <int> read length of one read.[default 100]
        
        -mlen <int> physical limit for the long DNA segment.[default 400000]
        
        -help Show this message.
