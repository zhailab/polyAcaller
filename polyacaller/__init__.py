__version__ = "0.2"

from .fast5readwrapper import Fast5ReadWrapper, extract_fast5_polyA, main

"""
write by Jinbu Jia, 2020.05.30

The input fast5 file is a file which has been basecalled by Guppy with
flipflop mode.

For extract the best polyA region:
python polyAcaller.py test_data/test.fast5 test.out 

For extract all potential polyA region:
python polyAcaller.py test_data/test.fast5 test.out 0 #default 1

It's better to limit the polyA search region in a specific region. 
Such as the region between adapter and genome mapping region. It's
better to set some padding region, such as:

                    genome mapping region (|)      unmapping region (-)  3' adapter (*)
read: ---***----|||||||||||||||||||||||||||||||-------------------------***************----
search region:                       ||||||||||-------------------------*****
                                     |--10nt--|                         |5nt|

For extract the best polyA region in a specific region:
python polyAcaller.py test_data/test.fast5 test.out 0 test_data/test_sample2.search_region.txt

test_data/test_sample2.search_region.txt (tab sperated):
read_id 	base 	search_start_base 	search_end_base
000478c6-4c63-4cb7-baca-ca4954a12ee6 	A 	266 	274

The next version will implemented adaper searching and genome mapping. But now
you can use minimap2 to perform genome mapping and use blastn or other tools to 
find adapter, and then set the search region based on the mapping information.
"""



        
