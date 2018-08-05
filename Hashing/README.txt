Currently (july 2018) hashes all jpg files in the current directory where the .exe is located. 

Creates 1 txt file per jpg in the current directory with the hamming distance to every other jpg file in the current directory. 

To change output: Change the returnstr on line 432.

To change hashing algorithm: Change the text on line 456 to either aHash, dHash, pHash or gHash.