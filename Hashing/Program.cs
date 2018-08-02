using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using System.Drawing;
using System.IO;
using System.Windows.Media.Imaging;
using System.Collections.Concurrent;
using System.Threading;

namespace Hashing
{
    class Program
    {
        // The hashing part (aka all below) is borrowed from https://github.com/juliasin/SimilarByHashes
        public List<int> set1 { get; set; }
        public List<int> set2 { get; set; }
        public double dist { get; set; }
        public string ImagePath { get; set; }
        public string DirectoryPath { get; set; }

        static public UInt64 aHash(Bitmap img)
        {
            //Bitmap bmp = BitmapImage2Bitmap(img);
            Bitmap bmp = img;
            bmp = ResizeImage(bmp, new System.Drawing.Size(8, 8));
            bmp = GreyScalling(bmp);
            double[,] k = Matrix(bmp);
            BitArray bits = SetOfBits1(k);
            UInt64 hash = Hash(bits);
            return hash;
        }

        static public UInt64 pHash(Bitmap img)
        {
            Bitmap bmp = img;
            bmp = ResizeImage(bmp, new System.Drawing.Size(32, 32));
            bmp = GreyScalling(bmp);
            double[,] k = Matrix(bmp);
            double[,] dct = DCT(k);
            double[,] reduce = Reduce(dct, 8);
            BitArray bits = SetOfBits1(reduce);
            UInt64 hash = Hash(bits);
            return hash;
        }

        static public UInt64 dHash(Bitmap img)
        {
            Bitmap bmp = img;
            bmp = ResizeImage(bmp, new System.Drawing.Size(9, 8));
            bmp = GreyScalling(bmp);
            double[,] k = Matrix(bmp);
            double[,] difmatr = differenceMatrix(k);
            BitArray bits = SetOfBitsDHash(difmatr);
            UInt64 hash = Hash(bits);
            return hash;
        }

        static public UInt64 gHash(Bitmap img)
        {
            Bitmap bmp = img;
            bmp = ResizeImage(bmp, new System.Drawing.Size(32, 32));
            bmp = GreyScalling(bmp);
            double[,] k = Matrix(bmp);
            BitArray bits = SetOfBitsGHash(k);
            UInt64 hash = Hash(bits);
            return hash;
        }

        static public BitArray SetOfBitsGHash(double[,] d)
        {
            int s = d.GetLength(0);
            double[] sumcols = new double[s];
            for (int j = 0; j < s; j++)
            {
                for (int i = 0; i < s; i++)
                {
                    sumcols[j] += d[i, j];
                }
            }
            double[] sumrows = new double[s];
            for (int i = 0; i < s; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    sumrows[i] += d[i, j];
                }
            }
            BitArray arr = new BitArray(64);
            int k = 63;
            for (int i = 0; i < sumcols.Length - 1; i++)
            {
                if (sumcols[i] > sumcols[i + 1]) arr[k] = true; else arr[k] = false;
                k--;
            }
            if (sumcols[sumcols.Length - 1] > sumcols[0]) arr[k] = true; else arr[k] = false;
            k--;
            for (int i = 0; i < sumrows.Length - 1; i++)
            {
                if (sumrows[i] > sumrows[i + 1]) arr[k] = true; else arr[k] = false;
                k--;
            }
            if (sumrows[sumrows.Length - 1] > sumrows[0]) arr[k] = true; else arr[k] = false;
            k--;

            return arr;
        }

        static public double[,] differenceMatrix(double[,] d)
        {
            double[,] difmatr = new double[8, 8];
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    difmatr[i, j] = d[i, j + 1] - d[i, j];
                }
            }
            return difmatr;
        }

        static public BitArray SetOfBitsDHash(double[,] y)
        {
            BitArray arr = new BitArray(8 * 8); int k = 0;
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    if (y[i, j] > 0) arr[k] = true; else arr[k] = false;
                    k++;
                }
            }
            return arr;
        }

        static public double[,] DCT(double[,] matrix)
        {
            int size = matrix.GetLength(0);
            double[,] resultmatr = new double[size, size];
            for (int u = 0; u < size; u++)
            {
                for (int v = 0; v < size; v++)
                {
                    double result = 0d;
                    for (int i = 0; i < size; i++)
                    {
                        for (int j = 0; j < size; j++)
                        {
                            result += (Alpha(i) * Alpha(j) * Math.Cos(((Math.PI * u) / (2 * size)) * (2 * i + 1)) * Math.Cos(((Math.PI * v) / (2 * size)) * (2 * j + 1)) * matrix[i, j]);
                        }
                    }
                    resultmatr[u, v] = Math.Round((result * 2 / Math.Sqrt(size * size)) * (1d / size + 1d / size));
                }
            }
            return resultmatr;
        }

        static private double Alpha(int u)
        {
            if (u == 0)
                return 1 / Math.Sqrt(2);
            else
                return 1;
        }

        static public double[,] Reduce(double[,] resultmatr, int newsize)
        {
            // double[,] tmp = new double[oldsize, oldsize];
            double[,] reducematrix = new double[newsize, newsize];
            for (int i = 0; i < newsize; i++)
            {
                for (int j = 0; j < newsize; j++)
                {
                    reducematrix[i, j] = resultmatr[i, j];
                }
            }
            return reducematrix;
        }

        static public Int64 Hamming(UInt64 x, UInt64 y)
        {
            Int64 dist = 0;
            UInt64 val = x ^ y;

            // Count the number of bits set
            while (val != 0)
            {
                // A bit is set, so increment the count and clear the bit
                dist++;
                val &= val - 1;
            }

            // Return the number of differing bits
            return dist;
        }

        static public BitArray SetOfBits1(double[,] cl)
        {
            int size = cl.GetLength(0);
            double avg = 0;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    avg += cl[i, j];
                }
            }
            avg = avg / (size * size);
            BitArray arr = new BitArray(size * size); int k = 0;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (cl[i, j] < avg) arr[k] = false; else arr[k] = true;
                    k++;
                }
            }
            return arr;
        }

        static public UInt64 Hash(BitArray arr)
        {
            var bytes = new byte[8];
            arr.CopyTo(bytes, 0);
            UInt64 hash = BitConverter.ToUInt64(bytes, 0);
            return hash;
        }

        static public Bitmap ResizeImage(Bitmap imgToResize, System.Drawing.Size size)
        {
            return new Bitmap(imgToResize, size);
        }

        static public Bitmap GreyScalling(Bitmap c)
        {
            Bitmap d = new Bitmap(c.Width, c.Height);
            for (int i = 0; i < c.Width; i++)
            {
                for (int x = 0; x < c.Height; x++)
                {
                    System.Drawing.Color oc = c.GetPixel(i, x);
                    int grayScale = (int)((oc.R * 0.3) + (oc.G * 0.59) + (oc.B * 0.11));
                    System.Drawing.Color nc = System.Drawing.Color.FromArgb(oc.A, grayScale, grayScale, grayScale);
                    d.SetPixel(i, x, nc);
                }
            }
            return d;
        }

        static public double[,] Matrix(Bitmap b)
        {
            System.Drawing.Color[,] clrs = new System.Drawing.Color[b.Height, b.Width];
            double[,] d = new double[b.Height, b.Width];
            for (int i = 0; i < b.Height; i++)
                for (int j = 0; j < b.Width; j++)
                {
                    clrs[i, j] = b.GetPixel(j, i);
                    d[i, j] = clrs[i, j].B;
                }
            return d;

        }

        // End of borrowed code.


        
        public string FilenameAnd_a_Hash(string file)
        {
            Bitmap bm = new Bitmap(file);
            UInt64 filehash = aHash(bm);
            bm.Dispose();
            //string filename = file + "__a_hash__" + filehash+".jpg";
            string filename = file + "__a_hash__" + filehash; // x:\folder\filename.jpg
            
            int removeCount = Directory.GetCurrentDirectory().Count(); // count of x:\folder\
            string actualFilename = file.Remove(0, removeCount+1); // 0 based starting count, result: filename.jpg

            string returnName = actualFilename+";"+filehash; // filename.jpg;HASH
            //string returnName = actualFilename.Remove(actualFilename.Count() - 4, 4);

            Console.WriteLine(returnName);
            return returnName;
        }
        
        public string FilenameAnd_d_Hash(string file)
        {
            Bitmap bm = new Bitmap(file);
            UInt64 filehash = dHash(bm);
            bm.Dispose();
            //string filename = file + "__d_hash__" + filehash + ".jpg";
            string filename = file + "__d_hash__" + filehash;
            string actualFilename = file.Remove(0, Directory.GetCurrentDirectory().Count()+ 1); // 0 based starting count, result: filename.jpg

            string returnName = actualFilename + ";" + filehash; // filename.jpg;HASH


            return returnName;
        }
        public string FilenameAnd_g_Hash(string file)
        {
            Bitmap bm = new Bitmap(file);
            UInt64 filehash = gHash(bm);
            bm.Dispose();
            //string filename = file + "__g_hash__" + filehash + ".jpg";
            string filename = file + "__g_hash__" + filehash;
            int removeCount = Directory.GetCurrentDirectory().Count(); // count of x:\folder\
            string actualFilename = file.Remove(0, removeCount + 1); // 0 based starting count, result: filename.jpg

            string returnName = actualFilename + ";" + filehash; // filename.jpg;HASH

            return returnName;
        }
        public string FilenameAnd_p_Hash(string file)
        {
            Bitmap bm = new Bitmap(file);
            UInt64 filehash = pHash(bm);
            bm.Dispose();
            //string filename = file + "__p_hash__" + filehash + ".jpg";
            string filename = file + "__p_hash__" + filehash;
            int removeCount = Directory.GetCurrentDirectory().Count(); // count of x:\folder\
            string actualFilename = file.Remove(0, removeCount + 1); // 0 based starting count, result: filename.jpg

            string returnName = actualFilename + ";" + filehash; // filename.jpg;HASH

            return returnName;
        }

        public string csvInput_aHash(string currentfile, string comparefile)
        {
            Bitmap bmCurrent = new Bitmap(currentfile);
            UInt64 filehashCurrent = aHash(bmCurrent);
            bmCurrent.Dispose();

            Bitmap bmCompare = new Bitmap(comparefile);
            UInt64 filehashCompare = aHash(bmCompare);
            bmCompare.Dispose();
            string actualCurrent = currentfile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg
            string actualcompare = comparefile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg


            Int64 hammingDistance = Hamming(filehashCurrent, filehashCompare);

            return actualCurrent+";"+actualcompare + ";" + hammingDistance;
        } // result: currentfile;comparefile;hammingDistance
        public string csvInput_dHash(string currentfile, string comparefile)
        {
            Bitmap bmCurrent = new Bitmap(currentfile);
            UInt64 filehashCurrent = dHash(bmCurrent);
            bmCurrent.Dispose();

            Bitmap bmCompare = new Bitmap(comparefile);
            UInt64 filehashCompare = dHash(bmCompare);
            bmCompare.Dispose();
            string actualCurrent = currentfile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg
            string actualcompare = comparefile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg


            Int64 hammingDistance = Hamming(filehashCurrent, filehashCompare);

            return actualCurrent + ";" + actualcompare + ";" + hammingDistance;
        }
        public string csvInput_gHash(string currentfile, string comparefile)
        {
            Bitmap bmCurrent = new Bitmap(currentfile);
            UInt64 filehashCurrent = gHash(bmCurrent);
            bmCurrent.Dispose();

            Bitmap bmCompare = new Bitmap(comparefile);
            UInt64 filehashCompare = gHash(bmCompare);
            bmCompare.Dispose();

            string actualCurrent = currentfile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg
            string actualcompare = comparefile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg

            
            Int64 hammingDistance = Hamming(filehashCurrent, filehashCompare);

            return actualCurrent + ";" + actualcompare + ";" + hammingDistance;
        }
        public string csvInput_pHash(string currentfile, string comparefile)
        {
            Bitmap bmCurrent = new Bitmap(currentfile);
            UInt64 filehashCurrent = pHash(bmCurrent);
            bmCurrent.Dispose();

            Bitmap bmCompare = new Bitmap(comparefile);
            UInt64 filehashCompare = pHash(bmCompare);
            bmCompare.Dispose();

            string actualCurrent = currentfile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg
            string actualcompare = comparefile.Remove(0, Directory.GetCurrentDirectory().Count() + 1); // 0 based starting count, result: filename.jpg


            Int64 hammingDistance = Hamming(filehashCurrent, filehashCompare);

            return actualCurrent + ";" + actualcompare + ";" + hammingDistance;
        }

        public string hashToCluster(string currentfile, string hashAlgo)
        {
            Bitmap bm = new Bitmap(currentfile);
            UInt64 filehash = 0;
            if (hashAlgo=="aHash")
            {
                filehash = aHash(bm);
            }
            else if (hashAlgo=="dHash")
            {
                filehash = dHash(bm);
            }
            else if (hashAlgo=="gHash")
            {
                filehash = gHash(bm);
            }
            else if (hashAlgo=="pHash")
            {
                filehash = pHash(bm);
            }
            else
            {
                Console.WriteLine("Set hash: aHash, dhash, gHash or pHash");
            }
                       
            bm.Dispose();

            return filehash.ToString();

        }

        static void Main(string[] args)
        {
            // 2 stages: 
            // stage 1: For all jpgs in the current folder,
            // Do hashing (all 3? To be tested)
            // Stage 1.1: make a file with currentName;hashName
            // Stage 2
            // Sort the hashes based on the hamming distance
            // 2.1 Select one image and do the sorting from that

            Program p = new Program();

            // stage 1: 
            string Currentdir = Directory.GetCurrentDirectory();
            string SearchPattern = "*.jpg";
            string[] files = Directory.GetFiles(Currentdir, SearchPattern);
            List<string> ahashlist = new List<string>();
            List<string> dhashlist = new List<string>();
            List<string> ghashlist = new List<string>();
            List<string> phashlist = new List<string>();

            ConcurrentBag<string> ahashcb = new ConcurrentBag<string>();
            ConcurrentBag<string> dhashcb = new ConcurrentBag<string>();
            ConcurrentBag<string> ghashcb = new ConcurrentBag<string>();
            ConcurrentBag<string> phashcb = new ConcurrentBag<string>();


            //string fileToWritea = "hashes_ahash.txt";
            //string fileToWrited = "hashes_dhash.txt";
            //string fileToWriteg = "hashes_ghash.txt";
            //string fileToWritep = "hashes_phash.txt";


            //TextWriter twa = new StreamWriter(fileToWritea);
            //TextWriter twd = new StreamWriter(fileToWrited);
            //TextWriter twg = new StreamWriter(fileToWriteg);
            //TextWriter twp = new StreamWriter(fileToWritep);

           

            DateTime now = DateTime.Now;
            Console.WriteLine(now);
            
            int totalLinesWritten = 0;
            string hashname = "pHash";
            TextWriter twh = new StreamWriter("hash_"+hashname+".txt");

            Parallel.ForEach(files, (currentFile) =>
            {
                //string actualCurrentFile = currentFile.Remove(0, Directory.GetCurrentDirectory().Count() + 1);
            

                ahashcb.Add(p.hashToCluster(currentFile, hashname));
                

            });

            // sequential through all the files, compare
            /*foreach (var currentFile in files)
            {
                string actualCurrentFile = currentFile.Remove(0, Directory.GetCurrentDirectory().Count() + 1);

                TextWriter actualTW = new StreamWriter(actualCurrentFile+".txt"); // Create a txt file for each jpg file in the current dir
                // Do I want one file per image or one file per image per hashing algorithm?
                // I think I want one file per image with 1 algorithm testing all 4 to find the best one. 
                // the way hashToCluster is written is to make it easy to switch between the algorithms. 

                Parallel.ForEach(files, (search) =>
                {
                    if(currentFile!=search)
                    {
                        //ahashcb.Add(p.csvInput_aHash(currentFile, search));
                        //dhashcb.Add(p.csvInput_dHash(currentFile, search));
                        //ghashcb.Add(p.csvInput_gHash(currentFile, search));
                        //phashcb.Add(p.csvInput_pHash(currentFile, search));
                        ahashcb.Add(p.hashToCluster(search, "aHash"));
                        Interlocked.Increment(ref filenumber);
                    }
                    else if (currentFile==search)
                    {

                    }
                    
                    
                    
                });
                */
            //Console.WriteLine(filenumber);
            //Console.WriteLine(ahashcb.Count());

            // TODO: Make the writing more effecient + in parallel.

            foreach (var item in ahashcb)
            {
                twh.WriteLine(item);
                totalLinesWritten++;
            }


            foreach (var item in ahashcb)
                {
                    //actualTW.WriteLine(item);
                    //totalLinesWritten++;
                }
                
                foreach (var item in dhashcb)
                {
                    //twd.WriteLine(item);
                }
                foreach (var item in ghashcb)
                {
                    //twg.WriteLine(item);
                }
                foreach (var item in phashcb)
                {
                    //twp.WriteLine(item);
                }

                totalLinesWritten = 0;
                while (ahashcb.Count > 0)

                {

                    string element;

                    if (ahashcb.TryTake(out element))

                    {

                    }

                }


            



            DateTime afterWork = DateTime.Now;
            Console.WriteLine("Finished at "+afterWork);
        }






        // Stage 2
        // Sort the hashes based on the hamming distance
        // 2.1 Select one image and do the sorting from that

        // The _*list_ now looks like this: filename.jpg;HASH  filename.jpg;HASH etc
        // I want something that ties those two together, so when I sort the hash, the filenames are sorted too
        // Key:Value thing an idea?
        // Also, how do I sort by hamming distance? I have a function ^^ that does something like it,
        // but do I use it on the entire list and sort by lowest difference (and how do I get lowest difference?)
        //

        // Ideally I want a database type thing where any sorting change results in moving of all the attached data
        // hash attached to filename
        // hamming distance attached to filename.
        // 

        //Idea: 
        // For each jpg, create a csv file with all the jpg-names;hamming distance from current file
        // this way we can do some excel sorting without resorting to something crazy like a 3-layer list. 






        // OLD: 
        //Stage 1.1: 
        // Create text file, keep it open
        //string fileToWritea = "hashesa.txt";
        //string fileToWrited = "hashesd.txt";
        //string fileToWriteg = "hashesg.txt";
        //string fileToWritep = "hashesp.txt";


        //TextWriter twa = new StreamWriter(fileToWritea);
        //TextWriter twd = new StreamWriter(fileToWrited);
        //TextWriter twg = new StreamWriter(fileToWriteg);
        //TextWriter twp = new StreamWriter(fileToWritep);


        // Line 1: 
        // Hash type

        // in the loop, write filename;hashName


        //string hasha = "";
        //string hashd = "";
        //string hashg = "";
        //string hashp = "";

        //twa.WriteLine("aHash");
        //twd.WriteLine("dHash");
        //twg.WriteLine("gHash");
        //twp.WriteLine("pHash");

        /*
        foreach (string file in files)
        {
            //aHash   
            hasha = p.FilenameAnd_a_Hash(file);
            twa.WriteLine(hasha);

        }
        //tw.WriteLine("dHash");
        foreach (string file in files)
        {
            // dhash
            hashd = p.FilenameAnd_d_Hash(file);
            twd.WriteLine(hashd);

        }
        //tw.WriteLine("gHash");
        foreach (string file in files)
        {
            //gHash
            hashg = p.FilenameAnd_g_Hash(file);
            twg.WriteLine(hashg);
        }
        //tw.WriteLine("pHash");
        foreach (string file in files)
        {
            // pHash
            hashp = p.FilenameAnd_p_Hash(file);
            twp.WriteLine(hashp);
        }
        */
        /*
        TextWriter.Synchronized(tw);

        tw.WriteLine("aHash");
        Parallel.ForEach(files, (currentFile) =>
        {
            tw.WriteLine(p.FilenameAnd_a_Hash(currentFile));
        });
        tw.WriteLine("dHash");
        Parallel.ForEach(files, (currentFile) =>
        {
            tw.WriteLine(p.FilenameAnd_d_Hash(currentFile));
        });
        tw.WriteLine("gHash");
        Parallel.ForEach(files, (currentFile) =>
        {
            tw.WriteLine(p.FilenameAnd_g_Hash(currentFile));
        });
        tw.WriteLine("pHash");
        Parallel.ForEach(files, (currentFile) =>
        {
            tw.WriteLine(p.FilenameAnd_p_Hash(currentFile));
        });
        */

        // close file. 
        //twa.Close();
        //twd.Close();
        //twg.Close();
        //twp.Close();



        //ahashlist.Add(p.FilenameAnd_a_Hash(currentFile)); // add    filename;aHash  to the list
        //dhashlist.Add(p.FilenameAnd_d_Hash(currentFile));
        //ghashlist.Add(p.FilenameAnd_g_Hash(currentFile));
        //phashlist.Add(p.FilenameAnd_p_Hash(currentFile));

    }
    }
