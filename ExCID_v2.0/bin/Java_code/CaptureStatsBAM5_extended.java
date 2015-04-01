import java.io.*;
import java.util.*;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Generates capture statistics from a BAM input file.
 * Generates 3 output files:
 * @author bainbrid
 *
 */
public class CaptureStatsBAM5_extended
{
	static final int prime_size = 1000;
	static int[] fivePrime = new int[prime_size]; //stores data about coverage upstream of the target 
	static int[] threePrime = new int[prime_size];//stores data about coverage downstream of the target 
	static int[] targetCov = new int[101]; //stores data about the coverage across a target
	static final int BUFFER = 100; //buffer region around a target
	static int tempcnt = 0; 
	static int totalReadsProduced = 0; //total reads contains in the bam
	static int totalReadsAligned = 0; //total reads aligned to a target region
	static long totalTargetCoverage = 0;  //total number of read bases aligned to the target
	static int totalReadsPaired = 0; //total number of reads with mate pairs (if any)
	static int totalPairedreadsWithmappedMates = 0; //total number of aligned reads which have mapped mates
	static int offtargetReadHitCount = 0; //total number of reads which do not align to a target region
	static int ontargetReadHitCount = 0; //total number of reads which align to a target region
	static int inbufferReadHitCount = 0; //total number of reads which align to the buffer region
	static int readsAlignedINCLUDINGDUPES = 0; //total reads aligned including duplicates
	static int duplicateReads = 0; //total number of duplicate reads
	static long totalAlignedBases = 0; //total number of aligned bases
	static int totalTargetedBases = 0; //total number of bases targeted
	static int totalBufferBases = 0; //total number of bases in the buffer region
	static int basesWithOneHitorMore = 0;   //total targeted bases with at least  1 coverage
	static int basesWith10HitsorMore = 0;  //total targeted bases with at least  10 coverage
	static int basesWith20HitsorMore = 0; //total targeted bases with at least  20 coverage
	static int basesWith40HitsorMore = 0; //total targeted bases with at least 40  coverage
        static int basesWith50HitsorMore = 0; //total targeted bases with at least 50  coverage
        static int basesWith100HitsorMore = 0; //total targeted bases with at least 100  coverage
        static int basesWith500HitsorMore = 0; //total targeted bases with at least 500  coverage
        static int basesWith1000HitsorMore = 0; //total targeted bases with at least 1000  coverage
	static int totalTargets = 0; //total taregted regions
	static int hitTargetCount = 0;  //total targets with at least 1 read aligned to them
	static int hitTarget_bufferonly_Count = 0; //total targets with no hits, except in buffer
	static int[] dupHisto = new int[9]; //depcrecated
	static int[] covHisto = new int[1001]; //coverage histogram
	static int nonTragetGoodHits = 0; //regions that have high coverage but are not in the target
	static String VERSION = "CapStatsV2.6 2015-03-25";
	static Hashtable<String, SAMFileWriter> fht ; 
	static final String dummy = "dummy";
	static boolean removeDupes = false; //do we not consider duplicate reads
	static boolean writeWGC = false; //write whole genome coverage statistics
	static double _percentage = 1.0; //number of read to take (to randomly dump some)
	static Random RAND = new Random(88651); //random number generator, good for removing a proportion of the reads
	static String[] targetChrs;
	static int[] targetStarts;
	static int[] targetStops;
	static double minmapscore = -1.0;
        static double minbasescore = -1.0;
	static int size =0 ; //Holds the size of the Chromosome being analyzed.
        static int[] coverage_forMedian = new int[1000]; //Used for calculating the median coverage.
	static int median_coverage = 0;
        static long READLENGTH = 0; //Gets READ Length
	
	/**
	 * ONLY WORKS ON A SINGLE FILES THAT HAS BEEN SORTED!!!!!!!!!!
	 * looks like it works in 4000M of ram.
	 * 
	 * usage: -o <output directory & file base name> -t <targetfile> - -i <BAM FILE> [-d] [-w] [-m] [-b] [-p]
	 * 
	 * i: BAM File
	 * t: target file
	 * o: output directory and base file name
	 * d: remove duplicates, and do not use them for statistics
	 * w: write whole genome coverage
	 * m: minimum mapscore (mapscore must be >= to this)
         * b: minimum base score (basescore must be >= to this)
	 * p: only take this proportion of reads into consideration (for scale back experiments)
	 * 
	 * 
	 * @throws Exception
	 * 
	 */
	public static void main(String[] args) throws Exception
	{
            if(args.length == 0) {usage();System.exit(0);};
            String[] validargs = {"r","t","o","d","i","w","p","m","b"};
            String[] arguments = new String[validargs.length];
            String warns = ParseOpts.parse(validargs, args, arguments);
            if(arguments[0] == null)
            {
                //deprecated.  no longer needed.
            }
            if(arguments[1] == null)
            {
                System.err.println("No target file specified!!! Exiting\n");
                usage();
                System.exit(2);
            }
            if(arguments[2] == null)
            {
                System.err.println("No output file specified!!! Exiting\n");
                usage();
                System.exit(3);
            }
            int index = arguments[2].lastIndexOf(File.separator);
            String fileName = arguments[2].substring(index + 1);
            if(fileName == null)
            {
                System.err.println("Please provide proper output file!!! Exiting\n");
                usage();
                System.exit(3);
            }
            if(arguments[3]!= null )
            {
                removeDupes = true;
            }
            if(arguments[4] == null)
            {
                System.err.println("No alignment file specified!!! Exiting\n");
                usage();
                System.exit(4);
            }
            if(arguments[5] != null)
            {
                writeWGC = true;
            }
            if(arguments[6] != null)
            {
                _percentage = Double.parseDouble(arguments[6]);
            }
            if(arguments[7] != null)
            {
                minmapscore = Double.parseDouble(arguments[7]);
            }
            if(arguments[8] != null)
            {
                minbasescore = Double.parseDouble(arguments[8]);
            }
            String targetFile = arguments[1];
            String alignmentFile = arguments[4];
            checkFile(alignmentFile);
            checkFile(targetFile);
            String outfile  = arguments[2];
            String covFile = outfile+".cov.fasta";
            System.out.println("Writing to: "+outfile);
            FileWriter wgcFasta = null;
            if(writeWGC)
            {
                wgcFasta =  new FileWriter(outfile+".wholeGenomeCov.fasta");
            }
            FileWriter covFasta = new FileWriter(covFile);
            FileWriter wig = new FileWriter(outfile+".missReads.wig");
            FileWriter missTar = new FileWriter(outfile+".missTargets.txt");
            fht = new Hashtable<String, SAMFileWriter>(50); 
            wig.write("track type=wiggle_0 name="+alignmentFile+"\n");
            loadTargets(targetFile);
            readBAM(alignmentFile,  wig, covFasta, missTar, wgcFasta);
            covFasta.close();
            if(writeWGC) wgcFasta.close();
            wig.close();
            writeReport(outfile);
            missTar.close();
	}
	
	/**
	 * load the target regions into memory and leave them there
	 * @param targetFile
	 */
	public static void loadTargets(String targetFile) throws Exception
	{
            BufferedReader br = new BufferedReader(new FileReader(targetFile));
            String line;
            int cnt = 0;
            int start = 0;
            int stop = 0;
            while((line = br.readLine())!=null)
            {
                String[] tokens = line.split("[ \t]+");
                if(tokens.length < 3) {continue;}
                //if(!tokens[0].substring(0,3).equalsIgnoreCase("chr"))continue;
                try
                {	
                    start = Integer.parseInt(tokens[1]);
                    stop = Integer.parseInt(tokens[2]);
                }
                catch(NumberFormatException e){continue;}
                cnt++;
            }
            targetChrs = new String[cnt];
            targetStarts = new int[cnt];
            targetStops = new int[cnt];
            cnt = 0;
            br.close();
            br = new BufferedReader(new FileReader(targetFile));
            while((line = br.readLine())!=null)
            {
                String[] tokens = line.split("[ \t]+");
                if(tokens.length < 3) {continue;}
                //if(!tokens[0].substring(0,3).equalsIgnoreCase("chr"))continue;
                try
                {	
                    start = Integer.parseInt(tokens[1]);
                    stop = Integer.parseInt(tokens[2]);
                }
                catch(NumberFormatException e){continue;}
                targetChrs[cnt] = tokens[0];
                targetStarts[cnt]=start;
                targetStops[cnt] = stop;
                cnt++;
            }
            totalTargets = cnt;     /* Should be used to count Targets Oct3_2014*///totalTargets = TR.length;
	}
	
	
	
	/**
	 * Major workhorse of the application.  Reads in the BAM file and processes each record.
	 * @param bamfile  The bam file to read in
	 * @param wig The output wig file to write offtarget high coverage regions
	 * @param covFasta the coverage fasta output file
	 * @param missTraget the miss target file
	 * @param wgCoverage whole genome coverage file
	 * @throws Exception
	 * 
	 */
	public static void readBAM(String bamfile,  FileWriter wig, FileWriter covFasta, FileWriter missTraget, FileWriter wgCoverage) throws Exception
	{
            SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
            SAMFileReader sfr  = new SAMFileReader(new File(bamfile));
            SAMFileHeader header = sfr.getFileHeader();
            Iterator<SAMRecord> iter = sfr.iterator();
            String lastchr = "SOMEVERYFAKETHINGGOESHERE";
            short[] COVERAGE = new short[0];
            char[] TR = new char[0];
            int cnt = 0;
            SAMRecord rec = null;
            while(iter.hasNext())
            {
                try
                {
                    rec = iter.next();
                    if(_percentage < 1.0)
                    {
                        double rd = RAND.nextDouble();
                        if(rd < _percentage) continue;
                    }

                    totalReadsProduced++;
                    if(rec.getMappingQuality() < minmapscore)
                    {
                        continue;
                    }
                    if(rec.getReadFailsVendorQualityCheckFlag()){
                        continue;
                    }
                                
                    if(rec.getNotPrimaryAlignmentFlag()){
                        continue;
                    }
                    
                    if(rec.getReadUnmappedFlag())
                    {
                        //System.out.println("Unmapped! region.  breaking!!!");
                        continue;
                    }
                    totalReadsAligned++;
                    
                    if(rec.getReadPairedFlag())
                    {
                        totalReadsPaired++;
                        if(!rec.getMateUnmappedFlag()){
                            totalPairedreadsWithmappedMates++;
                        }
                    }
                    
                    if(rec.getDuplicateReadFlag())
                    {
                        duplicateReads++;
                        if(removeDupes){continue;}
                    }
                    
                    /////////////////Added on Feb 11, 2015/////////////////////////////
                    READLENGTH = rec.getReadLength();
 
                    ///////////////////////////////////////////////////////////////////
                                
                    String currchr = rec.getReferenceName();
                    if(!currchr.equals(lastchr))
                    {
                        if(!lastchr.equals("SOMEVERYFAKETHINGGOESHERE"))
                        {
                            getTargetsAndWriteCoverage(lastchr,  COVERAGE, covFasta,  missTraget,  wgCoverage);
                            findWhereReadsHit(lastchr,  COVERAGE,wig);
                        }
                        lastchr = currchr;
                        System.out.println(currchr);
                        size = header.getSequence(currchr).getSequenceLength()+1;
                        COVERAGE = new short[size];
                        TR = getTargetPos(currchr, size);
                        if(TR == null || COVERAGE == null)
                        {
                            System.err.println("COVERAGE or TR is null! "+currchr);
                        }

                    }
                    processRecord(rec, TR, COVERAGE);
                }
                catch(Exception e)
                {
                    System.err.println("Error on record: "+cnt+"\n"+e.getMessage()+" "+Arrays.toString(e.getStackTrace()));
                    System.err.println(rec.toString()+" "+rec.getReferenceName()+" "+rec.getAlignmentStart()+" "+rec.getAlignmentEnd());
            //	throw e;
                }
                cnt++;
            }
            System.out.println("Done read bam");
            getTargetsAndWriteCoverage(lastchr,  COVERAGE, covFasta,  missTraget,  wgCoverage);
            findWhereReadsHit(lastchr,  COVERAGE,wig);
            COVERAGE = null;
	}
	
	/**
	 * Processes read record from a BAM.  Adding the alignment to coverage array.
	 * @param rec
	 * @param TR
	 * @param COVERAGE
	 */
	public static void processRecord(SAMRecord rec, char[] TR, short[] COVERAGE) throws Exception
	{
            boolean inbuffer = false;
            boolean ontarget = false;
            int start = rec.getAlignmentStart();
            int stop = rec.getAlignmentEnd();
            int referenceposition = 0;
            byte[] baseQual = rec.getBaseQualities();
            /*if(rec.getReadNegativeStrandFlag())
            {
                    start-=25;
            }
            else
            {
                    stop+=25;
            }*/

            for(int i = start; i <= stop; i++)
            {
                totalAlignedBases++;
                try
                {
                    if(TR[i] == 1)
                    {
                        ontarget = true;
                    }
                    else if(TR[i] == 2)
                    {
                        inbuffer = true;
                    }
                    
                    referenceposition = rec.getReferencePositionAtReadPosition(i-start+1);
                    if(referenceposition==0){continue;}
                    
                    if(minbasescore > 0){
                        if((double)baseQual[i-start] >= minbasescore) {
                            COVERAGE[referenceposition]++;
                        }
                    }else {
                        COVERAGE[referenceposition]++;
                    }
                }
                catch(Exception e)
                {
                    System.err.println("array size:"+COVERAGE.length);
                    System.err.println("start stop:"+start+"  "+stop);
                    System.err.println(e.getMessage()+" -- "+e.getLocalizedMessage()+" -- "+e.getCause()+" -- "+Arrays.toString(e.getStackTrace()));
                    //throw e;
                    break;
                }
            }
            if(ontarget)
            {
                ontargetReadHitCount++;
            }
            else if(inbuffer)
            {
                inbufferReadHitCount++;
            }
            else
            {
                offtargetReadHitCount++;
            }
	}
	
	/**
	 * makes sure a file exists
	 * @param s
	 * @throws Exception
	 */
	public static void checkFile(String s) throws Exception
	{
            File f = new File(s);
            if(!f.exists())
            {
                throw new FileNotFoundException("No such file as \""+s+"\".  File not found.");
            }
	}
	
	public static void usage()
	{
            String s= "Version: "+VERSION+"\nUsage: -o <output directory & file base name> -t <targetfile> -i <BAM FILE ...> [-d] [-w] [-m <value>] [-b <value>]\n\t* t: target file\n\t* o: output directory and base file name";
            s+="\n\t* d: remove duplicates, and do not use them for statistics\n\t* i: alignment file (multiple files are not allowed)\n\t* w: write whole genome coverage\n\t* m: minimum mapscore (mapscore must be >= to this)\n\t* b: minimum base quality\n";
            System.out.println(s);                                                                                                 	 
	}
	
	/**
	 * removes the "chr" portion of any chromosome name
	 * @param c
	 * @return
	 */
	public static String removechr(String c)
	{
            if(c.length() > 3 && c.substring(0, 3).equalsIgnoreCase("chr")){
                c = c.substring(3);
            }
            return c;
	}
	
	/**
	 * converts fractions into percentages with 2 decimal positions
	 * @param num
	 * @param dom
	 * @return
	 */
	public static double pc(int num, int dom)
	{
            double pc = (double)num/(double)dom;
            pc*=10000.0;pc+=0.5; int ipc = (int)pc; pc = (double)ipc/100;
            return pc;
	}
	
	/**
	 * Writes all the statistical information to an output file.
	 * @param fname output file name
	 * @throws Exception
	 */
	public static void writeReport(String fname) throws Exception
	{
            int nonduplicatereads = totalReadsAligned - duplicateReads;
            if(totalTargetedBases == 0)
            {
                System.err.println("Total targeted bases is zero.  This means that no read has aligned to a chromosome that contains a target. No target matches a chromosome in the BAM, or something else very weird.  Aborting.");
                System.exit(1);
            }
            if(totalReadsAligned == 0)
            {
                System.err.println("No reads aligned. Aborting.");
                System.exit(2);
            }
            if(nonduplicatereads == 0)
            {
                System.err.println("All reads are duplicates. Aborting.");
                System.exit(3);
            }
            if(totalTargets == 0)
            {
                //I don't think we should ever see this error, as its dealt with above.
                System.err.println("No target regions given.  Aborting.");
                System.exit(4);
            }


            int sum =0;
            for(int i = 0; i<coverage_forMedian.length; i++){
                if(sum >= (totalTargetedBases/2)){
                    median_coverage = i;
                    break;
                }else{
                    sum+=coverage_forMedian[i];
                }
            }


            FileWriter report = new FileWriter(fname+".CoverageReport.csv");
            report.write("Version: CapStatsV2.5 2014-05-08\n");
            report.write("BUFFER size:,"+BUFFER+"\n");
            report.write("Read Stats\n");
            report.write("Total Reads Produced:,"+totalReadsProduced+"\n");
            report.write("Total Yield Produced:,"+(READLENGTH * totalReadsProduced)+"\n");
            report.write("Total Unique Yield Produced:,"+(READLENGTH * (totalReadsAligned-duplicateReads))+"\n");
            report.write("Duplicate Reads:,"+duplicateReads+",("+pc(duplicateReads,totalReadsAligned)+"%)\n");
            report.write("Total Reads Aligned:,"+totalReadsAligned+",("+pc(totalReadsAligned,totalReadsProduced)+"%)");
            report.write(",reads paired:,"+totalReadsPaired);
            report.write(",reads paired with mapped mates:,"+totalPairedreadsWithmappedMates+"\n");
            report.write("Aligned Reads On-Buffer:,"+inbufferReadHitCount+",("+pc(inbufferReadHitCount,totalReadsAligned)+"%)\n");
            report.write("Aligned Reads On-Target:,"+ontargetReadHitCount+",("+pc(ontargetReadHitCount,totalReadsAligned)+"%)\n");
            report.write("Average Coverage:,-,("+((int)(totalTargetCoverage/totalTargetedBases))+")\n");
            report.write("Median Coverage:,-,("+median_coverage+")\n");
            int hittot = inbufferReadHitCount+ontargetReadHitCount;
            report.write("Reads that hit target or buffer:,"+hittot+",("+pc(hittot,totalReadsAligned)+"%)\n");
            report.write("Total Aligned Reads (expected):,"+totalReadsAligned+"\n");
            report.write("Total Aligned Reads (calculated):,"+(offtargetReadHitCount+inbufferReadHitCount+ontargetReadHitCount)+"\n");
            report.write("Target Stats\n");
            report.write("Targets Hit:,"+hitTargetCount+",("+pc(hitTargetCount,totalTargets)+"%)\n");
            report.write("Target Buffers Hit:,"+hitTarget_bufferonly_Count+",("+pc(hitTarget_bufferonly_Count,totalTargets)+"%)\n");
            report.write("Total Targets:,"+totalTargets+"\n");
            report.write("Non target regions with high coverage:,"+nonTragetGoodHits+"\n");
            report.write("Base Stats\n");
            report.write("Bases Targeted:,"+totalTargetedBases+"\n");
            report.write("Buffer Bases:,"+totalBufferBases+"\n");
            report.write("Bases with 1+ coverage:,"+basesWithOneHitorMore+",("+pc(basesWithOneHitorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 10+ coverage:,"+basesWith10HitsorMore+",("+pc(basesWith10HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 20+ coverage:,"+basesWith20HitsorMore+",("+pc(basesWith20HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 40+ coverage:,"+basesWith40HitsorMore+",("+pc(basesWith40HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 50+ coverage:,"+basesWith50HitsorMore+",("+pc(basesWith50HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 100+ coverage:,"+basesWith100HitsorMore+",("+pc(basesWith100HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 500+ coverage:,"+basesWith500HitsorMore+",("+pc(basesWith500HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Bases with 1000+ coverage:,"+basesWith1000HitsorMore+",("+pc(basesWith1000HitsorMore,totalTargetedBases)+"%)\n");
            report.write("Duplicate read distribution\n");
            report.write("1,2,3,4,5,6to10,11to20,20plus\n");
            int[] mult = {1,1,2,3,4,5,8,15,25};
            for(int i = 1 ; i < dupHisto.length; i++)
            {
                report.write(pc(dupHisto[i]*mult[i],totalReadsAligned)+"%,");
            }
            report.write("\n");
            report.write("Coverage Histogram (may look weird if target regions overlap...)\n");
            for(int i = 0; i < covHisto.length; i++){
                report.write(i+",");
            }
            report.write("\n");
            for(int i = 0; i < covHisto.length; i++){
                report.write(covHisto[i]+",");
            }
            report.write("\n");

            report.write("Target and region coverage plot\n");
            report.write("Position,5'count,3'count\n");
            for(int i = 20; i <= prime_size; i+=20)
            {
                report.write(i+","+fivePrime[fivePrime.length-(i-1)-1]+","+threePrime[i-1]+"\n");
            }
            report.write("%tar-Pos,count\n");
            for(int i = 0; i < 101; i+=2)
            {
                report.write(i+","+targetCov[i]+"\n");
            }
            report.close();
	}
	
	/**
	 * This method is destructive to the data structure, no further work can be done after this method has ran.  
	 * Works out whether reads are on or off target and how far off target they are
	 *  
	 * @param chromo -- current chromosome
	 * @param COVERAGE -- the coverage of the genome
	 * @param wig -- writes a wig file
	 * @throws Exception
	 */
	public static void findWhereReadsHit(String chromo,  short COVERAGE[], FileWriter wig) throws Exception
	{
            for(int j = 0; j < targetChrs.length; j++)
            {
                if(!removechr(targetChrs[j]).equals(removechr(chromo))) continue;
                int start = targetStarts[j];
                int end = targetStops[j];
                for(int i = start - 500; i < end +500; i++)
                {
                    if(i < 0 || i >= size) {continue;}
                    COVERAGE[i] = 0;
                }
            }
            
            for(int i = 0; i < COVERAGE.length; i++)
            {
                if(COVERAGE[i] > 20)
                {
                    int j = i;
                    nonTragetGoodHits++;
                    while(i < COVERAGE.length && COVERAGE[i] > 0)
                    {
                        i++;
                    }

                    while(j > 0 && COVERAGE[j] > 0)
                    {
                        j--;
                    }
                    wig.write("fixedStep  chrom="+chromo+"  start="+j+"  step=1\n");
                    for(int h = j; h < i; h++)
                    {
                        wig.write(COVERAGE[h]+"\n");
                    }
                }
            }
            wig.flush();
	}
	
	static boolean supertets = false;
	
	/**
	 * Gets the target regions from the target file, and writes over the coverage fasta files, as well as determines many of the coverage statistics.
	 * @param chromo  Current chromosome
	 * @param COVERAGE The array which contains the coverage of every base in the genome
	 * @param covFasta The filewriter for the target-specific coverage
	 * @param missTraget  Write a "wig" format file (good for ucsc) which shows you where all off-target regions with high coverage are
	 * @param wgCoverage A filewriter for the whole genome coverage... if null, this file won't be written
	 * @throws Exception
	 */
	public static void getTargetsAndWriteCoverage(String chromo,  short COVERAGE[], FileWriter covFasta, FileWriter missTraget, FileWriter wgCoverage) throws Exception
        {
            if(wgCoverage != null)
            {
                wgCoverage.write(">"+chromo);
                for(int i = 0; i < COVERAGE.length; i++)
                {
                    if(i%100==0) wgCoverage.write("\n");
                    wgCoverage.write(COVERAGE[i]+" ");
                }
                wgCoverage.write("\n");
            }
            for(int j = 0; j < targetChrs.length; j++)
            {
                if(!removechr(targetChrs[j]).equals(removechr(chromo))) {continue;}
                //totalTargets++;
                int start = targetStarts[j];
                int end = targetStops[j];
                int length = end - start+1;
                boolean collectTargetCov = length > 99 ; 

                if(supertets)
                {
                    System.out.println(targetChrs[j]+" "+start+" "+end);
                }
                if(collectTargetCov)
                {
                    for(int i = 0; i < prime_size; i++)
                    {
                        if((start - i) < 0 || (end+i) >= size){
                                ///System.err.println("The BED Target "+targetChrs[j]+" "+start+" "+end+" is going out of Bound!!!\n");
                                continue;
                        }
                        fivePrime[i]+=COVERAGE[start-i];
                        threePrime[i]+=COVERAGE[end+i];

                        ///fivePrime[i]+=COVERAGE[end-i+300];
                        //threePrime[i]+=COVERAGE[start+i-300];

                    }
                }

                if(supertets)
                {

                    for(int i = 0; i < 500; i++)
                    {
                        if((start-i) < 0) {continue;}
                        System.out.print( (start-i)+" ");
                    }
                    System.out.print("\n");
                    for(int i = 0; i < 500; i++)
                    {
                        if((end+i) >= size) {continue;}
                        System.out.print( (end+i)+" ");
                    }
                    System.out.print("\n");


                    supertets= false;
                }

                boolean targetHit = false;
                short[] pc = new short[101];
                short[] pc2 = new short[101];

                covFasta.write(">"+chromo+" "+start+" "+end+"\n");
                boolean spaceit = false;
                if(end - start > 10000) spaceit = true;
                for(int i = 0; i < length; i++)
                {
                    if((i+start) >= size) {continue;}
                    if(spaceit && i%100 == 0) covFasta.write("\n");
                    short cov = COVERAGE[i+start];
                    if(cov < 0)
                    {
                        cov = Short.MAX_VALUE;
                        System.err.println("Coverage less than 0!!!!!!!\n");
                    }
                    short temp_cov = cov;
                    if(temp_cov >= covHisto.length){
                        temp_cov = (short) (covHisto.length-1);
                    }
                    covHisto[temp_cov]++;
                    totalTargetCoverage+=cov;
                    if(cov > 0)
                    {
                        targetHit=true;
                        basesWithOneHitorMore++;
                    }
                    if(cov > 9){
                            basesWith10HitsorMore++;}
                    if(cov > 19){
                            basesWith20HitsorMore++;}
                    if(cov > 39){
                            basesWith40HitsorMore++;}
                    if(cov > 49){
                            basesWith50HitsorMore++;}
                    if(cov > 99){
                            basesWith100HitsorMore++;}
                    if(cov > 499){
                            basesWith500HitsorMore++;}
                    if(cov > 999){
                            basesWith1000HitsorMore++;}
                    
                    covFasta.write(cov+" ");

                    if(cov < coverage_forMedian.length){
                        coverage_forMedian[cov]++;
                    }else{
                        int[] tmp = new int[coverage_forMedian.length];
                        System.arraycopy(coverage_forMedian, 0, tmp, 0, coverage_forMedian.length);
                        coverage_forMedian = new int[cov+1];
                        System.arraycopy(tmp, 0, coverage_forMedian, 0, tmp.length);
                        coverage_forMedian[cov]++;
                    }

                    if(collectTargetCov)
                    {
                        int pcpos = (int)((double)i/(double)length*100+0.5);
                        pc[pcpos] += cov;
                        pc2[pcpos]++;
                    }
                }
                covFasta.write("\n");
                
                
                for(int index = 0; index < pc.length; index++)
                {
                    if(pc2[index] != 0)
                    {
                        int d = (int) (((double)pc[index]/(double)pc2[index])+0.5);
                        pc[index] = (short) d;	
                    }
                }
                
                for(int i = 0; i < 101; i++)
                {
                    targetCov[i]+=pc[i];
                }
                if(targetHit)
                {
                    hitTargetCount++;
                }
                else
                {
                    missTraget.write(targetChrs[j]+"\t"+targetStarts[j]+"\t"+targetStops[j]+"\n");
                    boolean hit = false;
                    for(int i = start - BUFFER; i < start && !hit; i++)
                    {
                        if(i < 0) {continue;}
                        if(COVERAGE[i] > 0){
                            hit=true;
                        }
                    }
                    for(int i = end; i < end+BUFFER && !hit; i++)
                    {
                        if(i >= size) {continue;}
                        if(COVERAGE[i] > 0){
                            hit=true;
                        }
                    }
                    if(hit){
                        hitTarget_bufferonly_Count++;
                    }
                }
                /*for(int i = 0; i < length; i++){
                    if((i+start) >= size) {continue;}
                    COVERAGE[i+start]=0;
                }*/   //Commented out on June 13th 2014, to include non flattened bedfile.
            }
	}
		
	/**
	 * 
	 * @param chromo the current chromosome to load
	 * @param size the size of the chromosome
	 * @return
	 * @throws Exception
	 */
	public static char[]  getTargetPos(String chromo, int size) throws Exception
	{
            char[] TR = new char[size];
            chromo = removechr(chromo);
            for(int j = 0; j < targetChrs.length; j++)
            {
                try{
                        if(!removechr(targetChrs[j]).equals(chromo))continue;
                        int start = targetStarts[j];
                        int end = targetStops[j];
                        for(int i = start; i <= end; i++)
                        {
                            if(i >= size) {
                                continue;
                            }else{
                                TR[i] = 1;
                            }
                        }
                        for(int i = start -BUFFER; i <start; i++)
                        {
                            if(i < 0) {
                                continue;
                            }else{
                                if(TR[i] == 0){
                                    TR[i] = 2;
                                }
                            }
                        }
                        if(end < (size-1)) {
                            for(int i = end+1; i < end+BUFFER; i++)
                            {
                                if(TR[i] == 0){
                                    TR[i] = 2;
                                }
                            }
                        }else{
                            continue;
                        }
                }
                catch(Exception e)
                {
                    System.err.println("HUGE ERROR TRYING TO PLACE THE TARGET SEQUENCE.  ARE THE TARGETS AND GENOME THE SAME VERSION????");
                    break;
                }
            }		
            for(int i =1; i < TR.length; i++)
            {
                if(TR[i] == 1){
                    totalTargetedBases++;
                }else if(TR[i] == 2){
                    totalBufferBases++;
                }
            }
            return TR;
        }
}
