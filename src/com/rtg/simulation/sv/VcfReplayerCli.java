/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.simulation.sv;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.PrereadNames;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.vcf.Adjacency;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Replays VCF file on given genomes.
 * Currently only handles simple homozygous breakpoints with known mates.
 * Eg. does not handle alt Calls containing '.' (unmated breakpoint) or ',' (heterozygous) yet.
 *
 */
public class VcfReplayerCli extends AbstractCli {

  private static final String REPLAY = "replay";

  /** Set this to true to print out debug messages */
  public static final boolean VERBOSE = false;

  /**
   * main method
   * @param args input arguments
   */
  public static void main(String[] args) {
    new VcfReplayerCli().mainExit(args);
  }

  void process(CFlags flags) throws IOException {
    final File template = (File) flags.getValue(CommonFlags.TEMPLATE_FLAG);
    final File output = (File) flags.getValue(CommonFlags.OUTPUT_FLAG);
    final File replayVcf = (File) flags.getValue(REPLAY);
    // we use a separate set for each chromosome.
    final Map<String, TreeSet<Adjacency>> setA = new HashMap<>();
    final Map<String, TreeSet<Adjacency>> setB = new HashMap<>();

    boolean gtWarned = false;
    try (VcfReader reader = VcfReader.openVcfReader(replayVcf)) {
      final int[] gts = new int[2];
      while (reader.hasNext()) {
        final VcfRecord rec = reader.next();
        final List<String> strAlts = rec.getAltCalls();
        final Adjacency[] alts = new Adjacency[strAlts.size()];
        for (int i = 0; i < strAlts.size(); ++i) {
          alts[i] = Adjacency.parseAdjacency(rec.getSequenceName(), rec.getOneBasedStart(), strAlts.get(i));
        }
        if (!rec.hasFormat(VcfUtils.FORMAT_GENOTYPE)) {
          if (!gtWarned) {
            Diagnostic.warning("Records without GT field assumed to be homozygous.");
            gtWarned = true;
          }
          gts[0] = 1;
          gts[1] = 1;
        } else {
          final List<String> gtList = rec.getFormat(VcfUtils.FORMAT_GENOTYPE);
          if (gtList.size() > 1) {
            throw new NoTalkbackSlimException("Replaying multi sample VCF files is unsupported");
          }
          final String[] gtSplit = gtList.get(0).split("/");
          for (int i = 0; i < gtSplit.length; ++i) {
            //if its . just go with first thing in the alt field
            gts[i] = gtSplit[i].equals(".") ? 1 : Integer.parseInt(gtSplit[i]);
          }
          if (gtSplit.length == 1) {
            gts[1] = gts[0];
          }
        }
        if (gts[0] > 0) {
          final Adjacency adjA = alts[gts[0] - 1];
          if (adjA != null) {
            final String key = adjA.thisChromosome();
            if (setA.get(key) == null) {
              setA.put(key, new TreeSet<Adjacency>());
            }
            setA.get(key).add(adjA);
          }
        }
        if (gts[1] > 0) {
          final Adjacency adjB = alts[gts[1] - 1];
          if (adjB != null) {
            final String key = adjB.thisChromosome();
            if (setB.get(key) == null) {
              setB.put(key, new TreeSet<Adjacency>());
            }
            setB.get(key).add(adjB);
          }
        }
      }
    }
    // TODO: Check we got complementary pairs of Adjacency objects
    // TODO: check that every named chromosome is one of the chromosomes in our reference.
    SdfUtils.validateHasNames(template);
    try (SequencesReader reader = SequencesReaderFactory.createMemorySequencesReader(template, true, LongRange.NONE)) {
      SdfUtils.validateNoDuplicates(reader, false);

      try (SdfWriter w1 = new SdfWriter(new File(output, "left"), com.rtg.util.Constants.MAX_FILE_SIZE, reader.getPrereadType(), false, true, true, reader.type());
           SdfWriter w2 = new SdfWriter(new File(output, "right"), com.rtg.util.Constants.MAX_FILE_SIZE, reader.getPrereadType(), false, true, true, reader.type())) {
        doWork(setA, reader, w1);
        doWork(setB, reader, w2);
      }
    }
  }

  private void doWork(Map<String, TreeSet<Adjacency>> set, SequencesReader reader, SdfWriter w1) throws IOException {
    final ArrayList<String> names = loadNames(reader);
    final Set<String> doneEnds = new HashSet<>();
    for (int chrNum = 0; chrNum < names.size() ; ++chrNum) {
      final String chrName = names.get(chrNum);
      generateChromosome(chrNum, chrName, true, set, names, reader, doneEnds, w1);
    }
    // now do any remaining reference sequences, starting from the end.
    for (int chrNum = 0; chrNum < names.size() ; ++chrNum) {
      final String chrName = names.get(chrNum);
      generateChromosome(chrNum, chrName, false, set, names, reader, doneEnds, w1);
    }

    // TODO: support Heterozygous Left/Right, ie. generate w2 as well.
    // (exactly like above, starting with a newly empty doneEnds set.)
    doneEnds.clear();
    // We may be able to do this by marking all the Adjacency objects we have 'used'
    // during the left pass, and removing them from the tree sets before the second pass?
    // But we need to change the TreeSets so they do not throw away duplicates (the homozygous case).
    // In fact, we would need to *add* duplicate Adjacency objects for the homozygous case.
    // We may be able to do this by deleting each Adjacency object as it is used?

//    for (int chrNum = 0; chrNum < names.size() ; ++chrNum) {
//      final String chrName = names.get(chrNum);
//      generateChromosome(chrNum, chrName, true, set, names, reader, doneEnds, w2);
//    }
//    // now do any remaining reference sequences, starting from the end.
//    for (int chrNum = 0; chrNum < names.size() ; ++chrNum) {
//      final String chrName = names.get(chrNum);
//      generateChromosome(chrNum, chrName, false, set, names, reader, doneEnds, w2);
//    }
  }

  
  protected void generateChromosome(int chrNum, String chromoName, boolean fwd, Map<String, TreeSet<Adjacency>> set, ArrayList<String> nameToNum, SequencesReader reader, final Set<String> doneEnds, SdfWriter writer) throws IOException {
    String chrName = chromoName;
    int currentChr = chrNum;
    boolean currentForward = fwd;
    int currentPos = fwd ? 1 : reader.length(chrNum);
    if (doneEnds.contains(chrName + ":" + currentPos)) {
      return;
    }
    if (VERBOSE) {
      System.out.println("\nSTARTING " + chrNum + " = " + chrName + ":" + currentPos);
    }
    writer.startSequence(chrName);
    while (true) {
      final byte[] chromosome = reader.read(currentChr);
      final int searchFrom = currentPos + (currentForward ? 1 : -1);
      final Adjacency adj = new Adjacency(chrName, searchFrom, currentForward);
      if (VERBOSE) {
        System.out.println("   looking " + (currentForward ? "fwd" : "bkwd") + " from " + chrName + ":" + searchFrom + " " + currentForward);
      }
      final TreeSet<Adjacency> chrSet = set.get(chrName);
      final Adjacency next = chrSet == null ? null
          : (currentForward ? chrSet.higher(adj) : chrSet.lower(adj));
      if (next != null && currentForward == next.thisIsForward()) {
        // TODO: mark this Adjacency object for removal after this pass, since we have used it.
        write(writer, chromosome, currentPos, next.thisEndPos() + (currentForward ? -1 : +1));
        if (VERBOSE) {
          System.out.println("jumping to breakpoint: " + next.mateChromosome() + ":" + next.mateStartPos() + next.mateIsForward());
        }
        final String bases = next.thisIsForward() ? next.thisBases() : DnaUtils.reverseComplement(next.thisBases());
        write(writer, DnaUtils.encodeString(bases), 1, next.thisBases().length());
        chrName = next.mateChromosome();
        currentPos = next.mateStartPos();
        currentForward = next.mateIsForward();
        currentChr = nameToNum.indexOf(chrName);
      } else {
        // no breakpoints to follow, so we stop at the end/start of this region/chromosome
        final int finish;
        if (next != null) {
          finish = next.thisEndPos() + (currentForward ? -1 : +1);
        } else {
          finish = currentForward ? chromosome.length : 1;
        }
        write(writer, chromosome, currentPos, finish);
        if (VERBOSE) {
          System.out.println("finished at " + chrName + ":" + finish);
        }
        doneEnds.add(chrName + ":" + finish);
        break;
      }
    }
    writer.endSequence();
  }

  /**
   * Writes the given region of <code>chromosome</code> into the SDF writer.
   * If <code>end &lt; start</code>, then the region is reversed.
   *
   * @param w writer
   * @param chromosome chromosome bytes
   * @param start 1 based start position inclusive
   * @param end 1 based end position inclusive
   * @throws IOException if error
   */
  private void write(SdfWriter w, byte[] chromosome, int start, int end) throws IOException {
    final byte[] b = new byte[Math.abs(end - start) + 1];
    if (start <= end) {
      //System.out.println("copy chr[" + (start - 1) + "..] into b.  len=" + b.length);
      System.arraycopy(chromosome, start - 1, b, 0, b.length);
    } else {
      // copy and reverse
      //System.out.println("copy rev chr[" + (start - 1 - (b.length-1)) + ".." + (start - 1) + "] into b.  len=" + b.length);
      for (int i = 0; i < b.length; ++i) {
        b[i] = DNA.complement(chromosome[start - 1 - i]);
      }
    }
    w.write(b, null, b.length);
  }

  private ArrayList<String> loadNames(SequencesReader reader) throws IOException {
    final PrereadNames names = new PrereadNames(reader.path(), LongRange.NONE);
    final ArrayList<String> list = new ArrayList<>();
    for (int i = 0; i < reader.numberSequences(); ++i) {
      list.add(names.name(i));
    }
    return list;
  }

  static void getFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlags.initOutputDirFlag(flags);
    flags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, "SDF", "template SDF");
    flags.registerRequired('v', REPLAY, File.class, "FILE", "VCF file to replay");
    flags.setValidator(new Validator() {

      @Override
      public boolean isValid(CFlags flags) {
        if (!CommonFlags.validateOutputDirectory(flags)
          || !CommonFlags.validateInputFile(flags, REPLAY)
          || !CommonFlags.validateTemplate(flags)) {
          return false;
        }
        return true;
      }
    });
  }

  @Override
  protected void initFlags() {
    getFlags(mFlags);
  }


  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    process(mFlags);
    return 0;
  }

  @Override
  public String moduleName() {
    return "vcfsvreplay";
  }

  @Override
  public String description() {
    return null;
  }

}
