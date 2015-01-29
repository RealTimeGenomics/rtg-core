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
package com.rtg.reader;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.MissingResourceException;
import java.util.Queue;
import java.util.Set;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.Protein;
import com.rtg.mode.SequenceType;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;
import com.rtg.sam.SamUtils;
import com.rtg.taxonomy.SequenceToTaxonIds;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Prints out statistics about a given preread.
 *
 */
public final class SdfStatistics extends AbstractCli {

  private static final String MODULE_NAME = "sdfstats";
  private static final String SDF_FLAG = "SDF";
  private static final String NS_FLAG = "unknowns";
  private static final String POSITIONS_FLAG = "position";
  private static final String QS_FLAG = "quality";
  private static final String SEX_FLAG = "sex";
  private static final String TAXONOMY_FLAG = "taxonomy";
  private static final String NAMES_AND_LENGTHS_FLAG = "lengths";
  private static final String SAM_FLAG = "Xsam-header";

  //private static final String SUMMARY_FLAG = "summary";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * construct a flags object
   * @param flags flags object to initialise
   */
  public void initFlags(final CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Print statistics that describe a directory of SDF formatted data.");
    CommonFlagCategories.setCategories(flags);

    try {
      final Flag inFlag = flags.registerRequired(File.class, SDF_FLAG, "SDF directories");
      inFlag.setMinCount(1);
      inFlag.setMaxCount(Integer.MAX_VALUE);
      inFlag.setCategory(INPUT_OUTPUT);

      flags.registerOptional('n', NS_FLAG, "display info about unknown bases (Ns)").setCategory(REPORTING);
      flags.registerOptional('p', POSITIONS_FLAG, "only display info about unknown bases (Ns) by read position").setCategory(REPORTING);

      flags.registerOptional('q', QS_FLAG, "display mean of quality").setCategory(REPORTING);
      flags.registerOptional(SEX_FLAG, Sex.class, "sex", "display reference sequence list for the given sex, if defined").setCategory(REPORTING).setMaxCount(Integer.MAX_VALUE);
      flags.registerOptional(TAXONOMY_FLAG, "display information about taxonomy").setCategory(REPORTING);

      flags.registerOptional(NAMES_AND_LENGTHS_FLAG, "print out the name and length of each sequence. (Not recommended for read sets)").setCategory(REPORTING);
      flags.registerOptional(SAM_FLAG, "print out a SAM format header corresponding to this SDF").setCategory(REPORTING);
      //flags.setDescription("prints statistics about sequences");

      flags.setValidator(VALIDATOR);

    } catch (final MissingResourceException e) {
      throw new SlimException(e);
    }
  }

  private static final Validator VALIDATOR = new Validator() {
      @Override
      public boolean isValid(final CFlags flags) {

        final Collection<Object> values = flags.getAnonymousValues(0);
        for (final Object o : values) {
          final File f = (File) o;
          if (!CommonFlags.validateSDF(f)) {
            return false;
          }
          if (flags.isSet(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, f)) {
            return false;
          }
          if (flags.isSet(TAXONOMY_FLAG)) {
            if (!new File(f, Taxonomy.TAXONOMY_FILE).isFile()) {
              flags.error("--" + TAXONOMY_FLAG + " was specified but " + f + " is missing a '" + Taxonomy.TAXONOMY_FILE + "'");
              return false;
            }
            if (!new File(f, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE).isFile()) {
              flags.error("--" + TAXONOMY_FLAG + " was specified but " + f + " is missing a '" + SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE + "'");
              return false;
            }
          }
        }
        return true;
      }
    };

  /**
   * Gathers and prints statistics about a group for slim data files.
   *
   * @param dir SDF to print statistics for
   * @param out where to write statistics to
   * @param nStats whether to produce statistics on N occurrences
   * @param pStats whether to produce only statistics regarding position based N occurrences.
   * @param quality whether to output the quality score per read histogram
   * @exception IOException if an error occurs.
   */
  public void performStatistics(File dir, final PrintStream out, final boolean nStats, final boolean pStats, final boolean quality) throws IOException {
    final AnnotatedSequencesReader reader;
    try {
      reader = SequencesReaderFactory.createDefaultSequencesReader(dir);
    } catch (final FileNotFoundException e) {
      if (dir.isDirectory()) {
        throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, dir.toString());
      } else if (dir.exists()) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified file, \"" + dir.getPath() + "\", is not an SDF.");
      } else {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified SDF, \"" + dir.getPath() + "\", does not exist.");
      }
    }
    try (AnnotatedSequencesReader sr = reader) {
      performStatistics(sr, dir, out, nStats, pStats, quality);
      printReadMe(sr, out);
    }
  }

  /**
   * Gather statistics about slim a slim data file and write to an output stream
   * @param reader sequences reader to print statistics for
   * @param dir directory for the sequences reader
   * @param out where to write statistics to
   * @param nStats whether to produce statistics on N occurrences
   * @param pStats whether to produce only statistics regarding position based N occurrences.
   * @param quality whether to output the quality score per read histogram
   * @exception IOException if an error occurs.
   */
  public static void performStatistics(AnnotatedSequencesReader reader, File dir, final PrintStream out, final boolean nStats, final boolean pStats, final boolean quality) throws IOException {
    out.append("Location           : ");
    out.append(dir.getAbsolutePath());
    out.append(StringUtils.LS);
    if (reader.commandLine() != null) {
      out.append("Parameters         : ");
      out.append(reader.commandLine());
      out.append(StringUtils.LS);
    }
    if (reader.comment() != null) {
      out.append("Comment            : ");
      out.append(reader.comment());
      out.append(StringUtils.LS);
    }
    if (reader.samReadGroup() != null) {
      out.append("SAM read group     : ");
      out.append(reader.samReadGroup());
      out.append(StringUtils.LS);
    }
    out.append("SDF Version        : ");
    out.append(Long.toString(reader.sdfVersion()));
    out.append(StringUtils.LS);

    if (!pStats) {
      out.append("Type               : ");
      out.append(reader.type().toString());
      out.append(StringUtils.LS);
      out.append("Source             : ");
      out.append(reader.getPrereadType().toString());
      out.append(StringUtils.LS);
      out.append("Paired arm         : ");
      out.append(reader.getArm().toString());
      out.append(StringUtils.LS);
      final SdfId sdfId = reader.getSdfId();
      if (sdfId.available()) {
        out.append("SDF-ID             : ");
        out.append(sdfId.toString());
        out.append(StringUtils.LS);
      }
      out.append("Number of sequences: ");
      out.append(Long.toString(reader.numberSequences()));
      out.append(StringUtils.LS);
      final long max = reader.maxLength();
      final long min = reader.minLength();
      if (max >= min) {
        out.append("Maximum length     : ");
        out.append(Long.toString(max));
        out.append(StringUtils.LS);
        out.append("Minimum length     : ");
        out.append(Long.toString(min));
        out.append(StringUtils.LS);
      }
      out.append("Sequence names     : ");
      if (reader.hasNames()) {
        out.append("yes");
      } else {
        out.append("no");
      }
      out.append(StringUtils.LS);
      final long[] counts = reader.residueCounts();
      long sum = 0;
      for (int i = 0 ; i < counts.length; i++) {
        out.append((reader.type() == SequenceType.DNA) ? DNA.values()[i].toString() : Protein.values()[i].toString());
        out.append("                  : ");
        out.append(Long.toString(counts[i]));
        out.append(StringUtils.LS);
        sum += counts[i];
      }
      out.append("Total residues     : ");
      out.append(Long.toString(sum));
      out.append(StringUtils.LS);
      if (quality) {
        printQualityHistogram(reader, out);
      }
      out.append("Residue qualities  : ");
      if (reader.hasQualityData() && reader.hasHistogram()) {
        out.append("yes");
      } else {
        out.append("no");
      }
      out.append(StringUtils.LS);
      if (nStats) {
        printNBlocks(reader, out);
      }
      out.append(StringUtils.LS);
    } else {
      printPositionBlock(reader, out);
    }
  }

  static void printSAMHeader(File dir, final Appendable out) throws IOException {
    out.append("SAM Header:").append(StringUtils.LS);
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
      final SAMFileHeader header = new SAMFileHeader();
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      final SAMProgramRecord pg = new SAMProgramRecord(Constants.APPLICATION_NAME);
      pg.setProgramVersion(Environment.getVersion());
      final int[] lengths = reader.sequenceLengths(0, reader.numberSequences());
      for (int i = 0; i < lengths.length; i++) {
        final SAMSequenceRecord record = new SAMSequenceRecord(reader.name(i), lengths[i]);
        header.addSequence(record);
      }
      if (reader.getSdfId().available()) {
        header.addComment(SamUtils.TEMPLATE_SDF_ATTRIBUTE + reader.getSdfId().toString());
      }
      final ByteArrayOutputStream bo = new ByteArrayOutputStream();
      new SAMFileWriterFactory().makeSAMWriter(header, true, bo).close();
      out.append(bo.toString());
    }
  }


  static void printReferenceSequences(File dir, Sex sex, final Appendable out) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
      final ReferenceGenome rg = new ReferenceGenome(reader, sex);
      out.append("Sequences for sex=").append(String.valueOf(sex)).append(":").append(StringUtils.LS);
      out.append(rg.toString()).append(StringUtils.LS);
      out.append(StringUtils.LS);
    }
    out.append(StringUtils.LS);
  }

  private static void printPositionBlock(final SequencesReader reader, final Appendable out) throws IOException {
    if (reader.hasHistogram()) {
      final long[] posHisto = reader.posHistogram();
      out.append("Histogram of N position frequencies");
      out.append(StringUtils.LS);
      out.append(printHistogram(posHisto, true));
    } else {
      out.append("Histogram of N position frequencies is not available for this SDF").append(StringUtils.LS);
    }
  }

  private static void printQualityHistogram(final SequencesReader reader, final Appendable out) throws IOException {

    if (reader.hasQualityData() && reader.hasHistogram()) {
      final double globalQualityAverage = reader.globalQualityAverage();
      out.append("Average quality    : ").append(Utils.realFormat(MathUtils.phred(globalQualityAverage), 1));
      out.append(StringUtils.LS);
      out.append("Average qual / pos : ");
      out.append(StringUtils.LS);
      out.append(printQualityHistogram(reader.positionQualityAverage()));
    } else {
      out.append("Quality statistics are not available for this SDF").append(StringUtils.LS);
    }
  }

  private static String printQualityHistogram(final double[] histo) {
    final StringBuilder str = new StringBuilder();
    // find the last valid non-zero value
    int lastIndex = -1;
    for (int i = histo.length - 1; i >= 0; i--) {
      if (histo[i] != 0) {
        lastIndex = i;
        break;
      }
    }

    final long addition = 1;

    for (int i = 0; i <= lastIndex; i++) {
      //if (histo[i] > 0) {
      final String number;
      if ((long) i == SdfWriter.MAX_HISTOGRAM) {
        number = Long.toString((long) i + addition) + "+";
      } else {
        number = Long.toString((long) i + addition);
      }
      for (int j = 0; j < 18 - number.length(); j++) {
        str.append(" ");
      }
      str.append(number);
      str.append(" : ");
      str.append(Utils.realFormat(MathUtils.phred(histo[i]), 1));
      str.append(StringUtils.LS);
      //}
    }
    return str.toString();
  }

  private static void printNBlocks(final SequencesReader reader, final Appendable out) throws IOException {
    final long nBlocks = reader.nBlockCount();
    final long longestBlock = reader.longestNBlock();

    if (reader.hasHistogram()) {
      out.append("Blocks of Ns       : ");
      out.append(Long.toString(nBlocks));
      out.append(StringUtils.LS);
      out.append("Longest block of Ns: ");
      out.append(Long.toString(longestBlock));
      out.append(StringUtils.LS);

      out.append("Histogram of N frequencies");
      out.append(StringUtils.LS);
      final long[] histo = reader.histogram();
      out.append(printHistogram(histo, false));
    } else {
      out.append("N counts are not available on this SDF").append(StringUtils.LS);
    }
  }

  static void printSequenceNameAndLength(File dir, PrintStream out) throws IOException {
    out.println("Sequence lengths: ");
    try (SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
      while (sr.nextSequence()) {
        if (sr.hasNames()) {
          out.println(sr.currentName() + "\t" + sr.currentLength());
        } else {
          out.println("sequence_" + sr.currentSequenceId() + "\t" + sr.currentLength());
        }
      }
    }
    out.println();
  }

  static void printReadMe(SequencesReader reader, PrintStream out) throws IOException {
    final String readme = reader.getReadMe();
    if (readme != null) {
      out.println("Additional Info:");
      out.println(readme);
      out.println();
    }
  }

  static void printTaxonomyStatistics(File dir, PrintStream out) throws IOException {
    final Taxonomy tax = new Taxonomy();
    try (FileInputStream fis = new FileInputStream(new File(dir, Taxonomy.TAXONOMY_FILE))) {
      tax.read(fis);
    }
    out.append("Taxonomy nodes     : ");
    out.append(Integer.toString(tax.size()));
    out.append(StringUtils.LS);

    final File sequenceTaxonFile = new File(dir, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE);
    final Map<String, Integer> sequenceToId = SequenceToTaxonIds.sequenceToIds(sequenceTaxonFile);
    final Set<Integer> uniqTaxIds = new HashSet<>();
    uniqTaxIds.addAll(sequenceToId.values());
    out.append("Sequence nodes     : ");
    out.append(Integer.toString(uniqTaxIds.size()));
    out.append(StringUtils.LS);
    out.append("Other nodes        : ");
    out.append(Integer.toString(tax.size() - uniqTaxIds.size()));
    out.append(StringUtils.LS);
  }


  private static String printHistogram(final long[] histo, final boolean addOne) {
    final StringBuilder str = new StringBuilder();
    // find the last valid non-zero value
    int lastIndex = -1;
    for (int i = histo.length - 1; i >= 0; i--) {
      if (histo[i] != 0) {
        lastIndex = i;
        break;
      }
    }

    final long addition = addOne ? 1 : 0;

    for (int i = 0; i <= lastIndex; i++) {
      //if (histo[i] > 0) {
      final String number;
      if ((long) i == SdfWriter.MAX_HISTOGRAM) {
        number = Long.toString((long) i + addition) + "+";
      } else {
        number = Long.toString((long) i + addition);
      }
      for (int j = 0; j < 18 - number.length(); j++) {
        str.append(" ");
      }
      str.append(number);
      str.append(" : ");
      str.append(Long.toString(histo[i]));
      str.append(StringUtils.LS);
      //}
    }
    return str.toString();
  }

  /**
   * Main program for building and searching. Use -h to get help.
   * @param args command line arguments.
   */
  public static void main(final String[] args) {
    new SdfStatistics().mainExit(args);
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final Collection<Object> values = mFlags.getAnonymousValues(0);

    final Queue<File> dirs = new LinkedList<>();
    for (final Object o : values) {
      final File f = (File) o;
      if (ReaderUtils.isPairedEndDirectory(f)) {
        dirs.add(ReaderUtils.getLeftEnd(f));
        dirs.add(ReaderUtils.getRightEnd(f));
      } else {
        dirs.add(f);
      }
    }
    if (dirs.size() == 0) {
      throw new NoTalkbackSlimException(ErrorType.NO_VALID_INPUTS);
    }
    final PrintStream ps = new PrintStream(out);
    try {
      for (File dir : dirs) {
        performStatistics(dir, ps, mFlags.isSet(NS_FLAG), mFlags.isSet(POSITIONS_FLAG), mFlags.isSet(QS_FLAG));
        if (mFlags.isSet(NAMES_AND_LENGTHS_FLAG)) {
          printSequenceNameAndLength(dir, ps);
        }
        if (mFlags.isSet(SEX_FLAG)) {
          for (Object o : mFlags.getValues(SEX_FLAG)) {
            final Sex s = (Sex) o;
            printReferenceSequences(dir, s, ps);
          }
        }
        if (mFlags.isSet(TAXONOMY_FLAG)) {
          printTaxonomyStatistics(dir, ps);
        }
        if (mFlags.isSet(SAM_FLAG)) {
          printSAMHeader(dir, ps);
        }
      }
    } finally {
      ps.flush();
    }
    return 0;
    // don't delete logs on exceptions...
  }

}
