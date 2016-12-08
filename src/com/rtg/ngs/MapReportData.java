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

package com.rtg.ngs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.rtg.sam.SamUtils;
import com.rtg.util.Environment;
import com.rtg.util.Histogram;
import com.rtg.util.License;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMRecord;

/**
 * Accumulates histogram data from alignments.
 *
 */
public final class MapReportData {

  private static final String CL = "#CL";

  /** Map reporter file format version */
  public static final String VERSION = "map report data v3";

  private static final String VERSION_HEADER = "#Version\t" + Environment.getVersion() + ", " + VERSION + StringUtils.LS;

  private static final String SEPARATOR = "\t";

  private final ArrayList<String> mCommandLines = new ArrayList<>();
  private final HashMap<DistributionType, Histogram> mHistograms = new HashMap<>();

  /**
   * The name of the map report text file
   */
  public static final String MAP_REPORT_FILE_NAME = "reportdata.tsv";

  /**
   * Possible distributions in map report
   */
  public enum DistributionType {

    /** Mapping Counts */
    MAPC("Mapping Counts By Status", License.isDeveloper(), License.isDeveloper()),
    /** Alignment Score - single end and unmated paired end left arm */
    AS("Alignment Score Distribution"),
    /** Alignment Score - unmated paired end right arm */
    AS2("Right Arm Alignment Score"),
    /** Alignment Score - mated paired end left arm */
    ASM("Mated Left Arm Alignment Score"),
    /** Alignment Score - mated paired end right arm */
    AS2M("Mated Right Arm Alignment Score"),
    /** Length of reads */
    RLEN("Read Length Distribution", License.isDeveloper(), License.isDeveloper()),
    /** Length of right arm reads */
    RLEN2("Right Arm Read Length", License.isDeveloper(), License.isDeveloper()),
    /** Length of left arm mated reads */
    RLENM("Mated Left Arm Read Length", License.isDeveloper(), License.isDeveloper()),
    /** Length of right arm mated reads */
    RLEN2M("Mated Right Arm Read Length", License.isDeveloper(), License.isDeveloper()),
    /** unmapped read lengths */
    RLENU("Read Length Distribution (Unmapped Only)", License.isDeveloper(), License.isDeveloper()),
    /** Number of hits */
    NH("Number Of Hits Distribution", License.isDeveloper(), License.isDeveloper()),
    /** Right Arm Number of hits */
    NH2("Right Arm Number Of Hits", License.isDeveloper(), License.isDeveloper()),
    /** Mated Left Arm Number of hits */
    NHM("Mated Left Arm Number Of Hits", License.isDeveloper(), License.isDeveloper()),
    /** Mated Right Arm Number of hits */
    NH2M("Mated Right Arm Number Of Hits", License.isDeveloper(), License.isDeveloper()),
    /** Number of mismatches */
    NM("Number Of Mismatches Distribution", License.isDeveloper(), License.isDeveloper()),
    /** Right Arm Number of mismatches */
    NM2("Right Arm Number Of Mismatches", License.isDeveloper(), License.isDeveloper()),
    /** Mated Left Arm Number of mismatches */
    NMM("Mated Left Arm Number Of Mismatches", License.isDeveloper(), License.isDeveloper()),
    /** Mated Right Arm Number of mismatches */
    NM2M("Mated Right Arm Number Of Mismatches", License.isDeveloper(), License.isDeveloper()),
    /** Mapping orientation */
    ORI("Mapping Orientation", true, License.isDeveloper()),
    /** Right Arm Mapping orientation */
    ORI2("Right Arm Mapping Orientation", true, License.isDeveloper()),
    /** Mated Left Arm Mapping orientation */
    ORIM("Mated Left Arm Mapping Orientation", true, License.isDeveloper()),
    /** Mated Right Arm Mapping orientation */
    ORI2M("Mated Right Arm Mapping Orientation", true, License.isDeveloper()),
    /** Mate orientation */
    MORI("Mate Orientation", License.isDeveloper(), License.isDeveloper()),
    /** Fragment length */
    FLEN("Fragment Length Distribution"),
    /** Mapping Quality (MAPQ) */
    MAPQ("Mapping Quality (MAPQ) Distribution");

    private final String mLongName;
    private final boolean mShowData;
    private final boolean mShowImage;

    /**
     *
     */
    DistributionType(String longName) {
      this(longName, true, true);
    }

    DistributionType(String longName, boolean showData, boolean showImage) {
      mLongName = longName;
      mShowData = showData;
      mShowImage = showImage;
    }

    /**
     * Get left/right mated/unmated version of this.
     * @param isRight want right version
     * @param isMated want mated version
     * @return a {@link DistributionType}
     */
    public DistributionType getType(boolean isRight, boolean isMated) {
      int dt = this.ordinal();
      if (isRight) {
        dt += 1;
      }
      if (isMated) {
        dt += 2;
      }
      return DistributionType.values()[dt];
    }

    /**
     * Returns long name for distribution type
     * @return a name
     */
    public String longName() {
      return mLongName;
    }

    /**
     * @return true if data should be shown for this report
     */
    public boolean showData() {
      return mShowData;
    }

    /**
     * @return true if image should be shown for this report
     */
    public boolean showImage() {
      return mShowImage;
    }
  }

  /**
   * Handles merging of multiple map report data sets.
   */
  public static final class Merger {

    private final List<MapReportData> mReportData = Collections.synchronizedList(new ArrayList<MapReportData>());

    /**
     * Creates a new map report data and adds it to the collective.
     * @return a MapReportData
     */
    public MapReportData createMapReportData() {
      final MapReportData mr = new MapReportData();
      mReportData.add(mr);
      return mr;
    }

    /**
     * Creates a MapReportData from contents of a stream which are the output of another MapReportData.
     * Adds it to the collective.
     * @param in input stream to read from
     * @return new MapReportData
     * @throws IOException if an I/O problem occurs
     */
    public MapReportData createMapReportData(InputStream in) throws IOException {
      final MapReportData mr = createMapReportData();

      // some magic for reading from file...
      final BufferedReader br = new BufferedReader(new InputStreamReader(in));
      String line;
      while ((line = br.readLine()) != null) {
        if (line.startsWith("#")) {
          // hmm any header stuff - like version number...
          if (line.startsWith("#Version") && !line.contains(VERSION)) {
            throw new NoTalkbackSlimException("Unsupported map statistics version: " + line);
          } else if (line.startsWith(CL)) {
            mr.mCommandLines.add(line.split("\t", 2)[1]);
          }
        } else {
          final String[] parts = line.split(SEPARATOR);
          mr.mHistograms.get(DistributionType.valueOf(parts[0])).increment(Integer.parseInt(parts[1]), Long.parseLong(parts[2]));
        }
      }
      return mr;
    }

    /**
     * Creates a MapReportData that has the combined statistics of all the other reporters.
     * @return a MapReportData
     */
    public MapReportData blendReportData() {
      final MapReportData mr = new MapReportData();
      for (MapReportData imr : mReportData) {
        mr.merge(imr);
      }
      return mr;
    }
  }

  /**
   * Constructor for single MapReportData
   */
  public MapReportData() {
    for (DistributionType d : DistributionType.values()) {
      mHistograms.put(d, new Histogram());
    }
  }


  /**
   * Base types of paired end distribution types
   */
  private static final DistributionType[] PAIRED_BASE_TYPES = {DistributionType.AS, DistributionType.RLEN, DistributionType.NH, DistributionType.NM, DistributionType.ORI};

  /**
   * Base types of paired end distribution types
   * @return the base types of the paired end distribution types
   */
  public static DistributionType[] getPairedBaseTypes() {
    return Arrays.copyOf(PAIRED_BASE_TYPES, PAIRED_BASE_TYPES.length);
  }

  /**
   * Check if this map report data contains paired end information
   * @return true if this map report data contains paired end information
   */
  public boolean isPairedEnd() {
    for (DistributionType dt : PAIRED_BASE_TYPES) {
      for (int i = 1; i < 4; ++i) {
        if (mHistograms.get(DistributionType.values()[dt.ordinal() + i]).getLength() > 0) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Get the histogram
   * @param type the type of the histogram to fetch
   * @return the histogram for the given type
   */
  public Histogram getHistogram(DistributionType type) {
    return mHistograms.get(type);
  }

  /**
   * Get the command line(s) associated with this report
   * @return a list of command lines
   */
  public List<String> getCommandLines() {
    return new ArrayList<>(mCommandLines);
  }

  private void increment(DistributionType distribution, int position) {
    mHistograms.get(distribution).increment(position);
  }

  private void increment(DistributionType distribution, Integer position) {
    if (position != null) {
      mHistograms.get(distribution).increment(position);
    }
  }

  /**
   * Extracts statistic counts from the SAM record.
   *
   * @param sam The SAM record for a read.
   */
  public void processRead(SAMRecord sam) {
    final boolean isPaired = sam.getReadPairedFlag();
    final boolean isRight = isPaired && sam.getSecondOfPairFlag();
    final boolean isMated = isPaired && sam.getProperPairFlag();

    // Read length
    final int readLength = sam.getReadLength();

    // Read Orientation
    if (!sam.getReadUnmappedFlag()) {
      // total mapped count is 0
      increment(DistributionType.MAPC, 0);

      // mapped read length
      increment(DistributionType.RLEN.getType(isRight, isMated), readLength);

      // mapping quality value
      increment(DistributionType.MAPQ, sam.getMappingQuality());

      increment(DistributionType.ORI.getType(isRight, isMated), sam.getReadNegativeStrandFlag() ? 1 : 0);

      if (isMated) {
        // total mated count is 1
        increment(DistributionType.MAPC, 1);
        final int insertSize = sam.getInferredInsertSize();
        if (insertSize > 0) {
          // Fragment length - only for +ve side of things
          increment(DistributionType.FLEN, insertSize);

          // Mate orientation, L is leftmost arm on reference, R is rightmost arm on reference
          // 0 => 00 = FF (L then R, both forward)
          // 1 => 01 = FR (L is forward, R is reverse)
          // 2 => 10 = RF (L is reverse, R is forward)
          // 3 => 11 = RR (L then R, both reverse)
          int ori = sam.getReadNegativeStrandFlag() ? 2 : 0;
          ori |= sam.getMateNegativeStrandFlag() ? 1 : 0;
          increment(DistributionType.MORI, ori);
        }
      } else if (isPaired) {
        // total unmated count is 2
        increment(DistributionType.MAPC, 2);
      }
      // Alignment score
      increment(DistributionType.AS.getType(isRight, isMated), sam.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE));

      // Number of hits
      increment(DistributionType.NH.getType(isRight, isMated), sam.getIntegerAttribute(SamUtils.ATTRIBUTE_NH));

      // Number of mismatches
      increment(DistributionType.NM.getType(isRight, isMated), sam.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES));
    } else {
      // total unmapped count is 3
      increment(DistributionType.MAPC, 3);
      // unmapped read length
      increment(DistributionType.RLENU, readLength);
    }
  }

  /**
   * Merge the contents of another MapReportData.
   * @param other the other MapReportData
   */
  void merge(MapReportData other) {
    mCommandLines.addAll(other.mCommandLines);
    for (DistributionType dist : DistributionType.values()) {
      mHistograms.get(dist).addHistogram(other.mHistograms.get(dist));
    }
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    try {
      writeHistograms(sb, false);
    } catch (IOException ioe) {
      // unpossible?
    }
    return sb.toString();
  }

  private void writeHeader(Appendable out) throws IOException {
    out.append(VERSION_HEADER);
    if (mCommandLines.size() == 0) {
      // only add command line if we are writing the file
      if (CommandLine.getCommandLine() != null) {
        mCommandLines.add(CommandLine.getCommandLine());
      } else {
        mCommandLines.add("Internal");
      }
    }

    for (String cl : mCommandLines) {
      out.append(CL)
        .append("\t")
        .append(cl)
        .append(StringUtils.LS);
    }
    out.append("#Hist")
      .append(SEPARATOR)
      .append("Bin")
      .append(SEPARATOR)
      .append("Count")
      .append(StringUtils.LS);
  }

  private void writeHistograms(Appendable out, boolean addSectionHeader) throws IOException {
    for (DistributionType dist : DistributionType.values()) {
      final Histogram histogram = mHistograms.get(dist);
      if (addSectionHeader) {
        out.append("#").append(dist.longName()).append(StringUtils.LS);
      }
      for (int i = 0; i < histogram.getLength(); ++i) {
        final long value = histogram.getValue(i);
        if (value > 0) {
          out.append(dist.toString())
          .append(SEPARATOR)
          .append(Integer.toString(i))
          .append(SEPARATOR)
          .append(Long.toString(value))
          .append(StringUtils.LS);
        }
      }
    }
  }


  /**
   * Writes data to a stream.
   * @param out a print stream to write to
   * @throws IOException if an I/O problem occurs
   */
  public void write(PrintStream out) throws IOException {
    writeHeader(out);
    //writeSummary(out);
    writeHistograms(out, true);
  }

  /**
   * Writes data to a file.
   * @param outFile the file to write to
   * @throws IOException if an I/O problem occurs
   */
  public void write(File outFile) throws IOException {
    try (PrintStream ps = new PrintStream(outFile)) {
      write(ps);
    }
  }

}
