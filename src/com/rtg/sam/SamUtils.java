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
package com.rtg.sam;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.alignment.CgGotohEditDistance;
import com.rtg.reader.FastaUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.RuntimeIOException;

/**
 *
 */
public final class SamUtils {

  /** Length of CG raw reads */
  public static final int CG_RAW_READ_LENGTH = 35;

  /** The position of the overlap in a CG read (left arm). */
  public static final int CG_OVERLAP_POSITION = 5;


  //Sam spec defined attributes

  /** Alignment score attribute name */
  public static final String ATTRIBUTE_ALIGNMENT_SCORE = "AS";

  /** Number of mismatches attribute name */
  public static final String ATTRIBUTE_NUM_MISMATCHES = "NM";

  /** The number of stored alignments in the SAM file that match the query in the current record. */
  public static final String ATTRIBUTE_IH = "IH";

  /** The number of reported alignments that match the query in the current record. */
  public static final String ATTRIBUTE_NH = "NH";

  /** Mismatching positions attribute name */
  public static final String ATTRIBUTE_MISMATCH_POSITIONS = "MD";


  // Deprecated in current Sam spec, retained for backwards compatibility

  /** CG: The name of the attribute storing the bases in the overlap region */
  public static final String ATTRIBUTE_CG_OVERLAP_BASES = "GS";

  /** CG: The name of the attribute storing the quality information in the overlap region */
  public static final String ATTRIBUTE_CG_OVERLAP_QUALITY = "GQ";

  /** CG: The name of the attribute storing the instructions on how to recreate the original read */
  public static final String ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS = "GC";


  //BELOW are our own attributes.

  /** Split Read */
  public static final String ATTRIBUTE_SPLIT_READ = "XX";

  /** Combo score */
  public static final String ATTRIBUTE_COMBO_SCORE = "XA";

  /** The unmapped etymology. */
  public static final String ATTRIBUTE_UNMAPPED_ETYMOLOGY = "XC";

  /** This SAM attribute is for the CG Gotoh "super cigar", to handle CG overlaps. */
  public static final String CG_SUPER_CIGAR = "XU";

  /**
   * This SAM attribute stores just the differences from the reference to the read.
   * It must always appear together with the super cigar attribute.
   */
  public static final String CG_READ_DELTA = "XR";

  /**
   * This SAM attribute stores just the extra quality bytes for the part of
   * the CG overlap region that was omitted from the flattened read.
   */
  public static final String CG_OVERLAP_QUALITY = "XQ";

  /** This SAM attribute stores the alignment score of the mate. */
  public static final String ATTRIBUTE_MATE_ALIGNMENT_SCORE = "XM";

  /** This SAM attribute stores the end position of the mate. */
  public static final String ATTRIBUTE_MATE_END = "XN";

  /** Used to show if read mapping accuracy determined a record to be correct or incorrect */
  public static final String ATTRIBUTE_READ_ACCURACY_STATUS = "XE";


  //Attributes defined by other tools

  /** BWA: Type: Unique/Repeat/N/Mate-sw */
  public static final String ATTRIBUTE_BWA_TYPE = "XT";

  /** BWA: Number of best hits */
  public static final String ATTRIBUTE_BWA_NUM_BEST_HITS = "X0";


  /** CIGAR value for SAME */
  public static final char CIGAR_SAME = '=';
  /** CIGAR value for MISMATCH */
  public static final char CIGAR_MISMATCH = 'X';
  /** CIGAR value for SAME OR MISMATCH */
  public static final char CIGAR_SAME_OR_MISMATCH = 'M';
  /** CIGAR value for DELETION FROM REFERENCE */
  public static final char CIGAR_DELETION_FROM_REF = 'D';
  /** CIGAR value for INSERTION INTO REFERENCE */
  public static final char CIGAR_INSERTION_INTO_REF = 'I';
  /** CIGAR value for GAP IN READ */
  public static final char CIGAR_GAP_IN_READ = 'N';
  /** CIGAR value for SOFT CLIP */
  public static final char CIGAR_SOFT_CLIP = 'S';
  /** CIGAR value for HARD CLIP */
  public static final char CIGAR_HARD_CLIP = 'H';
  /** Pseudo-CIGAR value for unmapped */
  public static final char CIGAR_UNMAPPED = 'u';

  private static final char[] CIGAR_CODES = {CIGAR_SAME_OR_MISMATCH, CIGAR_INSERTION_INTO_REF, CIGAR_DELETION_FROM_REF, CIGAR_GAP_IN_READ, CIGAR_SOFT_CLIP, 'H', 'P', CIGAR_SAME, CIGAR_MISMATCH};

  /** The attribute used in a comment line to indicate the SDF ID of the reads */
  public static final String READ_SDF_ATTRIBUTE = "READ-SDF-ID:";

  /** The attribute used in a comment line to indicate the SDF ID of the reference */
  public static final String TEMPLATE_SDF_ATTRIBUTE = "TEMPLATE-SDF-ID:";

  /** The attribute used in a comment line to indicate the gender used during mapping */
  public static final String GENDER_ATTRIBUTE = "MAPPING-GENDER:";

  /** The attribute used in a comment line to indicate the ID of this run */
  public static final String RUN_ID_ATTRIBUTE = "RUN-ID:";

  /** The filename extension used for SAM files */
  public static final String SAM_SUFFIX = ".sam";

  /** The filename extension used for BAM files */
  public static final String BAM_SUFFIX = ".bam";

  /** The filename extension used for BAM index files */
  public static final String BAI_SUFFIX = BamIndexer.BAM_INDEX_EXTENSION;

//  private static final boolean IGNORE_HEADER_INCOMPATIBILITY = GlobalFlags.isSet(GlobalFlags.SAM_IGNORE_INCOMPATIBLE_HEADERS_FLAG);

  private SamUtils() {
  }

  /**
   * Cats multiple gzipped SAM files together
   * @param destination destination for cat
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void samCat(final OutputStream destination, final File... inputFiles) throws IOException {
    samCat(inputFiles.length > 0 && FileUtils.isGzipFilename(inputFiles[0]), destination, false, inputFiles);
  }

  /**
   * Cats multiple gzipped SAM files together
   * @param inputIsGzipped true if the input files are gzipped
   * @param destination destination for cat
   * @param deleteIntermediate if true, delete intermediate files as soon as we are finished with them
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void samCat(boolean inputIsGzipped, OutputStream destination, boolean deleteIntermediate, File... inputFiles) throws IOException {
    final byte[] buff = new byte[FileUtils.BUFFERED_STREAM_SIZE];
    final OneShotTimer timer = new OneShotTimer("samCat");
    for (int i = 0; i < inputFiles.length; i++) {
      final long t0 = System.nanoTime();
      boolean dropHeader = i > 0;
      boolean scanNewLine = false;
      final File current = inputFiles[i];
      final long length = current.length();
      if (length > 0) {
        try (InputStream input = inputIsGzipped ? FileUtils.createGzipInputStream(current, false) : FileUtils.createFileInputStream(current, false)) {
          int len;
          while ((len = input.read(buff)) > 0) {
            int currentPos = 0;
            if (dropHeader) {
              for (; currentPos < len; currentPos++) {
                final char c = (char) buff[currentPos];
                if (scanNewLine && c == '\n') {
                  scanNewLine = false;
                } else if (!scanNewLine) {
                  if (c == '@') {
                    scanNewLine = true;
                  } else {
                    dropHeader = false;
                    break;
                  }
                }
              }
            }
            //if dropHeader == true here, then len should == currentPos
            assert dropHeader == (len == currentPos);
            destination.write(buff, currentPos, len - currentPos);
          }
        }
      }
      devLogSamConcat(t0, length, current.getAbsolutePath());
      if (deleteIntermediate) {
        if (!current.delete()) {
          Diagnostic.developerLog("cannot delete " + current.getPath());
        }
      }
    }
    timer.stopLog();
  }

  @JumbleIgnore
  private static void devLogSamConcat(long startTime, long bytesLength, String filePath) {
    final long diff = System.nanoTime() - startTime;
    Diagnostic.developerLog("sam concat file=" + filePath + " bytes=" + bytesLength
        + " time=" + (diff / 1000000) + "ms"
        + " bytes/sec=" + Utils.realFormat(bytesLength * 1.0e9 / diff, 2));
  }

  /**
   * @return The available Cigar codes
   */
  public static char[] getCigarCodes() {
    return CIGAR_CODES.clone();
  }

  /**
   * Convert null to empty string.
   *
   * @param s string
   * @return non-null string
   */
  public static String allowEmpty(final String s) {
    return s == null ? "" : s;
  }

  /**
   * Get the reads SDF-ID from the SAM header
   * @param header SAM header containing the SDF-ID
   * @return GUID for the reads SDF from the SAM header, 0 if no GUID is specified in the SAM header
   * @throws NoTalkbackSlimException if the SAM header contains a malformed id
   */
  public static SdfId getReadsGuid(SAMFileHeader header) {
    return getSdfGuid(header.getComments(), READ_SDF_ATTRIBUTE);
  }

  /**
   * Get the reference SDF-ID from the SAM header
   * @param header SAM header containing the SDF-ID
   * @return GUID for the reference SDF from the SAM header, 0 if no GUID is specified in the SAM header
   * @throws NoTalkbackSlimException if the SAM header contains a malformed id
   */
  public static SdfId getReferenceGuid(SAMFileHeader header) {
    return getSdfGuid(header.getComments(), TEMPLATE_SDF_ATTRIBUTE);
  }

  private static SdfId getSdfGuid(List<String> header, String sdfIdLabel) {
    for (final String comment : header) {
      if (comment.replaceAll("@CO\t", "").startsWith(sdfIdLabel)) {
        final String stringGuid = comment.substring(comment.indexOf(':') + 1);
        try {
          return new SdfId(stringGuid);
        } catch (final NumberFormatException e) {
          throw new NoTalkbackSlimException("Malformed " + sdfIdLabel + " attribute from SAM header : '" + stringGuid + "'.");
        }
      }
    }
    return new SdfId(0);
  }

  private static final HashSet<String> ALREADY_REPORTED_SAM = new HashSet<>();

  /**
   * Emits the Run id of a SAM header to diagnostic user log
   * @param header the SAM file header
   */
  public static synchronized void logRunId(SAMFileHeader header) {
    for (final String comment : header.getComments()) {
      if (comment.replaceAll("@CO\t", "").startsWith(RUN_ID_ATTRIBUTE)) {
        final String stringGuid = comment.substring(comment.indexOf(':') + 1);
        if (ALREADY_REPORTED_SAM.add(stringGuid)) {
          Diagnostic.userLog("Referenced SAM file with RUN-ID: " + stringGuid);
        }
      }
    }
  }

  /**
   * Update the Run-Id comment in the header with the current id of this run
   * @param header the header in which to replace the Run-Id comment
   */
  public static void updateRunId(SAMFileHeader header) {
    final ArrayList<String> newComments = new ArrayList<>();
    for (final String comment : header.getComments()) {
      if (comment.replaceAll("@CO\t", "").startsWith(RUN_ID_ATTRIBUTE)) {
        newComments.add(SamUtils.RUN_ID_ATTRIBUTE + CommandLine.getRunId().toString());
      } else {
        newComments.add(comment);
      }
    }
    header.setComments(newComments);
  }

  /**
   * Check if the given reads SDF id matches that specified in the SAM header
   * @param header SAM header containing the GUID
   * @param sdfId SDF if to check
   * @throws NoTalkbackSlimException if the id does not match
   */
  public static void checkReadsGuid(final SAMFileHeader header, final SdfId sdfId) {
    final SdfId myGuid = getReadsGuid(header);
    if (!myGuid.check(sdfId)) {
      throw new NoTalkbackSlimException("SDF-ID of given SDF does not match SDF used during mapping.");
    } else if (!myGuid.available()) {
      Diagnostic.warning("No READ-SDF-ID found in SAM header, unable to verify read-id correctness.");
    }
  }

  /**
   * Check if the given reference SDF id matches that specified in the SAM header
   * @param header SAM header containing the GUID
   * @param referenceSdfId reference SDF id to check against
   * @throws NoTalkbackSlimException if the id does not match
   */
  public static void checkReferenceGuid(final SAMFileHeader header, final SdfId referenceSdfId) {
    checkReferenceGuid(getReferenceGuid(header), referenceSdfId);
  }


  private static void checkReferenceGuid(SdfId samGuid, SdfId referenceGuid) {
    if (!samGuid.check(referenceGuid)) {
      throw new NoTalkbackSlimException("TEMPLATE-SDF-ID of given file does not match reference SDF-ID used during mapping.");
    } else if (!samGuid.available()) {
      Diagnostic.warning("No TEMPLATE-SDF-ID found in SAM header, unable to verify reference correctness.");
    }
  }

  /**
   * Gets the reference sequence names from a SAM header, in the order they are declared.
   * @param header the SAM header
   * @return a List of the sequence names
   */
  public static List<String> getSequenceNames(SAMFileHeader header) {
    final List<String> names = new ArrayList<>();
    for (final SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
      names.add(rec.getSequenceName());
    }
    return names;
  }

  /**
   * Gets a lookup from sequence names to sequence ids
   * @param dict the sequence dictionary
   * @return the lookup.
   */
  public static Map<String, Integer> getSequenceIdLookup(SAMSequenceDictionary dict) {
    final Map<String, Integer> lookup = new HashMap<>();
    for (SAMSequenceRecord rec : dict.getSequences()) {
      lookup.put(rec.getSequenceName(), rec.getSequenceIndex());
    }
    return lookup;
  }

  /**
   * Creates a suitable program record from command line arguments and adds it to the header,
   * chaining to other program records as required.
   * @param header header to apply to
   */
  public static void addProgramRecord(SAMFileHeader header) {
    final SAMProgramRecord pg = new SAMProgramRecord(Constants.APPLICATION_NAME);
    if (CommandLine.getCommandLine() != null) {
      pg.setCommandLine(CommandLine.getCommandLine());
    } else {
      pg.setCommandLine("Internal");
    }
    pg.setProgramVersion(Environment.getVersion());
    addProgramRecord(header, pg);
  }

  /**
   * Adds given program record to the header, chaining to other program records
   * as required
   * @param header header to apply to
   * @param record record to add
   */
  public static void addProgramRecord(SAMFileHeader header, SAMProgramRecord record) {
    if (record.getProgramName() == null) {
      record.setProgramName(record.getId());
    }
    final HashSet<String> idsUsed = new HashSet<>();

    final HashMap<String, SAMProgramRecord> toAppend = new HashMap<>();
    for (final SAMProgramRecord r : header.getProgramRecords()) {
      idsUsed.add(r.getId());
      toAppend.put(r.getId(), r);
    }
    for (final SAMProgramRecord r : header.getProgramRecords()) {
      if (r.getPreviousProgramGroupId() != null && toAppend.containsKey(r.getPreviousProgramGroupId())) {
        toAppend.remove(r.getPreviousProgramGroupId());
      }
    }
    if (toAppend.size() == 0) {
      header.addProgramRecord(record);
    } else {
      for (final SAMProgramRecord r : toAppend.values()) {
        String idStr;
        if (r.getId().matches("^[a-zA-Z]+-[0-9]+$")) {
          final String[] split = StringUtils.split(r.getId(), '-');
          int id = Integer.parseInt(split[1]);
          do {
            id++;
            idStr = split[0] + "-" + id;
          } while (idsUsed.contains(idStr));
        } else {
          int id = 1;
          idStr = r.getId() + "-" + id;
          while (idsUsed.contains(idStr)) {
            id++;
            idStr = r.getId() + "-" + id;
          }
        }
        final SAMProgramRecord real = new SAMProgramRecord(idStr);
        real.setCommandLine(record.getCommandLine());
        real.setPreviousProgramGroupId(r.getId());
        real.setProgramName(record.getProgramName());
        real.setProgramVersion(record.getProgramVersion());
        header.addProgramRecord(real);
      }
    }
  }

  /**
   * @param zippedSam zipped file to check
   * @return true if it starts with a <code>SAM</code> header (i.e. first character is @)
   * @throws IOException if an IO Error occurs
   */
  public static boolean looksLikeSam(File zippedSam) throws IOException {
    final boolean result;
    try (Reader fr = new InputStreamReader(FileUtils.createGzipInputStream(zippedSam, false))) {
      result = fr.read() == (int) '@';
    }
    return result;
  }

  /**
   * @param file the file to check.
   * @return true if this looks like a BAM file.
   * @throws IOException if an IO Error occurs
   */
  public static boolean isBAMFile(final File file) throws IOException {
    final boolean result;
    try (BufferedInputStream bis = new BufferedInputStream(new FileInputStream(file))) {
      if (!BlockCompressedInputStream.isValidFile(bis)) {
        return false;
      }
      final int buffSize = BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE;
      bis.mark(buffSize);
      final byte[] buffer = new byte[buffSize];
      final int len = IOUtils.readAmount(bis, buffer, 0, buffSize);
      bis.reset();
      final byte[] magicBuf = new byte[4];
      final int magicLength = IOUtils.readAmount(new BlockCompressedInputStream(new ByteArrayInputStream(buffer, 0, len)), magicBuf, 0, 4);
      //checks we read 4 bytes and they were "BAM\1" in ascii
      result = magicLength == 4 && Arrays.equals(new byte[]{(byte) 66, (byte) 65, (byte) 77, (byte) 1}, magicBuf);

    }
    return result;
  }


  /**
   * Get the NH or IH value
   * @param sam the sam record
   * @return the value contained in the NH attribute, or if that doesn't exist,
   * the IH attribute. If neither exist, returns null.
   */
  public static Integer getNHOrIH(SAMRecord sam) {
    final Integer nh = sam.getIntegerAttribute(ATTRIBUTE_NH);
    if (nh != null) {
      return nh;
    }
    return sam.getIntegerAttribute(ATTRIBUTE_IH);
  }


  /**
   * Expand a SAM record's quality field using Super Cigar quality overlap field.
   * @param qualityBuffer a 35-long byte array buffer for the expanded quality
   * @param qualityOverlap the overlap quality values
   * @param samQualities the flattened qualities from the SAM record
   * @param first true if first of pair
   * @param reverseCompliment true if reverse compliment
   * @param phredifyQualityOverlap true if the qualities in the overlap region should be converted to phred scale (i.e. have 33 subtracted)
   * @return the expanded quality array
   * @throws BadSuperCigarException if the <code>samQualities</code> length plus overlap length does not equal the expected value
   */
  public static byte[] expandCgSuperCigarQualities(byte[] samQualities, byte[] qualityBuffer, String qualityOverlap, boolean first, boolean reverseCompliment, boolean phredifyQualityOverlap) throws BadSuperCigarException {
    assert qualityBuffer.length == SamUtils.CG_RAW_READ_LENGTH;
    final int overlapQualityOffset = phredifyQualityOverlap ? FastaUtils.PHRED_LOWER_LIMIT_CHAR : 0;
    final int xqlength = qualityOverlap == null ? 0 : qualityOverlap.length();
    if (samQualities.length + xqlength != qualityBuffer.length) {
      throw new BadSuperCigarException("SAM record qualities plus " + SamUtils.CG_OVERLAP_QUALITY + " not expected length. Was: " + (samQualities.length + xqlength) + " expected: " + qualityBuffer.length);
    }
    if ((first && !reverseCompliment) || (reverseCompliment && !first)) {
      System.arraycopy(samQualities, 0, qualityBuffer, 0, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
      for (int i = 0; i < xqlength; i++) {
        qualityBuffer[CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + i] = (byte) (qualityOverlap.charAt(i) - overlapQualityOffset);
      }
      System.arraycopy(samQualities, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, qualityBuffer, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + xqlength, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
    } else {
      System.arraycopy(samQualities, 0, qualityBuffer, 0, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
      for (int i = 0; i < xqlength; i++) {
        qualityBuffer[samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + i] = (byte) (qualityOverlap.charAt(i) - overlapQualityOffset);
      }
      System.arraycopy(samQualities, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, qualityBuffer, qualityBuffer.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
    }
    if (reverseCompliment) {
      Utils.reverseInPlace(qualityBuffer);
    }
    return qualityBuffer;
  }

  /**
   * @param header combined sam header
   * @return creates a map of read group to sample id
   */
  public static Map<String, String> getReadGroupToSampleId(final SAMFileHeader header) {
    final HashMap<String, String> readGroupToSampleMap = new HashMap<>();
    for (final SAMReadGroupRecord rec : header.getReadGroups()) {
      //System.err.println("k=" + rec.getReadGroupId() + " v=" + rec.getSample());
      readGroupToSampleMap.put(rec.getReadGroupId(), rec.getSample());
    }
    return readGroupToSampleMap;
  }

  /**
   * Get a list of the sample names mentioned in the header in sorted order.
   * @param header combined sam header
   * @return creates an array of sample names
   */
  public static String[] getSampleNames(final SAMFileHeader header) {
    final Set<String> sampleNames = new TreeSet<>();
    for (final SAMReadGroupRecord rec : header.getReadGroups()) {
      if (rec.getSample() != null) {
        sampleNames.add(rec.getSample());
      }
    }
    return sampleNames.toArray(new String[sampleNames.size()]);
  }

  /**
   * Converts a CIGAR that may contain new style operators for match and
   * mismatch into the legacy style CIGAR.
   * @param cigar the CIGAR to convert
   * @return the legacy CIGAR
   */
  public static Cigar convertToLegacyCigar(Cigar cigar) {
    final Cigar cg = new Cigar();
    int count = 0;
    for (int i = 0; i < cigar.numCigarElements(); i++) {
      final CigarElement ce = cigar.getCigarElement(i);
      if (ce.getOperator().equals(CigarOperator.EQ) || ce.getOperator().equals(CigarOperator.X)) {
        count += ce.getLength();
      } else {
        if (count > 0) {
          cg.add(new CigarElement(count, CigarOperator.M));
        }
        cg.add(ce);
        count = 0;
      }
    }
    if (count > 0) {
      cg.add(new CigarElement(count, CigarOperator.M));
    }
    return cg;
  }

  /**
   * Updates the CIGAR in the supplied record to be legacy-compatible.
   * @param record the <code>SAMRecord</code> to update.
   */
  public static void convertToLegacyCigar(SAMRecord record) {
    record.setCigar(convertToLegacyCigar(record.getCigar()));
  }

  /**
   * Converts a SAM file name into one with an appropriate extension, depending on whether the file should be gzipped.
   * @param gzip true if the output file is destined to be GZIP compressed.
   * @param file the input file
   * @return the appropriately adjusted file
   */
  public static File getZippedSamFileName(final boolean gzip, final File file) {
    final String name = file.getName().toLowerCase(Locale.getDefault());
    final File samfile = name.endsWith(SAM_SUFFIX) || gzip && name.endsWith(SAM_SUFFIX + FileUtils.GZ_SUFFIX) ? file : new File(file.getPath() + SAM_SUFFIX);
    return FileUtils.getZippedFileName(gzip, samfile);
  }

  /**
   * Prune Illumina arm marker.
   * @param name read name
   * @param isPaired whether in paired end mode
   * @return pruned read name
   */
  public static String samReadName(String name, boolean isPaired) {
    if (isPaired && name != null && name.length() > 2
        && name.charAt(name.length() - 2) == '/'
        && (name.endsWith("/1") || name.endsWith("/2"))) {
      return name.substring(0, name.length() - 2);
    }
    return name;
  }

  /**
   * Method to check the equivalence of two SAM headers
   * @param fh a <code>SAMFileHeader</code> value
   * @param lh a <code>SAMFileHeader</code> value
   * @return true if the headers are compatible.
   */
  public static boolean checkHeaderDictionary(final SAMFileHeader fh, final SAMFileHeader lh) {
    if (fh.getSortOrder() != lh.getSortOrder()) {
      return false;
    }
    final List<SAMSequenceRecord> flist = fh.getSequenceDictionary().getSequences();
    final List<SAMSequenceRecord> llist = lh.getSequenceDictionary().getSequences();
    final Iterator<SAMSequenceRecord> fi = flist.iterator();
    final Iterator<SAMSequenceRecord> li = llist.iterator();
    while (fi.hasNext()) {
      if (!li.hasNext()) {
        return false;
      }
      final SAMSequenceRecord fsr = fi.next();
      final SAMSequenceRecord lsr = li.next();
      if (!fsr.getSequenceName().equals(lsr.getSequenceName()) || fsr.getSequenceLength() != lsr.getSequenceLength()) {
        return false;
      }
    }
    if (li.hasNext()) {
      return false;
    }
    return true;
  }

  /**
   * Open file, return header, close file
   * @param f SAM file
   * @return the header
   * @throws IOException if an I/O error occurs
   */
  public static SAMFileHeader getSingleHeader(File f) throws IOException {
    final SAMFileHeader result;
    try (SamReader sr = SamUtils.makeSamReader(f)) {
      result = sr.getFileHeader();
    }
    return result;
  }

  /**
   * Read header from stream and then close
   * @param is the InputStream
   * @return the header
   * @throws IOException if an I/O error occurs
   */
  public static SAMFileHeader getSingleHeader(InputStream is) throws IOException {
    final SAMFileHeader result;
    try (SamReader sr = SamUtils.makeSamReader(is, null)) {
      result = sr.getFileHeader();
    }
    return result;
  }

  /**
   * creates a header with the contents from the first file except the read group information is merged from all headers. Does some checking that headers are compatible
   *
   * @param files SAM files
   * @return the combined header
   * @throws IOException if an I/O error occurs
   */
  public static SAMFileHeader getUberHeader(Collection<File> files) throws IOException {
    return getUberHeader(files, false, null);
  }

  /**
   * creates a header with the contents from the first file except the read group information is merged from all headers. Does some checking that headers are compatible
   *
   * @param files SAM files
   * @param ignoreHeaderIncompatibility true if should not care about incompatible header
   * @param expectedSamples if non-null, check that headers contain sample information that overlaps the supplied names
   * @return the combined header
   * @throws IOException if an I/O error occurs
   */
  public static SAMFileHeader getUberHeader(Collection <File> files, boolean ignoreHeaderIncompatibility, String[] expectedSamples) throws IOException {
    if (files.size() == 0) {
      throw new IllegalArgumentException("File list is empty!");
    }
    final HashMap<String, SAMReadGroupRecord> readGroups = new HashMap<>();
    final HashMap<String, String> readGroupsSampleMap = new HashMap<>();
    SAMFileHeader first = null;
    File firstFile = null;
    final StringBuilder errorMessage = new StringBuilder();
    SdfId currentGuid = null;
    final HashSet<String> expectedSamplesSet = expectedSamples == null ? null : new HashSet<>(Arrays.asList(expectedSamples));
    boolean guidMismatch = false;
    for (final File file : files) {
      try (SamReader sfr = SamUtils.makeSamReader(file)) {
        if (first == null) {
          first = sfr.getFileHeader();
          firstFile = file;
          currentGuid = getReferenceGuid(first);
        } else if (!checkHeaderDictionary(first, sfr.getFileHeader())) {
          Diagnostic.warning(WarningType.SAM_INCOMPATIBLE_HEADERS, firstFile.getPath(), file.getPath());
          if (!ignoreHeaderIncompatibility) {
            throw new NoTalkbackSlimException(ErrorType.SAM_INCOMPATIBLE_HEADER_ERROR, "1");
          }
        }
        final SdfId fileGuid = getReferenceGuid(sfr.getFileHeader());
        if (!currentGuid.check(fileGuid)) {
          guidMismatch = true;
        } else if (!currentGuid.available()) {
          currentGuid = fileGuid;
        }
        if (sfr.getFileHeader().getReadGroups().size() != 0) {
          for (final SAMReadGroupRecord r : sfr.getFileHeader().getReadGroups()) {
            final String sample = r.getSample();
            if (!readGroups.containsKey(r.getReadGroupId())) {
              readGroups.put(r.getReadGroupId(), r);
              readGroupsSampleMap.put(r.getReadGroupId(), sample);
            } else {
              //check that the sample isn't different for the same read group id
              if (sample != null && !r.getSample().equals(readGroupsSampleMap.get(r.getReadGroupId()))) {
                Diagnostic.warning(file.getPath() + " contained read group with id \"" + r.getId() + "\" and sample \"" + sample + "\" but this read group has already been associated with sample \"" + readGroupsSampleMap.get(r.getReadGroupId()) + "\"");
                if (!ignoreHeaderIncompatibility) {
                  throw new NoTalkbackSlimException(ErrorType.SAM_INCOMPATIBLE_HEADER_ERROR, "1");
                }
              }
            }
            if (expectedSamplesSet != null) {
              if (sample != null) {
                if (!expectedSamplesSet.contains(sample)) {
                  errorMessage.append("Unexpected read group sample name: ").append(sample).append('.').append(StringUtils.LS);
                }
              } else {
                errorMessage.append("Input file \"").append(file.getPath()).append("\" contains a read group with no sample tag: ").append(r.getId()).append('.').append(StringUtils.LS);
              }
            }
          }
        } else if (expectedSamples != null) {
          errorMessage.append("Input file \"").append(file.getPath()).append("\" does not contain read group information.").append(StringUtils.LS);
        }
      }
    }
    if (guidMismatch) {
      Diagnostic.warning("Input SAM files contain mismatching template GUIDs");
    }
    if (errorMessage.length() > 0) {
      throw new NoTalkbackSlimException(errorMessage.toString().trim());
    }

    final SAMFileHeader header = first;
    if (readGroups.size() > 0) {
      final List<SAMReadGroupRecord> recList = new ArrayList<>();
      for (final SAMReadGroupRecord r : readGroups.values()) {
        recList.add(r);
      }
      header.setReadGroups(recList);
    }
    if (currentGuid.available() && !getReferenceGuid(header).available()) {
      header.addComment(TEMPLATE_SDF_ATTRIBUTE + currentGuid.toString());
    }
    return header;
  }

  /**
   * Check the compatibility of SAM files against a reference SDF.  Throws {@link NoTalkbackSlimException} if strict mode is set, warnings otherwise.
   * @param reference reference SDF sequences reader
   * @param uberHeader the combined header from a set of SAM files to check against
   * @param strict whether to throw exception for sequence dictionary mismatches
   * @throws IOException if error occurs while reading files
   * @throws NoTalkbackSlimException if SAM header sequence dictionaries mismatch, and in strict mode
   */
  public static void checkUberHeaderAgainstReference(SequencesReader reference, SAMFileHeader uberHeader, boolean strict) throws IOException, NoTalkbackSlimException {
    final SdfId refId = reference.getSdfId();

    final long numSequences = reference.numberSequences();
    final HashSet<String> refNames = new HashSet<>();
    for (long i = 0; i < numSequences; i++) {
      refNames.add(reference.name(i));
    }
    final SdfId samGuid = getReferenceGuid(uberHeader);
    if (!samGuid.check(refId)) {
      Diagnostic.warning("TEMPLATE-SDF-ID does not match reference SDF-ID.");
    }

    final StringBuilder sb = new StringBuilder();
    boolean error = false;
    final HashSet<String> notSeen = new HashSet<>(refNames);
    for (SAMSequenceRecord samRecord : uberHeader.getSequenceDictionary().getSequences()) {
      if (!refNames.contains(samRecord.getSequenceName())) {
        sb.append("Sequence '").append(samRecord.getSequenceName()).append("' not in reference SDF").append(StringUtils.LS);
        error = true;
      } else {
        notSeen.remove(samRecord.getSequenceName());
      }
    }
    if (error && strict) {
      throw new NoTalkbackSlimException(sb.toString());
    }
    for (String name : notSeen) {
      sb.append("Sequence '").append(name).append("' not in SAM header").append(StringUtils.LS);
    }
    if (sb.length() > 0) {
      Diagnostic.warning(sb.toString());
    }
  }

  private static SamReaderFactory getSamReaderFactory() {
    return SamReaderFactory.make()
      .validationStringency(ValidationStringency.SILENT);
  }

  /**
   * Entry point for specifically creating a SamReader given a pre-positioned stream, header, and known type
   * @param stream the stream to read from. Must already be performing decompression if required.
   * @param headerOverride the pre-determined SAM header
   * @param assumeType the type of input to assume.
   * @return the SamReader
   * @throws IOException if an I/O problem occurs opening the file
   */
  public static SamReader makeSamReader(InputStream stream, SAMFileHeader headerOverride, SamReader.Type assumeType) throws IOException {
    if (assumeType == null) {
      throw new NullPointerException();
    }
    try {
      return getSamReaderFactory()
        .open(SamInputResource.of(stream).header(headerOverride).assumeType(assumeType));
    } catch (final RuntimeIOException e) {
      throw (IOException) e.getCause();
    }
  }

  /**
   * Entry point for specifically creating a SamReader given a provided stream and header, but let
   * htsjdk decide the underlying format (including working out whether the input is compressed).
   * @param stream the stream to read from
   * @param headerOverride the pre-determined SAM header (or null to use the header from the stream)
   * @return the SamReader
   * @throws IOException if an I/O problem occurs opening the file
   */
  public static SamReader makeSamReader(InputStream stream, SAMFileHeader headerOverride) throws IOException {
    try {
      return getSamReaderFactory().open(SamInputResource.of(stream).header(headerOverride));
    } catch (final RuntimeIOException e) {
      throw (IOException) e.getCause();
    }
  }

  /**
   * Entry point for creating SamReaders using our preferences
   * @param stream the stream to read from
   * @return the SamReader
   * @throws IOException if an I/O problem occurs opening the file
   */
  public static SamReader makeSamReader(InputStream stream) throws IOException {
    return makeSamReader(stream, null);
  }

  /**
   * Entry point for creating SamReaders using our preferences
   * @param file the file to open
   * @return the SamReader
   * @throws IOException if an I/O problem occurs opening the file
   */
  public static SamReader makeSamReader(File file) throws IOException {
    try {
      return getSamReaderFactory()
        .open(file);
    } catch (final RuntimeIOException e) {
      throw (IOException) e.getCause();
    }
  }
}
