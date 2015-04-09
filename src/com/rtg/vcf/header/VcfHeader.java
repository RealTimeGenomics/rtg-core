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

package com.rtg.vcf.header;

import static com.rtg.util.StringUtils.TAB;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Class holds VCF header lines and sample names
 */
public class VcfHeader {

  /** current VCF version */
  public static final String VERSION = "4.1";

  /** current VCF version */
  public static final String VERSION_VALUE = "VCFv" + VERSION;

  /** Comment character for VCF files */
  public static final char COMMENT_CHAR = '#';

  /** string used to indicate meta lines */
  public static final String META_STRING = "##";

  /** Start of reference sequence lines*/
  public static final String CONTIG_STRING = META_STRING + "contig";
  /** Start of alt lines*/
  public static final String ALT_STRING = META_STRING + "ALT";
  /** Start of info lines*/
  public static final String INFO_STRING = META_STRING + "INFO";
  /** Start of filter lines */
  public static final String FILTER_STRING = META_STRING + "FILTER";
  /** Start of format lines */
  public static final String FORMAT_STRING = META_STRING + "FORMAT";
  /** Start of sample lines */
  public static final String SAMPLE_STRING = META_STRING + "SAMPLE";
  /** Start of pedigree lines */
  public static final String PEDIGREE_STRING = META_STRING + "PEDIGREE";

  /** file format line prefix */
  public static final String VERSION_LINE_PREFIX = META_STRING + "fileformat";
  /** full version string */
  public static final String VERSION_LINE = VERSION_LINE_PREFIX + "=" + VERSION_VALUE;

  /** header line for VCF files */
  public static final String HEADER_BASE = "" + COMMENT_CHAR
      + "CHROM" + TAB
      + "POS" + TAB
      + "ID" + TAB
      + "REF" + TAB
      + "ALT" + TAB
      + "QUAL" + TAB
      + "FILTER" + TAB
      + "INFO";

  /** header portion for format column */
  public static final String FORMAT_HEADER_STRING = "FORMAT";

  /** header line for VCF files with samples */
  public static final String HEADER_LINE = HEADER_BASE
          + TAB + FORMAT_HEADER_STRING;

  /** Minimal string that can be used as a VCF header */
  public static final String MINIMAL_HEADER = VERSION_LINE + '\n' + HEADER_LINE;

  private static final String[] HEADER_COLUMNS = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};

  private String mVersionLine;
  private final List<String> mGenericMetaInformationLines;
  private final List<ContigField> mContigLines;
  private final List<AltField> mAltLines;
  private final List<FilterField> mFilterLines;
  private final List<InfoField> mInfoLines;
  private final List<FormatField> mFormatLines;
  private final List<SampleField> mSampleLines;
  private final List<PedigreeField> mPedigreeLines;
  private final List<String> mSampleNames;
  private final HashMap<String, Integer> mNameToColumn;

  /**
   * create a new VCF header
   */
  public VcfHeader() {
    mGenericMetaInformationLines = new ArrayList<>();
    mSampleNames = new ArrayList<>();
    mContigLines = new ArrayList<>();
    mAltLines = new ArrayList<>();
    mFilterLines = new ArrayList<>();
    mInfoLines = new ArrayList<>();
    mFormatLines = new ArrayList<>();
    mSampleLines = new ArrayList<>();
    mPedigreeLines = new ArrayList<>();
    mNameToColumn = new HashMap<>();
  }

  /**
   * Create a copy of this header
   * @return the copy
   */
  public VcfHeader copy() {
    final VcfHeader copy = new VcfHeader();
    copy.mVersionLine = mVersionLine;
    copy.mGenericMetaInformationLines.addAll(mGenericMetaInformationLines);
    copy.mSampleNames.addAll(mSampleNames);
    copy.mContigLines.addAll(mContigLines);
    copy.mAltLines.addAll(mAltLines);
    copy.mFilterLines.addAll(mFilterLines);
    copy.mInfoLines.addAll(mInfoLines);
    copy.mFormatLines.addAll(mFormatLines);
    copy.mSampleLines.addAll(mSampleLines);
    copy.mPedigreeLines.addAll(mPedigreeLines);
    copy.mNameToColumn.putAll(mNameToColumn);
    return copy;
  }

  /**
   * Add the common header fields used in typical new VCF files
   */
  public void addCommonHeader() {
    addLine(VERSION_LINE);
    final Calendar cal = Calendar.getInstance();
    final SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd");
    addLine(META_STRING + "fileDate=" + sdf.format(cal.getTime()));
    addLine(META_STRING + "source=" + Environment.getVersion());
    addRunInfo();
  }

  /**
   * Add header fields giving run id and command line args, for appending to an existing VCF file
   */
  public void addRunInfo() {
    addLine(META_STRING + "CL=" + CommandLine.getCommandLine());
    addLine(META_STRING + "RUN-ID=" + CommandLine.getRunId());
  }

  /**
   * Add a reference line corresponding to the supplied reader
   * @param reference the reference SequencesReader
   */
  public void addReference(SequencesReader reference) {
    final SdfId sdfId = reference.getSdfId();
    if (sdfId != null && sdfId.available()) {
      addLine(VcfHeader.META_STRING + "TEMPLATE-SDF-ID=" + sdfId);
    }
    addLine(VcfHeader.META_STRING + "reference=" + reference.path());
  }

  /**
   * Add contig lines corresponding to the sequences present in a SAM header.
   * @param header the SAM header.
   */
  public void addContigFields(SAMFileHeader header) {
    final SAMSequenceDictionary dic =  header.getSequenceDictionary();
    for (final SAMSequenceRecord seq : dic.getSequences()) {
      final String astag = (seq.getAttribute(SAMSequenceRecord.ASSEMBLY_TAG) == null) ? "" : (",as=" + seq.getAttribute(SAMSequenceRecord.ASSEMBLY_TAG));
      final String md5tag = (seq.getAttribute(SAMSequenceRecord.MD5_TAG) == null) ? "" : (",md5=" + seq.getAttribute(SAMSequenceRecord.MD5_TAG));
      final String speciestag = (seq.getAttribute(SAMSequenceRecord.SPECIES_TAG) == null) ? "" : (",species=\"" + seq.getAttribute(SAMSequenceRecord.SPECIES_TAG) + "\"");
      addLine(META_STRING + "contig=<ID=\"" + seq.getSequenceName() + "\",length=" + seq.getSequenceLength() + astag + md5tag + speciestag + ">");
    }
  }

  /**
   * Add contig lines corresponding to the sequences present in a sequences reader.
   * @param reader the sequences reader.
   * @throws IOException if an IO error occurs whilst reading
   */
  public void addContigFields(SequencesReader reader) throws IOException {
    for (long i = 0; i < reader.numberSequences(); i++) {
      addLine(META_STRING + "contig=<ID=\"" + reader.name(i) + "\",length=" + reader.length(i) + ">");
    }
  }

  private <T extends IdField<T>> void addIdField(List<T> dest, T field) {
    for (T f : dest) {
      if (f.getId().equals(field.getId())) {
        throw new IllegalArgumentException("VCF header contains multiple field declarations with the same ID=" + field.getId() + StringUtils.LS
          + f.toString() + StringUtils.LS
          + field.toString());
      }
    }
    dest.add(field);
  }

  /**
   * Add an alt field
   * @param field the new alt field
   */
  public void addAltField(AltField field) {
    addIdField(mAltLines, field);
  }

  /**
   * Add a filter field
   * @param id the filter field identifier
   * @param description the field description
   */
  public void addFilterField(String id, String description) {
    addFilterField(new FilterField(id, description));
  }

  /**
   * Add a filter field
   * @param field the filter field
   */
  public void addFilterField(FilterField field) {
    addIdField(mFilterLines, field);
  }

  /**
   * Add an info field
   * @param id the info field identifier
   * @param type the type of value
   * @param number the specifier for the number of occurrences
   * @param description the field description
   */
  public void addInfoField(String id, MetaType type, VcfNumber number, String description) {
    addInfoField(new InfoField(id, type, number, description));
  }

  /**
   * Add an info field
   * @param field the new info field
   */
  public void addInfoField(InfoField field) {
    addIdField(mInfoLines, field);
  }

  /**
   * Add a format field
   * @param id the format field identifier
   * @param type the type of value
   * @param number the specifier for the number of occurrences
   * @param description the field description
   */
  public void addFormatField(String id, MetaType type, VcfNumber number, String description) {
    addFormatField(new FormatField(id, type, number, description));
  }

  /**
   * Add a format field
   * @param field the new format field
   */
  public void addFormatField(FormatField field) {
    addIdField(mFormatLines, field);
  }

  /**
   * Ensure that the header contains the specified info field (or one that is compatible)
   * @param field the new format field
   */
  public void ensureContains(FilterField field) {
    for (FilterField f : getFilterLines()) {
      if (f.getId().equals(field.getId())) {
        return; // Field already present
      }
    }
    addFilterField(field);
  }

  /**
   * Ensure that the header contains the specified info field (or one that is compatible)
   * @param field the new format field
   */
  public void ensureContains(InfoField field) {
    for (InfoField f : getInfoLines()) {
      if (f.getId().equals(field.getId())) {
        if (f.getType() == field.getType() && f.getNumber().equals(field.getNumber())) {
          return; // Field already present
        } else {
          throw new NoTalkbackSlimException("A VCF INFO field " + field.getId() + " which is incompatible is already present in the VCF header.");
        }
      }
    }
    addInfoField(field);
  }

  /**
   * Ensure that the header contains the specified format field (or one that is compatible)
   * @param field the new format field
   */
  public void ensureContains(FormatField field) {
    for (FormatField f : getFormatLines()) {
      if (f.getId().equals(field.getId())) {
        if (f.getType() == field.getType() && f.getNumber().equals(field.getNumber())) {
          return; // Field already present
        } else {
          throw new NoTalkbackSlimException("A VCF FORMAT field " + field.getId() + " which is incompatible is already present in the VCF header.");
        }
      }
    }
    addFormatField(field);
  }

  /**
   * @return get the meta line containing file format version
   */
  public String getVersionLine() {
    return mVersionLine;
  }

  /**
   * @param ver set the file format version value
   */
  public void setVersionValue(String ver) {
    mVersionLine = VERSION_LINE_PREFIX + "=" + ver;
  }

  /**
   * @return get just the value of the file format version line
   */
  public String getVersionValue() {
    if (mVersionLine == null) {
      return null;
    }
    final String[] split = mVersionLine.split("=", 2);
    if (split.length < 2) {
      //unpossible
      throw new NoTalkbackSlimException("VCF version line does not contain a value");
    }
    return split[1];
  }

  /**
   * @return meta information lines
   */
  public List<String> getGenericMetaInformationLines() {
    return mGenericMetaInformationLines;
  }
  public List<ContigField> getContigLines() {
    return mContigLines;
  }
  public List<AltField> getAltLines() {
    return mAltLines;
  }
  public List<FilterField> getFilterLines() {
    return mFilterLines;
  }
  public List<InfoField> getInfoLines() {
    return mInfoLines;
  }
  public List<FormatField> getFormatLines() {
    return mFormatLines;
  }
  public List<SampleField> getSampleLines() {
    return mSampleLines;
  }
  public List<PedigreeField> getPedigreeLines() {
    return mPedigreeLines;
  }

  /**
   * @param line meta information line
   * @return this, for call chaining
   */
  public VcfHeader addMetaInformationLine(String line) {
    if (isMetaLine(line)) {
      if (isContigLine(line)) {
        addIdField(mContigLines, parseContigLine(line));
      } else if (isAltLine(line)) {
        addAltField(parseAltLine(line));
      } else if (isFilterLine(line)) {
        addFilterField(parseFilterLine(line));
      } else if (isInfoLine(line)) {
        addInfoField(parseInfoLine(line));
      } else if (isFormatLine(line)) {
        addFormatField(parseFormatLine(line));
      } else if (isSampleLine(line)) {
        addIdField(mSampleLines, parseSampleLine(line));
      } else if (isPedigreeLine(line)) {
        mPedigreeLines.add(parsePedigreeLine(line));
      } else {
        if (isVersionLine(line)) {
          if (mVersionLine != null) {
            throw new NoTalkbackSlimException("More than one version line found");
          }
          mVersionLine = line;
        } else {
          mGenericMetaInformationLines.add(line);
        }
      }
    } else {
      throw new NoTalkbackSlimException("Not a meta information line: " + line);
    }
    return this;
  }

  /**
   * parse and add an arbitrary header/meta line
   * @param line line to add
   * @return this for chaining
   */
  public VcfHeader addLine(String line) {
    if (isMetaLine(line)) {
      return addMetaInformationLine(line);
    } else if (line.startsWith("#")) {
      return addColumnHeaderLine(line);
    } else {
      throw new NoTalkbackSlimException("Not a header line: " + line);
    }
  }

  /**
   * Add sample names contained in given header line
   * @param line column header line
   * @return this, for chain calling
   * @throws IllegalArgumentException if header is not correctly formed
   */
  public VcfHeader addColumnHeaderLine(String line) {
    final String[] split = line.split("\t");
    if (split.length < 8) {
      throw new NoTalkbackSlimException("VCF header line missing required columns");
    }
    if (split.length == 9) {
     throw new NoTalkbackSlimException("VCF header line contains format field but no sample fields");
    }
    for (int i = 0; i < 8; i++) {
      if (!split[i].equals(HEADER_COLUMNS[i])) {
        throw new NoTalkbackSlimException("Incorrect VCF header column " + (i + 1) + " expected \"" + HEADER_COLUMNS[i] + "\" was \"" + split[i] + "\"");
      }
    }
    if (split.length > 9) {
      for (int i = 9; i < split.length; i++) {
        addSampleName(split[i]);
      }
    }
    return this;
  }


  /**
   * convert <code>contig</code> line into <code>ContigField</code> object
   * @param line the line
   * @return the object
   */
  public static ContigField parseContigLine(String line) {
    return new ContigField(line);
  }

  /**
   * convert info line into <code>InfoField</code> object
   * @param line the line
   * @return the object
   */
  public static InfoField parseInfoLine(String line) {
    return new InfoField(line);
  }
  /**
   * convert alt line into <code>AltField</code> object
   * @param line the line
   * @return the object
   */
  public static AltField parseAltLine(String line) {
    return new AltField(line);
  }
  /**
   * convert filter line into <code>FilterField</code> object
   * @param line the line
   * @return the object
   */
  public static FilterField parseFilterLine(String line) {
    return new FilterField(line);
  }

  /**
   * convert format line into <code>FormatField</code> object
   * @param line the line
   * @return the object
   */
  public static FormatField parseFormatLine(String line) {
    return new FormatField(line);
  }

  /**
   * convert format line into <code>SampleField</code> object
   * @param line the line
   * @return the object
   */
  public static SampleField parseSampleLine(String line) {
    return new SampleField(line);
  }

  /**
   * convert format line into <code>PedigreeField</code> object
   * @param line the line
   * @return the object
   */
  public static PedigreeField parsePedigreeLine(String line) {
    return new PedigreeField(line);
  }

  /**
   * @param line line of <code>VCF</code> file
   * @return true if it is a file format version meta line
   */
  public static boolean isVersionLine(String line) {
    return line.startsWith(VERSION_LINE_PREFIX + "=");
  }

  /**
   * @param line line of <code>VCF</code> file
   * @return if line is meta information
   */
  public static boolean isMetaLine(String line) {
    return line.startsWith(META_STRING);
  }

  /**
   * @param line line of <code>VCF</code> file
   * @return if line is reference sequence line
   */
  public static boolean isContigLine(String line) {
    return line.startsWith(CONTIG_STRING);
  }

  /**
   * @param line line of <code>VCF</code> file
   * @return if line is information line
   */
  public static boolean isInfoLine(String line) {
    return line.startsWith(INFO_STRING);
  }
  /**
   * @param line line of <code>VCF</code> file
   * @return if line is alt line
   */
  public static boolean isAltLine(String line) {
    return line.startsWith(ALT_STRING);
  }
  /**
   * @param line line of <code>VCF</code> file
   * @return if line is filter line
   */
  public static boolean isFilterLine(String line) {
    return line.startsWith(FILTER_STRING);
  }
  /**
   * @param line line of <code>VCF</code> file
   * @return if line is format line
   */
  public static boolean isFormatLine(String line) {
    return line.startsWith(FORMAT_STRING);
  }
  /**
   * @param line line of <code>VCF</code> file
   * @return if line is sample line
   */
  public static boolean isSampleLine(String line) {
    return line.startsWith(SAMPLE_STRING);
  }
  /**
   * @param line line of <code>VCF</code> file
   * @return if line is pedigree line
   */
  public static boolean isPedigreeLine(String line) {
    return line.startsWith(PEDIGREE_STRING);
  }

  /**
   * add sample name
   * @param name name to be added
   * @return this, for call chaining
   */
  public VcfHeader addSampleName(String name) {
    if (mSampleNames.contains(name)) {
      throw new NoTalkbackSlimException("Duplicate sample name \"" + name + "\" in VCF header");
    }
    mNameToColumn.put(name, mSampleNames.size());
    mSampleNames.add(name);
    return this;
  }

  /**
   * Get the column index of the specified sample
   * @param name the sample name
   * @return the index of the specified sample, or null if the sample is unknown
   */
  public Integer getSampleIndex(String name) {
    return mNameToColumn.get(name);
  }

  /**
   * Remove sample names from this record
   * @param names hash set of names to remove
   */
  public void removeSamples(HashSet<String> names) {
    final ArrayList<String> samplesclone = new ArrayList<>(mSampleNames);
    for (final String sampleName : samplesclone) {
      if (names.contains(sampleName)) {
        mSampleNames.remove(sampleName);
      }
    }
    //now reindex the mNameToColumn thing.
    mNameToColumn.clear();
    for (int i = 0; i < mSampleNames.size(); i++) {
      mNameToColumn.put(mSampleNames.get(i), i);
    }

    //and clear out any sample lines
    final Iterator<SampleField> iterator = mSampleLines.iterator();
    while (iterator.hasNext()) {
      final SampleField sample = iterator.next();
      if (names.contains(sample.getId())) {
        iterator.remove();
      }
    }
  }

  /**
   * Remove all samples from this header.
   */
  public void removeAllSamples() {
    mNameToColumn.clear();
    mSampleLines.clear();
    mSampleNames.clear();
  }

  /**
   * @return sample names
   */
  public List<String> getSampleNames() {
    return mSampleNames;
  }

  /**
   * @return number of samples in current header
   */
  public int getNumberOfSamples() {
    return mSampleNames.size();
  }

  /**
   * @return header string
   */
  public String getColumnHeader() {
    final StringBuilder sb;
    if (getNumberOfSamples() > 0) {
      sb = new StringBuilder(HEADER_LINE);
      for (int i = 0; i < getNumberOfSamples(); i++) {
        sb.append(TAB).append(mSampleNames.get(i));
      }
    } else {
      sb = new StringBuilder(HEADER_BASE);
    }
    return sb.toString();
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    if (mVersionLine != null) {
      sb.append(mVersionLine).append('\n');
    }
    for (String line : getGenericMetaInformationLines()) {
      sb.append(line).append('\n');
    }
    for (ContigField mContigLine : mContigLines) {
      sb.append(mContigLine).append('\n');
    }
    for (AltField mAltLine : mAltLines) {
      sb.append(mAltLine).append('\n');
    }
    for (FilterField mFilterLine : mFilterLines) {
      sb.append(mFilterLine).append('\n');
    }
    for (InfoField mInfoLine : mInfoLines) {
      sb.append(mInfoLine).append('\n');
    }
    for (FormatField mFormatLine : mFormatLines) {
      sb.append(mFormatLine).append('\n');
    }
    for (SampleField mSampleLine : mSampleLines) {
      sb.append(mSampleLine).append('\n');
    }
    for (PedigreeField mPedigreeLine : mPedigreeLines) {
      sb.append(mPedigreeLine).append('\n');
    }
    sb.append(getColumnHeader()).append("\n");
    return sb.toString();
  }

  private static final Pattern KEY_VALUE_PATTERN = Pattern.compile("^([^=]+)=([^,\"]*)\\s*,?");
  private static final Pattern KEY_VALUE_WITH_QUOTES_PATTERN = Pattern.compile("^([^=]+)=\"([^\"\\\\]*(\\\\.[^\"\\\\]*)*)\"\\s*,?");

  static LinkedHashMap<String, String> parseMetaLine(String line, Pattern linePattern) {
    final Matcher m = linePattern.matcher(line.trim());
    if (!m.matches()) {
      throw new NoTalkbackSlimException("Could not parse header line: " + line);
    }
    final LinkedHashMap<String, String> ret = new LinkedHashMap<>(4);
    //        V----------------------------------V
    //##INFO=<ID=a,Number=b,Type=c,Description="d">
    String inner = m.group(1).trim();
    while (inner.length() > 0) {
      Matcher keyValMatcher = KEY_VALUE_WITH_QUOTES_PATTERN.matcher(inner);
      if (keyValMatcher.find()) {
        final String key = keyValMatcher.group(1).trim();
        String val = keyValMatcher.group(2);
        if (val != null) {
          val = StringUtils.removeBackslashEscapes(val.trim());
        }
        ret.put(key, val);
      } else {
        keyValMatcher = KEY_VALUE_PATTERN.matcher(inner);
        if (keyValMatcher.find()) {
          final String key = keyValMatcher.group(1).trim();
          ret.put(key, keyValMatcher.group(2).trim());
        } else {
          throw new NoTalkbackSlimException("Could not parse header line: " + line);
        }
      }
      inner = inner.substring(keyValMatcher.end());
    }
    return ret;
  }

  static void checkRequiredMetaKeys(Map<String, String> provided, String line, String... required) {
    for (String key : required) {
      if (!provided.containsKey(key)) {
        throw new NoTalkbackSlimException("VCF header missing " + key + " declaration on line: " + line);
      }
    }
  }

}
