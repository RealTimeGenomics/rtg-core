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
package com.rtg.visualization;

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

import com.rtg.bed.BedRecord;
import com.rtg.launcher.AbstractCli;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reader.SourceTemplateReadWriter;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamUtils;
import com.rtg.simulation.MutatedSampleOffsets;
import com.rtg.simulation.SimulatedReadNameParser;
import com.rtg.simulation.SimulatedReadNameParserFactory;
import com.rtg.util.Pair;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * This class take the read data and lays it out in <code>ASCII</code> mode
 *
 */
public final class Aview extends AbstractCli {

  private static final String LOCATION_LABEL = "LOC";
  private static final String REFERENCE_LABEL = "REF";
  private static final String READ_LABEL = "READ";
  private static final String UNMAPPED_LABEL = "UNMP";

  private static final int FAILED_BG = DisplayHelper.RED;
  private static final int GEN_BG = DisplayHelper.GREEN;
  private static final int CALL_BG = DisplayHelper.YELLOW;
  private static final int BED_BG = DisplayHelper.BLUE;


  private PrintStream mOut = null;
  private AviewParams mParams = null;

  private SamAssistance mCgLegacyAssistant;
  private SamAssistance mCgAssistant;
  private SamAssistance mSimpleAssistant;
  private DisplayHelper mDisplayHelper;

  private AviewModel mModel = null;

  private List<SequencesReader> mSdfs = null;
  private Map<Pair<SdfId, PrereadArm>, SequencesReader> mGuidToSdfs = null;
  private Map<SdfId, SdfId[]> mSimulatorParents = null;
  private Map<SdfId, SdfId> mSimulatorOriginalReferences = null;

  // Some view variables
  private int mScreenSelectionStart = 0;
  private int mScreenSelectionEnd = 0;

  // Some state used during parsing
  private Map<String, SimulatedReadNameParser> mParsers = null;
  private boolean mParserAvailable = false;
  private SequencesReader mCurrentReadSdf = null;
  private SdfId[] mCurrentReadSimulatorParents = null;
  private SdfId mCurrentReadSimulatorOriginalReference = null;
  private boolean mAlternateBackgrounds = false;
  private int mAlternateCount = 0;

  private boolean mWarnedTemplate = false;


  @Override
  public String moduleName() {
    return "aview";
  }

  @Override
  public String description() {
    return "ASCII read mapping and variant viewer";
  }

  @Override
  protected void initFlags() {
    AviewParams.initFlags(mFlags);
  }

  private void warnInvalidTemplate() {
    if (!mWarnedTemplate) {
      printOnScreen(mDisplayHelper.decorateBackground("WARNING: Supplied template does not match the template used during mapping", DisplayHelper.RED));
      mWarnedTemplate = true;
    }
  }

  private void process() throws IOException, BadSuperCigarException {
    mDisplayHelper = mParams.displayHelper();
    openSdfs(mParams);
    try {
      String ref = DnaUtils.bytesToSequenceIncCG(mModel.template());
      final SdfId refGuid = ReaderUtils.getSdfId(mParams.referenceFile());
      final int zeroBasedStart = mModel.zeroBasedStart();
      final int zeroBasedEnd = mModel.zeroBasedEnd();

      printOnScreen(mDisplayHelper.header());
      final String region = mParams.sequenceName() + ":" + mParams.start() + ((mParams.end() - 1) == mParams.start() ? "" : "-" + mParams.end());
      printOnScreen(mDisplayHelper.decorateLabel(LOCATION_LABEL) + region);
      printCoordIndicators(ref);

      ref = expandReference(ref);
      printOnScreen(mDisplayHelper.decorateLabel(REFERENCE_LABEL), mDisplayHelper.decorateWithHighlight(ref, null, 0, mParams.colorBases()), "");

      mAlternateBackgrounds = true;
      final boolean[] highlightMask = new boolean[ref.length() + 1]; // Set true at each display position that needs markup during read display
      final int highlightBg = printTracks(zeroBasedStart, zeroBasedEnd, highlightMask);
      final String prettyRef = mDisplayHelper.decorateWithHighlight(ref, highlightMask, highlightBg, mParams.colorBases());
      int printHeaderCount = 0;
      for (final SAMRecord r : mModel.records()) {
        final SdfId mappingRefGuid = SamUtils.getReferenceGuid(r.getHeader());
        if (!refGuid.check(mappingRefGuid)) {
          warnInvalidTemplate();
        }
        if ((mParams.headerLineRepeat() > 0) && (printHeaderCount == mParams.headerLineRepeat())) {
          printOnScreen(mDisplayHelper.decorateLabel(REFERENCE_LABEL), prettyRef, "");
          printHeaderCount = 0;
        }
        ++printHeaderCount;
        //System.err.println(ref + "\n" + Arrays.toString(inserts) + " " + r.getAlignmentStart());

        final int screenReadStart = r.getAlignmentStart() - zeroBasedStart;
        final int readStart = screenReadStart + getInsertsBetween(0, screenReadStart) - 1;
        //System.err.println("sRS=" + screenReadStart + " r.gAS=" + r.getAlignmentStart() + " zBS=" + zeroBasedStart);

        setCurrentReadSdf(r);
        final String readName = getReadName(r);
        boolean evaluated = false;
        boolean correct = false;
        long diff = 0;
        setParsers(r, readName);
        if (mParserAvailable) {
          if (readName == null) { // Given parser could not process the read
            throw new NoTalkbackSlimException("Mixing renamed and unrenamed simulated reads but no reads SDF was provided");
          }
          final SimulatedReadNameParser parser = getParser(r);
          if (parser != null) {
            parser.setReadInfo(readName, r.getCigar().getReferenceLength());
            // Only evaluate if the template GUID matches that from which the read was generated
            if (mCurrentReadSimulatorOriginalReference != null && mParams.baselineFile() != null && mCurrentReadSimulatorOriginalReference.check(refGuid)
                || mCurrentReadSimulatorParents != null && mCurrentReadSimulatorParents[parser.templateSet()].check(refGuid)) {
              evaluated = true;
              if (parser.templateName().equals(r.getReferenceName())) {
                diff = Math.abs(parser.templatePosition() - r.getAlignmentStart());
                if (diff <= mParams.mappingTolerance()) {
                  correct = true;
                }
              }
            }
          }
        }
        final String[] rawRead;
        try {
          if (mParams.unflattenCgi() && r.hasAttribute(SamUtils.CG_SUPER_CIGAR)) {
            rawRead = mCgAssistant.samToReads(r, ref, mModel.template(), readStart, mParams.displayDots(), mParams.showSoftClippedBases());
          } else if (mParams.unflattenCgi() && r.hasAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS)) {
            rawRead = mCgLegacyAssistant.samToReads(r, ref, mModel.template(), readStart, mParams.displayDots(), mParams.showSoftClippedBases());
          } else if (r.getReadUnmappedFlag()) {
            final int bases = Math.min(ref.length() - readStart, r.getReadLength());
            rawRead = new String[]{mDisplayHelper.getSpaces(readStart) + r.getReadString().substring(0, bases)};
          } else {
            // System.err.println(mModel.template());
            rawRead = mSimpleAssistant.samToReads(r, ref, mModel.template(), readStart, mParams.displayDots(), mParams.showSoftClippedBases());
          }
        } catch (IllegalStateException e) {
          printOnScreen(mDisplayHelper.decorateBackground("WARNING: Cannot display alignment (" + e.getMessage() + "): " + r.getSAMString(), DisplayHelper.RED));
          continue;
        }

        boolean first = true;
        for (final String str : rawRead) {
          final String readLabel = first ? mDisplayHelper.decorateLabel(r.getReadUnmappedFlag() ? UNMAPPED_LABEL : READ_LABEL) : mDisplayHelper.getSpaces(DisplayHelper.LABEL_LENGTH + 1);
          final String descr = first ? getReadDesc(r, readName, evaluated, correct, diff) : "";
          String readStr = str + mDisplayHelper.getSpaces(ref.length() - str.length());
          if ((r.getReadPairedFlag() && r.getProperPairFlag())
              || (!r.getReadPairedFlag() && !r.getReadUnmappedFlag())) {
            readStr = mDisplayHelper.decorateBold(readStr);
          }
          readStr = mDisplayHelper.decorateWithHighlight(readStr, highlightMask, highlightBg, mParams.colorBases());
          printOnScreen(readLabel, readStr, descr);
          first = false;
        }
      }
      if (mParams.unmappedFiles() != null) {
        processUnmapped(mParams.unmappedFiles(), mParams.sequenceName(), zeroBasedStart, zeroBasedEnd);
      }
      printOnScreen(mDisplayHelper.footer());
    } finally {
      closeSdfs();
    }
  }

  private String getReadDesc(SAMRecord r, String readName, boolean evaluated, boolean correct, long diff) {
    final String nh = SamUtils.getNHOrIH(r) != null ? " NH:" + SamUtils.getNHOrIH(r) : "";
    final String as = r.hasAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE) ? " AS:" + r.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE) : "";
    final String mated;
    final char direction = r.getReadNegativeStrandFlag() ? '<' : '>';
    char armCh = ' ';
    if (r.getReadPairedFlag()) {
      mated = r.getReadUnmappedFlag() ? " UM" : r.getProperPairFlag() ? " MA" : " UN";
      armCh = r.getFirstOfPairFlag() ? '1' : '2';  // Avoid 'F', confusable with Forward
    } else {
      mated = r.getReadUnmappedFlag() ? " UM" : " SE";
    }
    final StringBuilder extraSb = new StringBuilder();
    if (evaluated) {
      if (correct) {
        extraSb.append(mDisplayHelper.decorateForeground("  =", DisplayHelper.GREEN));
      } else {
        final String diffStr = (diff < 100) ? (diff < 10 ? " " : "") + diff : ">>";
        extraSb.append(mDisplayHelper.decorateForeground(mDisplayHelper.escape(diffStr) + "X", DisplayHelper.RED));
      }
    }
    extraSb.append(" ").append(mDisplayHelper.escape(Character.toString(direction) + Character.toString(armCh) + mated + nh + as));
    if (mParams.printCigars()) {
      final String superCigar = r.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
      extraSb.append(" ").append(superCigar == null ? r.getCigarString() : superCigar);
    }
    if (mParams.printReadName()) {
      extraSb.append(" ").append(readName != null ? readName : r.getReadName());
    }
    if (mParams.printReadGroup()) {
      final String rg = r.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE);
      if (rg != null) {
        extraSb.append(" RG:").append(rg);
      }
    }
    if (mParams.printSample()) {
      final String sample = r.getReadGroup() == null ? null : r.getReadGroup().getSample();
      if (sample != null) {
        extraSb.append(" SM:").append(sample);
      }
    }
    if (mParams.printMapQ()) {
      extraSb.append(" ").append(r.getMappingQuality());
    }
    if (mParams.printMatePosition()) {
      extraSb.append(" ").append(r.getMateAlignmentStart());
    }
    return extraSb.toString();
  }

  private SimulatedReadNameParser getParser(final SAMRecord r) {
    final SimulatedReadNameParser parser;
    if (mCurrentReadSimulatorOriginalReference != null && mParams.baselineFile() != null) {
      final String sample = r.getReadGroup() == null ? null : r.getReadGroup().getSample();
      parser = mParsers.get(sample);
    } else {
      parser = mParsers.get(null);
    }
    return parser;
  }

  private void setParsers(SAMRecord record, final String readName) throws IOException {
    if (mParsers == null && readName != null) {
      mParsers = new HashMap<>();
      if (mCurrentReadSimulatorOriginalReference != null && mParams.baselineFile() != null) {
        String[] samples = mParams.wantedSamples();
        if (samples == null && (record.getReadGroup() != null)) { // None explicitly specified, use the sample given in the record
          samples = new String[] {record.getReadGroup().getSample() };
        }
        mParserAvailable = false;
        if (samples == null) {
          Diagnostic.warning("Unable to evaluate correctness of simulated mapping: no sample name present in read group.");
        } else {
          final MutatedSampleOffsets[] offsets = MutatedSampleOffsets.getOffsets(mParams.baselineFile(), new RegionRestriction(mParams.sequenceName(), 0, mParams.end()), samples);
          for (int i = 0; i < offsets.length; ++i) {
            if (offsets[i] == null) {
              Diagnostic.warning("Unable to evaluate correctness of sample \"" + samples[i] + "\" due to overlapping variants in generated VCF.");
            }
          }
          final SimulatedReadNameParser[] parsers = SimulatedReadNameParserFactory.getParsers(readName, offsets);
          if (parsers.length == samples.length) {
            for (int i = 0; i < parsers.length; ++i) {
              if (parsers[i] != null) {
                mParsers.put(samples[i], parsers[i]);
                mParserAvailable = true;
              }
            }
          }
        }
      } else {
        final SimulatedReadNameParser parser = SimulatedReadNameParserFactory.getParser(readName);
        mParserAvailable = parser != null;
        if (mParserAvailable) {
          mParsers.put(null, parser);
        }
      }
    }
  }

  private int printTracks(final int zeroBasedStart, final int zeroBasedEnd, final boolean[] highlightMask) {
    final int projectTrack = mParams.projectTrackId();
    int trackId = 1;
    int highlightBg = CALL_BG;
    for (final AviewModel.AviewTrack track : mModel.tracks()) {
      if (track instanceof AviewModel.SnpSet) {
        final AviewModel.SnpSet called = (AviewModel.SnpSet) track;
        for (final Map.Entry<String, ArrayList<AviewVariant>> sampleVars : called.sampleSnps().entrySet()) {
          final boolean project = projectTrack <= 0 || trackId == projectTrack;
          final boolean[] hl = project ? highlightMask : null;
          final boolean isBaseline = track instanceof AviewModel.BaselineSet;
          if (project && projectTrack > 0) {
            highlightBg = isBaseline ? GEN_BG : CALL_BG;
          }
          final String sampleName = sampleVars.getKey();
          printVariants(sampleVars.getValue(), zeroBasedStart, zeroBasedEnd, isBaseline, called.type(), "(" + trackId + ") SM:" + sampleName + " " + called.name(), hl);
          ++trackId;
        }
      } else if (track instanceof AviewModel.BedSet) {
        final AviewModel.BedSet bed = (AviewModel.BedSet) track;
        final boolean project = projectTrack <= 0 || trackId == projectTrack;
        final boolean[] hl = project ? highlightMask : null;
        if (project && projectTrack > 0) {
          highlightBg = BED_BG;
        }
        printBeds(bed.bedRecords(), zeroBasedStart, zeroBasedEnd, bed.type(), "(" + trackId + ") " + bed.name(), hl);
        ++trackId;
      }
    }
    return highlightBg; // Return the background color used for highlighting
  }

  private void printCoordIndicators(final String ref) {
    final String label = mDisplayHelper.getSpaces(DisplayHelper.LABEL_LENGTH + 1);
    final int selectionStart = mParams.start() - mModel.oneBasedStart();
    mScreenSelectionStart = getInsertsBetween(0, selectionStart) + selectionStart;

    final String startSpaces = mDisplayHelper.getSpaces(mScreenSelectionStart);
    final int difference = Math.abs(mParams.end() - mParams.start());
    final String indicatorStr;
    if (difference > 1) {
      final int selectionLength = Math.min(mParams.end()  - mParams.start() - 1, mModel.inserts().length - selectionStart); //if they've asked for more region than available, just give what we've got.
      final int newSelectionLength = getInsertsBetween(selectionStart, selectionStart + selectionLength) + selectionLength;
      mScreenSelectionEnd = mScreenSelectionStart + newSelectionLength + 1;

      final String midSpaces = mDisplayHelper.decorateUnderline(mDisplayHelper.getSpaces(newSelectionLength - 1));
      indicatorStr = startSpaces + '|' + midSpaces + '|';
    } else {
      mScreenSelectionEnd = mScreenSelectionStart + 1;
      indicatorStr = startSpaces + '|';
    }

    final int digits = (int) (Math.log(ref.length()) / Math.log(10)) + 1;
    final StringBuffer[] posStrs = new StringBuffer[digits];
    final char[] lastCol = new char[digits];
    final String initPos = String.valueOf(mModel.oneBasedStart() - 1);
    for (int i = 0; i < digits; ++i) {
      posStrs[i] = new StringBuffer();
      final int j = initPos.length() - i - 1;
      lastCol[i] = (j >= 0) ? initPos.charAt(j) : ' ';
    }
    for (int i = mModel.oneBasedStart(); i < mModel.oneBasedStart() + ref.length(); ++i) {
      final String pos = String.valueOf(i);
      for (int d = 0; d < digits; ++d) {
        final int j = pos.length() - d - 1;
        final char c = (j >= 0) ? pos.charAt(j) : ' ';
        posStrs[d].append(c != lastCol[d] ? c : ' ');
        lastCol[d] = c;
      }
    }
    for (int i = digits - 1; i >= 0; --i) {
      printOnScreen(label, mDisplayHelper.decorateBackground(expandReference(posStrs[i].toString()), DisplayHelper.WHITE_PLUS), "");
    }
    printOnScreen(label, indicatorStr, "");
  }

  /**
   *
   * @param ref the original reference
   * @return a string with underscores at the insertion positions
   */
  protected String expandReference(String ref) {
    return expandReference(mModel.inserts(), ref, mDisplayHelper);
  }
  static String expandReference(final int[] inserts, String ref, DisplayHelper displayHelper) {
    String localref = ref;
    int n = 0;
    for (int i = 0; i < inserts.length; ++i) {
      if (inserts[i] > 0) {
        localref = localref.substring(0, i + n) + displayHelper.getInserts(inserts[i]) + localref.substring(i + n);
        n += inserts[i];
      }
    }
    //System.out.println("inserts: " + Arrays.toString(inserts));
    //System.out.println("reference: " + ref);
    //System.out.println("result: " + localref);
    return localref;
  }

  private void closeSdfs() throws IOException {
    for (final SequencesReader r : mSdfs) {
      r.close();
    }
  }

  private void openSdfs(final AviewParams params) throws IOException {
    mSdfs = new ArrayList<>();
    mSimulatorParents = new HashMap<>();
    mSimulatorOriginalReferences = new HashMap<>();
    mGuidToSdfs = new HashMap<>();
    if (params.readFiles() != null) {
      for (final File f : params.readFiles()) {
        SequencesReader s;
        if (ReaderUtils.isPairedEndDirectory(f)) {
          s = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(ReaderUtils.getLeftEnd(f));
          mSdfs.add(s);
          mGuidToSdfs.put(new Pair<>(s.getSdfId(), s.getArm()), s);

          s = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(ReaderUtils.getRightEnd(f));
          mGuidToSdfs.put(new Pair<>(s.getSdfId(), s.getArm()), s);
          mSdfs.add(s);
        } else {
          s = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(f);
          mGuidToSdfs.put(new Pair<>(s.getSdfId(), s.getArm()), s);
          mSdfs.add(s);
        }
        SdfId[] templateMap = SourceTemplateReadWriter.readTemplateMap(s.path());
        if (templateMap == null) {
          templateMap = new SdfId[] {
              new SdfId(0L)
          };
        }
        mSimulatorParents.put(s.getSdfId(), templateMap);
        mSimulatorOriginalReferences.put(s.getSdfId(), SourceTemplateReadWriter.readMutationMap(s.path()));
      }
    }
  }

  private SequencesReader getReadSdf(SdfId guid, PrereadArm arm) {
    return mGuidToSdfs.get(new Pair<>(guid, arm));
  }

  private void setCurrentReadSdf(SAMRecord r) {
    final SdfId guid = SamUtils.getReadsGuid(r.getHeader());
    if (guid.available()) {
      final PrereadArm arm = !r.getReadPairedFlag() ? PrereadArm.UNKNOWN : r.getFirstOfPairFlag() ? PrereadArm.LEFT : PrereadArm.RIGHT;
      mCurrentReadSdf = getReadSdf(guid, arm);
    } else {
      mCurrentReadSdf = null;
    }
    final SdfId currentReadSdfGuid = mCurrentReadSdf == null ? new SdfId(0) : mCurrentReadSdf.getSdfId();
    mCurrentReadSimulatorParents = !currentReadSdfGuid.available() ? null : mSimulatorParents.get(currentReadSdfGuid);
    mCurrentReadSimulatorOriginalReference = !currentReadSdfGuid.available() ? null : mSimulatorOriginalReferences.get(currentReadSdfGuid);
  }

  private String getReadName(SAMRecord r) throws IOException {
    try {
      final long readId = Long.parseLong(r.getReadName());
      // id looks like an actual ID, so look up it's name from the SDF. We can also
      if (mCurrentReadSdf != null) {
        return mCurrentReadSdf.fullName(readId);
      }
      return null;
    } catch (final NumberFormatException nfe) {
      return r.getReadName(); // name given in SAM record wasn't a long, assume it's already textual. We should probably fail here, since we can't validate sdf ids fully
    }
  }

  private void processUnmapped(File[] files, String sequenceName, final int correctStart, final int correctEnd) throws IOException {
    final Set<String> unmappedLines = new HashSet<>();
    for (final File file : files) {
      final SamReader unmappedReader = SamUtils.makeSamReader(file);
      for (final SAMRecord r : unmappedReader) {
        if (r.getReadUnmappedFlag()) {
          setCurrentReadSdf(r);
          final String readName = getReadName(r);
          setParsers(r, readName);
          if (mParserAvailable) {
            final SimulatedReadNameParser parser = getParser(r);
            if (parser != null) {
              parser.setReadInfo(readName, r.getCigar().getReferenceLength());
              if (parser.templateName().equals(sequenceName)) {
                final long pos = parser.templatePosition();
                if (pos >= correctStart + 1 && pos + r.getReadLength() <= correctEnd - 1) {
                  // This doesn't do anything special with CG read structure
                  final String read = SamHelper.getCorrectComplement(parser.forwardFrame(), r.getReadString());
                  final int end = (int) (pos - correctStart - 1);
                  final String seq = mDisplayHelper.getSpaces(end + getInsertsBetween(0, end)) + read + mDisplayHelper.getSpaces(mModel.inserts().length - (end + read.length()) + 3);
                  if (!unmappedLines.contains(seq + readName)) {
                    unmappedLines.add(seq + readName);
                    printOnScreen(mDisplayHelper.decorateLabel(UNMAPPED_LABEL), seq, readName);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  private void printBeds(final ArrayList<BedRecord> records, final int correctStart, final int correctEnd, final String label, String endLabel, boolean[] highlightMask) {
    if (records != null) {
      Integer charPosition = 0; // Current character position (i.e. including inserts)
      int refPosition = correctStart; // Reference position we are currently up to, zero-based. (correctStart, correctEnd are zero-based)
      int residualInsertLength = 0;
      final StringBuilder line = new StringBuilder();
      for (final BedRecord record : records) {
        //System.err.println("snp for SM:" + sample + " at:" + line.getLocation() + " refPosition:" + refPosition + " charPosition:" + charPosition + " residual:" + residualInsertLength + " ref:" + line.getReference() + " pred:" + line.getPrediction());
        if (record.getEnd() > correctStart && record.getStart() < correctEnd) {

          char startChar = '(';
          final char endChar = ')';

          // First output spaces to bring us to the correct start position
          int newRefPosition = record.getStart();
          int refLength = record.getLength();
          if (record.getStart() < correctStart) { // If we are a region coming in from the left, truncate it at correctStart
            newRefPosition = correctStart;
            startChar = '*';
            refLength = refLength - (correctStart - record.getStart());
          }
          final int numInserts = getInsertsBetween(refPosition - correctStart, newRefPosition - correctStart); // From one past ref position to new position
          final int numspaces = newRefPosition - refPosition + numInserts + residualInsertLength;
          if (numspaces < 0) {
            // Overlapping calls
            printOnScreen(mDisplayHelper.decorateBackground("WARNING: Overlapping BED regions between position " + (newRefPosition + 1) + " and " + (refPosition + 1) + ". Skipping last.", DisplayHelper.RED));
            continue;
          }
          final String spaces = mDisplayHelper.getSpaces(numspaces);
          line.append(spaces);
          charPosition += spaces.length();
          refPosition = newRefPosition;
          // At this point the output line is complete to the position of the call. charPostition and refPosition are also correct.


          newRefPosition = Math.min(refPosition + refLength, correctEnd);
          // Count displayed inserts between start and end of displayed call, including insert at the terminal coordinate
          final int displayedRefInserts = getInsertsBetween(refPosition - correctStart, Math.min(newRefPosition, correctEnd) - correctStart); // From one past ref position to new position. Also grab adjacent insert
          final int refDisplayLength = refLength + displayedRefInserts;

          String name;
          // If the call is shorter than the actual number of reference bases, illustrate this by extending to the number of reference bases
          if (refDisplayLength <= 1) {
            name = "O";
          } else {
            final StringBuilder sb = new StringBuilder().append(startChar);
            final String regionName = record.getAnnotations().length == 0 ? "" : record.getAnnotations()[0];
            sb.append(regionName.length() > (refDisplayLength - 2) ? regionName.substring(0, refDisplayLength - 2) : regionName);
            while (sb.length() < refDisplayLength - 1) {
              sb.append('-');
            }
            sb.append(endChar);
            name = sb.toString();
          }

          final int maxPermittedCallLength = correctEnd - refPosition + displayedRefInserts;
          if (name.length() > maxPermittedCallLength) { // Handle the case where a long call would exceed the length of the display by truncating
            name = name.substring(0, maxPermittedCallLength - 1) + "*";
            residualInsertLength = 0;
          } else {
            residualInsertLength = refDisplayLength - name.length();
          }

          // Now output the region itself.
          line.append(mDisplayHelper.decorate(name, DisplayHelper.WHITE, BED_BG));
          if (highlightMask != null) {
            for (int i = 0; i < name.length() && (charPosition + i) < highlightMask.length; ++i) {
              highlightMask[charPosition + i] = true;
            }
          }
          charPosition += name.length();

          refPosition = newRefPosition;
        }
      }
      String newEndLabel = "";
      if (endLabel != null) {
        final int numInserts = getInsertsBetween(refPosition - correctStart, correctEnd - correctStart); // From one past ref position to new position
        final int numspaces = correctEnd - refPosition + numInserts + residualInsertLength;
        if (numspaces > 0) {
          final String spaces = mDisplayHelper.getSpaces(numspaces);
          line.append(spaces);
        }
        newEndLabel = " " + endLabel;
      }
      printOnScreen(mDisplayHelper.decorateLabel(label), line.toString(), newEndLabel);
    }
  }


  private void printVariants(final ArrayList<AviewVariant> snpArray, final int correctStart, final int correctEnd, boolean baseline, final String label, String endLabel, boolean[] highlightMask) {
    if (snpArray != null) {
      Integer charPosition = 0; // Current character position (i.e. including inserts)
      int refPosition = correctStart; // Reference position we are currently up to, zero-based. (correctStart, correctEnd are zero-based)
      int residualInsertLength = 0;
      final StringBuilder strand1 = new StringBuilder();
      final StringBuilder strand2 = new StringBuilder();
      for (final AviewVariant var : snpArray) {
        //System.err.println("snp for SM:" + sample + " at:" + line.getLocation() + " refPosition:" + refPosition + " charPosition:" + charPosition + " residual:" + residualInsertLength + " ref:" + line.getReference() + " pred:" + line.getPrediction());
        final int vStart = var.getPosition() - 1;
        if ((vStart + var.referenceLength()) >= correctStart && vStart < correctEnd) {

          // First output spaces to bring us to the correct start position
          if (vStart < correctStart) { // If we have a region coming in from the left, ignore it but warn if it would intersect our region of interest
            if (vStart + var.referenceLength() > mModel.zeroBasedStart()) {
              printOnScreen(mDisplayHelper.decorateBackground("WARNING: Cannot display variant starting outside region at " + (vStart + 1) + ". Skipping.", DisplayHelper.RED));
              continue;
            }
          }
          int newRefPosition = vStart;
          final int numInserts = getInsertsBetween(refPosition - correctStart, newRefPosition - correctStart); // From one past ref position to new position
          final int numspaces = newRefPosition - refPosition + numInserts + residualInsertLength;
          if (numspaces < 0) {
            // Overlapping calls
            printOnScreen(mDisplayHelper.decorateBackground("WARNING: Overlapping calls between position " + (newRefPosition + 1) + " and " + (refPosition + 1) + ". Skipping last.", DisplayHelper.RED));
            continue;
          }
          final String spaces = mDisplayHelper.getSpaces(numspaces);
          strand1.append(spaces);
          strand2.append(spaces);
          charPosition += spaces.length();
          refPosition = newRefPosition;
          // At this point each strand is complete to the position of the call. charPostition and refPosition are also correct.


          final boolean hetCall = var.ntAlleleB() != null;
          String firstCall = DnaUtils.bytesToSequenceIncCG(var.ntAlleleA());
          String secondCall = hetCall ? DnaUtils.bytesToSequenceIncCG(var.ntAlleleB()) : firstCall;
          final int lengthDiff = Math.abs(firstCall.length() - secondCall.length());
          // Make displayed representation of each side of the call the same length
          if (hetCall) {
            //if (!line.getReference().contains("i") && firstCall.length() > 0 && secondCall.length() > 0 && lengthDiff > 0) {
            if (lengthDiff > 0) {
              if (firstCall.length() < secondCall.length()) {
                firstCall += mDisplayHelper.getInserts(lengthDiff);
              } else {
                secondCall += mDisplayHelper.getInserts(lengthDiff);
              }
            }
          }

          assert firstCall.length() == secondCall.length();

          final int refLength = var.referenceLength();
          newRefPosition = Math.min(refPosition + refLength, correctEnd);
          // Count displayed inserts between start and end of displayed call
          final int displayedRefInserts = getInsertsBetween(refPosition - correctStart,  Math.min(newRefPosition, correctEnd) - correctStart); // From one past ref position to new position (also grab subsequent insert for use if need be)
          final int refDisplayLength = refLength + displayedRefInserts;

          // If the call is shorter than the actual number of reference bases, illustrate this by extending to the number of reference bases
          if (refDisplayLength > firstCall.length()) {
            final String extra = mDisplayHelper.getInserts(refDisplayLength - firstCall.length());
            firstCall += extra;
            secondCall += extra;
          }

          final int maxPermittedCallLength = correctEnd - refPosition + displayedRefInserts;
          if (firstCall.length() > maxPermittedCallLength) { // Handle the case where a long call would exceed the length of the display by truncating
            firstCall = firstCall.substring(0, maxPermittedCallLength - 1) + "*";
            secondCall = secondCall.substring(0, maxPermittedCallLength - 1) + "*";
            residualInsertLength = 0;
          } else {
            residualInsertLength = refDisplayLength - firstCall.length();
          }

          // Now output the call itself.
          final boolean failed = var.isFiltered();
          final int bgCol = failed ? FAILED_BG : baseline ? GEN_BG : CALL_BG;
          strand1.append(mDisplayHelper.decorateBackground(firstCall, bgCol));
          strand2.append(mDisplayHelper.decorateBackground(secondCall, bgCol));
          if (highlightMask != null) {
            for (int i = 0; i < firstCall.length(); ++i) {
              highlightMask[charPosition + i] = true;
            }
          }
          charPosition += firstCall.length();

          refPosition = newRefPosition;
        }
      }
      String newEndLabel = "";
      if (endLabel != null) {
        final int numInserts = getInsertsBetween(refPosition - correctStart, correctEnd - correctStart); // From one past ref position to new position
        final int numspaces = correctEnd - refPosition + numInserts + residualInsertLength;
        if (numspaces > 0) {
          final String spaces = mDisplayHelper.getSpaces(numspaces);
          strand1.append(spaces);
          strand2.append(spaces);
        }
        newEndLabel = " " + endLabel;
      }

      printOnScreen(mDisplayHelper.decorateLabel(label), mDisplayHelper.decorateWithHighlight(strand1.toString(), null, 0, mParams.colorBases()), newEndLabel);
      printOnScreen(mDisplayHelper.decorateLabel(label), mDisplayHelper.decorateWithHighlight(strand2.toString(), null, 0, mParams.colorBases()), newEndLabel);
    }
  }

  private void printOnScreen(final String prefix, final String sequence, String suffix) {
    mOut.print(prefix);
    String clipped = clipSequence(sequence);
    if (mAlternateBackgrounds) {
      if (mDisplayHelper.supportsNesting() && (mAlternateCount++ & 1) == 1) {
        clipped = mDisplayHelper.decorateBackground(clipped, DisplayHelper.WHITE_PLUS);
      }
    }
    mOut.print(clipped);
    mOut.println(suffix);
  }

  private void printOnScreen(final String string) {
    if (string != null) {
      mOut.println(string);
    }
  }

  private String clipSequence(String sequence) {
    final int padding = mParams.regionPadding();
    if (padding > 0) {
      final int clipStart = mScreenSelectionStart - padding;
      final int clipEnd = mScreenSelectionEnd + padding;
      return mDisplayHelper.clipSequence(sequence, clipStart, clipEnd);
    } else {
      return sequence;
    }
  }

  private int getInsertsBetween(final int start, final int end) {
    final int[] inserts = mModel.inserts();
    int n = 0;
    for (int i = start; i < end; ++i) {
      n += inserts[i];
    }
    return n;
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    mOut = new PrintStream(out);
    try {
      mParams = AviewParams.makeParams(mFlags);
      initSamAssistants();
      mModel = new AviewModel(mParams);
      try {
        process();
      } catch (final BadSuperCigarException e) {
        throw new IOException(e);
      }
      return 0;
    } finally {
      mOut.flush();
    }
  }

  private void initSamAssistants() {
    mCgLegacyAssistant = new SamAssistanceCgLegacy();
    mCgAssistant = new SamAssistanceCg();
    mSimpleAssistant = new SamAssistanceSimple();
  }

}
