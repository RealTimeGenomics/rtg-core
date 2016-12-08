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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.bed.BedUtils;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.sam.SamUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.vcf.VcfUtils;

import htsjdk.samtools.SAMRecord;

/**
 */
final class AviewModel {

  private static final String GENERATED_LABEL = "GEN"; // Depending on how you feel, interpret as any one of: Generated, or Genotype, or Genome
  private static final String CALLED_LABEL = "CALL";
  private static final String BED_LABEL = "BED";

  interface AviewTrack {

    /**
     * A label for the type of track, appears on the left hand side of the display
     * @return the track label
     */
    String type();

    /**
     * The name of the track, appears on the right hand side of the display
     * @return the track name
     */
    String name();
  }

  static class SnpSet implements AviewTrack {
    private final String mType;
    private final String mSource;
    // Key = sample name, value is all variants for the sample
    private final LinkedHashMap<String, ArrayList<AviewVariant>> mSampleSnps;
    protected SnpSet(String type, String source, LinkedHashMap<String, ArrayList<AviewVariant>> snps) {
      mType = type;
      mSource = source;
      mSampleSnps = snps;
    }
    @Override
    public String type() {
      return mType;
    }
    @Override
    public String name() {
      return mSource;
    }
    public LinkedHashMap<String, ArrayList<AviewVariant>> sampleSnps() {
      return mSampleSnps;
    }
  }
  static class BaselineSet extends SnpSet {
    BaselineSet(String source, LinkedHashMap<String, ArrayList<AviewVariant>> snps) {
      super(GENERATED_LABEL, source, snps);
    }
  }
  static class CallSet extends SnpSet {
    CallSet(String source, LinkedHashMap<String, ArrayList<AviewVariant>> snps) {
      super(CALLED_LABEL, source, snps);
    }
  }
  static class BedSet implements AviewTrack {
    private final String mSource;
    private final ArrayList<BedRecord> mBedRecords;
    BedSet(String source, ArrayList<BedRecord> records) {
      mSource = source;
      mBedRecords = records;
    }
    @Override
    public String type() {
      return BED_LABEL;
    }
    @Override
    public String name() {
      return mSource;
    }
    public ArrayList<BedRecord> bedRecords() {
      return mBedRecords;
    }
  }

  private final byte[] mTemplate;

  private final ArrayList<AviewTrack> mTracks = new ArrayList<>();
  private final ArrayList<SAMRecord> mRecords;

  private final int[] mInserts; // Contains count of inserts at each reference position. Index is relative to mOneBasedStart
  private final int mOneBasedStart;

  private final int mZeroBasedStart;
  private final int mZeroBasedEnd; // Exclusive, dammit


  AviewModel(AviewParams p) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(p.referenceFile())) {
      final Map<String, Long> names = ReaderUtils.getSequenceNameMap(reader);
      final Long sequenceId = names.get(p.sequenceName());
      if (sequenceId == null) {
        throw new NoTalkbackSlimException("Given sequence name not found : " + p.sequenceName());
      }
      final long seqId = sequenceId;
      final int len = reader.length(seqId);
      if (p.start() > len) {
        throw new NoTalkbackSlimException("start is greater than template length. start = " + p.start() + ", template length = " + len);
      }
      if (p.end() <= p.start()) {
        throw new NoTalkbackSlimException("start is greater than end. start = " + p.start() + ", end = " + p.end());
      }
      mRecords = SamHelper.loadAlignments(p, reader);
      if (p.sortReads()) {
        SamHelper.sortAlignments(mRecords);
      }
      if (p.sortReadGroup()) {
        SamHelper.sortAlignmentsWithReadGroups(mRecords);
      }

      // Scan alignments to auto-expand the range if needed
      RegionRestriction region = expandRegion(mRecords, p.region());

      // Also expand if need be to ensure that we get features within the padded zone
      if (p.regionPadding() > 0) {
        final RegionRestriction paddedInitial = new RegionRestriction(region.getSequenceName(),
          Math.max(0, p.region().getStart() - p.regionPadding()),
          Math.min(len, p.region().getEnd() + p.regionPadding()));
        region = union(region, paddedInitial, 0, len);
      }

      // Grab any variants/beds overlapping the new region
      loadTracks(p, region);

      // Scan snps to auto-expand the range if needed (just to allow complete display of baseline variants near the edges
      for (final AviewTrack track : mTracks) {
        if (track instanceof SnpSet) {
          region = expandRegion((SnpSet) track, region);
        }
      }

      mOneBasedStart = Math.max(1, region.getStart() + 1);
      mZeroBasedStart = Math.max(0, region.getStart());
      mZeroBasedEnd = Math.min(region.getEnd(), len); // Converting from one-based inclusive to zero-based exclusive
      final int tlen = mZeroBasedEnd - mZeroBasedStart;

      mTemplate = new byte[tlen];
      reader.read(seqId, mTemplate, mZeroBasedStart, tlen);

      // Locate positions of inserts based on alignments
      mInserts = new int[tlen];
      for (final SAMRecord r : mRecords) {
        final String superCigar = r.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
        final String cigar = superCigar == null ? r.getCigarString() : superCigar;
        if (!CigarHelper.isMissing(cigar)) {
          CigarHelper.locateInserts(cigar, r.getAlignmentStart() - mOneBasedStart, mInserts);
        }
      }

      // Locate positions of inserts based on variants
      calculateSnpInserts();
    }
  }

  private RegionRestriction union(SequenceNameLocus a, SequenceNameLocus b, int min, int max) {
    return new RegionRestriction(a.getSequenceName(),
      Math.max(min, Math.min(a.getStart(), b.getStart())),
      Math.min(max, Math.max(a.getEnd(), b.getEnd())));
  }

  // Scan records to auto-expand the range if needed
  private static RegionRestriction expandRegion(ArrayList<SAMRecord> records, RegionRestriction region) {
    int minStart = region.getStart() + 1; // to one-based
    int maxEnd = region.getEnd();
    for (final SAMRecord r : records) {
//      System.err.println("r.start=" + r.getAlignmentStart() + " r.cigar=" + r.getCigarString() + " r.alignmentEnd=" + r.getAlignmentEnd() + " r.readLen=" + r.getReadLength());
      if (r.getAlignmentEnd() > maxEnd) {
        maxEnd = r.getAlignmentEnd();
      }
      if (r.getAlignmentStart() < minStart) {
        minStart = r.getAlignmentStart();
      }
    }
    return new RegionRestriction(region.getSequenceName(), minStart - 1, maxEnd);
  }

  // Scan variants to auto-expand the range if needed
  private static RegionRestriction expandRegion(SnpSet snpset, RegionRestriction region) {
    int minStart = region.getStart() + 1; // to one-based
    int maxEnd = region.getEnd();
    if (snpset != null) {
      for (final Map.Entry<String, ArrayList<AviewVariant>> sampleSnps : snpset.sampleSnps().entrySet()) {
        for (final AviewVariant v : sampleSnps.getValue()) {
          if (v.getPosition() < minStart) {
            minStart = v.getPosition();
          }
          final int variantEnd = v.getPosition() + maxPredictionLength(v) - 1; // end pos here is inclusive
          if (variantEnd > maxEnd) {
            maxEnd = variantEnd;
          }
        }
      }
    }
    return new RegionRestriction(region.getSequenceName(), minStart - 1, maxEnd);
  }

  private void loadTracks(AviewParams p, RegionRestriction region) throws IOException {
    // Baseline is always first, if present
    if (p.baselineFile() != null) {
      final Pair<File, String> named = parseNamedFileString(p.baselineFile().getPath());
      final File file = named.getA();
      final String name = named.getB();
      final BaselineSet baseline = new BaselineSet(name, new LinkedHashMap<String, ArrayList<AviewVariant>>());
      mTracks.add(baseline);
      VariantHelper.loadSnpRange(baseline.sampleSnps(), file, region, p.wantedSamples());
    }

    for (final File file : p.trackFiles()) {
      final Pair<File, String> named = parseNamedFileString(file.getPath());
      final File trackFile = named.getA();
      final String trackName = named.getB();
      loadTrack(trackFile, trackName, p, region);
    }
  }

  private void loadTrack(File trackFile, String trackName, AviewParams p, RegionRestriction region) throws IOException {
    final AviewTrack track;
    if (VcfUtils.isVcfExtension(trackFile)) {
      track = loadCallSet(trackFile, trackName, p, region);
    } else if (BedUtils.isBedExtension(trackFile)) {
      track = loadBedSet(trackFile, trackName, region);
    } else {
      throw new NoTalkbackSlimException("Unrecognized track file type for " + trackFile);
    }
    mTracks.add(track);
  }

  private SnpSet loadCallSet(File vcfFile, String vcfName, AviewParams p, RegionRestriction region) throws IOException {
    final SnpSet called = new CallSet(vcfName, new LinkedHashMap<String, ArrayList<AviewVariant>>());
    VariantHelper.loadSnpRange(called.sampleSnps(), vcfFile, region, p.wantedSamples());
    return called;
  }

  private BedSet loadBedSet(File bedFile, String bedName, RegionRestriction region) throws IOException {
    final BedSet bedset = new BedSet(bedName, new ArrayList<BedRecord>());
    final File index = new File(bedFile + TabixIndexer.TABIX_EXTENSION);
    try (BedReader r = BedReader.openBedReader(index.exists() ? region : null, bedFile, 0)) {
      while (r.hasNext()) {
        final BedRecord current = r.next();
        // Explicitly check the position in case this is from an untabixed source.
        if (!region.overlaps(current)) {
          continue;
        }
        //System.err.println("Loaded bed record: " + current.toString());
        bedset.bedRecords().add(current);
      }
    }
    return bedset;
  }

  private static Pair<File, String> parseNamedFileString(String fileName) {
    final String[] split = StringUtils.split(fileName, '=', 2);
    final File file = new File(split[0]);
    final String name;
    if (split.length > 1 && split[1].length() > 0) {
      name = split[1];
    } else {
      name = file.getName();
    }
    return new Pair<>(file, name);
  }

  byte[] template() {
    return mTemplate;
  }

  List<AviewTrack> tracks() {
    return mTracks;
  }

  ArrayList<SAMRecord> records() {
    return mRecords;
  }

  int[] inserts() {
    return mInserts;
  }

  private void calculateSnpInserts() {
    for (final AviewTrack track : mTracks) {
      if (track instanceof SnpSet) {
        for (final Map.Entry<String, ArrayList<AviewVariant>> sampleSnps : ((SnpSet) track).sampleSnps().entrySet()) {
          for (final AviewVariant line : sampleSnps.getValue()) {
            processVariantLine(line);
          }
        }
      }
    }
  }

  private static int maxPredictionLength(AviewVariant v) {
    final int maxPredictionLength;
    if (v.ntAlleleB() == null) {
      maxPredictionLength = v.ntAlleleA().length;
    } else {
      maxPredictionLength = Math.max(v.ntAlleleA().length, v.ntAlleleB().length);
    }
    return maxPredictionLength;
  }

  private void processVariantLine(final AviewVariant line) {
    final int position = line.getPosition() - mOneBasedStart;
    final int referenceLength = line.referenceLength();
    final int callDelta = maxPredictionLength(line) - referenceLength;
    if (callDelta > 0) { // The call represents an insertion of some type
      // Ensure the current position has at least that many inserts:
      //mInserts[position] = Math.max(mInserts[position], callDelta);

      // An alternative, may be better, which looks to see if there
      // are already inserts further down (depends on call reference
      // length) that will account for the inserts needed by this
      // call:
      // Check sum of inserts already stored over the length of reference covered
      final int endposition = Math.min(mInserts.length, position + referenceLength + 1);
      int insertTot = 0;
      for (int i = position; i < endposition; ++i) {
        insertTot += mInserts[i];
      }
      final int insertsNeeded = callDelta - insertTot;
      if (insertsNeeded > 0) { // Existing inserts do not cover what are needed
        mInserts[position] += insertsNeeded;
      }
    }
  }

  /** @return start of region (1 based inclusive) */
  int oneBasedStart() {
    return mOneBasedStart;
  }

  /** @return start of region (0 based inclusive) */
  int zeroBasedStart() {
    return mZeroBasedStart;
  }

  /** @return end of region (0 based exclusive) */
  int zeroBasedEnd() {
    return mZeroBasedEnd;
  }
}
